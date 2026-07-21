"""Source-derived OPD and isotropy-subgroup catalog."""

from __future__ import annotations

from dataclasses import asdict
from fractions import Fraction
from functools import lru_cache
import math
import re
from typing import Any

from ISODISTORT.Assembled.Backend.exactmath import integer_determinant3
from ISODISTORT.Assembled.Backend.isotropy.engine.get_isotropy import (
    generate_dynamic_point_isotropy_row,
    generate_dynamic_point_isotropy_row_magnetic,
    generate_dynamic_isotropy_rows,
    generate_dynamic_isotropy_rows_magnetic,
)
from ISODISTORT.Assembled.Backend.isotropy.engine.dynamic_isotropy_file import (
    DynamicIsotropyRow,
    sort_dynamic_rows_for_file,
)
from ISODISTORT.Assembled.Backend.isotropy.engine.source_data import SourceData
from ISODISTORT.Assembled.Backend.isotropy.engine.id_subgroup_magnetic import (
    id_subgroup_magnetic_identify_with_generator_block,
)
from ISODISTORT.Assembled.Backend.isotropy.lattice_conventions import (
    present_opd_basis_rows,
)
from ISODISTORT.Assembled.Backend.reciprocal.catalog import (
    _apply_fraction_k_matrix,
    _fmt_frac,
    _isodistort_k_records,
    _isodistort_star_vectors,
    _k_slot_info,
    _parse_k_record,
    _reciprocal_pml_to_cinter_matrix,
    _strip,
    evaluate_k_component,
    substitute_k_vector,
)
from ISODISTORT.Assembled.Backend.source.magnetic import data as magnetic_data
from ISODISTORT.Assembled.Backend.source.tables import SOURCE, SourceTables, source_tables


ISODISTORT_K_PARAMS = ("a", "b", "g")


def _fmt_complex(value: complex) -> str:
    re_part, im_part = round(value.real, 8), round(value.imag, 8)
    if abs(im_part) < 1e-8:
        if abs(re_part - round(re_part)) < 1e-8:
            return str(int(round(re_part)))
        return f"{re_part:g}"
    real_text = "" if abs(re_part) < 1e-8 else f"{re_part:g}"
    sign = "+" if im_part >= 0 else "-"
    imag_mag = abs(im_part)
    imag_text = "i" if abs(imag_mag - 1) < 1e-8 else f"{imag_mag:g}i"
    if not real_text:
        return imag_text if sign == "+" else f"-{imag_text}"
    return f"{real_text}{sign}{imag_text}"


def _opd_term(coef: str, param: str, first: bool) -> str:
    if coef == "0":
        return ""
    negative = coef.startswith("-")
    mag = coef[1:] if negative else coef
    body = param if mag == "1" else f"{mag}{param}"
    if negative:
        return f"-{body}"
    return body if first else f"+{body}"


def _opd_from_rows(rows: list[list[str]]) -> str:
    if not rows:
        return ""
    params = list("abcdefghijklmnopqrstuvwxyz")
    width = max(len(row) for row in rows)
    comps: list[str] = []
    for col in range(width):
        expr = ""
        for row_index, row in enumerate(rows):
            coef = row[col] if col < len(row) else "0"
            expr += _opd_term(coef, params[row_index], not expr)
        comps.append(expr or "0")
    return "(" + ",".join(comps) + ")"


def _fmt_opd_display_value(value: float, tol: float = 1e-8) -> str:
    if abs(value) < tol:
        return "0"
    rounded = round(value)
    if abs(value - rounded) < tol:
        return str(int(rounded))
    return f"{value:.3f}"


def _direction(dec: SourceTables, direction_id: int) -> dict[str, Any]:
    im = dec.images
    did = int(direction_id)
    if did <= 0 or did > len(im["image_dir_label"]):
        return {"id": did, "label": "?", "opd": "", "vectors": []}
    i = did - 1
    dim = int(im["image_dir_dimension"][i])
    subduction = int(im["image_dir_subduction"][i])
    ptr = int(im["image_dir_pointer"][i])
    rows: list[list[str]] = []
    if ptr > 0 and dim > 0:
        start = ptr - 1
        next_ptr = None
        for raw_ptr in im["image_dir_pointer"][i + 1:]:
            if int(raw_ptr) > 0:
                next_ptr = int(raw_ptr) - 1
                break
        end = next_ptr if next_ptr is not None else start + max(dim * max(subduction, 1), dim)
        raw = im["image_dir"][start:end]
        rows = [
            [_fmt_complex(complex(dec.const[int(code)])) for code in raw[j:j + dim]]
            for j in range(0, len(raw), dim)
            if raw[j:j + dim]
        ]
    return {
        "id": did,
        "label": str(im["image_dir_label"][i]).strip(),
        "dimension": dim,
        "subduction": subduction,
        "opd": _opd_from_rows(rows),
        "vectors": rows,
    }


def _matrix_to_opd(dec: SourceTables, rows: list[list[int]]) -> str:
    return _opd_from_rows([[_fmt_complex(complex(dec.const[int(code)])) for code in row] for row in rows])


def _stored_orderparam(dec: SourceTables, row_id: int) -> dict[str, Any]:
    iso = dec.isotropy
    idx = int(row_id) - 1
    dim = int(iso["isotropy_orderparam_dim"][idx])
    free = int(iso["isotropy_orderparam_freeparam"][idx])
    ptr = int(iso["isotropy_orderparam_pointer"][idx]) - 1
    raw = iso["isotropy_orderparam"][ptr:ptr + dim * free]
    rows = [raw[r * free:(r + 1) * free] for r in range(dim)]
    return {"dim": dim, "free": free, "pointer": ptr + 1, "opd": _matrix_to_opd(dec, rows)}


def _source_orderparam(dec: SourceTables, row_id: int, full_dim: int) -> dict[str, Any]:
    iso = dec.isotropy
    idx = int(row_id) - 1
    row_count = int(iso["isotropy_orderparam_dim"][idx])
    free = int(iso["isotropy_orderparam_freeparam"][idx])
    ptr = int(iso["isotropy_orderparam_pointer"][idx]) - 1
    raw = iso["isotropy_orderparam"][ptr:ptr + row_count * int(full_dim)]
    numeric_rows = [
        [
            float(dec.const[int(code)])
            for code in raw[r * int(full_dim):(r + 1) * int(full_dim)]
        ]
        for r in range(row_count)
    ]
    rows = [
        [
            _fmt_opd_display_value(value)
            for value in numeric_rows[r]
        ]
        for r in range(row_count)
    ]
    display_rows = rows[:free]
    return {
        "row_count": row_count,
        "free": free,
        "full_dim": int(full_dim),
        "pointer": ptr + 1,
        "matrix": rows,
        "display_rows": display_rows,
        "numeric_rows": numeric_rows,
        "opd": _opd_from_rows(display_rows),
    }


def _magnetic_source_orderparam(dec: SourceTables, row: Any) -> dict[str, Any]:
    row_count = int(row.orderparam_dim)
    free = int(row.orderparam_freeparam)
    raw = list(row.orderparam)
    numeric_rows = [
        [
            float(dec.const[int(code)])
            for code in raw[param * row_count:(param + 1) * row_count]
        ]
        for param in range(free)
    ]
    display_rows = [
        [
            _fmt_opd_display_value(value)
            for value in numeric_rows[param]
        ]
        for param in range(free)
    ]
    return {
        "row_count": row_count,
        "free": free,
        "full_dim": row_count,
        "matrix": display_rows,
        "display_rows": display_rows,
        "numeric_rows": numeric_rows,
        "opd": _opd_from_rows(display_rows),
    }


def _fixed_magnetic_subgroup_selection(
    *,
    parent_sg: int,
    gid: int,
    row: Any,
    source: dict[str, Any],
) -> tuple[dict[str, Any] | None, tuple[tuple[int, int, int, int, int], ...]]:
    """Reconstruct fixed-row BNS correspondence when it matches the Source table."""

    free = int(source["free"])
    full_dim = int(source["full_dim"])
    stride = [0.0] * (free * 48)
    for row_index, values in enumerate(source["numeric_rows"][:free]):
        for col, value in enumerate(values[:full_dim]):
            stride[row_index * 48 + col] = float(value)
    try:
        decoder = _subgroup_core_data()
        basis, operations = decoder.orderparam_to_subgroup_magnetic(
            int(gid), (), tuple(stride), free
        )
        result = id_subgroup_magnetic_identify_with_generator_block(
            decoder,
            int(parent_sg),
            basis,
            operations,
            0,
        )
    except (IndexError, KeyError, TypeError, ValueError, ZeroDivisionError):
        return None, ()
    if (
        result is None
        or result.selection is None
        or int(result.magnetic_group) != int(row.subgroup_magnetic_group)
        or tuple(int(value) for value in result.basis) != tuple(int(value) for value in row.basis)
        or tuple(int(value) for value in result.origin) != tuple(int(value) for value in row.origin)
    ):
        return None, ()
    return (
        asdict(result.selection),
        tuple(tuple(int(value) for value in record) for record in operations),
    )


@lru_cache(maxsize=1)
def _subgroup_core_data() -> SourceData:
    return SourceData(SOURCE, tables=source_tables().iso)


def _origin_text(raw: tuple[int, int, int, int]) -> str:
    den = int(raw[3] or 1)
    values = []
    for value in raw[:3]:
        fraction = Fraction(int(value), den)
        if fraction.denominator == 1:
            values.append(str(fraction.numerator))
        else:
            values.append(f"{fraction.numerator}/{fraction.denominator}")
    return "(" + ",".join(values) + ")"


def _basis_value(value: int | Fraction, denominator: int = 1) -> int | float:
    if isinstance(value, Fraction):
        fraction = value
    else:
        fraction = Fraction(int(value), int(denominator or 1))
    return int(fraction) if fraction.denominator == 1 else float(fraction)


def _basis_values(rows: list[list[int | Fraction]], denominator: int = 1) -> list[list[int | float]]:
    return [[_basis_value(value, denominator) for value in row] for row in rows]


def _basis_text(rows: list[list[int | Fraction]], denominator: int = 1) -> str:
    den = int(denominator or 1)
    return ",".join(
        "(" + ",".join(_fmt_frac(value if isinstance(value, Fraction) else Fraction(int(value), den)) for value in row) + ")"
        for row in rows
    )


def _identity_basis_rows(rows: list[list[int | Fraction]]) -> bool:
    if len(rows) != 3 or any(len(row) != 3 for row in rows):
        return False
    for row in range(3):
        for col in range(3):
            expected = Fraction(1) if row == col else Fraction(0)
            value = rows[row][col]
            if not isinstance(value, Fraction):
                value = Fraction(int(value), 1)
            if value != expected:
                return False
    return True


def _display_origin_text(
    parent_sg: int,
    subgroup: int,
    basis_rows: list[list[int | Fraction]],
    origin: tuple[int, int, int, int],
) -> str:
    if int(parent_sg) == int(subgroup) and tuple(int(value) for value in origin[:3]) == (0, 0, 0):
        return "(0,0,0)"
    if int(parent_sg) == int(subgroup) and _identity_basis_rows(basis_rows):
        return "(0,0,0)"
    return _origin_text(origin)


def _dynamic_kparam(
    dec: SourceTables,
    *,
    gid: int,
    k_params: dict[str, str] | None,
) -> tuple[int, int, int, int] | None:
    if not k_params:
        return None
    for key in ("_source_kparam", "source_kparam", "kparam"):
        raw_kparam = k_params.get(key)
        if raw_kparam is None:
            continue
        values = [int(value) for value in re.findall(r"-?\d+", str(raw_kparam))]
        if len(values) >= 4 and values[3] != 0:
            return (values[0], values[1], values[2], values[3])
    source_kparam = _source_kparam_from_visible(dec, gid=gid, k_params=k_params)
    if source_kparam is not None:
        return source_kparam
    sg = int(dec.little["little_irr_space_group"][int(gid) - 1])
    kslot = int(dec.little["little_irr_k"][int(gid) - 1])
    lattice = int(dec.space["ispace_lattice"][sg - 1])
    lattice_slot = (lattice - 1) * 27 + kslot - 1
    dim = int(dec.little["little_k_dim"][lattice_slot])
    values: list[Fraction] = [Fraction(0, 1), Fraction(0, 1), Fraction(0, 1)]
    # Source get_isotropy_ expects values by param-slot order.  The Web form
    # names fields by the public k-vector expression (a,b,g), which is not the
    # same thing as the nonzero coordinate axis in little_k's internal basis.
    # Use the displayed ISODISTORT k-vector as the authoritative bridge.
    display_vector = str(_k_slot_info(dec, sg, kslot).get("isodistort_kvector") or "")
    display_keys: list[str] = []
    for key in re.findall(r"(?<![A-Za-z])([abg])(?![A-Za-z])", display_vector):
        if key not in display_keys:
            display_keys.append(key)
    for param_index in range(min(dim, 3)):
        key = display_keys[param_index] if param_index < len(display_keys) else ISODISTORT_K_PARAMS[param_index]
        raw_value = k_params.get(key)
        if raw_value is None:
            raw_value = k_params.get(ISODISTORT_K_PARAMS[param_index])
        raw = raw_value
        if raw is None or str(raw).strip() == "":
            continue
        values[param_index] = Fraction(str(raw).strip())
    denominator = 1
    for value in values:
        denominator = math.lcm(denominator, value.denominator)
    return (
        int(values[0] * denominator),
        int(values[1] * denominator),
        int(values[2] * denominator),
        int(denominator),
    )


def _source_kparam_from_visible(
    dec: SourceTables,
    *,
    gid: int,
    k_params: dict[str, str] | None,
) -> tuple[int, int, int, int] | None:
    """Convert ISODISTORT-visible k params to Source/iso internal kparam.

    ISODISTORT displays parametric k coordinates in its conventional
    presentation.  Source/iso dynamic isotropy generation expects the
    Miller-Love/internal coordinates carried by ``little_k_star_conv2ml``.
    Web stores that result as hidden ``kparam1``; local code must derive it
    from Source data instead of reading the hidden field.
    """

    if not k_params:
        return None
    sg = int(dec.little["little_irr_space_group"][int(gid) - 1])
    kslot = int(dec.little["little_irr_k"][int(gid) - 1])
    lattice = int(dec.space["ispace_lattice"][sg - 1])
    lattice_slot = (lattice - 1) * 27 + kslot - 1
    dim = int(dec.little["little_k_dim"][lattice_slot])
    if dim <= 0:
        return None
    display_vector = str(_k_slot_info(dec, sg, kslot).get("isodistort_kvector") or "")
    match = re.fullmatch(r"\((.*)\)", display_vector.strip())
    if not match:
        return None
    visible: list[Fraction] = []
    for component in match.group(1).split(","):
        value = evaluate_k_component(component, k_params)
        if value is None:
            return None
        visible.append(value)
    if len(visible) != 3:
        return None
    sg_slot = (sg - 1) * 27 + kslot - 1
    pointer = int(dec.little.get("little_k_star_conv2ml_pointer", [0])[sg_slot])
    if pointer <= 0:
        return None
    raw = [int(x) for x in dec.little["little_k_star"][16 * (pointer - 1):16 * pointer]]
    if len(raw) != 16:
        return None
    denominator = int(raw[15] or 1)
    converted = [
        Fraction(raw[12 + col], denominator)
        + sum(visible[row] * Fraction(raw[4 * row + col], denominator) for row in range(3))
        for col in range(3)
    ]
    den = 1
    for value in converted:
        den = math.lcm(den, value.denominator)
    return (
        int(converted[0] * den),
        int(converted[1] * den),
        int(converted[2] * den),
        int(den),
    )


def _opd_subgroup_setting_id(
    parent_sg: int,
    ordinary_subgroup: int,
    parent_setting_id: int | None,
) -> int | None:
    """Use one inter-setting on both sides of a same-SG OPD embedding.

    DISPLAY ISOTROPY keeps the selected parent setting from appearing as an
    extra subgroup-origin shift when the subgroup has the same ordinary space
    group number.  A proper subgroup remains in its default setting.
    """

    if parent_setting_id is None or int(parent_sg) != int(ordinary_subgroup):
        return None
    return int(parent_setting_id)


def _dynamic_source_matrix(row: DynamicIsotropyRow) -> list[list[str]]:
    return [list(matrix_row) for matrix_row in row.matrix]


def _display_opd(opd: str) -> str:
    text = str(opd).strip()
    match = re.fullmatch(r"\((.*)\)", text)
    if not match:
        return text
    parts = [part.strip() for part in match.group(1).split(",")]
    if len(parts) >= 4 and len(parts) % 2 == 0:
        chunks = [",".join(parts[index:index + 2]) for index in range(0, len(parts), 2)]
        text = "(" + ";".join(chunks) + ")"

    def round_coeff(match: re.Match[str]) -> str:
        value = float(match.group(0))
        if abs(value - round(value)) < 1e-8:
            rounded = str(int(round(value)))
        else:
            rounded = f"{value:.3f}"
        return "0" if rounded == "-0" else rounded

    return re.sub(r"-?\d+\.\d+(?=[A-Za-z])", round_coeff, text)


def _dynamic_isotropy_row(
    dec: SourceTables,
    *,
    parent_sg: int,
    gid: int,
    index: int,
    row: DynamicIsotropyRow,
    parent_setting_id: int | None,
    parent_cell: tuple[float, float, float, float, float, float] | None,
    source_kparam: tuple[int, int, int, int] | None = None,
    k_params: dict[str, str] | None = None,
    display_k_params: dict[str, str] | None = None,
) -> dict[str, Any]:
    subgroup = int(row.subgroup_number)
    basis_values = tuple(int(value) for value in row.basis_values)
    origin_values = tuple(int(value) for value in row.origin_values)
    displayed = dec.subgroup_change_setting_cinter(
        int(parent_sg),
        subgroup,
        basis_values,
        origin_values,
        parent_setting_id=parent_setting_id,
        subgroup_setting_id=_opd_subgroup_setting_id(
            int(parent_sg), subgroup, parent_setting_id
        ),
    )
    basis_rows = _basis_fraction_rows(displayed["basis"], int(displayed["basis_denominator"]))
    basis_rows = present_opd_basis_rows(
        subgroup,
        basis_rows,
        parent_cell,
        data=_subgroup_core_data(),
        parametric=True,
    )
    # Web OPD `s` follows the native dynamic .iso row basis determinant.  The
    # displayed basis may be transformed to the parent setting for text output,
    # but that presentation transform must not redefine the subgroup size.
    det = abs(integer_determinant3(basis_values))
    source_matrix = _dynamic_source_matrix(row)
    display_opd = _display_opd(row.opd)
    full_dim = int(dec.little["little_irr_full_dim"][int(gid) - 1])
    k_active = _active_k_vectors_for_gid(
        dec,
        parent_sg,
        int(gid),
        full_dim,
        source_matrix,
        display_k_params if display_k_params is not None else k_params,
    )
    return {
        "row_id": int(index),
        "dynamic": True,
        "dynamic_gid": int(gid),
        "canonical": int(index) == 1,
        "opd_label": str(row.direction),
        "direction_id": 0,
        "direction_label": str(row.direction),
        "direction_opd": row.opd,
        "stored_opd": row.opd,
        "display_opd": display_opd,
        "source_opd": row.opd,
        "source_matrix": source_matrix,
        "source_display_rows": source_matrix,
        "source_kparam": list(source_kparam) if source_kparam is not None else None,
        "op_rows": int(row.free),
        "free": int(row.free),
        "arms": det,
        "s": det,
        "i": _subgroup_index(dec, int(parent_sg), subgroup, det),
        "k_active_vectors": k_active,
        "k_active": ",".join(k_active),
        "subgroup": {"number": subgroup, "symbol": dec.default_setting_space_symbol(subgroup)},
        "basis": _basis_values(basis_rows),
        "basis_text": _basis_text(basis_rows),
        "basis_denominator": int(displayed["basis_denominator"]),
        "source_basis_values": list(basis_values),
        "source_origin_values": list(origin_values),
        "det": det,
        "origin": _display_origin_text(parent_sg, subgroup, basis_rows, displayed["origin"]),
        "origin_values": list(row.origin_values),
        "ferroic": None,
    }


def _dynamic_opd_rows(
    dec: SourceTables,
    *,
    parent_sg: int,
    gid: int,
    k_params: dict[str, str] | None = None,
    display_k_params: dict[str, str] | None = None,
    selected_orderparam: int | str | None = None,
    parent_setting_id: int | None = None,
    parent_cell: tuple[float, float, float, float, float, float] | None = None,
) -> list[dict[str, Any]]:
    kparam = _dynamic_kparam(dec, gid=gid, k_params=k_params)
    if kparam is None:
        return []
    point_occurrence = _dynamic_point_occurrence(selected_orderparam)
    selected_row = (
        generate_dynamic_point_isotropy_row(
            _subgroup_core_data(),
            gid=int(gid),
            kparam=kparam,
            occurrence=point_occurrence,
        )
        if point_occurrence is not None
        else None
    )
    if selected_row is not None:
        # P rows precede every higher-free row in final dynamic-file order.
        # Preserve the complete-file row id even though this opt-in surface
        # returns only the explicitly selected row.
        indexed_rows = ((point_occurrence, selected_row),)
    else:
        rows = sort_dynamic_rows_for_file(
            generate_dynamic_isotropy_rows(_subgroup_core_data(), gid=int(gid), kparam=kparam)
        )
        indexed_rows = tuple(enumerate(rows, start=1))
    out: list[dict[str, Any]] = []
    for index, row in indexed_rows:
        iso_row = _dynamic_isotropy_row(
            dec,
            parent_sg=parent_sg,
            gid=int(gid),
            index=index,
            row=row,
            parent_setting_id=parent_setting_id,
            parent_cell=parent_cell,
            source_kparam=kparam,
            k_params=k_params,
            display_k_params=display_k_params,
        )
        direction = {
            "id": 0,
            "label": str(row.direction),
            "dimension": int(row.free),
            "subduction": 0,
            "opd": _display_opd(row.opd),
            "vectors": source_matrix if (source_matrix := _dynamic_source_matrix(row)) else [],
            "dynamic": True,
        }
        out.append({"direction": direction, "isotropy": iso_row})
    return out


def _dynamic_point_occurrence(orderparam: int | str | None) -> int | None:
    """Return the exact Source occurrence for an explicit ``P<n>`` selector."""

    if not isinstance(orderparam, str):
        return None
    match = re.fullmatch(r"P([1-9][0-9]*)", orderparam.strip())
    return int(match.group(1)) if match is not None else None


def _basis_fraction_rows(rows: list[list[int]], denominator: int) -> list[list[Fraction]]:
    den = int(denominator or 1)
    return [[Fraction(int(value), den) for value in row] for row in rows]


def _first_little_gid_for_old_id(dec: SourceTables, parent_sg: int, old_id: int) -> int | None:
    for gid, row_old_id in enumerate(dec.little["little_irr_old"], start=1):
        if int(row_old_id) != int(old_id):
            continue
        if int(dec.little["little_irr_space_group"][gid - 1]) == int(parent_sg):
            return gid
    return None


def _compact_vector_text(value: str) -> str:
    return re.sub(r"\s+", "", str(value))


def _is_nonzero_opd_value(value: Any) -> bool:
    try:
        return abs(float(value)) > 1e-9
    except (TypeError, ValueError):
        return str(value).strip() not in {"", "0", "0.0", "-0", "-0.0"}


def _display_records_for_star_arm(
    dec: SourceTables,
    *,
    parent_sg: int,
    lattice: int,
    symbol: str,
    label: str,
    source_records: tuple[tuple[Fraction, Fraction, Fraction], ...],
) -> tuple[tuple[Fraction, Fraction, Fraction], ...]:
    records = _isodistort_k_records(source_records, lattice, symbol, label)
    if int(lattice) in {9, 10, 13, 14}:
        records = _apply_fraction_k_matrix(source_records, _reciprocal_pml_to_cinter_matrix(dec, int(parent_sg)))
    return records


def _visible_k_values_for_gid(
    dec: SourceTables,
    *,
    parent_sg: int,
    gid: int,
    k_params: dict[str, str],
) -> tuple[Fraction, Fraction, Fraction] | None:
    kslot = int(dec.little["little_irr_k"][int(gid) - 1])
    display_vector = str(_k_slot_info(dec, int(parent_sg), kslot).get("isodistort_kvector") or "")
    match = re.fullmatch(r"\((.*)\)", display_vector.strip())
    if not match:
        return None
    values: list[Fraction] = []
    for component in match.group(1).split(","):
        value = evaluate_k_component(component, k_params)
        if value is None:
            return None
        values.append(value)
    if len(values) != 3:
        return None
    return values[0], values[1], values[2]


def _solve_star_parameter_values(
    records: tuple[tuple[Fraction, Fraction, Fraction], ...],
    visible: tuple[Fraction, Fraction, Fraction],
) -> list[Fraction] | None:
    if not records:
        return None
    dim = len(records) - 1
    if dim <= 0:
        return []
    equations = [
        [records[param_index + 1][axis] for param_index in range(dim)]
        + [visible[axis] - records[0][axis]]
        for axis in range(3)
    ]
    row = 0
    pivots: list[int] = []
    for col in range(dim):
        pivot = next((candidate for candidate in range(row, 3) if equations[candidate][col]), None)
        if pivot is None:
            continue
        equations[row], equations[pivot] = equations[pivot], equations[row]
        scale = equations[row][col]
        equations[row] = [value / scale for value in equations[row]]
        for other in range(3):
            if other == row:
                continue
            factor = equations[other][col]
            if factor:
                equations[other] = [
                    equations[other][index] - factor * equations[row][index]
                    for index in range(dim + 1)
                ]
        pivots.append(col)
        row += 1
        if row == 3:
            break
    for equation in equations:
        if not any(equation[:dim]) and equation[dim]:
            return None
    if len(pivots) != dim:
        return None
    values = [Fraction(0) for _ in range(dim)]
    for equation in equations:
        pivot = next((index for index, value in enumerate(equation[:dim]) if value), None)
        if pivot is not None:
            values[pivot] = equation[dim]
    return values


def _numeric_isodistort_star_vectors(
    dec: SourceTables,
    parent_sg: int,
    gid: int,
    k_params: dict[str, str] | None,
) -> list[str]:
    if not k_params:
        return []
    kslot = int(dec.little["little_irr_k"][int(gid) - 1])
    lattice = int(dec.space["ispace_lattice"][int(parent_sg) - 1])
    lattice_slot = (lattice - 1) * 27 + kslot - 1
    dim = int(dec.little["little_k_dim"][lattice_slot])
    sg_slot = (int(parent_sg) - 1) * 27 + kslot - 1
    if dim <= 0:
        return []
    count = int(dec.little["little_k_star_count"][sg_slot])
    ml_ptr = int(dec.little["little_k_star_ml_pointer"][sg_slot])
    if count <= 0 or ml_ptr <= 0:
        return []
    label = _strip(dec.little["little_k_label"][lattice_slot])
    symbol = dec.space_symbol(parent_sg).strip()
    raw = dec.little["little_k_star"]
    first_item = [int(x) for x in raw[16 * (ml_ptr - 1):16 * ml_ptr]]
    if len(first_item) != 16:
        return []
    first_source_records = _parse_k_record(first_item)[:dim + 1]
    first_records = _display_records_for_star_arm(
        dec,
        parent_sg=parent_sg,
        lattice=lattice,
        symbol=symbol,
        label=label,
        source_records=first_source_records,
    )
    visible = _visible_k_values_for_gid(dec, parent_sg=parent_sg, gid=gid, k_params=k_params)
    if int(lattice) == 3 and label == "GP" and visible is not None and count == 2:
        second_item = [int(x) for x in raw[16 * ml_ptr:16 * (ml_ptr + 1)]]
        second_source = _parse_k_record(second_item)[:dim + 1] if len(second_item) == 16 else ()
        first_sign = second_source[1][0] if len(second_source) > 1 else Fraction(-1)
        second = (visible[0], -visible[1], visible[2]) if first_sign > 0 else (-visible[0], visible[1], -visible[2])
        return [
            "(" + ",".join(_fmt_frac(value) for value in visible) + ")",
            "(" + ",".join(_fmt_frac(value) for value in second) + ")",
        ]
    parameter_values = _solve_star_parameter_values(first_records, visible) if visible is not None else None
    if parameter_values is None:
        return []
    out: list[str] = []
    for offset in range(count):
        item = [int(x) for x in raw[16 * (ml_ptr - 1 + offset):16 * (ml_ptr + offset)]]
        if len(item) != 16:
            break
        source_records = _parse_k_record(item)[:dim + 1]
        values: list[str] = []
        records = _display_records_for_star_arm(
            dec,
            parent_sg=parent_sg,
            lattice=lattice,
            symbol=symbol,
            label=label,
            source_records=source_records,
        )
        if len(records) - 1 != len(parameter_values):
            return []
        for axis in range(3):
            value = records[0][axis]
            for param_index, vector in enumerate(records[1:]):
                value += vector[axis] * parameter_values[param_index]
            values.append(_fmt_frac(value))
        out.append("(" + ",".join(values) + ")")
    return out


def _active_k_vectors(
    dec: SourceTables,
    parent_sg: int,
    old_id: int,
    source: dict[str, Any],
    k_params: dict[str, str] | None = None,
) -> list[str]:
    gid = _first_little_gid_for_old_id(dec, parent_sg, old_id)
    if gid is None:
        return []
    full_dim = int(source.get("full_dim") or 0)
    rows = source.get("display_rows") or source.get("matrix") or []
    return _active_k_vectors_for_gid(dec, parent_sg, gid, full_dim, rows, k_params)


def _active_k_vectors_for_gid(
    dec: SourceTables,
    parent_sg: int,
    gid: int,
    full_dim: int,
    rows: Any,
    k_params: dict[str, str] | None = None,
) -> list[str]:
    kslot = int(dec.little["little_irr_k"][gid - 1])
    lattice = int(dec.space["ispace_lattice"][int(parent_sg) - 1])
    lattice_slot = (lattice - 1) * 27 + kslot - 1
    dim = int(dec.little["little_k_dim"][lattice_slot])
    star = _numeric_isodistort_star_vectors(dec, parent_sg, gid, k_params)
    if not star:
        star = [
            _compact_vector_text(substitute_k_vector(vector, k_params))
            for vector in _isodistort_star_vectors(dec, parent_sg, kslot, dim)
        ]
    if not star:
        return []
    if full_dim <= 0 or full_dim % len(star) != 0:
        return []
    arm_width = full_dim // len(star)
    active: list[str] = []
    for arm_index, vector in enumerate(star):
        start = arm_index * arm_width
        stop = start + arm_width
        if any(any(_is_nonzero_opd_value(value) for value in row[start:stop]) for row in rows):
            active.append(vector)
    return active


def _subgroup_index(dec: SourceTables, parent_sg: int, subgroup: int, size: int) -> int | None:
    parent_order = dec.space_group_point_group_order(parent_sg)
    subgroup_order = dec.space_group_point_group_order(subgroup)
    value = int(size) * int(parent_order)
    if subgroup_order <= 0 or value % subgroup_order:
        return None
    return value // subgroup_order


def _image_directions(dec: SourceTables, image_id: int) -> list[dict[str, Any]]:
    im = dec.images
    index = int(image_id) - 1
    if index < 0 or index >= len(im["image_label"]):
        return []
    ptr = int(im["image_subgroup_pointer"][index]) - 1
    count = int(im["image_subgroup_count"][index])
    out: list[dict[str, Any]] = []
    for offset in range(count):
        entry = ptr + offset
        direction = _direction(dec, int(im["image_subgroup"][entry]))
        out.append({
            **direction,
            "maximal": bool(im["image_subgroup_maximal"][entry]),
            "minimal": bool(im["image_subgroup_min"][entry]),
            "rg": bool(im["image_subgroup_rg"][entry]),
        })
    return out


def _isotropy_rows(
    dec: SourceTables,
    old_id: int,
    *,
    k_params: dict[str, str] | None = None,
    display_k_params: dict[str, str] | None = None,
    parent_setting_id: int | None = None,
    parent_cell: tuple[float, float, float, float, float, float] | None = None,
) -> list[dict[str, Any]]:
    if old_id <= 0:
        return []
    parent_sg = int(dec.irreps["irrep_space_group"][int(old_id) - 1])
    full_dim = int(dec.image_record(old_id)["dimension"])
    iso = dec.isotropy
    out: list[dict[str, Any]] = []
    for row_id in dec.isotropy_row_ids_for_old_irrep(old_id):
        idx = row_id - 1
        basis = tuple(int(x) for x in iso["isotropy_basis"][idx * 9:(idx + 1) * 9])
        origin = tuple(int(x) for x in iso["isotropy_origin"][idx * 4:(idx + 1) * 4])
        stored = _stored_orderparam(dec, row_id)
        source = _source_orderparam(dec, row_id, full_dim)
        direction_id = int(iso["isotropy_direction"][idx])
        direction = _direction(dec, direction_id)
        subgroup = int(iso["isotropy_subgroup"][idx])
        displayed = dec.subgroup_change_setting_cinter(  # type: ignore[arg-type]
            parent_sg,
            subgroup,
            basis,
            origin,
            parent_setting_id=parent_setting_id,
            subgroup_setting_id=_opd_subgroup_setting_id(
                parent_sg, subgroup, parent_setting_id
            ),
        )
        size = int(displayed["size"])
        # ISODISTORT prints s from the native data_isotropy supercell basis,
        # before the display-setting reduction used for basis/origin text.
        supercell_size = abs(integer_determinant3(basis))
        subgroup_index = _subgroup_index(dec, parent_sg, subgroup, supercell_size)
        k_active = _active_k_vectors(
            dec,
            parent_sg,
            old_id,
            source,
            display_k_params if display_k_params is not None else k_params,
        )
        basis_rows = _basis_fraction_rows(displayed["basis"], int(displayed["basis_denominator"]))
        basis_rows = present_opd_basis_rows(
            subgroup,
            basis_rows,
            parent_cell,
            data=_subgroup_core_data(),
            parametric=bool(k_params),
        )
        out.append({
            "row_id": int(row_id),
            "canonical": False,
            "opd_label": str(iso["isotropy_orderparam_label"][idx]).strip(),
            "direction_id": direction_id,
            "direction_label": direction["label"],
            "direction_opd": direction["opd"],
            "stored_opd": stored["opd"],
            "display_opd": _display_opd(source["opd"]),
            "source_opd": source["opd"],
            "source_matrix": source["matrix"],
            "source_display_rows": source["display_rows"],
            "source_numeric_rows": source["numeric_rows"][: int(source["free"])],
            "op_rows": source["row_count"],
            "free": source["free"],
            "arms": size,
            "s": supercell_size,
            "i": subgroup_index,
            "k_active_vectors": k_active,
            "k_active": ",".join(k_active),
            "subgroup": {"number": subgroup, "symbol": dec.default_setting_space_symbol(subgroup)},
            "basis": _basis_values(basis_rows),
            "basis_text": _basis_text(basis_rows),
            "basis_denominator": int(displayed["basis_denominator"]),
            "source_basis_values": list(basis),
            "source_origin_values": list(origin),
            "det": size,
            "origin": _display_origin_text(parent_sg, subgroup, basis_rows, displayed["origin"]),
            "ferroic": {
                "raw": int(iso["isotropy_ferroic"][idx]),
                "ferroelectric": int(iso["isotropy_ferroelectric"][idx]),
                "ferroelastic": int(iso["isotropy_ferroelastic"][idx]),
            },
        })
    canonical_id = dec.canonical_isotropy_row_id(out)
    for row in out:
        row["canonical"] = int(row["row_id"]) == canonical_id
    return out


def _magnetic_subgroup_index_for_groups(parent_magnetic_group: int, subgroup_magnetic_group: int, size: int) -> int | None:
    table = magnetic_data().table
    parent_pg = int(table["mag_point_group"][int(parent_magnetic_group) - 1])
    subgroup_pg = int(table["mag_point_group"][int(subgroup_magnetic_group) - 1])
    parent_order = int(table["mag_point_group_order"][parent_pg - 1])
    subgroup_order = int(table["mag_point_group_order"][subgroup_pg - 1])
    value = int(size) * parent_order
    if subgroup_order <= 0 or value % subgroup_order:
        return None
    return value // subgroup_order


def _magnetic_subgroup_index(row: Any, size: int) -> int | None:
    return _magnetic_subgroup_index_for_groups(
        int(row.parent_magnetic_group),
        int(row.subgroup_magnetic_group),
        int(size),
    )


def _dynamic_magnetic_isotropy_row(
    dec: SourceTables,
    *,
    parent_sg: int,
    gid: int,
    index: int,
    row: DynamicIsotropyRow,
    parent_setting_id: int | None,
    parent_cell: tuple[float, float, float, float, float, float] | None,
    source_kparam: tuple[int, int, int, int] | None = None,
    k_params: dict[str, str] | None = None,
    display_k_params: dict[str, str] | None = None,
) -> dict[str, Any]:
    table = magnetic_data().table
    parent_magnetic_group = _subgroup_core_data().magnetic_orderparam_group_for_irrep(int(gid))
    subgroup_magnetic_group = int(row.subgroup_number)
    ordinary_subgroup = int(table["mag_space_group"][subgroup_magnetic_group - 1])
    basis_values = tuple(int(value) for value in row.basis_values)
    origin_values = tuple(int(value) for value in row.origin_values)
    displayed = dec.subgroup_change_setting_cinter(
        int(parent_sg),
        ordinary_subgroup,
        basis_values,
        origin_values,
        parent_setting_id=parent_setting_id,
        subgroup_setting_id=_opd_subgroup_setting_id(
            int(parent_sg), ordinary_subgroup, parent_setting_id
        ),
    )
    basis_rows = _basis_fraction_rows(displayed["basis"], int(displayed["basis_denominator"]))
    basis_rows = present_opd_basis_rows(
        ordinary_subgroup,
        basis_rows,
        parent_cell,
        data=_subgroup_core_data(),
        parametric=True,
    )
    det = abs(integer_determinant3(basis_values))
    source_matrix = _dynamic_source_matrix(row)
    display_opd = _display_opd(row.opd)
    full_dim = int(dec.little["little_irr_full_dim"][int(gid) - 1])
    k_active = _active_k_vectors_for_gid(
        dec,
        parent_sg,
        int(gid),
        full_dim,
        source_matrix,
        display_k_params if display_k_params is not None else k_params,
    )
    selection = row.magnetic_subgroup_selection
    return {
        "row_id": int(index),
        "dynamic": True,
        "dynamic_gid": int(gid),
        "canonical": int(index) == 1,
        "magnetic": True,
        "opd_label": str(row.direction),
        "direction_id": 0,
        "direction_label": str(row.direction),
        "direction_opd": row.opd,
        "stored_opd": row.opd,
        "display_opd": display_opd,
        "source_opd": row.opd,
        "source_matrix": source_matrix,
        "source_display_rows": source_matrix,
        "source_kparam": list(source_kparam) if source_kparam is not None else None,
        "op_rows": int(row.free),
        "free": int(row.free),
        "arms": det,
        "s": det,
        "i": _magnetic_subgroup_index_for_groups(parent_magnetic_group, subgroup_magnetic_group, det),
        "k_active_vectors": k_active,
        "k_active": ",".join(k_active),
        "subgroup": {
            "number": subgroup_magnetic_group,
            "symbol": str(table["mag_bns_label"][subgroup_magnetic_group - 1]).strip(),
            "ordinary_number": ordinary_subgroup,
            "display_label": str(table["mag_nlabel"][subgroup_magnetic_group - 1]).strip(),
        },
        "basis": _basis_values(basis_rows),
        "basis_text": _basis_text(basis_rows),
        "basis_denominator": int(displayed["basis_denominator"]),
        "source_basis_values": list(basis_values),
        "source_origin_values": list(origin_values),
        **(
            {"magnetic_subgroup_selection": asdict(selection)}
            if selection is not None
            else {}
        ),
        **(
            {
                "source_operation_records": [
                    list(record) for record in row.magnetic_operation_records
                ]
            }
            if row.magnetic_operation_records
            else {}
        ),
        "det": det,
        "origin": _display_origin_text(parent_sg, ordinary_subgroup, basis_rows, displayed["origin"]),
        "origin_values": list(origin_values),
        "ferroic": None,
        "magnetic_parent_group": int(parent_magnetic_group),
        "magnetic_subgroup_old": None,
    }


def magnetic_opd_rows(
    old_id: int,
    *,
    parent_sg: int | None = None,
    gid: int | None = None,
    k_params: dict[str, str] | None = None,
    display_k_params: dict[str, str] | None = None,
    selected_orderparam: int | str | None = None,
    parent_setting_id: int | None = None,
    parent_cell: tuple[float, float, float, float, float, float] | None = None,
) -> list[dict[str, Any]]:
    dec = source_tables()
    if int(old_id) <= 0:
        if parent_sg is None or gid is None:
            return []
        return _dynamic_magnetic_opd_rows(
            dec,
            parent_sg=int(parent_sg),
            gid=int(gid),
            k_params=k_params,
            display_k_params=display_k_params,
            selected_orderparam=selected_orderparam,
            parent_setting_id=parent_setting_id,
            parent_cell=parent_cell,
        )
    table = magnetic_data().table
    out: list[dict[str, Any]] = []
    for row in magnetic_data().magnetic_isotropy_rows_for_irrep(int(old_id)):
        basis = tuple(int(value) for value in row.basis)
        origin = tuple(int(value) for value in row.origin)
        native_size = abs(integer_determinant3(basis))
        ordinary_parent = int(table["mag_space_group"][int(row.parent_magnetic_group) - 1])
        ordinary_subgroup = int(table["mag_space_group"][int(row.subgroup_magnetic_group) - 1])
        displayed = dec.subgroup_change_setting_cinter(
            ordinary_parent,
            ordinary_subgroup,
            basis,
            origin,
            parent_setting_id=parent_setting_id,
            subgroup_setting_id=_opd_subgroup_setting_id(
                ordinary_parent, ordinary_subgroup, parent_setting_id
            ),
        )
        basis_rows = _basis_fraction_rows(displayed["basis"], int(displayed["basis_denominator"]))
        basis_rows = present_opd_basis_rows(
            ordinary_subgroup,
            basis_rows,
            parent_cell,
            data=_subgroup_core_data(),
            parametric=False,
        )
        source = _magnetic_source_orderparam(dec, row)
        selection_payload: dict[str, Any] | None = None
        source_operation_records: tuple[tuple[int, int, int, int, int], ...] = ()
        if gid is not None:
            selection_payload, source_operation_records = _fixed_magnetic_subgroup_selection(
                parent_sg=int(parent_sg or ordinary_parent),
                gid=int(gid),
                row=row,
                source=source,
            )
        k_active = (
            _active_k_vectors(
                dec,
                int(parent_sg or ordinary_parent),
                int(old_id),
                source,
                display_k_params,
            )
            if gid is not None
            else []
        )
        iso_row = {
            "row_id": int(row.row_id),
            "canonical": False,
            "magnetic": True,
            "opd_label": row.orderparam_label,
            "direction_id": int(row.row_id),
            "direction_label": row.orderparam_label,
            "direction_opd": source["opd"],
            "stored_opd": source["opd"],
            "display_opd": source["opd"],
            "source_opd": source["opd"],
            "source_matrix": source["matrix"],
            "source_display_rows": source["display_rows"],
            "source_numeric_rows": source["numeric_rows"],
            "op_rows": source["row_count"],
            "free": source["free"],
            "arms": int(displayed["size"]),
            "s": native_size,
            "i": _magnetic_subgroup_index(row, native_size),
            "k_active_vectors": k_active,
            "k_active": ",".join(k_active),
            "subgroup": {
                "number": int(row.subgroup_magnetic_group),
                "symbol": row.subgroup_bns_label,
                "ordinary_number": ordinary_subgroup,
                "display_label": row.subgroup_number_label,
            },
            "basis": _basis_values(basis_rows),
            "basis_text": _basis_text(basis_rows),
            "basis_denominator": int(displayed["basis_denominator"]),
            "source_basis_values": list(basis),
            "source_origin_values": list(origin),
            **(
                {"magnetic_subgroup_selection": selection_payload}
                if selection_payload is not None
                else {}
            ),
            **(
                {
                    "source_operation_records": [
                        list(record) for record in source_operation_records
                    ]
                }
                if source_operation_records
                else {}
            ),
            "det": int(displayed["size"]),
            "origin": _display_origin_text(ordinary_parent, ordinary_subgroup, basis_rows, displayed["origin"]),
            "magnetic_parent_group": int(row.parent_magnetic_group),
            "magnetic_subgroup_old": int(row.subgroup_old),
        }
        out.append({
            "direction": {
                "id": int(row.row_id),
                "label": row.orderparam_label,
                "opd": source["opd"],
                "dimension": int(row.orderparam_dim),
                "subduction": int(row.orderparam_freeparam),
                "vectors": source["matrix"],
            },
            "isotropy": iso_row,
        })
    if out:
        out[0]["isotropy"]["canonical"] = True
    return out


def _dynamic_magnetic_opd_rows(
    dec: SourceTables,
    *,
    parent_sg: int,
    gid: int,
    k_params: dict[str, str] | None = None,
    display_k_params: dict[str, str] | None = None,
    selected_orderparam: int | str | None = None,
    parent_setting_id: int | None = None,
    parent_cell: tuple[float, float, float, float, float, float] | None = None,
) -> list[dict[str, Any]]:
    kparam = _dynamic_kparam(dec, gid=gid, k_params=k_params)
    if kparam is None:
        return []
    point_occurrence = _dynamic_point_occurrence(selected_orderparam)
    selected_row = (
        generate_dynamic_point_isotropy_row_magnetic(
            _subgroup_core_data(),
            gid=int(gid),
            kparam=kparam,
            occurrence=point_occurrence,
        )
        if point_occurrence is not None
        else None
    )
    if selected_row is not None:
        indexed_rows = ((point_occurrence, selected_row),)
    else:
        rows = sort_dynamic_rows_for_file(
            generate_dynamic_isotropy_rows_magnetic(
                _subgroup_core_data(),
                gid=int(gid),
                kparam=kparam,
            )
        )
        indexed_rows = tuple(enumerate(rows, start=1))
    out: list[dict[str, Any]] = []
    for index, row in indexed_rows:
        iso_row = _dynamic_magnetic_isotropy_row(
            dec,
            parent_sg=parent_sg,
            gid=int(gid),
            index=index,
            row=row,
            parent_setting_id=parent_setting_id,
            parent_cell=parent_cell,
            source_kparam=kparam,
            k_params=k_params,
            display_k_params=display_k_params,
        )
        direction = {
            "id": 0,
            "label": str(row.direction),
            "dimension": int(row.free),
            "subduction": 0,
            "opd": _display_opd(row.opd),
            "vectors": source_matrix if (source_matrix := _dynamic_source_matrix(row)) else [],
            "dynamic": True,
            "magnetic": True,
        }
        out.append({"direction": direction, "isotropy": iso_row})
    return out


def opd_rows(
    old_id: int,
    *,
    parent_sg: int | None = None,
    gid: int | None = None,
    k_params: dict[str, str] | None = None,
    display_k_params: dict[str, str] | None = None,
    selected_orderparam: int | str | None = None,
    parent_setting_id: int | None = None,
    parent_cell: tuple[float, float, float, float, float, float] | None = None,
) -> list[dict[str, Any]]:
    dec = source_tables()
    if int(old_id) <= 0:
        if parent_sg is None or gid is None:
            return []
        return _dynamic_opd_rows(
            dec,
            parent_sg=int(parent_sg),
            gid=int(gid),
            k_params=k_params,
            display_k_params=display_k_params,
            selected_orderparam=selected_orderparam,
            parent_setting_id=parent_setting_id,
            parent_cell=parent_cell,
        )
    if old_id <= 0:
        return []
    image = dec.image_record(old_id)
    directions = _image_directions(dec, int(image["id"]))
    by_direction: dict[int, list[dict[str, Any]]] = {}
    for row in _isotropy_rows(
        dec,
        old_id,
        k_params=k_params,
        display_k_params=display_k_params,
        parent_setting_id=parent_setting_id,
        parent_cell=parent_cell,
    ):
        by_direction.setdefault(int(row["direction_id"]), []).append(row)
    rows: list[dict[str, Any]] = []
    seen: set[int] = set()
    for direction in directions:
        did = int(direction["id"])
        seen.add(did)
        concrete = by_direction.get(did) or [None]
        for iso_row in concrete:
            rows.append({"direction": direction, "isotropy": iso_row})
    for did in sorted(set(by_direction) - seen):
        direction = _direction(dec, did)
        for iso_row in by_direction[did]:
            rows.append({"direction": direction, "isotropy": iso_row})
    return rows
