"""Mode request-context helpers."""

from __future__ import annotations

from collections.abc import Sequence
from fractions import Fraction
from numbers import Integral, Real
import re
from typing import Any
import numpy as np
from ISODISTORT.Assembled.Backend.fraction_expression import evaluate_fraction_expression
from ISODISTORT.Assembled.Backend.reciprocal import catalog as reciprocal_catalog
from ISODISTORT.Assembled.Backend.modes.engine.decoder import ModeDataDecoder
from ISODISTORT.Assembled.Backend.modes.engine.input import Case

from ISODISTORT.Assembled.Backend.modes.common import (
    _isotropy_from_opd_row,
    _k_params,
)

def _solve_linear_rational(
    rows: list[list[Fraction]],
    rhs: list[Fraction],
    width: int,
) -> tuple[Fraction, ...] | None:
    matrix = [list(row) + [rhs[index]] for index, row in enumerate(rows)]
    pivot_row = 0
    pivots: list[int] = []
    for col in range(width):
        pivot = next((row for row in range(pivot_row, len(matrix)) if matrix[row][col] != 0), None)
        if pivot is None:
            continue
        matrix[pivot_row], matrix[pivot] = matrix[pivot], matrix[pivot_row]
        factor = matrix[pivot_row][col]
        matrix[pivot_row] = [value / factor for value in matrix[pivot_row]]
        for row_index, row in enumerate(matrix):
            if row_index == pivot_row or row[col] == 0:
                continue
            factor = row[col]
            matrix[row_index] = [row[i] - factor * matrix[pivot_row][i] for i in range(width + 1)]
        pivots.append(col)
        pivot_row += 1
    for row in matrix:
        if all(row[col] == 0 for col in range(width)) and row[width] != 0:
            return None
    if len(pivots) < width:
        return None
    solution = [Fraction(0) for _ in range(width)]
    for row_index, col in enumerate(pivots):
        if col < width:
            solution[col] = matrix[row_index][width]
    return tuple(solution)



def _source_kparam_record(selected_opd: dict[str, Any] | None) -> tuple[int, int, int, int] | None:
    iso = _isotropy_from_opd_row(selected_opd)
    if not iso:
        return None
    raw = iso.get("source_kparam")
    if (
        isinstance(raw, (str, bytes, bytearray))
        or not isinstance(raw, Sequence)
        or len(raw) != 4
        or any(
            isinstance(value, bool) or not isinstance(value, Integral)
            for value in raw
        )
        or raw[3] == 0
    ):
        return None
    return tuple(int(value) for value in raw)  # type: ignore[return-value]



def _selected_dynamic_gid(
    selected_irrep: dict[str, Any],
    selected_opd: dict[str, Any] | None,
) -> int | None:
    for raw in (selected_irrep.get("gid"), selected_irrep.get("dynamic_gid")):
        if not isinstance(raw, bool) and isinstance(raw, Integral) and raw > 0:
            return int(raw)
        try:
            if int(raw) > 0:
                return None
        except OverflowError:
            return None
        except (TypeError, ValueError):
            pass
    iso = _isotropy_from_opd_row(selected_opd)
    if iso:
        raw = iso.get("dynamic_gid")
        if not isinstance(raw, bool) and isinstance(raw, Integral) and raw > 0:
            return int(raw)
    return None



def _finite_real_matrix(value: object) -> np.ndarray | None:
    if (
        isinstance(value, (str, bytes, bytearray))
        or not isinstance(value, Sequence)
        or not value
    ):
        return None
    rows: list[list[float]] = []
    width: int | None = None
    for raw_row in value:
        if (
            isinstance(raw_row, (str, bytes, bytearray))
            or not isinstance(raw_row, Sequence)
            or not raw_row
        ):
            return None
        if width is None:
            width = len(raw_row)
        elif len(raw_row) != width:
            return None
        if any(isinstance(item, bool) or not isinstance(item, Real) for item in raw_row):
            return None
        try:
            rows.append([float(item) for item in raw_row])
        except (OverflowError, TypeError, ValueError):
            return None
    matrix = np.asarray(rows, dtype=float)
    if matrix.ndim != 2 or not np.isfinite(matrix).all():
        return None
    return matrix


def _exact_typed_transpose(
    rows: Sequence[Sequence[Real]],
    transpose: Sequence[Sequence[Real]],
) -> bool:
    return all(
        type(transpose[column][row]) is type(rows[row][column])
        and transpose[column][row] == rows[row][column]
        for row in range(len(rows))
        for column in range(len(rows[row]))
    )


def _complete_fixed_secondary_type3_source_rows(
    decoder: ModeDataDecoder,
    *,
    coupled: bool,
    spec: dict[str, Any],
) -> bool:
    if coupled is not False or spec.get("primary") is not False:
        return False
    gid = spec.get("gid")
    old_id = spec.get("old_id")
    if any(
        isinstance(value, bool)
        or not isinstance(value, Integral)
        or value <= 0
        for value in (gid, old_id)
    ):
        return False
    try:
        little = decoder.little_record_by_gid(int(gid))
    except (IndexError, KeyError, TypeError, ValueError):
        return False
    if (
        getattr(little, "gid", None) != gid
        or getattr(little, "old_id", None) != old_id
        or getattr(little, "irrep_type", None) != 3
    ):
        return False

    source_rows = _finite_real_matrix(spec.get("source_numeric_rows"))
    direction = _finite_real_matrix(spec.get("direction_matrix"))
    if source_rows is None or direction is None:
        return False
    row_count, coordinate_count = source_rows.shape
    raw_source_rows = spec["source_numeric_rows"]
    raw_direction = spec["direction_matrix"]
    if (
        row_count < 4
        or row_count % 2 != 0
        or coordinate_count < row_count
        or direction.shape != (coordinate_count, row_count)
        or not _exact_typed_transpose(raw_source_rows, raw_direction)
    ):
        return False
    try:
        if (
            np.linalg.matrix_rank(source_rows, tol=1e-10) != row_count
            or np.linalg.matrix_rank(direction, tol=1e-10) != row_count
        ):
            return False
    except np.linalg.LinAlgError:
        return False

    source_free_count = spec.get("source_free_count")
    return source_free_count is None or bool(
        not isinstance(source_free_count, bool)
        and isinstance(source_free_count, Integral)
        and source_free_count == row_count
    )


def _should_emit_all_source_opd_rows(
    decoder: ModeDataDecoder,
    *,
    coupled: bool,
    spec: dict[str, Any],
) -> bool:
    """Use every faithful row when complete mode details expose independent OPDs."""

    rows = list(spec.get("source_numeric_rows") or [])
    if not rows:
        return False
    try:
        source_rows = np.asarray(rows, dtype=float)
        source_free_count = int(spec.get("source_free_count") or 0)
        independent_selected_rows = bool(
            spec.get("primary")
            and source_free_count == len(rows)
            and source_rows.ndim == 2
            and np.linalg.matrix_rank(source_rows, tol=1e-10) == len(rows)
        )
    except (TypeError, ValueError, np.linalg.LinAlgError):
        independent_selected_rows = False
    if independent_selected_rows:
        return True
    old_id = spec.get("old_id")
    if isinstance(old_id, bool) or not isinstance(old_id, Integral) or old_id < 0:
        return False
    if old_id > 0:
        if _complete_fixed_secondary_type3_source_rows(
            decoder,
            coupled=coupled,
            spec=spec,
        ):
            return True
        return bool(
            coupled
            and int(spec.get("little_type") or 0) == 3
            and len(rows) > 1
        )
    return bool(not spec.get("primary") and len(rows) > 1)



def _source_record_to_case_k_params(
    decoder: ModeDataDecoder,
    *,
    gid: object,
    source_kparam: object,
) -> tuple[Fraction, ...] | None:
    if isinstance(gid, bool) or not isinstance(gid, Integral) or gid <= 0:
        return None
    if (
        isinstance(source_kparam, (str, bytes, bytearray))
        or not isinstance(source_kparam, Sequence)
        or len(source_kparam) != 4
        or any(
            isinstance(value, bool) or not isinstance(value, Integral)
            for value in source_kparam
        )
    ):
        return None
    denominator = int(source_kparam[3])
    if denominator == 0:
        return None
    gid = int(gid)
    target = tuple(Fraction(int(source_kparam[axis]), denominator) for axis in range(3))
    sg = int(decoder.iso.little["little_irr_space_group"][gid - 1])
    kslot = int(decoder.iso.little["little_irr_k"][gid - 1])
    lattice = int(decoder.iso.space["ispace_lattice"][sg - 1])
    lattice_slot = (lattice - 1) * 27 + kslot - 1
    dim = int(decoder.iso.little["little_k_dim"][lattice_slot])
    if dim <= 0:
        return ()
    sg_slot = (sg - 1) * 27 + kslot - 1
    pointer = int(decoder.iso.little.get("little_k_star_conv2ml_pointer", [0])[sg_slot])
    if pointer <= 0:
        return None
    raw = [int(value) for value in decoder.iso.little["little_k_star"][16 * (pointer - 1):16 * pointer]]
    if len(raw) != 16:
        return None
    denominator = int(raw[15] or 1)
    converted_rows = [
        [Fraction(raw[4 * visible_axis + converted_axis], denominator) for visible_axis in range(3)]
        for converted_axis in range(3)
    ]
    try:
        coordinate_map = reciprocal_catalog.k_coordinate_map_for_slot(sg, kslot)
    except ValueError:
        return None
    if len(coordinate_map.parameter_names) != dim:
        return None
    converted_zero = [
        Fraction(raw[12 + converted_axis], denominator)
        + sum(
            converted_rows[converted_axis][axis] * coordinate_map.origin[axis]
            for axis in range(3)
        )
        for converted_axis in range(3)
    ]
    converted_columns = [
        tuple(
            sum(
                converted_rows[converted_axis][axis]
                * coordinate_map.columns[param_index][axis]
                for axis in range(3)
            )
            for converted_axis in range(3)
        )
        for param_index in range(dim)
    ]
    equations = [
        [converted_columns[param_index][axis] for param_index in range(dim)]
        for axis in range(3)
    ]
    solution = _solve_linear_rational(
        equations,
        [target[axis] - converted_zero[axis] for axis in range(3)],
        dim,
    )
    return tuple(solution) if solution is not None else None


def _exact_source_k_identity(
    gid: object,
    source_kparam: object,
) -> tuple[int, tuple[int, int, int, int]] | None:
    if isinstance(gid, bool) or not isinstance(gid, Integral) or gid <= 0:
        return None
    if (
        isinstance(source_kparam, (str, bytes, bytearray))
        or not isinstance(source_kparam, Sequence)
        or len(source_kparam) != 4
        or any(
            isinstance(value, bool) or not isinstance(value, Integral)
            for value in source_kparam
        )
    ):
        return None
    values = tuple(int(value) for value in source_kparam)
    if values[3] == 0:
        return None
    return int(gid), (values[0], values[1], values[2], values[3])


def _source_display_case(
    decoder: ModeDataDecoder,
    case: Case,
    *,
    gid: object,
    source_kparam: object,
) -> Case | None:
    """Prove that ``case`` is the exact Source record named by a render spec."""

    identity = _exact_source_k_identity(gid, source_kparam)
    if identity is None:
        return None
    source_gid, source_record = identity
    try:
        source_sg = int(decoder.iso.little["little_irr_space_group"][source_gid - 1])
        kslot = int(decoder.iso.little["little_irr_k"][source_gid - 1])
        lattice = int(decoder.iso.space["ispace_lattice"][source_sg - 1])
        source_label = str(
            decoder.iso.little["little_k_label"][(lattice - 1) * 27 + kslot - 1]
        ).strip()
        source_params = _source_parameter_record_to_case_k_params(
            decoder,
            gid=source_gid,
            source_kparam=source_record,
        )
    except (IndexError, KeyError, TypeError, ValueError, ZeroDivisionError):
        return None
    if (
        source_sg != int(case.sg)
        or source_label != str(case.k_label)
        or source_params is None
        or tuple(case.k_params) != source_params
    ):
        return None
    return case


def _source_parameter_record_to_case_k_params(
    decoder: ModeDataDecoder,
    *,
    gid: int,
    source_kparam: tuple[int, int, int, int] | None,
) -> tuple[Fraction, ...] | None:
    """Return the internal parameters substituted into ``little_k_star``."""

    if source_kparam is None or int(source_kparam[3]) == 0:
        return None
    gid = int(gid)
    sg = int(decoder.iso.little["little_irr_space_group"][gid - 1])
    kslot = int(decoder.iso.little["little_irr_k"][gid - 1])
    lattice = int(decoder.iso.space["ispace_lattice"][sg - 1])
    dim = int(decoder.iso.little["little_k_dim"][(lattice - 1) * 27 + kslot - 1])
    denominator = int(source_kparam[3])
    return tuple(Fraction(int(source_kparam[index]), denominator) for index in range(max(0, dim)))



def _pml_vector_to_case_k_params(
    decoder: ModeDataDecoder,
    *,
    gid: int,
    reciprocal_vector_pml: tuple[Fraction, Fraction, Fraction] | list[Fraction],
) -> tuple[Fraction, ...] | None:
    """Express a Source PML reciprocal point in the visible k-family parameters."""

    gid = int(gid)
    sg = int(decoder.iso.little["little_irr_space_group"][gid - 1])
    kslot = int(decoder.iso.little["little_irr_k"][gid - 1])
    lattice = int(decoder.iso.space["ispace_lattice"][sg - 1])
    lattice_slot = (lattice - 1) * 27 + kslot - 1
    dim = int(decoder.iso.little["little_k_dim"][lattice_slot])
    if dim <= 0:
        return ()
    try:
        coordinate_map = reciprocal_catalog.k_coordinate_map_for_slot(sg, kslot)
    except ValueError:
        return None
    if len(coordinate_map.parameter_names) != dim:
        return None
    equations = [
        [coordinate_map.columns[param_index][axis] for param_index in range(dim)]
        for axis in range(3)
    ]
    target = decoder.reciprocal_setting_change_vector(
        sg,
        "pml",
        "cinter",
        tuple(Fraction(value) for value in reciprocal_vector_pml),
    )
    candidates: list[tuple[Fraction, ...]] = []
    for sx in (-1, 0, 1):
        for sy in (-1, 0, 1):
            for sz in (-1, 0, 1):
                shifted = (target[0] + sx, target[1] + sy, target[2] + sz)
                solution = _solve_linear_rational(
                    equations,
                    [
                        shifted[axis] - coordinate_map.origin[axis]
                        for axis in range(3)
                    ],
                    dim,
                )
                if solution is not None and all(Fraction(0) <= value < Fraction(1) for value in solution):
                    candidates.append(tuple(solution))
    return min(candidates) if candidates else None



def _case_k_params(
    decoder: ModeDataDecoder,
    *,
    selected_k: dict[str, Any],
    selected_irrep: dict[str, Any],
    selected_opd: dict[str, Any] | None,
    spec: dict[str, Any],
    k_params: dict[str, str] | None,
) -> tuple[Fraction, ...]:
    if not spec.get("coupled") and spec.get("case_k_params") is not None:
        return tuple(Fraction(str(value)) for value in spec.get("case_k_params") or ())
    spec_gid = spec.get("gid")
    spec_source_kparam = None
    if spec.get("magnetic") is False:
        spec_source_kparam = spec.get("carrier_source_kparam")
    if spec_source_kparam is None:
        spec_source_kparam = spec.get("source_kparam")
    if spec_gid is not None and spec_source_kparam is not None:
        try:
            source_record = tuple(int(value) for value in spec_source_kparam)
        except (TypeError, ValueError):
            source_record = ()
        if len(source_record) == 4:
            solved = _source_parameter_record_to_case_k_params(
                decoder,
                gid=int(spec_gid),
                source_kparam=source_record,  # type: ignore[arg-type]
            )
            if solved is not None:
                return solved
    if spec.get("case_k_params") is not None:
        return tuple(Fraction(str(value)) for value in spec.get("case_k_params") or ())
    if str(spec.get("k_label")) != str(selected_k.get("label")):
        return ()
    if spec.get("primary"):
        gid = _selected_dynamic_gid(selected_irrep, selected_opd)
        solved = (
            _source_parameter_record_to_case_k_params(
                decoder,
                gid=gid,
                source_kparam=_source_kparam_record(selected_opd),
            )
            if gid is not None
            else None
        )
        if solved is not None:
            return solved
    return _k_params(k_params)



def _evaluated_kvector_text(selected_k: dict[str, Any], k_params: dict[str, str] | None = None) -> str:
    text = str(selected_k.get("isodistort_kvector") or selected_k.get("kvector") or "?").strip().strip("()")
    params = {
        **(selected_k.get("display_parameters") or selected_k.get("parameters") or {}),
        **(k_params or {}),
    }
    components = [item.strip() for item in text.split(",")]
    if len(components) != 3 or not isinstance(params, dict):
        return text

    def evaluate(expression: str) -> Fraction:
        expression = re.sub(r"(?<![A-Za-z0-9_])(\d+(?:/\d+)?)([abg])\b", r"\1*\2", expression)
        message = f"unsupported k-vector expression {expression!r}"
        return evaluate_fraction_expression(
            expression,
            params,
            message,
            message,
        )

    values: list[str] = []
    for component in components:
        try:
            values.append(str(evaluate(component)))
        except (SyntaxError, ValueError, ZeroDivisionError):
            values.append(component)
    return ",".join(values)



def _opd_direction_text(
    selected_opd: dict[str, Any] | None,
    selected_irrep: dict[str, Any] | None = None,
    group_width: int | None = None,
) -> str:
    if not isinstance(selected_opd, dict):
        return ""
    direction = selected_opd.get("direction") or {}
    isotropy = selected_opd.get("isotropy") or {}
    text = str(isotropy.get("display_opd") or direction.get("opd") or "")
    little_dim = int((selected_irrep or {}).get("little_dim") or group_width or 0)
    if little_dim <= 0:
        return text
    body = text.strip().strip("()")
    parts = [part.strip() for part in body.replace(";", ",").split(",")]
    if len(parts) <= little_dim:
        return "(" + ",".join(parts) + ")"
    if len(parts) % little_dim:
        return text
    groups = [",".join(parts[index:index + little_dim]) for index in range(0, len(parts), little_dim)]
    return "(" + ";".join(groups) + ")"



def _direction_matrix_text(matrix: Any, tol: float = 1e-7, group_width: int | None = None) -> str:
    rows = [list(row) for row in (matrix or [])]
    if not rows:
        return ""
    free = len(rows[0])
    variables = [chr(ord("a") + index) for index in range(free)]
    components: list[str] = []
    for row in rows:
        terms: list[str] = []
        for coefficient, variable in zip(row, variables, strict=False):
            value = float(coefficient)
            if abs(value) <= tol:
                continue
            sign = "-" if value < 0 else "+"
            magnitude = abs(value)
            body = variable if abs(magnitude - 1.0) <= tol else f"{magnitude:.3f}".rstrip("0").rstrip(".") + variable
            terms.append(("-" if not terms and sign == "-" else "" if not terms else sign) + body)
        components.append("".join(terms) or "0")
    width = int(group_width or 0)
    if width > 0 and len(components) > width and len(components) % width == 0:
        return "(" + ";".join(",".join(components[index:index + width]) for index in range(0, len(components), width)) + ")"
    return "(" + ",".join(components) + ")"



def _little_star_arm_width(decoder: ModeDataDecoder, gid: int) -> int | None:
    record = decoder.little_record_by_gid(int(gid))
    sg = int(decoder.iso.little["little_irr_space_group"][int(gid) - 1])
    kslot = int(decoder.iso.little["little_irr_k"][int(gid) - 1])
    star_count = int(decoder.iso.little["little_k_star_count"][(sg - 1) * 27 + kslot - 1])
    if star_count <= 0 or int(record.full_dim) % star_count:
        return None
    return int(record.full_dim) // star_count


def _spec_opd_direction_text(
    decoder: ModeDataDecoder,
    *,
    spec: dict[str, Any],
    selected_direction: str,
) -> str:
    """Render secondary and coupled Source OPD coordinates in mode labels."""

    direction_matrix = spec.get("direction_matrix")
    if not direction_matrix or not (
        not spec.get("primary") or spec.get("coupled")
    ):
        return selected_direction
    return _direction_matrix_text(
        direction_matrix,
        group_width=(
            _little_star_arm_width(decoder, int(spec.get("gid")))
            if spec.get("gid") is not None
            else None
        ),
    )


def _spec_display_kvector(decoder: ModeDataDecoder, case: Case, spec: dict[str, Any], fallback: str) -> str:
    reciprocal_vector = spec.get("reciprocal_vector_pml")
    gid = spec.get("gid")
    raw_display_params = spec.get("display_case_k_params")
    selected_display_params = (
        tuple(raw_display_params)
        if spec.get("coupled") is True
        and spec.get("primary") is True
        and isinstance(raw_display_params, Sequence)
        and not isinstance(raw_display_params, (str, bytes, bytearray))
        and all(isinstance(value, Fraction) for value in raw_display_params)
        else None
    )
    source_case = (
        None
        if selected_display_params is not None
        else _source_display_case(
            decoder,
            case,
            gid=gid,
            source_kparam=spec.get("source_kparam"),
        )
    )
    if spec.get("primary") and not spec.get("coupled") and source_case is None:
        return fallback
    if reciprocal_vector is not None and gid is not None:
        try:
            sg = int(decoder.iso.little["little_irr_space_group"][int(gid) - 1])
            kslot = int(decoder.iso.little["little_irr_k"][int(gid) - 1])
            lattice = int(decoder.iso.space["ispace_lattice"][sg - 1])
            kdim = int(decoder.iso.little["little_k_dim"][(lattice - 1) * 27 + kslot - 1])
            if kdim <= 0:
                if spec.get("primary") and not spec.get("coupled"):
                    return fallback
                raise ValueError("fixed k point uses the Source display representative")
            if source_case is not None:
                visible = decoder.display_k_vector_from_case(source_case)
                return ",".join(str(Fraction(value)) for value in visible)
            coordinate_map = reciprocal_catalog.k_coordinate_map_for_slot(sg, kslot)
            requested = selected_display_params
            if requested is None:
                requested = spec.get("case_k_params")
            if requested is None:
                requested = _source_record_to_case_k_params(
                    decoder,
                    gid=gid,
                    source_kparam=spec.get("source_kparam"),
                )
            if (
                requested is not None
                and len(requested) >= len(coordinate_map.parameter_names)
            ):
                visible = coordinate_map.evaluate(
                    tuple(
                        Fraction(value)
                        for value in requested[: len(coordinate_map.parameter_names)]
                    )
                )
            else:
                raw_visible = decoder.display_k_vector_from_case(case)
                solved = _solve_linear_rational(
                    [
                        [
                            coordinate_map.columns[param_index][axis]
                            for param_index in range(len(coordinate_map.columns))
                        ]
                        for axis in range(3)
                    ],
                    [
                        Fraction(raw_visible[axis]) - coordinate_map.origin[axis]
                        for axis in range(3)
                    ],
                    len(coordinate_map.columns),
                )
                if solved is None:
                    raise ValueError("could not recover displayed k-family parameters")
                reduced = tuple(
                    value if Fraction(-1, 2) <= value < 1 else value % 1
                    for value in solved
                )
                visible = coordinate_map.evaluate(reduced)
            return ",".join(str(Fraction(value)) for value in visible)
        except (IndexError, KeyError, TypeError, ValueError, ZeroDivisionError):
            pass
    try:
        return ",".join(str(value) for value in decoder.display_k_vector_from_case(case))
    except Exception:
        return fallback
