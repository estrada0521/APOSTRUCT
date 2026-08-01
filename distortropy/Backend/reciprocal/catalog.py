"""Source-derived reciprocal-space catalog and display coordinates."""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from functools import lru_cache
import math
import re
from typing import Any, Sequence

from distortropy.Backend.exactmath import fraction_matrix_inverse3
from distortropy.Backend.source.iso_data import (
    reciprocal_to_cinter_matrix_from_table,
)
from distortropy.Backend.source.tables import SourceTables, source_tables


GREEK = {"GM": "Γ", "DT": "Δ", "LD": "Λ", "SM": "Σ"}
PARAMS = ("", "α", "β", "γ")
DISPLAY_K_PARAMS = ("a", "b", "g")
TYPE_LABEL = {1: "real", 2: "type 2", 3: "type 3"}


@dataclass(frozen=True, slots=True)
class KCoordinateMap:
    """Exact affine coordinates behind one displayed K-family expression."""

    parameter_names: tuple[str, ...]
    origin: tuple[Fraction, Fraction, Fraction]
    columns: tuple[tuple[Fraction, Fraction, Fraction], ...]

    def evaluate(
        self,
        parameters: Sequence[Fraction],
    ) -> tuple[Fraction, Fraction, Fraction]:
        if len(parameters) != len(self.columns):
            raise ValueError("K-family parameter count does not match its coordinate map")
        return tuple(
            self.origin[axis]
            + sum(
                Fraction(parameters[index]) * column[axis]
                for index, column in enumerate(self.columns)
            )
            for axis in range(3)
        )  # type: ignore[return-value]


def _strip(value: Any) -> str:
    return str(value).strip()


def _frac(num: int, den: int) -> Fraction:
    return Fraction(int(num), int(den or 1))


def _fmt_frac(value: Fraction) -> str:
    return str(value.numerator) if value.denominator == 1 else f"{value.numerator}/{value.denominator}"


def _term(value: Fraction, symbol: str) -> str:
    if not symbol:
        return _fmt_frac(value)
    if value == 1:
        return symbol
    if value == -1:
        return "-" + symbol
    return f"{_fmt_frac(value)}{symbol}"


def _fmt_component(records: tuple[tuple[Fraction, Fraction, Fraction], ...], axis: int) -> str:
    parts: list[str] = []
    for index, vector in enumerate(records):
        value = vector[axis]
        if value == 0:
            continue
        text = _term(value, PARAMS[index])
        parts.append(text if not parts or text.startswith("-") else f"+{text}")
    return "".join(parts) if parts else "0"


def _fmt_vector(records: tuple[tuple[Fraction, Fraction, Fraction], ...]) -> str:
    return "(" + ", ".join(_fmt_component(records, axis) for axis in range(3)) + ")"


def evaluate_k_component(text: str, k_params: dict[str, str] | None) -> Fraction | None:
    """Evaluate one Source k-vector component from exact request parameters."""

    if not k_params:
        return None
    total = Fraction(0)
    terms = re.findall(r"[+-]?[^+-]+", str(text).replace(" ", ""))
    if not terms:
        return None
    for term in terms:
        if not term:
            continue
        sign = -1 if term.startswith("-") else 1
        body = term[1:] if term[0:1] in {"+", "-"} else term
        match = re.fullmatch(r"(?:(\d+(?:/\d+)?)?)?([abg])", body)
        if match:
            raw_value = k_params.get(match.group(2))
            if raw_value is None and k_params.get("_single") is not None:
                raw_value = k_params.get("_single")
            if raw_value is None or str(raw_value).strip() == "":
                return None
            coefficient = Fraction(match.group(1)) if match.group(1) else Fraction(1)
            total += sign * coefficient * Fraction(str(raw_value).strip())
            continue
        try:
            total += sign * Fraction(body)
        except (ValueError, ZeroDivisionError):
            return None
    return total


def substitute_k_vector(vector: str, k_params: dict[str, str] | None) -> str:
    """Substitute exact request parameters into one displayed k vector."""

    if not k_params:
        return vector
    if "_single" not in k_params:
        public_values = [
            str(k_params[key]).strip()
            for key in DISPLAY_K_PARAMS
            if k_params.get(key) is not None and str(k_params.get(key)).strip() != ""
        ]
        if len(set(public_values)) == 1:
            k_params = {**k_params, "_single": public_values[0]}
    text = str(vector).strip()
    match = re.fullmatch(r"\((.*)\)", text)
    if not match:
        return text
    values: list[str] = []
    for component in match.group(1).split(","):
        value = evaluate_k_component(component, k_params)
        if value is None:
            return text
        values.append(_fmt_frac(value))
    return "(" + ",".join(values) + ")"


def _permute_records(
    records: tuple[tuple[Fraction, Fraction, Fraction], ...],
    perm: tuple[int, int, int],
) -> tuple[tuple[Fraction, Fraction, Fraction], ...]:
    return tuple(tuple(record[index] for index in perm) for record in records)  # type: ignore[return-value]


def _apply_k_matrix(
    records: tuple[tuple[Fraction, Fraction, Fraction], ...],
    matrix: tuple[tuple[int, int, int], tuple[int, int, int], tuple[int, int, int]],
) -> tuple[tuple[Fraction, Fraction, Fraction], ...]:
    out: list[tuple[Fraction, Fraction, Fraction]] = []
    for index, record in enumerate(records):
        transformed = tuple(
            sum(Fraction(matrix[row][col]) * record[col] for col in range(3))
            for row in range(3)
        )
        if index > 0:
            nonzero = [abs(value.numerator) for value in transformed if value and value.denominator == 1]
            if nonzero:
                scale = nonzero[0]
                for value in nonzero[1:]:
                    scale = _gcd(scale, value)
                if scale > 1:
                    transformed = tuple(value / scale for value in transformed)
        out.append(transformed)  # type: ignore[arg-type]
    return tuple(out)


def _apply_fraction_k_matrix(
    records: tuple[tuple[Fraction, Fraction, Fraction], ...],
    matrix: tuple[tuple[Fraction, Fraction, Fraction], ...],
) -> tuple[tuple[Fraction, Fraction, Fraction], ...]:
    return tuple(
        tuple(sum(record[row] * matrix[row][col] for row in range(3)) for col in range(3))
        for record in records
    )  # type: ignore[return-value]


def _rref_k_directions(
    rows: tuple[tuple[Fraction, Fraction, Fraction], ...],
) -> tuple[tuple[tuple[Fraction, Fraction, Fraction], ...], tuple[int, ...]]:
    work = [list(row) for row in rows if any(row)]
    pivots: list[int] = []
    out_row = 0
    for col in range(3):
        pivot = next((row for row in range(out_row, len(work)) if work[row][col]), None)
        if pivot is None:
            continue
        work[out_row], work[pivot] = work[pivot], work[out_row]
        scale = work[out_row][col]
        work[out_row] = [value / scale for value in work[out_row]]
        for row in range(len(work)):
            if row == out_row:
                continue
            factor = work[row][col]
            if factor:
                work[row] = [work[row][i] - factor * work[out_row][i] for i in range(3)]
        pivots.append(col)
        out_row += 1
        if out_row == len(work):
            break
    basis = tuple(tuple(row) for row in work[:out_row])
    return basis, tuple(pivots)  # type: ignore[return-value]


def _canonicalize_k_display_records(
    records: tuple[tuple[Fraction, Fraction, Fraction], ...],
) -> tuple[tuple[tuple[Fraction, Fraction, Fraction], ...], dict[int, str]]:
    if not records:
        return records, {}
    constant = list(records[0])
    basis, pivots = _rref_k_directions(records[1:])
    for index, pivot in enumerate(pivots):
        shift = constant[pivot]
        if not shift:
            continue
        constant = [constant[axis] - shift * basis[index][axis] for axis in range(3)]
    param_names = {
        index + 1: DISPLAY_K_PARAMS[pivot]
        for index, pivot in enumerate(pivots)
    }
    return (tuple(constant), *basis), param_names  # type: ignore[return-value]


def _display_k_records(
    records: tuple[tuple[Fraction, Fraction, Fraction], ...],
    lattice: int,
    symbol: str,
    label: str,
) -> tuple[tuple[Fraction, Fraction, Fraction], ...]:
    # Primitive monoclinic little-k data are held in the Miller-Love axis
    # convention, while the public surface displays Kovalev
    # b-unique coordinates.  In the decoded tables this is the cyclic
    # reciprocal-axis permutation (h,k,l)_web=(k,l,h)_ML.
    if int(lattice) == 2:
        return _permute_records(records, (1, 2, 0))
    if int(lattice) == 3:
        transformed = _apply_k_matrix(records, ((0, 1, -1), (0, 1, 0), (1, 0, 0)))
        if label == "V" and transformed and transformed[0] == (Fraction(-1, 2), Fraction(0), Fraction(0)):
            # Centered monoclinic V-star arms are displayed in the same
            # representative strip as the primary V arm.
            transformed = ((Fraction(-1, 2), Fraction(1, 2), Fraction(0)), *transformed[1:])
        if label in {"Y", "M"} and transformed and transformed[0][0] == 0 and transformed[0][1] == Fraction(1, 2):
            constant = (transformed[0][0], Fraction(1), transformed[0][2])
            transformed = (constant, *transformed[1:])
        return transformed
    if int(lattice) == 5 and symbol.startswith("C"):
        return _apply_k_matrix(records, ((1, 1, 0), (-1, 1, 0), (0, 0, 1)))
    if int(lattice) == 5 and symbol.startswith("A"):
        transformed = _apply_k_matrix(records, ((0, 0, 1), (-1, 1, 0), (1, 1, 0)))
        if label in {"Y", "R", "S", "T", "D", "H"} and transformed and transformed[0][2] > 0:
            constant = (transformed[0][0], transformed[0][1], -transformed[0][2])
            transformed = (constant, *transformed[1:])
        return transformed
    if int(lattice) in {6, 9}:
        return _apply_k_matrix(records, ((0, 1, 1), (1, 0, 1), (1, 1, 0)))
    if int(lattice) == 7:
        return _apply_k_matrix(records, ((-1, 1, 1), (1, -1, 1), (1, 1, -1)))
    return records


def _fmt_display_k_component(
    records: tuple[tuple[Fraction, Fraction, Fraction], ...],
    axis: int,
    param_names: dict[int, str],
) -> str:
    if not records:
        return "0"
    constant = records[0][axis]
    param_parts: list[str] = []
    for index, vector in enumerate(records[1:], start=1):
        value = vector[axis]
        if value == 0:
            continue
        symbol = param_names.get(index, DISPLAY_K_PARAMS[min(index - 1, 2)])
        param_parts.append(_term(value, symbol))
    if param_parts:
        if constant.denominator == 1 and constant > 1:
            constant = constant % 1
        if constant.denominator == 1 and constant > 0 and param_parts[0].startswith("-"):
            parts = list(param_parts)
            parts.append(_fmt_frac(constant))
            return "".join(part if index == 0 or part.startswith("-") else f"+{part}" for index, part in enumerate(parts))

    parts: list[str] = []
    if constant:
        parts.append(_fmt_frac(constant))
    for part in param_parts:
        parts.append(part if not parts or part.startswith("-") else f"+{part}")
    return "".join(parts) if parts else "0"


def _fmt_display_k_vector(
    records: tuple[tuple[Fraction, Fraction, Fraction], ...],
    lattice: int,
    symbol: str,
    label: str,
    dec: SourceTables | None = None,
    sg: int | None = None,
) -> str:
    display_records = _display_k_records(records, lattice, symbol, label)
    param_names: dict[int, str] = {}
    if dec is not None and sg is not None and int(lattice) in {9, 10, 13, 14}:
        display_records = _apply_fraction_k_matrix(
            records,
            reciprocal_to_cinter_matrix_from_table(dec.space, int(sg), "pml"),
        )
        display_records, param_names = _canonicalize_k_display_records(display_records)
    if display_records and all(value == 0 for value in display_records[0]) and len(display_records) > 1:
        if label == "GP":
            return "(a,b,g)"
        if int(lattice) == 3 and label == "B":
            return "(a,0,g)"
        params = display_records[1:]
        if len(params) == 2:
            if all(record[2] == 0 for record in params):
                return "(a,b,0)"
            if all(record[0] == 0 for record in params):
                return "(0,b,g)"
            if all(record[1] == 0 for record in params):
                return "(a,0,g)"
    if display_records and len(display_records) == 3:
        constant = display_records[0]
        params = display_records[1:]
        if constant[0] == 0 and constant[1] == 0 and constant[2] == Fraction(1, 2):
            if all(record[2] == 0 for record in params):
                return "(a,b,1/2)"
        if constant[0] == Fraction(1, 2) and constant[1] == 0 and constant[2] == 0:
            if all(record[0] == 0 for record in params):
                return "(1/2,b,g)"
    if not param_names:
        for index, vector in enumerate(display_records[1:], start=1):
            for axis, value in enumerate(vector):
                if value:
                    param_names[index] = DISPLAY_K_PARAMS[axis]
                    break
    return "(" + ",".join(_fmt_display_k_component(display_records, axis, param_names) for axis in range(3)) + ")"


def _parse_k_record(raw: list[int]) -> tuple[tuple[Fraction, Fraction, Fraction], ...]:
    vectors: list[tuple[Fraction, Fraction, Fraction]] = []
    for start in range(0, 16, 4):
        x, y, z, den = (int(v) for v in raw[start:start + 4])
        if den == 0:
            vectors.append((Fraction(0), Fraction(0), Fraction(0)))
        else:
            vectors.append((_frac(x, den), _frac(y, den), _frac(z, den)))
    return tuple(vectors)


def _coordinate_map(expression: str) -> KCoordinateMap:
    match = re.fullmatch(r"\((.*)\)", str(expression).strip())
    if match is None:
        raise ValueError("K-family display expression must contain three components")
    components = tuple(part.strip() for part in match.group(1).split(","))
    if len(components) != 3:
        raise ValueError("K-family display expression must contain three components")
    present = set(
        re.findall(r"(?<![A-Za-z])([abg])(?![A-Za-z])", match.group(1))
    )
    names = tuple(name for name in DISPLAY_K_PARAMS if name in present)

    def evaluate(values: tuple[Fraction, ...]) -> tuple[Fraction, Fraction, Fraction]:
        params = {
            name: str(value)
            for name, value in zip(names, values, strict=True)
        }
        if not params:
            params = {"_single": "0"}
        result = tuple(evaluate_k_component(component, params) for component in components)
        if any(value is None for value in result):
            raise ValueError("K-family display expression is not affine rational")
        return result  # type: ignore[return-value]

    zero = tuple(Fraction(0) for _ in names)
    origin = evaluate(zero)
    columns = []
    for index in range(len(names)):
        values = list(zero)
        values[index] = Fraction(1)
        point = evaluate(tuple(values))
        columns.append(
            tuple(point[axis] - origin[axis] for axis in range(3))
        )
    return KCoordinateMap(names, origin, tuple(columns))  # type: ignore[arg-type]


def _gcd(left: int, right: int) -> int:
    while right:
        left, right = right, left % right
    return abs(left)


def _k_slot_info(dec: SourceTables, sg: int, kslot: int) -> dict[str, Any]:
    lattice = int(dec.space["ispace_lattice"][sg - 1])
    symbol = dec.space_symbol(sg).strip()
    lattice_slot = (lattice - 1) * 27 + kslot - 1
    sg_slot = (sg - 1) * 27 + kslot - 1
    label = _strip(dec.little["little_k_label"][lattice_slot])
    raw_kov_id = int(dec.little["little_k_kov"][sg_slot])
    kov_id: int | None = raw_kov_id
    if raw_kov_id == 0 and label != "GP":
        kov_id = None
    dim = int(dec.little["little_k_dim"][lattice_slot])
    raw = [int(x) for x in dec.little["little_k"][lattice_slot * 16:(lattice_slot + 1) * 16]]
    basis = _parse_k_record(raw)
    return {
        "kslot": kslot,
        "kid": kov_id,
        "label": label,
        "display": GREEK.get(label, label),
        "dimension": dim,
        "basis_records": [[_fmt_frac(x) for x in row] for row in basis],
        "kvector": _fmt_vector(basis[:dim + 1]),
        "display_kvector": _fmt_display_k_vector(basis[:dim + 1], lattice, symbol, label, dec, sg),
        "star_size": int(dec.little["little_k_star_count"][sg_slot]),
        "minus_k": int(dec.little["little_k_star_minusk"][sg_slot]),
        "little_order": int(dec.little["little_ops_count"][sg_slot]),
    }


def _star_vectors(dec: SourceTables, sg: int, kslot: int, dim: int) -> list[str]:
    sg_slot = (sg - 1) * 27 + kslot - 1
    count = int(dec.little["little_k_star_count"][sg_slot])
    if count <= 0:
        return []
    ml_ptr = int(dec.little["little_k_star_ml_pointer"][sg_slot])
    if ml_ptr <= 0:
        return []
    if dim == 0:
        raw = dec.little["little_k_star_special"]
        arms = []
        for offset in range(count):
            item = [int(x) for x in raw[4 * (ml_ptr - 1 + offset):4 * (ml_ptr + offset)]]
            if len(item) != 4:
                break
            den = item[3] or 1
            arms.append(_fmt_vector(((_frac(item[0], den), _frac(item[1], den), _frac(item[2], den)),)))
        return arms
    raw = dec.little["little_k_star"]
    arms = []
    for offset in range(count):
        item = [int(x) for x in raw[16 * (ml_ptr - 1 + offset):16 * (ml_ptr + offset)]]
        if len(item) != 16:
            break
        arms.append(_fmt_vector(_parse_k_record(item)[:dim + 1]))
    return arms


def _display_star_vectors(dec: SourceTables, sg: int, kslot: int, dim: int) -> list[str]:
    sg_slot = (sg - 1) * 27 + kslot - 1
    count = int(dec.little["little_k_star_count"][sg_slot])
    if count <= 0:
        return []
    ml_ptr = int(dec.little["little_k_star_ml_pointer"][sg_slot])
    if ml_ptr <= 0:
        return []
    lattice = int(dec.space["ispace_lattice"][sg - 1])
    lattice_slot = (lattice - 1) * 27 + kslot - 1
    symbol = dec.space_symbol(sg).strip()
    label = _strip(dec.little["little_k_label"][lattice_slot])
    if dim == 0:
        raw = dec.little["little_k_star_special"]
        arms = []
        for offset in range(count):
            item = [int(x) for x in raw[4 * (ml_ptr - 1 + offset):4 * (ml_ptr + offset)]]
            if len(item) != 4:
                break
            den = item[3] or 1
            record = ((_frac(item[0], den), _frac(item[1], den), _frac(item[2], den)),)
            arms.append(_fmt_display_k_vector(record, lattice, symbol, label, dec, sg))
        if int(lattice) == 3 and label == "L" and len(arms) == 2:
            arms[1] = "(-1/2,1/2,1/2)"
        return arms
    raw = dec.little["little_k_star"]
    arms = []
    for offset in range(count):
        item = [int(x) for x in raw[16 * (ml_ptr - 1 + offset):16 * (ml_ptr + offset)]]
        if len(item) != 16:
            break
        arms.append(_fmt_display_k_vector(_parse_k_record(item)[:dim + 1], lattice, symbol, label, dec, sg))
    return arms


def _k_label(label: str) -> str:
    match = re.match(r"^[A-Za-z]+", label)
    return match.group(0) if match else label


def _little_irreps(dec: SourceTables, sg: int, kslot: int) -> list[dict[str, Any]]:
    out: list[dict[str, Any]] = []
    for gid, row_sg in enumerate(dec.little["little_irr_space_group"], start=1):
        if int(row_sg) != sg:
            continue
        if int(dec.little["little_irr_k"][gid - 1]) != kslot:
            continue
        little = dec.little_record_by_gid(gid)
        nmod = int(dec.little["little_k_star_count"][(sg - 1) * 27 + kslot - 1])
        small_dim = little.full_dim // nmod if nmod and little.full_dim % nmod == 0 else None
        image = dec.image_record(little.old_id) if little.old_id > 0 else None
        out.append({
            "gid": little.gid,
            "old_id": little.old_id,
            "symbol": little.label,
            "kov_label": _strip(dec.little["little_irr_full_kov_label"][gid - 1]),
            "k_label": _k_label(little.label),
            "full_dim": little.full_dim,
            "little_dim": small_dim,
            "type": little.irrep_type,
            "type_label": TYPE_LABEL.get(little.irrep_type, f"type {little.irrep_type}"),
            "lif": little.lif,
            "real_pointer": little.real_pointer,
            "real2_pointer": little.real2_pointer,
            "image": image,
        })
    return out


def _has_extra_k_slot(dec: SourceTables, sg: int, kslot: int) -> bool:
    lattice = int(dec.space["ispace_lattice"][sg - 1])
    lattice_slot = (lattice - 1) * 27 + kslot - 1
    sg_slot = (sg - 1) * 27 + kslot - 1
    label = _strip(dec.little["little_k_label"][lattice_slot])
    if not label:
        return False
    if int(dec.little["little_k_star_count"][sg_slot]) > 0:
        return True
    return any(
        int(row_sg) == sg and int(dec.little["little_irr_k"][gid - 1]) == kslot
        for gid, row_sg in enumerate(dec.little["little_irr_space_group"], start=1)
    )


@lru_cache(maxsize=256)
def k_coordinate_map_for_slot(sg_number: int, kslot: int) -> KCoordinateMap:
    """Return the exact displayed affine map for one Source K slot."""

    sg = int(sg_number)
    slot = int(kslot)
    dec = source_tables()
    info = _k_slot_info(dec, sg, slot)
    return _coordinate_map(str(info["display_kvector"]))


def source_kparam_from_visible(
    dec: SourceTables,
    *,
    sg: int,
    kslot: int,
    k_params: dict[str, str] | None,
) -> tuple[int, int, int, int] | None:
    """Convert public K parameters to Source's internal values."""

    if not k_params:
        return None
    lattice = int(dec.space["ispace_lattice"][int(sg) - 1])
    lattice_slot = (lattice - 1) * 27 + int(kslot) - 1
    if int(dec.little["little_k_dim"][lattice_slot]) <= 0:
        return None
    display_vector = str(_k_slot_info(dec, int(sg), int(kslot))["display_kvector"])
    match = re.fullmatch(r"\((.*)\)", display_vector.strip())
    if match is None:
        return None
    visible = tuple(
        evaluate_k_component(component, k_params)
        for component in match.group(1).split(",")
    )
    if len(visible) != 3 or any(value is None for value in visible):
        return None
    sg_slot = (int(sg) - 1) * 27 + int(kslot) - 1
    pointer = int(dec.little.get("little_k_star_conv2ml_pointer", [0])[sg_slot])
    if pointer <= 0:
        return None
    raw = [
        int(value)
        for value in dec.little["little_k_star"][16 * (pointer - 1):16 * pointer]
    ]
    if len(raw) != 16:
        return None
    denominator = int(raw[15] or 1)
    converted = tuple(
        Fraction(raw[12 + col], denominator)
        + sum(
            Fraction(visible[row]) * Fraction(raw[4 * row + col], denominator)
            for row in range(3)
        )
        for col in range(3)
    )
    common_denominator = 1
    for value in converted:
        common_denominator = math.lcm(common_denominator, value.denominator)
    return (
        int(converted[0] * common_denominator),
        int(converted[1] * common_denominator),
        int(converted[2] * common_denominator),
        common_denominator,
    )


def source_k_little_order(
    dec: SourceTables,
    *,
    sg: int,
    kslot: int,
    parameters: Sequence[Fraction],
) -> int:
    """Return the exact little-group order at one Source K-family point."""

    lattice = int(dec.space["ispace_lattice"][int(sg) - 1])
    lattice_slot = (lattice - 1) * 27 + int(kslot) - 1
    dimension = int(dec.little["little_k_dim"][lattice_slot])
    if len(parameters) != dimension:
        raise ValueError("K-family parameter count does not match its Source slot")
    raw = [
        int(value)
        for value in dec.little["little_k"][lattice_slot * 16:(lattice_slot + 1) * 16]
    ]
    basis = _parse_k_record(raw)
    vector = tuple(
        basis[0][axis]
        + sum(
            Fraction(parameters[index]) * basis[index + 1][axis]
            for index in range(dimension)
        )
        for axis in range(3)
    )
    units = (
        (Fraction(1), Fraction(0), Fraction(0)),
        (Fraction(0), Fraction(1), Fraction(0)),
        (Fraction(0), Fraction(0), Fraction(1)),
    )
    point_operations: list[int] = []
    for record in dec.generate_space_group_records(int(sg)):
        point_op = int(record[4])
        if point_op not in point_operations:
            point_operations.append(point_op)
    order = 0
    for point_op in point_operations:
        columns = [dec.vrot_fraction(int(sg), point_op, unit) for unit in units]
        direct = tuple(
            tuple(columns[col][row] for col in range(3))
            for row in range(3)
        )
        inverse = fraction_matrix_inverse3(direct)
        rotated = tuple(
            sum(inverse[col][row] * vector[col] for col in range(3))
            for row in range(3)
        )
        if all((rotated[axis] - vector[axis]).denominator == 1 for axis in range(3)):
            order += 1
    return order


def ensure_source_k_not_special(
    dec: SourceTables,
    *,
    sg: int,
    kslot: int,
    label: str,
    parameters: Sequence[Fraction],
) -> None:
    """Mirror Source ``is_k_special`` by detecting a larger little group."""

    lattice = int(dec.space["ispace_lattice"][int(sg) - 1])
    lattice_slot = (lattice - 1) * 27 + int(kslot) - 1
    expected = int(
        dec.little["little_ops_count"][(int(sg) - 1) * 27 + int(kslot) - 1]
    )
    if int(dec.little["little_k_dim"][lattice_slot]) <= 0:
        return
    if source_k_little_order(
        dec,
        sg=int(sg),
        kslot=int(kslot),
        parameters=parameters,
    ) > expected:
        raise ValueError(
            f"{label} parameters are at a special K point; select the corresponding fixed K label"
        )


def kpoints(sg_number: int) -> dict[str, Any]:
    dec = source_tables()
    sg = int(sg_number)
    lattice = int(dec.space["ispace_lattice"][sg - 1])
    count = int(dec.little["little_k_count"][lattice - 1])
    point_group = int(dec.space["ispace_point_group"][sg - 1])
    out = []
    for kslot in range(1, 28):
        is_extra = kslot > count
        if is_extra and not _has_extra_k_slot(dec, sg, kslot):
            continue
        info = _k_slot_info(dec, sg, kslot)
        irreps = _little_irreps(dec, sg, kslot)
        row = {
            **info,
            "sg_specific": is_extra,
            "star": _star_vectors(dec, sg, kslot, int(info["dimension"])),
            "n_irreps": len(irreps),
            "irreps": irreps,
        }
        if not is_extra:
            out.append(row)
            continue
        insert_at = len(out)
        for index, item in enumerate(out):
            if int(item["dimension"]) > int(row["dimension"]):
                insert_at = index
                break
        out.insert(insert_at, row)
    return {
        "source": "Source/data_*",
        "space_group": {
            "number": sg,
            "symbol": dec.space_symbol(sg),
            "ml_symbol": _strip(dec.space["space_label_ml"][sg - 1]),
            "crystal_class": _strip(dec.space["space_label"][sg - 1]),
            "lattice": lattice,
            "lattice_type": _strip(dec.space["lattice_label"][lattice - 1]),
            "point_group": point_group,
            "point_group_order": int(dec.space["ipoint_group_order"][point_group - 1]),
        },
        "kpoints": out,
    }
