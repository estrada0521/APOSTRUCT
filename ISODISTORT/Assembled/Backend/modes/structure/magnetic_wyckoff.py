"""Source-only magnetic Wyckoff orbit grouping for complete-mode structures."""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from functools import lru_cache
import itertools
import math
from typing import Sequence

from sympy import Matrix

from ISODISTORT.Assembled.Backend.exactmath import (
    fraction_matrix_inverse3 as _matrix_inverse,
    fraction_matrix_multiply3 as _matrix_multiply,
    fraction_row_multiply3 as _row_multiply,
)
from ISODISTORT.Assembled.Backend.source.magnetic import data as magnetic_data
from ISODISTORT.Assembled.Backend.source.tables import source_tables


FractionPoint = tuple[Fraction, Fraction, Fraction]
FractionRecord = tuple[int, int, int, int]
MagneticOperationRecord = tuple[int, int, int, int, int]
WyckoffFormula = tuple[FractionPoint, FractionPoint, FractionPoint, FractionPoint]
FractionMatrix = tuple[
    tuple[Fraction, Fraction, Fraction],
    tuple[Fraction, Fraction, Fraction],
    tuple[Fraction, Fraction, Fraction],
]


def _input_fraction(value: Fraction | int | float) -> Fraction:
    """Recover the small rational represented by a presentation float."""

    if isinstance(value, Fraction):
        return value
    fraction = Fraction(value)
    if isinstance(value, float) or fraction.denominator > 1_000_000_000:
        decimal = Fraction(str(round(float(value), 9)))
        crystallographic = decimal.limit_denominator(384)
        if abs(crystallographic - decimal) <= Fraction(1, 200_000_000):
            return crystallographic
        return decimal
    if fraction.denominator > 384:
        crystallographic = fraction.limit_denominator(384)
        if abs(crystallographic - fraction) <= Fraction(1, 200_000_000):
            return crystallographic
    return fraction


@dataclass(frozen=True)
class MagneticGroupSetting:
    magnetic_group: int
    bns_number: str
    magnetic_type: int
    ordinary_space_group: int
    reference_space_group: int
    reference_basis: tuple[int, ...]
    reference_origin: FractionRecord


@dataclass(frozen=True)
class MagneticWyckoffRow:
    source_ordinal: int
    label: str
    formula: WyckoffFormula
    formula_records: tuple[
        FractionRecord, FractionRecord, FractionRecord, FractionRecord
    ]


@dataclass(frozen=True)
class MagneticWyckoffBranch:
    """One Source-ordered equivalent-formula branch."""

    source_ordinal: int
    operation_index: int
    operation_record: MagneticOperationRecord
    formula: WyckoffFormula
    formula_records: tuple[
        FractionRecord, FractionRecord, FractionRecord, FractionRecord
    ]


@dataclass(frozen=True)
class MagneticWyckoffIdentification:
    """Branch-aware magnetic Wyckoff identification result."""

    row: MagneticWyckoffRow
    branch: MagneticWyckoffBranch
    parameters: FractionPoint
    representative: FractionPoint


@dataclass(frozen=True)
class MagneticOrbitGroup:
    source_ordinal: int
    ordinary_orbit_indices: tuple[int, ...]
    multiplicity: int
    wyckoff_label: str
    standard_representative: FractionPoint


def _fraction_record(values: Sequence[Fraction]) -> FractionRecord:
    denominator = 1
    for value in values:
        denominator = math.lcm(denominator, Fraction(value).denominator)
    numerators = [int(Fraction(value) * denominator) for value in values]
    common = math.gcd(
        math.gcd(abs(numerators[0]), abs(numerators[1])), abs(numerators[2])
    )
    common = math.gcd(common, denominator)
    if common:
        numerators = [value // common for value in numerators]
        denominator //= common
    return numerators[0], numerators[1], numerators[2], denominator


def _fraction_values(record: Sequence[int]) -> FractionPoint:
    denominator = int(record[3])
    if denominator == 0:
        return Fraction(0), Fraction(0), Fraction(0)
    return tuple(Fraction(int(record[axis]), denominator) for axis in range(3))  # type: ignore[return-value]


def _fold_point(point: Sequence[Fraction]) -> FractionPoint:
    return tuple(Fraction(value) % 1 for value in point)  # type: ignore[return-value]


def _periodic_point_close(
    left: Sequence[Fraction],
    right: Sequence[Fraction],
    tolerance: Fraction,
) -> bool:
    for axis in range(3):
        delta = (Fraction(left[axis]) - Fraction(right[axis])) % 1
        if min(delta, 1 - delta) > tolerance:
            return False
    return True


@lru_cache(maxsize=None)
def _cinter_centering_translations(sg: int) -> tuple[FractionPoint, ...]:
    count = _centering_count(int(sg))
    matrix, _origin = _pml_to_cinter_affine(int(sg))
    out: list[FractionPoint] = []
    for values in itertools.product(range(count), repeat=3):
        image = _fold_point(
            _row_multiply(tuple(Fraction(value) for value in values), matrix)
        )
        if image not in out:
            out.append(image)
        if len(out) == count:
            break
    if len(out) != count:
        raise ValueError(
            f"centering translation count {len(out)} != {count} for SG{sg}"
        )
    return tuple(out)


def _inverse3(matrix: Sequence[int]) -> tuple[tuple[Fraction, Fraction, Fraction], ...]:
    source = tuple(
        tuple(Fraction(int(matrix[3 * row + column])) for column in range(3))
        for row in range(3)
    )
    try:
        return _matrix_inverse(source)
    except ValueError:
        raise ValueError(f"singular magnetic reference basis: {tuple(matrix[:9])}")


def _pml_to_cinter_affine(sg: int) -> tuple[FractionMatrix, FractionPoint]:
    data = source_tables()
    return (
        tuple(
            tuple(Fraction(value) for value in row)
            for row in data.pml_to_cinter_matrix(int(sg))
        ),
        tuple(Fraction(value) for value in data.cml_to_cinter_origin(int(sg))),
    )  # type: ignore[return-value]


def _pml_point_operation_matrix(sg: int, point_op: int) -> FractionMatrix:
    data = source_tables()
    lattice = int(data.space["ispace_lattice"][int(sg) - 1])
    code = int(
        data.space["ipoint_op_psettings"][(lattice - 1) * 72 + int(point_op) - 1]
    )
    if code == 0:
        raise KeyError(f"missing PML point operation SG{sg} op{point_op}")
    entries: list[Fraction] = []
    for _ in range(9):
        entries.append(Fraction(code % 3 - 1))
        code //= 3
    return tuple(tuple(entries[row * 3 + col] for col in range(3)) for row in range(3))  # type: ignore[return-value]


def _cinter_point_operation(matrix: FractionMatrix) -> int:
    data = source_tables()
    target = tuple(value for row in matrix for value in row)
    if any(value.denominator != 1 for value in target):
        raise ValueError(f"non-integral cinter point operation: {matrix}")
    integer_target = tuple(int(value) for value in target)
    for point_op in range(1, 73):
        start = (point_op - 1) * 9
        if (
            tuple(int(value) for value in data.space["ipoint_op"][start : start + 9])
            == integer_target
        ):
            return point_op
    raise KeyError(f"cinter point operation not found: {integer_target}")


def _pml_point_operation(sg: int, matrix: FractionMatrix) -> int:
    target = tuple(value for row in matrix for value in row)
    for point_op in range(1, 73):
        try:
            candidate = _pml_point_operation_matrix(int(sg), point_op)
        except KeyError:
            continue
        if tuple(value for row in candidate for value in row) == target:
            return point_op
    raise KeyError(f"PML point operation not found for SG{sg}: {target}")


def _magnetic_point_operation(ordinary_point_op: int, antiunitary: bool) -> int:
    table = magnetic_data().table
    for index, mapped in enumerate(table["mag_point_op_mag2nonmag"]):
        if int(mapped) == int(ordinary_point_op) and bool(
            table["mag_point_op_r"][index]
        ) == bool(antiunitary):
            return index + 1
    raise KeyError(
        f"magnetic point operation not found: ordinary={ordinary_point_op}, r={antiunitary}"
    )


def _pml_record_to_cinter(
    setting: MagneticGroupSetting,
    record: MagneticOperationRecord,
) -> MagneticOperationRecord:
    table = magnetic_data().table
    sg = int(setting.ordinary_space_group)
    ordinary_pml = int(table["mag_point_op_mag2nonmag"][int(record[4]) - 1])
    pml_rotation = _pml_point_operation_matrix(sg, ordinary_pml)
    pml_to_cinter, origin = _pml_to_cinter_affine(sg)
    cinter_to_pml = _matrix_inverse(pml_to_cinter)
    cinter_rotation = _matrix_multiply(
        _matrix_multiply(cinter_to_pml, pml_rotation), pml_to_cinter
    )
    transformed_translation = _row_multiply(_fraction_values(record[:4]), pml_to_cinter)
    rotated_origin = _row_multiply(origin, cinter_rotation)
    cinter_translation = tuple(
        transformed_translation[axis] - rotated_origin[axis] + origin[axis]
        for axis in range(3)
    )
    ordinary_cinter = _cinter_point_operation(cinter_rotation)
    magnetic_cinter = _magnetic_point_operation(
        ordinary_cinter,
        bool(table["mag_point_op_r"][int(record[4]) - 1]),
    )
    x, y, z, denominator = _fraction_record(_fold_point(cinter_translation))
    return x, y, z, denominator, magnetic_cinter


def _cinter_record_to_pml(
    setting: MagneticGroupSetting,
    record: MagneticOperationRecord,
) -> MagneticOperationRecord:
    """Transform a cinter magnetic operation record back to PML coordinates."""

    data = source_tables()
    table = magnetic_data().table
    sg = int(setting.ordinary_space_group)
    ordinary_cinter = int(table["mag_point_op_mag2nonmag"][int(record[4]) - 1])
    start = (ordinary_cinter - 1) * 9
    cinter_rotation: FractionMatrix = tuple(
        tuple(
            Fraction(value)
            for value in data.space["ipoint_op"][start + row * 3 : start + row * 3 + 3]
        )
        for row in range(3)
    )  # type: ignore[assignment]
    pml_to_cinter, origin = _pml_to_cinter_affine(sg)
    cinter_to_pml = _matrix_inverse(pml_to_cinter)
    pml_rotation = _matrix_multiply(
        _matrix_multiply(pml_to_cinter, cinter_rotation),
        cinter_to_pml,
    )
    cinter_translation = _fraction_values(record[:4])
    rotated_origin = _row_multiply(origin, cinter_rotation)
    pml_translation = _row_multiply(
        tuple(
            cinter_translation[axis] + rotated_origin[axis] - origin[axis]
            for axis in range(3)
        ),
        cinter_to_pml,
    )
    ordinary_pml = _pml_point_operation(sg, pml_rotation)
    magnetic_pml = _magnetic_point_operation(
        ordinary_pml,
        bool(table["mag_point_op_r"][int(record[4]) - 1]),
    )
    x, y, z, denominator = _fraction_record(pml_translation)
    return x, y, z, denominator, magnetic_pml


@lru_cache(maxsize=None)
def magnetic_group_setting(magnetic_group: int) -> MagneticGroupSetting:
    table = magnetic_data().table
    group = int(magnetic_group)
    if group < 1 or group > len(table["mag_nlabel"]):
        raise IndexError(f"magnetic group out of range: {group}")
    index = group - 1
    magnetic_type = int(table["mag_type"][index])
    ordinary_sg = int(table["mag_space_group"][index])
    reference_sg = ordinary_sg
    basis = (1, 0, 0, 0, 1, 0, 0, 0, 1)
    origin: FractionRecord = (0, 0, 0, 1)
    if magnetic_type == 4:
        pointer = int(table["mag_type4_pointer"][index])
        if pointer <= 0:
            raise ValueError(f"type-4 magnetic group {group} has no setting pointer")
        reference_sg = int(table["mag_type4_supergroup"][pointer - 1])
        basis = tuple(
            int(value)
            for value in table["mag_type4_supergroup_basis"][
                (pointer - 1) * 9 : pointer * 9
            ]
        )
        origin = tuple(
            int(value)
            for value in table["mag_type4_supergroup_origin"][
                (pointer - 1) * 4 : pointer * 4
            ]
        )  # type: ignore[assignment]
    return MagneticGroupSetting(
        magnetic_group=group,
        bns_number=str(table["mag_nlabel"][index]).strip(),
        magnetic_type=magnetic_type,
        ordinary_space_group=ordinary_sg,
        reference_space_group=reference_sg,
        reference_basis=basis,
        reference_origin=origin,
    )


def _compose_operations(
    setting: MagneticGroupSetting,
    left: MagneticOperationRecord,
    right: MagneticOperationRecord,
) -> MagneticOperationRecord:
    data = source_tables()
    table = magnetic_data().table
    ordinary_left = int(table["mag_point_op_mag2nonmag"][int(left[4]) - 1])
    rotated_right = data.vrot_fraction(
        int(setting.ordinary_space_group), ordinary_left, _fraction_values(right[:4])
    )
    left_translation = _fraction_values(left[:4])
    translation = tuple(
        left_translation[axis] + rotated_right[axis] for axis in range(3)
    )
    point_op = int(
        table["mag_point_op_mlt"][(int(right[4]) - 1) * 144 + int(left[4]) - 1]
    )
    x, y, z, denominator = _fraction_record(_fold_point(translation))
    return x, y, z, denominator, point_op


@lru_cache(maxsize=None)
def _generate_magnetic_records_cached(
    magnetic_group: int,
) -> tuple[MagneticOperationRecord, ...]:
    table = magnetic_data().table
    group_setting = magnetic_group_setting(int(magnetic_group))
    group = int(magnetic_group)
    count = int(table["mag_gen_count"][group - 1])
    raw = table["mag_gen"][(group - 1) * 20 : group * 20]
    records = [
        tuple(int(value) for value in raw[index * 5 : index * 5 + 5])
        for index in range(count)
    ]
    by_point_op = {int(record[4]) for record in records}
    changed = True
    while changed:
        changed = False
        for left in tuple(records):
            for right in tuple(records):
                product = _compose_operations(group_setting, left, right)
                if int(product[4]) in by_point_op:
                    continue
                by_point_op.add(int(product[4]))
                records.append(product)
                changed = True
    return tuple(sorted(records, key=lambda record: int(record[4])))


@lru_cache(maxsize=None)
def _generate_binary_magnetic_records_cached(
    magnetic_group: int,
) -> tuple[MagneticOperationRecord, ...]:
    group_setting = magnetic_group_setting(int(magnetic_group))
    return tuple(
        _cinter_record_to_pml(
            group_setting,
            _pml_record_to_cinter(group_setting, record),
        )
        for record in _generate_magnetic_records_cached(int(magnetic_group))
    )


@lru_cache(maxsize=None)
def generate_magnetic_space_group_records(
    magnetic_group: int,
    *,
    setting: str = "binary",
) -> tuple[MagneticOperationRecord, ...]:
    """Return Source PML operation records or an explicit working setting."""

    requested = str(setting).strip().lower()
    if requested in {"binary", "pml"}:
        return _generate_binary_magnetic_records_cached(int(magnetic_group))
    closure = _generate_magnetic_records_cached(int(magnetic_group))
    if requested == "closure":
        return closure
    if requested != "cinter":
        raise KeyError(f"unsupported magnetic operation setting: {setting!r}")
    group_setting = magnetic_group_setting(int(magnetic_group))
    return tuple(_pml_record_to_cinter(group_setting, record) for record in closure)


def _apply_operation(
    setting: MagneticGroupSetting,
    record: MagneticOperationRecord,
    point: Sequence[Fraction | int],
    *,
    record_setting: str = "cinter",
) -> FractionPoint:
    data = source_tables()
    table = magnetic_data().table
    requested = str(record_setting).strip().lower()
    if requested == "pml":
        record = _pml_record_to_cinter(setting, record)
    elif requested != "cinter":
        raise KeyError(f"unsupported magnetic operation setting: {record_setting!r}")
    ordinary_point_op = int(table["mag_point_op_mag2nonmag"][int(record[4]) - 1])
    start = (ordinary_point_op - 1) * 9
    rotation = tuple(int(value) for value in data.space["ipoint_op"][start : start + 9])
    values = tuple(Fraction(value) for value in point)
    rotated = tuple(
        sum(values[row] * rotation[row * 3 + axis] for row in range(3))
        for axis in range(3)
    )
    translation = _fraction_values(record[:4])
    return _fold_point(tuple(rotated[axis] + translation[axis] for axis in range(3)))


def _transform_reference_formula(
    formula: WyckoffFormula,
    setting: MagneticGroupSetting,
) -> WyckoffFormula:
    if setting.magnetic_type != 4:
        return formula
    inverse = _inverse3(setting.reference_basis)
    origin = _fraction_values(setting.reference_origin)
    base = _row_multiply(
        tuple(formula[0][axis] - origin[axis] for axis in range(3)), inverse
    )
    parameters = tuple(_row_multiply(vector, inverse) for vector in formula[1:])
    return (base, *parameters)  # type: ignore[return-value]


def _nice_wyckoff_pml_formula(
    formula: WyckoffFormula,
    *,
    sg: int,
) -> WyckoffFormula:
    """Canonicalize a spatial formula with ``nice_wyckoff_`` rules."""

    pml_to_cinter, origin = _pml_to_cinter_affine(int(sg))
    cinter_to_pml = _matrix_inverse(pml_to_cinter)
    cinter_parameters = [_row_multiply(vector, pml_to_cinter) for vector in formula[1:]]
    active_slots = [
        index for index, vector in enumerate(cinter_parameters) if any(vector)
    ]
    parameter_matrix = Matrix(cinter_parameters)
    reduced, pivots = parameter_matrix.rref()
    canonical_cinter = [
        tuple(Fraction(reduced[row, col]) for col in range(3))
        for row in range(parameter_matrix.rows)
        if any(reduced[row, col] != 0 for col in range(3))
    ]
    if len(canonical_cinter) != len(active_slots):
        raise ValueError(
            f"dependent magnetic Wyckoff parameter slots for SG{sg}: {formula}"
        )

    transformed_base = _row_multiply(formula[0], pml_to_cinter)
    cinter_base = _fold_point(
        tuple(transformed_base[axis] + origin[axis] for axis in range(3))
    )
    for vector, pivot in zip(canonical_cinter, pivots):
        coefficient = cinter_base[int(pivot)]
        cinter_base = tuple(
            cinter_base[axis] - coefficient * vector[axis] for axis in range(3)
        )  # type: ignore[assignment]
    pml_base = _row_multiply(
        tuple(cinter_base[axis] - origin[axis] for axis in range(3)),
        cinter_to_pml,
    )

    canonical: list[FractionPoint] = [
        (Fraction(0), Fraction(0), Fraction(0)) for _ in range(3)
    ]
    for slot, vector in zip(active_slots, canonical_cinter):
        canonical[slot] = _row_multiply(vector, cinter_to_pml)
    return (pml_base, *canonical)  # type: ignore[return-value]


def _formula_pml_to_cinter(formula: WyckoffFormula, sg: int) -> WyckoffFormula:
    matrix, origin = _pml_to_cinter_affine(int(sg))
    transformed_base = _row_multiply(formula[0], matrix)
    base = tuple(transformed_base[axis] + origin[axis] for axis in range(3))
    parameters = tuple(_row_multiply(vector, matrix) for vector in formula[1:])
    return (base, *parameters)  # type: ignore[return-value]


@lru_cache(maxsize=None)
def _magnetic_wyckoff_rows_cached(
    magnetic_group: int,
    output_setting: str,
) -> tuple[MagneticWyckoffRow, ...]:
    data = source_tables()
    setting = magnetic_group_setting(int(magnetic_group))
    out: list[MagneticWyckoffRow] = []
    for ordinal, row in enumerate(
        data.wyckoff_rows(setting.reference_space_group), start=1
    ):
        formula = _nice_wyckoff_pml_formula(
            _transform_reference_formula(data.wyckoff_fraction_vectors(row), setting),
            sg=setting.ordinary_space_group,
        )
        if output_setting == "cinter":
            formula = _formula_pml_to_cinter(formula, setting.ordinary_space_group)
        formula_records = (
            _fraction_record(formula[0]),
            *tuple(
                (0, 0, 0, 0) if not any(vector) else _fraction_record(vector)
                for vector in formula[1:]
            ),
        )
        out.append(
            MagneticWyckoffRow(
                source_ordinal=ordinal,
                label=str(row.label),
                formula=formula,
                formula_records=formula_records,
            )
        )
    return tuple(out)


def magnetic_wyckoff_rows(
    magnetic_group: int,
    *,
    setting: str = "cinter",
) -> tuple[MagneticWyckoffRow, ...]:
    requested = str(setting).strip().lower()
    if requested in {"bns", "pml", "binary"}:
        requested = "bns"
    elif requested != "cinter":
        raise KeyError(f"unsupported magnetic Wyckoff setting: {setting!r}")
    return _magnetic_wyckoff_rows_cached(int(magnetic_group), requested)


def _reduce_bns_formula_base(base: FractionPoint, *, sg: int) -> FractionPoint:
    """Reduce a PML formula to cinter for ``get_wyckoff_all_magnetic_``."""

    pml_to_cinter, origin = _pml_to_cinter_affine(int(sg))
    cinter_to_pml = _matrix_inverse(pml_to_cinter)
    transformed = _row_multiply(base, pml_to_cinter)
    reduced = _fold_point(tuple(transformed[axis] + origin[axis] for axis in range(3)))
    return _row_multiply(
        tuple(reduced[axis] - origin[axis] for axis in range(3)),
        cinter_to_pml,
    )


def _same_equivalent_formula(left: WyckoffFormula, right: WyckoffFormula) -> bool:
    difference = tuple(left[0][axis] - right[0][axis] for axis in range(3))
    return all(value.denominator == 1 for value in difference) and left[1:] == right[1:]


@lru_cache(maxsize=None)
def _magnetic_wyckoff_formula_branches_cached(
    magnetic_group: int,
    row_ordinal: int,
    output_setting: str,
) -> tuple[MagneticWyckoffBranch, ...]:
    """Evaluate the spatial branch loop in ``get_wyckoff_all_magnetic_``."""

    group = int(magnetic_group)
    row_index = int(row_ordinal)
    setting = magnetic_group_setting(group)
    data = source_tables()
    table = magnetic_data().table
    rows = magnetic_wyckoff_rows(group, setting="bns")
    if row_index < 1 or row_index > len(rows):
        raise KeyError(
            f"magnetic Wyckoff row {row_index} out of range for group {group}"
        )
    row = rows[row_index - 1]
    records = generate_magnetic_space_group_records(group, setting="binary")
    if not records:
        raise ValueError(f"magnetic group {group} has no generated operations")

    branches: list[tuple[int, MagneticOperationRecord, WyckoffFormula]] = [
        (1, records[0], row.formula)
    ]
    for operation_index, record in enumerate(records[1:], start=2):
        ordinary_point_op = int(table["mag_point_op_mag2nonmag"][int(record[4]) - 1])
        translation = _fraction_values(record[:4])
        rotated_base = data.vrot_fraction(
            int(setting.ordinary_space_group), ordinary_point_op, row.formula[0]
        )
        base = tuple(rotated_base[axis] + translation[axis] for axis in range(3))
        parameters = tuple(
            data.vrot_fraction(
                int(setting.ordinary_space_group), ordinary_point_op, vector
            )
            for vector in row.formula[1:]
        )
        candidate: WyckoffFormula = (base, *parameters)  # type: ignore[assignment]
        if any(
            _same_equivalent_formula(candidate, existing)
            for _index, _record, existing in branches
        ):
            continue
        branches.append(
            (
                operation_index,
                record,
                (
                    _reduce_bns_formula_base(
                        base, sg=int(setting.ordinary_space_group)
                    ),
                    *parameters,
                ),  # type: ignore[arg-type]
            )
        )

    out: list[MagneticWyckoffBranch] = []
    for source_ordinal, (operation_index, record, formula) in enumerate(
        branches, start=1
    ):
        if output_setting == "cinter":
            formula = _formula_pml_to_cinter(formula, int(setting.ordinary_space_group))
        out.append(
            MagneticWyckoffBranch(
                source_ordinal=source_ordinal,
                operation_index=operation_index,
                operation_record=record,
                formula=formula,
                formula_records=tuple(_fraction_record(vector) for vector in formula),  # type: ignore[arg-type]
            )
        )
    return tuple(out)


def magnetic_wyckoff_formula_branches(
    magnetic_group: int,
    row_ordinal: int,
    *,
    setting: str = "bns",
) -> tuple[MagneticWyckoffBranch, ...]:
    """Return Source-ordered equivalent formulas for one magnetic row."""

    requested = str(setting).strip().lower()
    if requested in {"bns", "pml", "binary"}:
        requested = "bns"
    elif requested != "cinter":
        raise KeyError(f"unsupported magnetic Wyckoff branch setting: {setting!r}")
    return _magnetic_wyckoff_formula_branches_cached(
        int(magnetic_group), int(row_ordinal), requested
    )


@lru_cache(maxsize=4096)
def _formula_solver(
    parameters: tuple[FractionPoint, FractionPoint, FractionPoint],
) -> tuple[FractionMatrix, FractionMatrix, tuple[int, ...]]:
    """Return RREF and row operations for a Wyckoff parameter matrix."""

    reduced = [
        [Fraction(parameters[column][row]) for column in range(3)] for row in range(3)
    ]
    transform = [
        [Fraction(int(row == column)) for column in range(3)] for row in range(3)
    ]
    pivots: list[int] = []
    pivot_row = 0
    for column in range(3):
        source_row = next(
            (row for row in range(pivot_row, 3) if reduced[row][column] != 0),
            None,
        )
        if source_row is None:
            continue
        if source_row != pivot_row:
            reduced[pivot_row], reduced[source_row] = (
                reduced[source_row],
                reduced[pivot_row],
            )
            transform[pivot_row], transform[source_row] = (
                transform[source_row],
                transform[pivot_row],
            )
        scale = reduced[pivot_row][column]
        reduced[pivot_row] = [value / scale for value in reduced[pivot_row]]
        transform[pivot_row] = [value / scale for value in transform[pivot_row]]
        for row in range(3):
            if row == pivot_row:
                continue
            factor = reduced[row][column]
            if factor == 0:
                continue
            reduced[row] = [
                reduced[row][index] - factor * reduced[pivot_row][index]
                for index in range(3)
            ]
            transform[row] = [
                transform[row][index] - factor * transform[pivot_row][index]
                for index in range(3)
            ]
        pivots.append(column)
        pivot_row += 1
        if pivot_row == 3:
            break
    return (
        tuple(tuple(row) for row in reduced),
        tuple(tuple(row) for row in transform),
        tuple(pivots),
    )  # type: ignore[return-value]


def _formula_solution(
    formula: WyckoffFormula,
    point: FractionPoint,
) -> tuple[Fraction, Fraction, Fraction] | None:
    # A rank-zero Source formula is exactly one point on the periodic cell.
    if not any(value for parameter in formula[1:] for value in parameter):
        return (
            (Fraction(0), Fraction(0), Fraction(0))
            if _fold_point(formula[0]) == _fold_point(point)
            else None
        )
    reduced, transform, pivots = _formula_solver(formula[1:])
    shift_ranges = tuple(
        range(
            math.floor(formula[0][axis] - point[axis]) - 1,
            math.ceil(formula[0][axis] - point[axis]) + 2,
        )
        for axis in range(3)
    )
    for shift in itertools.product(*shift_ranges):
        target = tuple(
            point[axis] + shift[axis] - formula[0][axis] for axis in range(3)
        )
        transformed = tuple(
            sum(transform[row][axis] * target[axis] for axis in range(3))
            for row in range(3)
        )
        if any(
            all(reduced[row][column] == 0 for column in range(3))
            and transformed[row] != 0
            for row in range(3)
        ):
            continue
        values = [Fraction(0), Fraction(0), Fraction(0)]
        for row, column in enumerate(pivots):
            values[column] = transformed[row]
        return tuple(values)  # type: ignore[return-value]
    return None


def identify_magnetic_wyckoff_branch(
    magnetic_group: int,
    point: Sequence[Fraction | int],
    *,
    setting: str = "cinter",
) -> MagneticWyckoffIdentification:
    """Apply ``id_wyckoff2_magnetic_`` row/branch first-success traversal."""

    requested = str(setting).strip().lower()
    if requested in {"bns", "pml", "binary"}:
        requested = "bns"
    elif requested != "cinter":
        raise KeyError(f"unsupported magnetic identification setting: {setting!r}")
    folded = _fold_point(tuple(_input_fraction(value) for value in point))
    centerings = (
        _cinter_centering_translations(
            magnetic_group_setting(int(magnetic_group)).ordinary_space_group
        )
        if requested == "cinter"
        else ((Fraction(0), Fraction(0), Fraction(0)),)
    )
    rows = magnetic_wyckoff_rows(int(magnetic_group), setting=requested)
    for row in rows:
        for branch in magnetic_wyckoff_formula_branches(
            int(magnetic_group), row.source_ordinal, setting=requested
        ):
            for centering in centerings:
                formula = (
                    tuple(
                        branch.formula[0][axis] + centering[axis] for axis in range(3)
                    ),
                    *branch.formula[1:],
                )
                solution = _formula_solution(formula, folded)  # type: ignore[arg-type]
                if solution is None:
                    continue
                representative = tuple(
                    formula[0][axis]
                    + sum(
                        solution[index] * formula[index + 1][axis] for index in range(3)
                    )
                    for axis in range(3)
                )
                return MagneticWyckoffIdentification(
                    row=row,
                    branch=branch,
                    parameters=solution,
                    representative=_fold_point(representative),
                )
    raise KeyError(
        f"no magnetic Wyckoff formula for group {magnetic_group}, point={folded}"
    )


def identify_magnetic_wyckoff(
    magnetic_group: int,
    point: Sequence[Fraction | int],
) -> tuple[MagneticWyckoffRow, FractionPoint]:
    """Identify one cinter point, preserving the compact public return."""

    result = identify_magnetic_wyckoff_branch(int(magnetic_group), point)
    return result.row, result.representative


def magnetic_orbit_points(
    magnetic_group: int,
    point: Sequence[Fraction | int],
    *,
    magnetic_records: Sequence[MagneticOperationRecord] | None = None,
    record_setting: str = "cinter",
) -> tuple[FractionPoint, ...]:
    """Return one full centered orbit in Source magnetic-operation order."""

    setting = magnetic_group_setting(int(magnetic_group))
    records = (
        tuple(tuple(int(value) for value in record) for record in magnetic_records)
        if magnetic_records is not None
        else generate_magnetic_space_group_records(
            int(magnetic_group), setting=record_setting
        )
    )
    out: list[FractionPoint] = []
    for record in records:
        base_image = _apply_operation(
            setting, record, point, record_setting=record_setting
        )
        for centering in _cinter_centering_translations(setting.ordinary_space_group):
            image = _fold_point(
                tuple(base_image[axis] + centering[axis] for axis in range(3))
            )
            if image not in out:
                out.append(image)
    return tuple(out)


def _magnetic_unitary_orbit(
    setting: MagneticGroupSetting,
    point: FractionPoint,
    records: Sequence[MagneticOperationRecord],
    *,
    record_setting: str,
) -> tuple[FractionPoint, ...]:
    """Expand one orbit with the magnetic group's BNS unitary coset."""

    table = magnetic_data().table
    out: list[FractionPoint] = []
    centerings = _cinter_centering_translations(setting.ordinary_space_group)
    for record in records:
        if bool(table["mag_point_op_r"][int(record[4]) - 1]):
            continue
        base_image = _apply_operation(
            setting, record, point, record_setting=record_setting
        )
        for centering in centerings:
            image = _fold_point(
                tuple(base_image[axis] + centering[axis] for axis in range(3))
            )
            if image not in out:
                out.append(image)
    return tuple(out)


def _centering_count(sg: int) -> int:
    symbol = str(source_tables().space["space_label_bc"][int(sg) - 1]).strip()
    return {"P": 1, "A": 2, "B": 2, "C": 2, "I": 2, "F": 4, "R": 3}.get(symbol[:1], 1)


def group_ordinary_orbits_magnetic(
    magnetic_group: int,
    ordinary_orbits: Sequence[
        Sequence[Fraction | int] | Sequence[Sequence[Fraction | int]]
    ],
    *,
    magnetic_records: Sequence[MagneticOperationRecord] | None = None,
    record_setting: str = "cinter",
    presentation_tolerance: float | None = None,
) -> tuple[MagneticOrbitGroup, ...]:
    setting = magnetic_group_setting(int(magnetic_group))
    records = (
        tuple(tuple(int(value) for value in record) for record in magnetic_records)
        if magnetic_records is not None
        else generate_magnetic_space_group_records(
            int(magnetic_group), setting=record_setting
        )
    )
    normalized: list[tuple[FractionPoint, ...]] = []
    for item in ordinary_orbits:
        if len(item) == 3 and all(
            not isinstance(value, (list, tuple)) for value in item
        ):
            point = _fold_point(tuple(_input_fraction(value) for value in item))  # type: ignore[arg-type]
            normalized.append(
                _magnetic_unitary_orbit(
                    setting, point, records, record_setting=record_setting
                )
            )
        else:
            points = tuple(
                _fold_point(tuple(_input_fraction(value) for value in point))
                for point in item  # type: ignore[assignment]
            )
            normalized.append(tuple(dict.fromkeys(points)))

    point_to_orbit: dict[FractionPoint, int] = {}
    for index, orbit in enumerate(normalized):
        for point in orbit:
            previous = point_to_orbit.setdefault(point, index)
            if previous != index:
                raise ValueError(
                    f"ordinary input orbits overlap: {previous} and {index}, point={point}"
                )

    tolerance = (
        None
        if presentation_tolerance is None
        else Fraction(str(presentation_tolerance))
    )
    if tolerance is not None and tolerance <= 0:
        raise ValueError("presentation tolerance must be positive")

    def point_owner(point: FractionPoint) -> int | None:
        exact = point_to_orbit.get(point)
        if tolerance is None:
            return exact
        owners = {
            owner
            for candidate, owner in point_to_orbit.items()
            if _periodic_point_close(point, candidate, tolerance)
        }
        return next(iter(owners)) if len(owners) == 1 else None

    if tolerance is not None:
        for point, owner in point_to_orbit.items():
            if point_owner(point) != owner:
                raise ValueError(
                    "ordinary input orbits overlap within presentation tolerance: "
                    f"point={point}, tolerance={tolerance}"
                )

    adjacency: list[set[int]] = [{index} for index in range(len(normalized))]
    for index, orbit in enumerate(normalized):
        representative = orbit[0]
        for record in records:
            base_image = _apply_operation(
                setting, record, representative, record_setting=record_setting
            )
            for centering in _cinter_centering_translations(
                setting.ordinary_space_group
            ):
                image = _fold_point(
                    tuple(base_image[axis] + centering[axis] for axis in range(3))
                )
                target = point_owner(image)
                if target is None:
                    raise KeyError(
                        f"magnetic image is absent from ordinary coverage: "
                        f"group={magnetic_group}, orbit={index}, centering={centering}, image={image}"
                    )
                adjacency[index].add(target)
                adjacency[target].add(index)

    visited: set[int] = set()
    groups: list[MagneticOrbitGroup] = []
    for start in range(len(normalized)):
        if start in visited:
            continue
        component: set[int] = set()
        pending = [start]
        while pending:
            current = pending.pop()
            if current in component:
                continue
            component.add(current)
            pending.extend(adjacency[current] - component)
        visited.update(component)
        ordered_indices = tuple(sorted(component))
        representative = normalized[ordered_indices[0]][0]
        orbit_points: list[FractionPoint] = []
        for record in records:
            base_image = _apply_operation(
                setting, record, representative, record_setting=record_setting
            )
            for centering in _cinter_centering_translations(
                setting.ordinary_space_group
            ):
                image = _fold_point(
                    tuple(base_image[axis] + centering[axis] for axis in range(3))
                )
                if image not in orbit_points:
                    orbit_points.append(image)
        row, standard = identify_magnetic_wyckoff(int(magnetic_group), representative)
        groups.append(
            MagneticOrbitGroup(
                source_ordinal=len(groups) + 1,
                ordinary_orbit_indices=ordered_indices,
                multiplicity=len(orbit_points),
                wyckoff_label=str(row.label),
                standard_representative=standard,
            )
        )
    return tuple(groups)
