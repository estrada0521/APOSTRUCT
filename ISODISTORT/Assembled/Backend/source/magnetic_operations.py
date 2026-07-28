"""Exact Source magnetic operation records in named coordinate settings."""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from functools import lru_cache
import math
from typing import Sequence

from ISODISTORT.Assembled.Backend.exactmath import (
    fraction_matrix_inverse3,
    fraction_matrix_multiply3,
    fraction_row_multiply3,
)
from ISODISTORT.Assembled.Backend.source.magnetic import data as magnetic_data
from ISODISTORT.Assembled.Backend.source.tables import source_tables


FractionPoint = tuple[Fraction, Fraction, Fraction]
FractionRecord = tuple[int, int, int, int]
MagneticOperationRecord = tuple[int, int, int, int, int]
FractionMatrix = tuple[
    tuple[Fraction, Fraction, Fraction],
    tuple[Fraction, Fraction, Fraction],
    tuple[Fraction, Fraction, Fraction],
]


@dataclass(frozen=True)
class MagneticGroupSetting:
    magnetic_group: int
    bns_number: str
    magnetic_type: int
    ordinary_space_group: int
    reference_space_group: int
    reference_basis: tuple[int, ...]
    reference_origin: FractionRecord


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
    return tuple(
        tuple(entries[row * 3 + col] for col in range(3))
        for row in range(3)
    )  # type: ignore[return-value]


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
    cinter_to_pml = fraction_matrix_inverse3(pml_to_cinter)
    cinter_rotation = fraction_matrix_multiply3(
        fraction_matrix_multiply3(cinter_to_pml, pml_rotation), pml_to_cinter
    )
    transformed_translation = fraction_row_multiply3(
        _fraction_values(record[:4]), pml_to_cinter
    )
    rotated_origin = fraction_row_multiply3(origin, cinter_rotation)
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
    data = source_tables()
    table = magnetic_data().table
    sg = int(setting.ordinary_space_group)
    ordinary_cinter = int(table["mag_point_op_mag2nonmag"][int(record[4]) - 1])
    start = (ordinary_cinter - 1) * 9
    cinter_rotation: FractionMatrix = tuple(
        tuple(
            Fraction(value)
            for value in data.space["ipoint_op"][
                start + row * 3 : start + row * 3 + 3
            ]
        )
        for row in range(3)
    )  # type: ignore[assignment]
    pml_to_cinter, origin = _pml_to_cinter_affine(sg)
    cinter_to_pml = fraction_matrix_inverse3(pml_to_cinter)
    pml_rotation = fraction_matrix_multiply3(
        fraction_matrix_multiply3(pml_to_cinter, cinter_rotation),
        cinter_to_pml,
    )
    cinter_translation = _fraction_values(record[:4])
    rotated_origin = fraction_row_multiply3(origin, cinter_rotation)
    pml_translation = fraction_row_multiply3(
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
    """Return the Source magnetic group in an explicit coordinate setting."""

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
