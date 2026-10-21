"""Embed Source-ordered child operations in the parent setting."""

from __future__ import annotations

from fractions import Fraction
from typing import Any, Sequence

from verification.backend.exactmath import (
    fraction_matrix_inverse3 as _matrix_inverse,
    fraction_matrix_multiply3 as _matrix_multiply,
    fraction_record3,
    fraction_row_multiply3 as _row_multiply,
)
from verification.backend.source.magnetic import data as magnetic_data
from verification.backend import native
from verification.backend.source.magnetic_operations import _magnetic_point_operation


OperationRecord = tuple[int, int, int, int, int]
FractionVector = tuple[Fraction, Fraction, Fraction]
FractionMatrix = tuple[FractionVector, FractionVector, FractionVector]


def _point_matrix(decoder: Any, sg: int, point_op: int) -> FractionMatrix:
    units: FractionMatrix = (
        (Fraction(1), Fraction(0), Fraction(0)),
        (Fraction(0), Fraction(1), Fraction(0)),
        (Fraction(0), Fraction(0), Fraction(1)),
    )
    return tuple(
        tuple(decoder.vrot_fraction(int(sg), int(point_op), unit))
        for unit in units
    )  # type: ignore[return-value]


def _record_affine(decoder: Any, sg: int, record: OperationRecord) -> tuple[FractionMatrix, FractionVector]:
    denominator = int(record[3])
    if denominator == 0:
        raise ValueError(f"zero operation denominator: {record}")
    return (
        _point_matrix(decoder, int(sg), int(record[4])),
        tuple(Fraction(int(record[axis]), denominator) for axis in range(3)),  # type: ignore[arg-type]
    )


def _fraction_record(values: FractionVector, point_op: int) -> OperationRecord:
    return (*fraction_record3(values), int(point_op))






def embedded_child_operation_records(
    decoder: Any,
    *,
    parent_sg: int,
    child_sg: int,
    subgroup_basis: Sequence[int],
    subgroup_origin: Sequence[int],
) -> tuple[OperationRecord, ...]:
    """Map Source-ordered child operations into the raw parent PML embedding."""

    basis_values = tuple(int(value) for value in subgroup_basis)
    if len(basis_values) != 9:
        raise ValueError(f"expected 9 subgroup basis values, got {len(basis_values)}")
    origin_values = tuple(int(value) for value in subgroup_origin)
    if len(origin_values) != 4 or origin_values[3] == 0:
        raise ValueError(f"invalid subgroup origin: {subgroup_origin}")
    basis: FractionMatrix = tuple(
        tuple(Fraction(basis_values[row * 3 + col]) for col in range(3))
        for row in range(3)
    )  # type: ignore[assignment]
    inverse = _matrix_inverse(basis)
    origin: FractionVector = tuple(Fraction(origin_values[axis], origin_values[3]) for axis in range(3))  # type: ignore[assignment]
    parent_point_ops = {
        _point_matrix(decoder, int(parent_sg), int(record[4])): int(record[4])
        for record in decoder.generate_space_group_records(int(parent_sg))
    }
    out: list[OperationRecord] = []
    for record in decoder.generate_space_group_records(int(child_sg)):
        child_matrix, child_translation = _record_affine(decoder, int(child_sg), record)
        parent_matrix = _matrix_multiply(inverse, _matrix_multiply(child_matrix, basis))
        point_op = parent_point_ops.get(parent_matrix)
        if point_op is None:
            raise KeyError(
                f"child SG{child_sg} operation {record[4]} does not embed in parent SG{parent_sg}"
            )
        parent_translation = tuple(
            origin[axis]
            - _row_multiply(origin, parent_matrix)[axis]
            + _row_multiply(child_translation, basis)[axis]
            for axis in range(3)
        )
        out.append(_fraction_record(parent_translation, point_op))
    return tuple(out)


def embedded_magnetic_child_operation_records(
    decoder: Any,
    *,
    parent_sg: int,
    child_magnetic_group: int,
    subgroup_basis: Sequence[int],
    subgroup_origin: Sequence[int],
) -> tuple[OperationRecord, ...]:
    """Map a canonical BNS child group into a raw parent PML embedding."""

    basis_values = tuple(int(value) for value in subgroup_basis)
    if len(basis_values) != 9:
        raise ValueError(f"expected 9 subgroup basis values, got {len(basis_values)}")
    origin_values = tuple(int(value) for value in subgroup_origin)
    if len(origin_values) != 4 or origin_values[3] == 0:
        raise ValueError(f"invalid subgroup origin: {subgroup_origin}")

    table = magnetic_data().table
    magnetic_group = int(child_magnetic_group)
    child_sg = native.magnetic_ordinary_group(magnetic_group)
    basis: FractionMatrix = tuple(
        tuple(Fraction(basis_values[row * 3 + col]) for col in range(3))
        for row in range(3)
    )  # type: ignore[assignment]
    inverse = _matrix_inverse(basis)
    origin: FractionVector = tuple(
        Fraction(origin_values[axis], origin_values[3]) for axis in range(3)
    )  # type: ignore[assignment]
    parent_point_ops = {
        _point_matrix(decoder, int(parent_sg), int(record[4])): int(record[4])
        for record in decoder.generate_space_group_records(int(parent_sg))
    }

    out: list[OperationRecord] = []
    for record in native.magnetic_records(magnetic_group, "binary"):
        magnetic_point_op = int(record[4])
        ordinary_point_op = int(
            table["mag_point_op_mag2nonmag"][magnetic_point_op - 1]
        )
        ordinary_record = tuple(int(value) for value in record[:4]) + (
            ordinary_point_op,
        )
        child_matrix, child_translation = _record_affine(
            decoder,
            child_sg,
            ordinary_record,
        )
        parent_matrix = _matrix_multiply(
            inverse,
            _matrix_multiply(child_matrix, basis),
        )
        parent_point_op = parent_point_ops.get(parent_matrix)
        if parent_point_op is None:
            raise KeyError(
                f"magnetic child operation {magnetic_point_op} does not embed "
                f"in parent SG{parent_sg}"
            )
        parent_translation = tuple(
            origin[axis]
            - _row_multiply(origin, parent_matrix)[axis]
            + _row_multiply(child_translation, basis)[axis]
            for axis in range(3)
        )
        embedded_magnetic_point_op = _magnetic_point_operation(
            parent_point_op,
            bool(table["mag_point_op_r"][magnetic_point_op - 1]),
        )
        out.append(
            _fraction_record(parent_translation, embedded_magnetic_point_op)
        )
    return tuple(out)
