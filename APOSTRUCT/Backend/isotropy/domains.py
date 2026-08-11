"""Source-faithful ``map_subgroup2_`` / ``generate_domains2_`` ordering."""

from __future__ import annotations

from fractions import Fraction
import math

from APOSTRUCT.Backend.exactmath import (
    fraction_matrix_inverse3,
    fraction_row_multiply3,
)
from APOSTRUCT.Backend.isotropy.engine.lattice import get_new_fractionals
from APOSTRUCT.Backend.isotropy.engine.source_data import SourceData


Matrix4 = tuple[tuple[Fraction, ...], ...]
Operation = tuple[int, int, int, int, int]


def _mat4_mul(left: Matrix4, right: Matrix4) -> Matrix4:
    return tuple(tuple(sum(left[row][k] * right[k][col] for k in range(4)) for col in range(4)) for row in range(4))


def _inverse_affine(matrix: Matrix4) -> Matrix4:
    spatial = tuple(tuple(row[:3]) for row in matrix[:3])
    try:
        inverse = fraction_matrix_inverse3(spatial)
    except ValueError:
        raise ValueError("singular subgroup basis")
    out = [[Fraction(0) for _ in range(4)] for __ in range(4)]
    for row in range(3):
        for col in range(3):
            out[row][col] = inverse[row][col]
    translation = matrix[3][:3]
    transformed_translation = fraction_row_multiply3(translation, inverse)
    for col, value in enumerate(transformed_translation):
        out[3][col] = -value
    out[3][3] = Fraction(1)
    return tuple(tuple(row) for row in out)


def _setting_matrix(basis: tuple[int, ...], origin: tuple[int, int, int, int]) -> Matrix4:
    out = [[Fraction(0) for _ in range(4)] for __ in range(4)]
    for row in range(3):
        for col in range(3):
            out[row][col] = Fraction(int(basis[row * 3 + col]))
    for col in range(3):
        out[3][col] = Fraction(int(origin[col]), int(origin[3] or 1))
    out[3][3] = Fraction(1)
    return tuple(tuple(row) for row in out)


def _pml_point_matrix(data: SourceData, sg: int, point_op: int) -> tuple[int, ...]:
    lattice = int(data.space["ispace_lattice"][int(sg) - 1])
    code = int(data.space["ipoint_op_psettings"][(lattice - 1) * 72 + int(point_op) - 1])
    if code == 0:
        raise KeyError(f"missing pml point operation SG{sg}/op{point_op}")
    values: list[int] = []
    for _row in range(3):
        for _col in range(3):
            values.append(code % 3 - 1)
            code //= 3
    return tuple(values)


def _operation_matrix(data: SourceData, sg: int, record: Operation) -> Matrix4:
    out = [[Fraction(0) for _ in range(4)] for __ in range(4)]
    point = _pml_point_matrix(data, sg, int(record[4]))
    for row in range(3):
        for col in range(3):
            out[row][col] = Fraction(point[row * 3 + col])
    for col in range(3):
        out[3][col] = Fraction(int(record[col]), int(record[3]))
    out[3][3] = Fraction(1)
    return tuple(tuple(row) for row in out)


def _fraction_mod1(values: tuple[Fraction, Fraction, Fraction]) -> tuple[int, int, int, int]:
    denominator = 1
    for value in values:
        denominator = math.lcm(denominator, value.denominator)
    numerators = [int(value * denominator) for value in values]
    divisor = denominator
    for numerator in numerators:
        divisor = math.gcd(divisor, abs(numerator))
    if divisor:
        denominator //= divisor
        numerators = [value // divisor for value in numerators]
    numerators = [value % denominator for value in numerators]
    return numerators[0], numerators[1], numerators[2], denominator


def _matrix_to_parent_record(data: SourceData, sg: int, matrix: Matrix4) -> Operation:
    point = tuple(int(matrix[row][col]) for row in range(3) for col in range(3))
    if any(matrix[row][col].denominator != 1 for row in range(3) for col in range(3)):
        raise KeyError(f"non-integral mapped point operation SG{sg}")
    parent_records = data.generate_space_group_records(int(sg))
    point_by_matrix = {_pml_point_matrix(data, sg, int(record[4])): int(record[4]) for record in parent_records}
    point_op = point_by_matrix.get(point)
    if point_op is None:
        raise KeyError(f"mapped point operation not present in SG{sg}: {point}")
    fraction = _fraction_mod1((matrix[3][0], matrix[3][1], matrix[3][2]))
    for record in parent_records:
        if int(record[4]) != point_op:
            continue
        if all((Fraction(fraction[axis], fraction[3]) - Fraction(record[axis], record[3])) % 1 == 0 for axis in range(3)):
            return record
    raise KeyError(f"mapped operation not represented in SG{sg}: {fraction}/op{point_op}")


def mapped_subgroup_records(
    data: SourceData,
    *,
    parent_sg: int,
    child_sg: int,
    basis: tuple[int, ...],
    origin: tuple[int, int, int, int],
) -> tuple[Operation, ...]:
    setting = _setting_matrix(basis, origin)
    inverse = _inverse_affine(setting)
    return tuple(
        _matrix_to_parent_record(
            data,
            parent_sg,
            _mat4_mul(_mat4_mul(inverse, _operation_matrix(data, child_sg, child)), setting),
        )
        for child in data.generate_space_group_records(int(child_sg))
    )


def domain_operation_records(
    data: SourceData,
    *,
    parent_sg: int,
    child_sg: int,
    basis: tuple[int, ...],
    origin: tuple[int, int, int, int],
) -> tuple[Operation, ...]:
    """Return domain operations in ``generate_domains2_`` scan order."""

    mapped = mapped_subgroup_records(
        data,
        parent_sg=int(parent_sg),
        child_sg=int(child_sg),
        basis=basis,
        origin=origin,
    )
    return domain_operation_records_from_stabilizer(
        data,
        parent_sg=int(parent_sg),
        basis=basis,
        subgroup_operations=mapped,
    )


def domain_operation_records_from_stabilizer(
    data: SourceData,
    *,
    parent_sg: int,
    basis: tuple[int, ...],
    subgroup_operations: tuple[Operation, ...],
) -> tuple[Operation, ...]:
    """Order domains from the OPD stabilizer before child presentation."""

    parent_records = tuple(data.generate_space_group_records(int(parent_sg)))
    subgroup_points = {int(record[4]) for record in subgroup_operations}
    identity = (1, 0, 0, 0, 1, 0, 0, 0, 1)
    identity_index = next(index for index, record in enumerate(parent_records, start=1) if _pml_point_matrix(data, parent_sg, int(record[4])) == identity)
    marked = [False] * (len(parent_records) + 1)
    cosets = [identity_index]
    marked[identity_index] = True
    for index, record in enumerate(parent_records, start=1):
        if int(record[4]) in subgroup_points:
            marked[index] = True
    point_to_index = {int(record[4]): index for index, record in enumerate(parent_records, start=1)}
    multiplication = data.space["ipoint_op_mlt"]

    def point_product(left: int, right: int) -> int:
        # generate_domains2_ indexes the Fortran table as (right, left).
        return int(multiplication[(right - 1) * 72 + left - 1])

    target = len(parent_records) // len(subgroup_points)
    while len(cosets) < target:
        chosen = 0
        for index in range(len(parent_records), 0, -1):
            if not marked[index]:
                chosen = index
        if chosen == 0:
            raise RuntimeError("generate_domains2 point-coset scan did not close")
        cosets.append(chosen)
        chosen_point = int(parent_records[chosen - 1][4])
        for subgroup_point in subgroup_points:
            marked[point_to_index[point_product(chosen_point, subgroup_point)]] = True

    lattice = int(data.space["ispace_lattice"][int(parent_sg) - 1])
    out: list[Operation] = []
    for fraction in get_new_fractionals(basis):
        for coset in cosets:
            parent_record = parent_records[coset - 1]
            rotated = data._rotate_kernel_fraction_by_space_operation(  # noqa: SLF001
                fraction,
                parent_record,
                lattice=lattice,
            )
            out.append(data._vadd_fraction_operation(rotated, parent_record))  # noqa: SLF001
    return tuple(out)
