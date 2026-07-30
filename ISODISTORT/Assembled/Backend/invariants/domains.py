"""Invariant-facing domain operations backed by the Assembled domain kernel."""

from __future__ import annotations

from fractions import Fraction
import math

from ISODISTORT.Assembled.Backend.isotropy.domains import (
    _mat4_mul,
    _operation_matrix,
    _pml_point_matrix,
    mapped_subgroup_records,
)
from ISODISTORT.Assembled.Backend.isotropy.engine.lattice import get_new_fractionals


def _row_data(source, row_id: int):
    index = int(row_id) - 1
    table = source.isotropy
    return (
        int(table["isotropy_subgroup"][index]),
        tuple(int(value) for value in table["isotropy_basis"][index * 9:(index + 1) * 9]),
        tuple(int(value) for value in table["isotropy_origin"][index * 4:(index + 1) * 4]),
    )


def _det3(values) -> int:
    entries = tuple(int(value) for value in values)
    if len(entries) != 9:
        raise ValueError("domain basis must contain nine integers")
    return (
        entries[0] * (entries[4] * entries[8] - entries[5] * entries[7])
        - entries[1] * (entries[3] * entries[8] - entries[5] * entries[6])
        + entries[2] * (entries[3] * entries[7] - entries[4] * entries[6])
    )


def domain_count_from_isotropy_row(source, source_data, sg: int, row_id: int) -> int:
    child_sg, basis, _origin = _row_data(source, row_id)
    return domain_count_from_subgroup(source_data, sg=sg, child_sg=child_sg, basis=basis)


def domain_count_from_subgroup(source_data, *, sg: int, child_sg: int, basis, origin=(0, 0, 0, 1)) -> int:
    del origin
    parent_point_group = int(source_data.space["ispace_point_group"][int(sg) - 1])
    child_point_group = int(source_data.space["ispace_point_group"][int(child_sg) - 1])
    parent_order = int(source_data.space["ipoint_group_order"][parent_point_group - 1])
    child_order = int(source_data.space["ipoint_group_order"][child_point_group - 1])
    return 0 if child_order == 0 else parent_order * abs(_det3(basis)) // child_order


def domain_operation_record_from_isotropy_row(source, source_data, *, sg: int, row_id: int, domain: int):
    child_sg, basis, origin = _row_data(source, row_id)
    return domain_operation_record_from_subgroup(
        source_data,
        sg=sg,
        child_sg=child_sg,
        basis=basis,
        origin=origin,
        domain=domain,
    )


def domain_operation_record_from_subgroup(source_data, *, sg: int, child_sg: int, basis, origin, domain: int):
    if int(domain) <= 1:
        return None
    parent_records = tuple(source_data.generate_space_group_records(int(sg)))
    mapped = mapped_subgroup_records(
        source_data,
        parent_sg=int(sg),
        child_sg=int(child_sg),
        basis=tuple(int(value) for value in basis),
        origin=tuple(int(value) for value in origin),
    )
    subgroup_points = {int(record[4]) for record in mapped}
    if not subgroup_points:
        return None
    identity = (1, 0, 0, 0, 1, 0, 0, 0, 1)
    identity_index = next(
        (
            index
            for index, record in enumerate(parent_records, start=1)
            if _pml_point_matrix(source_data, int(sg), int(record[4])) == identity
        ),
        None,
    )
    if identity_index is None:
        return None
    marked = [False] * (len(parent_records) + 1)
    cosets = [identity_index]
    marked[identity_index] = True
    for index, record in enumerate(parent_records, start=1):
        if int(record[4]) in subgroup_points:
            marked[index] = True
    point_to_index = {
        int(record[4]): index for index, record in enumerate(parent_records, start=1)
    }
    multiplication = source_data.space["ipoint_op_mlt"]

    def point_product(left: int, right: int) -> int:
        return int(multiplication[(left - 1) * 72 + right - 1])

    target = len(parent_records) // len(subgroup_points)
    while len(cosets) < target:
        chosen = 0
        for index in range(len(parent_records), 0, -1):
            if not marked[index]:
                chosen = index
        if chosen == 0:
            return None
        cosets.append(chosen)
        chosen_point = int(parent_records[chosen - 1][4])
        for subgroup_point in subgroup_points:
            marked[point_to_index[point_product(chosen_point, subgroup_point)]] = True

    fractions = get_new_fractionals(tuple(int(value) for value in basis))
    maximum = len(cosets) * len(fractions)
    if not 1 <= int(domain) <= maximum:
        return None
    fraction_index, coset_index = divmod(int(domain) - 1, len(cosets))
    parent_record = parent_records[cosets[coset_index] - 1]
    rotated = source_data._rotate_kernel_fraction_by_space_operation(  # noqa: SLF001
        fractions[fraction_index],
        parent_record,
        lattice=int(source_data.space["ispace_lattice"][int(sg) - 1]),
    )
    return source_data._vadd_fraction_operation(rotated, parent_record)  # noqa: SLF001


def _domain_transform_matrix(basis, origin):
    matrix = [[Fraction(0) for _ in range(4)] for __ in range(4)]
    for row in range(3):
        for column in range(3):
            matrix[row][column] = Fraction(int(basis[row * 3 + column]))
    denominator = int(origin[3])
    if denominator == 0:
        raise ValueError("zero subgroup origin denominator")
    for column in range(3):
        matrix[3][column] = Fraction(int(origin[column]), denominator)
    matrix[3][3] = Fraction(1)
    return tuple(tuple(row) for row in matrix)


def _pml_operation_matrix(source_data, sg: int, record):
    return _operation_matrix(source_data, int(sg), tuple(int(value) for value in record))


def _fraction_record_mod1(values):
    denominator = 1
    for value in values:
        denominator = math.lcm(denominator, Fraction(value).denominator)
    numerators = [int(Fraction(value) * denominator) for value in values]
    divisor = denominator
    for numerator in numerators:
        divisor = math.gcd(divisor, abs(numerator))
    if divisor:
        denominator //= divisor
        numerators = [value // divisor for value in numerators]
    return tuple(value % denominator for value in numerators) + (denominator,)


__all__ = [
    "_domain_transform_matrix",
    "_fraction_record_mod1",
    "_mat4_mul",
    "_pml_operation_matrix",
    "domain_count_from_isotropy_row",
    "domain_count_from_subgroup",
    "domain_operation_record_from_isotropy_row",
    "domain_operation_record_from_subgroup",
]
