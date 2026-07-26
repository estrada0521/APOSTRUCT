"""Source-only Wyckoff splitting helper for undistorted superstructures.

This implements the ISO ``get_new_wyckoff_`` path used by
``SHOW WYCKOFF SUBGROUP``.  It intentionally depends only on Source-backed
``SourceTables`` tables, not on web captures or diagnostic fixtures.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
import itertools
import math
from numbers import Integral
import re
from functools import lru_cache
from typing import Any

from ISODISTORT.Assembled.Backend.exactmath import (
    fraction_matrix_inverse3,
    fraction_row_multiply3 as _row_multiply,
    integer_determinant3,
)
from ISODISTORT.Assembled.Backend.fraction_expression import (
    evaluate_fraction_expression,
    split_coordinate_expression3,
)
from ISODISTORT.Assembled.Backend.lattice_quotient import (
    integer_inverse_denominator,
    integral_row_images_source_order,
)
from ISODISTORT.Assembled.Backend.source.tables import SourceTables, WyckoffRow


FractionRecord = tuple[int, int, int, int]
WyckoffBranch = tuple[FractionRecord, FractionRecord, FractionRecord, FractionRecord]
Formula15 = tuple[int, int, int, int, int, int, int, int, int, int, int, int, int, int, int]


@dataclass(frozen=True)
class WyckoffSplitSourceOccurrence:
    flat_index: int
    kernel_fraction_ordinal: int
    parent_branch_ordinal: int
    kernel_fraction: FractionRecord
    parent_coset_record: tuple[int, int, int, int, int]
    parent_affine_branch: WyckoffBranch


@dataclass(frozen=True)
class WyckoffSplitRowProvenance:
    source_group_ordinal: int
    child_row_offset: int
    formula15: Formula15
    representative_flat_index: int
    orbit_flat_indices: tuple[int, ...]
    occurrences: tuple[WyckoffSplitSourceOccurrence, ...]


def _gcd_many(values: list[int]) -> int:
    nonzero = [abs(value) for value in values if value]
    if not nonzero:
        return 1
    result = nonzero[0]
    for value in nonzero[1:]:
        result = math.gcd(result, value)
    return result


def _reduce_fraction_record(record: FractionRecord) -> FractionRecord:
    x, y, z, den = (int(value) for value in record)
    if den == 0:
        raise ValueError("zero denominator")
    if den < 0:
        x, y, z, den = -x, -y, -z, -den
    factor = _gcd_many([x, y, z, den])
    return x // factor, y // factor, z // factor, den // factor


def _matinv_denominator(matrix: tuple[int, ...]) -> int:
    values = tuple(int(value) for value in matrix[:9])
    if len(values) < 9:
        raise IndexError("tuple index out of range")
    det = integer_determinant3(values)
    if det == 0:
        raise ValueError("singular matrix")
    return integer_inverse_denominator(values)


@lru_cache(maxsize=256)
def _get_new_fractionals(matrix: tuple[int, ...]) -> tuple[FractionRecord, ...]:
    """Return Source-ordered quotient translations for an integer basis.

    The shared exact quotient kernel generates the same ``(x, y, z)``
    solutions and restores the Source lexicographic scan order without
    enumerating the full modular cube.
    """

    return integral_row_images_source_order(matrix, _matinv_denominator(matrix))


def _fraction_record(values: tuple[Fraction, Fraction, Fraction]) -> FractionRecord:
    den = 1
    for value in values:
        den = math.lcm(den, value.denominator)
    nums = [int(value * den) for value in values]
    factor = _gcd_many(nums + [den])
    if factor > 1:
        nums = [num // factor for num in nums]
        den //= factor
    if den < 0:
        nums = [-num for num in nums]
        den = -den
    return nums[0], nums[1], nums[2], den


def _fraction_values(record: FractionRecord) -> tuple[Fraction, Fraction, Fraction]:
    return tuple(Fraction(int(record[axis]), int(record[3])) for axis in range(3))  # type: ignore[return-value]


def _operation_fraction_part(record: tuple[int, int, int, int, int]) -> FractionRecord:
    return int(record[0]), int(record[1]), int(record[2]), int(record[3])


def _fraction_add(left: FractionRecord, right: FractionRecord) -> FractionRecord:
    left_values = _fraction_values(left)
    right_values = _fraction_values(right)
    return _fraction_record(tuple(left_values[axis] + right_values[axis] for axis in range(3)))  # type: ignore[arg-type]


def _mod_record(record: FractionRecord) -> FractionRecord:
    values = tuple(value % 1 for value in _fraction_values(record))
    return _fraction_record(values)  # type: ignore[arg-type]


def _is_zero_record(record: FractionRecord) -> bool:
    return int(record[0]) == 0 and int(record[1]) == 0 and int(record[2]) == 0


def _matrix_inverse_3_fraction(
    basis: tuple[int, ...],
) -> tuple[tuple[Fraction, Fraction, Fraction], ...]:
    matrix = tuple(
        tuple(Fraction(int(basis[3 * row + column])) for column in range(3))
        for row in range(3)
    )
    return fraction_matrix_inverse3(matrix)


def _basis_tuple(basis: list[list[float]] | list[list[int]] | tuple[tuple[int, int, int], ...]) -> tuple[int, ...]:
    flat: list[int] = []
    for row in basis:
        if len(row) != 3:
            raise ValueError(f"expected 3-column basis row, got {row!r}")
        for value in row:
            fraction = Fraction(str(value))
            if fraction.denominator != 1:
                raise ValueError(f"non-integral subgroup basis entry: {value!r}")
            flat.append(int(fraction))
    if len(flat) != 9:
        raise ValueError(f"expected 3x3 subgroup basis, got {basis!r}")
    return tuple(flat)


def _origin_record(origin: list[float] | tuple[float, ...] | None) -> FractionRecord:
    if origin is None:
        return (0, 0, 0, 1)
    if len(origin) == 4:
        return _reduce_fraction_record(tuple(int(value) for value in origin))  # type: ignore[arg-type]
    values = tuple(Fraction(str(value)) for value in origin)
    if len(values) != 3:
        raise ValueError(f"expected 3-vector origin, got {origin!r}")
    return _fraction_record(values)  # type: ignore[arg-type]


def get_wyckoff_all(
    data: SourceTables,
    sg: int,
    wyckoff_index: int,
) -> tuple[WyckoffBranch, ...]:
    rows = data.wyckoff_rows(int(sg))
    index = int(wyckoff_index)
    if index < 1 or index > len(rows):
        raise IndexError(f"SG{sg} Wyckoff index out of range: {wyckoff_index}")
    row = rows[index - 1]
    raw_vectors = tuple(_fraction_record(vector) for vector in data.wyckoff_fraction_vectors(row))
    if len(raw_vectors) != 4:
        raise ValueError(f"unexpected Wyckoff vector count for SG{sg}/{row.label}: {len(raw_vectors)}")

    site_records = data.wyc_pg_elements_records(int(sg), row)
    cosets = data.wyc_pg_cosets_records(int(sg), site_records)
    out: list[WyckoffBranch] = [raw_vectors]  # type: ignore[list-item]
    for coset in cosets[1:]:
        rotated: list[FractionRecord] = []
        for vector in data.wyckoff_fraction_vectors(row):
            rotated.append(_fraction_record(data.vrot_fraction(int(sg), int(coset[4]), vector)))
        rotated[0] = _fraction_add(_operation_fraction_part(coset), rotated[0])
        out.append(tuple(rotated))  # type: ignore[arg-type]
    return tuple(out)


def transform_wyckoff_branches_to_child(
    data: SourceTables,
    parent_sg: int,
    wyckoff_index: int,
    subgroup_basis: tuple[int, ...],
    subgroup_origin: FractionRecord,
) -> tuple[tuple[WyckoffBranch, ...], ...]:
    basis = tuple(int(value) for value in subgroup_basis[:9])
    inverse = _matrix_inverse_3_fraction(basis)
    origin_child = _row_multiply(_fraction_values(subgroup_origin), inverse)
    branches = get_wyckoff_all(data, int(parent_sg), int(wyckoff_index))
    out: list[tuple[WyckoffBranch, ...]] = []
    for kernel_fraction in _get_new_fractionals(basis):
        kernel_values = _fraction_values(kernel_fraction)
        transformed_for_fraction: list[WyckoffBranch] = []
        for branch in branches:
            branch_vectors: list[FractionRecord] = []
            base_values = tuple(_fraction_values(branch[0])[axis] + kernel_values[axis] for axis in range(3))
            child_base = _row_multiply(base_values, inverse)
            child_base = tuple(child_base[axis] - origin_child[axis] for axis in range(3))
            branch_vectors.append(_mod_record(_fraction_record(child_base)))
            for parameter in branch[1:]:
                child_param = _row_multiply(_fraction_values(parameter), inverse)
                branch_vectors.append(_fraction_record(child_param))
            transformed_for_fraction.append(tuple(branch_vectors))  # type: ignore[arg-type]
        out.append(tuple(transformed_for_fraction))
    return tuple(out)


def _apply_child_operation(
    data: SourceTables,
    child_sg: int,
    operation: tuple[int, int, int, int, int],
    branch: WyckoffBranch,
) -> WyckoffBranch:
    out: list[FractionRecord] = []
    for index, vector in enumerate(branch):
        rotated = _fraction_record(data.vrot_fraction(int(child_sg), int(operation[4]), _fraction_values(vector)))
        if index == 0:
            rotated = _mod_record(_fraction_add(_operation_fraction_part(operation), rotated))
        out.append(rotated)
    return tuple(out)  # type: ignore[return-value]


def group_child_orbits(
    data: SourceTables,
    child_sg: int,
    child_branches_by_fraction: tuple[tuple[WyckoffBranch, ...], ...],
) -> tuple[dict[str, object], ...]:
    flat: list[WyckoffBranch] = [
        branch for fraction in child_branches_by_fraction for branch in fraction
    ]
    index_by_branch = {branch: index for index, branch in enumerate(flat)}
    child_ops = data.generate_space_group_records(int(child_sg))
    point_group_codes = [int(value) for value in data.space["ipoint_group_code"]]
    point_op_codes = [int(value) for value in data.space["ipoint_op_code"]]
    visited_group = [0] * len(flat)
    groups: list[dict[str, object]] = []
    for index, branch in enumerate(flat):
        if visited_group[index]:
            continue
        group_id = len(groups) + 1
        visited_group[index] = group_id
        active_parameters = tuple(i + 1 for i, vector in enumerate(branch[1:]) if not _is_zero_record(vector))
        stabilizer_point_ops: list[int] = []
        orbit_indices = {index}
        for operation in child_ops:
            transformed = _apply_child_operation(data, int(child_sg), operation, branch)
            target = index_by_branch.get(transformed)
            if target is None:
                raise KeyError(f"child operation image not found for SG{child_sg}: {operation} {branch} -> {transformed}")
            orbit_indices.add(target)
            if visited_group[target] not in (0, group_id):
                raise ValueError(f"orbit collision for SG{child_sg}: target={target} group={visited_group[target]} new={group_id}")
            visited_group[target] = group_id
            if target == index:
                stabilizer_point_ops.append(int(operation[4]))
        code = sum(point_op_codes[op - 1] for op in stabilizer_point_ops)
        try:
            site_pg = point_group_codes.index(code) + 1
        except ValueError:
            site_pg = None
        groups.append(
            {
                "group_id": group_id,
                "representative_index": index,
                "orbit_indices": tuple(sorted(orbit_indices)),
                "active_parameters": active_parameters,
                "stabilizer_point_ops": tuple(stabilizer_point_ops),
                "point_group_code": code,
                "site_pg": site_pg,
            }
        )
    return tuple(groups)


def _rref_fraction(
    values: list[list[Fraction]],
    coefficient_columns: int,
) -> tuple[list[list[Fraction]], tuple[int, ...]]:
    matrix = [list(row) for row in values]
    pivot_columns: list[int] = []
    pivot_row = 0
    for column in range(coefficient_columns):
        pivot = next(
            (row for row in range(pivot_row, 3) if matrix[row][column] != 0),
            None,
        )
        if pivot is None:
            continue
        matrix[pivot_row], matrix[pivot] = matrix[pivot], matrix[pivot_row]
        scale = matrix[pivot_row][column]
        matrix[pivot_row] = [value / scale for value in matrix[pivot_row]]
        for row in range(3):
            if row == pivot_row or matrix[row][column] == 0:
                continue
            factor = matrix[row][column]
            matrix[row] = [
                matrix[row][index] - factor * matrix[pivot_row][index]
                for index in range(len(matrix[row]))
            ]
        pivot_columns.append(column)
        pivot_row += 1
        if pivot_row == 3:
            break
    return matrix, tuple(pivot_columns)


@lru_cache(maxsize=None)
def _rank(columns: tuple[tuple[Fraction, Fraction, Fraction], ...]) -> int:
    if not columns:
        return 0
    _matrix, pivots = _rref_fraction(
        [[Fraction(column[row]) for column in columns] for row in range(3)],
        len(columns),
    )
    return len(pivots)


def _in_span(
    vector: tuple[Fraction, Fraction, Fraction],
    columns: list[tuple[Fraction, Fraction, Fraction]],
) -> bool:
    key = tuple(columns)
    return _rank(key) == _rank(key + (vector,))


def _affine_subset_of_child(
    base: tuple[Fraction, Fraction, Fraction],
    params: list[tuple[Fraction, Fraction, Fraction]],
    child_base: tuple[Fraction, Fraction, Fraction],
    child_params: list[tuple[Fraction, Fraction, Fraction]],
) -> bool:
    if any(not _in_span(param, child_params) for param in params):
        return False
    delta = tuple(base[axis] - child_base[axis] for axis in range(3))
    for offset in itertools.product(range(-2, 3), repeat=3):
        shifted = tuple(delta[axis] - Fraction(offset[axis]) for axis in range(3))
        if _in_span(shifted, child_params):
            return True
    return False


def _factor_column(values: list[int]) -> list[int]:
    divisor = 0
    for value in values:
        divisor = math.gcd(divisor, abs(value))
    if divisor <= 1:
        return values
    return [value // divisor for value in values]


def _rowop2_reduce(rows: list[list[int]]) -> list[list[int]]:
    if not rows:
        return rows
    nrows = len(rows)
    ncols = 3
    skipped = 0
    col = 0
    while col < ncols:
        pivot_row = skipped + col
        if pivot_row >= nrows:
            break
        while rows[pivot_row][col] == 0:
            swap_col = None
            for candidate_col in range(col + 1, ncols):
                if rows[pivot_row][candidate_col] != 0:
                    swap_col = candidate_col
                    break
            if swap_col is not None:
                for row in rows:
                    row[col], row[swap_col] = row[swap_col], row[col]
                break
            skipped += 1
            pivot_row = skipped + col
            if pivot_row >= nrows:
                break
        if pivot_row >= nrows:
            break
        pivot = rows[pivot_row][col]
        for target_col in range(ncols):
            if target_col == col or rows[pivot_row][target_col] == 0:
                continue
            target = rows[pivot_row][target_col]
            for row_index, row in enumerate(rows):
                if row_index != pivot_row:
                    row[target_col] = pivot * row[target_col] - row[col] * target
            rows[pivot_row][target_col] = 0
            column = _factor_column([row[target_col] for row in rows])
            for row, value in zip(rows, column, strict=True):
                row[target_col] = value
        col += 1
    for target_col in range(ncols):
        column = _factor_column([row[target_col] for row in rows])
        for row, value in zip(rows, column, strict=True):
            row[target_col] = value
    return rows


def _rowop_accepts_itry(
    child_branch: WyckoffBranch,
    group_branch: WyckoffBranch,
    parent_active_parameters: tuple[int, ...],
    itry: tuple[int, int, int],
) -> bool:
    child_active = tuple(index + 1 for index, vector in enumerate(child_branch[1:]) if not _is_zero_record(vector))
    child_vectors = tuple(tuple(int(value) for value in vector) for vector in child_branch)
    group_vectors = tuple(tuple(int(value) for value in vector) for vector in group_branch)
    denominator_product = 1
    for vector in child_vectors + group_vectors:
        denominator_product *= int(vector[3])

    child_base = list(child_vectors[0])
    for axis in range(3):
        child_base[axis] -= child_base[3] * itry[axis]
    rows: list[list[int]] = []
    for axis_index in child_active:
        vector = child_vectors[axis_index]
        rows.append([(denominator_product * int(vector[axis])) // int(vector[3]) for axis in range(3)])
    rows.append([
        (denominator_product * int(group_vectors[0][axis])) // int(group_vectors[0][3])
        - (denominator_product * int(child_base[axis])) // int(child_base[3])
        for axis in range(3)
    ])
    for axis_index in parent_active_parameters:
        vector = group_vectors[axis_index]
        rows.append([(denominator_product * int(vector[axis])) // int(vector[3]) for axis in range(3)])
    reduced = _rowop2_reduce([row[:] for row in rows])
    return all(row[col] == 0 for row in reduced for col in range(len(child_active), 3))


def child_row_matches_for_group(
    data: SourceTables,
    child_sg: int,
    child_branches_by_fraction: tuple[tuple[WyckoffBranch, ...], ...],
    group: dict[str, object],
) -> tuple[tuple[int, WyckoffBranch, WyckoffBranch], ...]:
    flat: list[WyckoffBranch] = [
        branch for fraction in child_branches_by_fraction for branch in fraction
    ]
    out: list[tuple[int, WyckoffBranch, WyckoffBranch]] = []
    group_indices = tuple(int(index) for index in group["orbit_indices"])  # type: ignore[index]
    parent_active_parameters = tuple(int(value) for value in group["active_parameters"])  # type: ignore[index]
    rowop_out: list[tuple[int, WyckoffBranch, WyckoffBranch]] = []
    for row in data.wyckoff_rows(int(child_sg)):
        if int(row.site_pg) != int(group["site_pg"]):
            continue
        child_branch = tuple(_fraction_record(vector) for vector in data.wyckoff_fraction_vectors(row))
        child_base = _fraction_values(child_branch[0])
        child_params = [_fraction_values(vector) for vector in child_branch[1:] if not _is_zero_record(vector)]
        # Affine membership is independent of the Source ``itry`` scan.
        affine_group_branches: list[WyckoffBranch] = []
        for index in group_indices:
            group_branch = flat[index]
            base = _fraction_values(group_branch[0])
            params = [_fraction_values(vector) for vector in group_branch[1:] if not _is_zero_record(vector)]
            if _affine_subset_of_child(base, params, child_base, child_params):
                affine_group_branches.append(group_branch)
        for itry in itertools.product((0, 1, -1), repeat=3):
            for group_branch in affine_group_branches:
                if _rowop_accepts_itry(child_branch, group_branch, parent_active_parameters, itry):
                    rowop_out.append((int(row.offset0) + 1, group_branch, child_branch))
    if rowop_out:
        return tuple(rowop_out)

    for row in data.wyckoff_rows(int(child_sg)):
        if int(row.site_pg) != int(group["site_pg"]):
            continue
        child_branch = tuple(_fraction_record(vector) for vector in data.wyckoff_fraction_vectors(row))
        child_base = _fraction_values(child_branch[0])
        child_params = [_fraction_values(vector) for vector in child_branch[1:] if not _is_zero_record(vector)]
        for index in group_indices:
            group_branch = flat[index]
            base = _fraction_values(group_branch[0])
            params = [_fraction_values(vector) for vector in group_branch[1:] if not _is_zero_record(vector)]
            if _affine_subset_of_child(base, params, child_base, child_params):
                out.append((int(row.offset0) + 1, group_branch, child_branch))
    return tuple(out)


def _solve_child_parameter_column(
    child_vectors: tuple[tuple[Fraction, Fraction, Fraction], ...],
    target: tuple[Fraction, Fraction, Fraction],
    *,
    allow_integer_shift: bool,
) -> tuple[Fraction, Fraction, Fraction] | None:
    matrix = tuple(
        tuple(Fraction(child_vectors[col + 1][row]) for col in range(3))
        for row in range(3)
    )
    if allow_integer_shift:
        primary_offsets = list(itertools.product((0, 1, -1), repeat=3))
        fallback_offsets = [
            offset
            for offset in itertools.product(range(-2, 3), repeat=3)
            if offset not in set(primary_offsets)
        ]
        offsets = primary_offsets + fallback_offsets
    else:
        offsets = [(0, 0, 0)]
    for offset in offsets:
        shifted = tuple(Fraction(target[row]) - int(offset[row]) for row in range(3))
        augmented, pivot_columns = _rref_fraction(
            [list(matrix[row]) + [shifted[row]] for row in range(3)],
            3,
        )
        if any(all(row[column] == 0 for column in range(3)) and row[3] != 0 for row in augmented):
            continue
        values = [Fraction(0), Fraction(0), Fraction(0)]
        for row, column in enumerate(pivot_columns):
            values[column] = augmented[row][3]
        vector = values[0], values[1], values[2]
        if all(
            sum(matrix[row][column] * vector[column] for column in range(3)) == shifted[row]
            for row in range(3)
        ):
            return vector
    return None


def _formula_block(
    constant: Fraction,
    coefficients: tuple[Fraction, Fraction, Fraction],
    *,
    force_emit: bool,
) -> tuple[int, int, int, int, int]:
    values = (constant % 1, coefficients[0], coefficients[1], coefficients[2])
    if not force_emit and all(value == 0 for value in values):
        return 0, 0, 0, 0, 0
    denominator = 1
    for value in values:
        denominator = math.lcm(denominator, value.denominator)
    return tuple(int(value * denominator) for value in values) + (denominator,)  # type: ignore[return-value]


def formula15_for_child_row(
    data: SourceTables,
    child_sg: int,
    child_row_offset: int,
    branch: WyckoffBranch,
    child_branch: WyckoffBranch | None = None,
) -> Formula15:
    rows = data.wyckoff_rows(int(child_sg))
    if int(child_row_offset) < 1 or int(child_row_offset) > len(rows):
        raise IndexError(f"child Wyckoff offset out of range for SG{child_sg}: {child_row_offset}")
    child_vectors = (
        tuple(_fraction_values(record) for record in child_branch)
        if child_branch is not None
        else data.wyckoff_fraction_vectors(rows[int(child_row_offset) - 1])
    )
    branch_vectors = tuple(_fraction_values(record) for record in branch)

    solved_columns: list[tuple[Fraction, Fraction, Fraction]] = []
    for column in range(4):
        target = tuple(
            branch_vectors[column][axis] - (child_vectors[0][axis] if column == 0 else Fraction(0))
            for axis in range(3)
        )
        solved = _solve_child_parameter_column(
            child_vectors,
            target,
            allow_integer_shift=column == 0,
        )
        if solved is None:
            raise ValueError(f"could not solve child formula for SG{child_sg} row offset {child_row_offset}")
        solved_columns.append(solved)

    active_axes = tuple(any(child_vectors[axis + 1]) for axis in range(3))
    blocks: list[int] = []
    for axis in range(3):
        blocks.extend(
            _formula_block(
                solved_columns[0][axis],
                (solved_columns[1][axis], solved_columns[2][axis], solved_columns[3][axis]),
                force_emit=active_axes[axis],
            )
        )
    return tuple(blocks)  # type: ignore[return-value]


@lru_cache(maxsize=256)
def _get_new_wyckoff_child_rows_formula15_provenance_cached(
    data: SourceTables,
    parent_sg: int,
    child_sg: int,
    wyckoff_index: int,
    subgroup_basis: tuple[int, ...],
    subgroup_origin: FractionRecord,
) -> tuple[WyckoffSplitRowProvenance, ...]:
    child_branches_by_fraction = transform_wyckoff_branches_to_child(
        data,
        int(parent_sg),
        int(wyckoff_index),
        subgroup_basis,
        subgroup_origin,
    )
    groups = group_child_orbits(
        data, int(child_sg), child_branches_by_fraction
    )
    parent_branches = get_wyckoff_all(data, int(parent_sg), int(wyckoff_index))
    parent_row = data.wyckoff_rows(int(parent_sg))[int(wyckoff_index) - 1]
    parent_cosets = data.wyc_pg_cosets_records(
        int(parent_sg),
        data.wyc_pg_elements_records(int(parent_sg), parent_row),
    )
    kernel_fractions = _get_new_fractionals(tuple(int(value) for value in subgroup_basis[:9]))
    if len(parent_branches) != len(parent_cosets):
        raise ValueError("parent Wyckoff branch/coset count mismatch")
    branch_count = len(parent_branches)
    out: list[WyckoffSplitRowProvenance] = []
    for source_group_ordinal, group in enumerate(groups):
        matches = child_row_matches_for_group(
            data, int(child_sg), child_branches_by_fraction, group
        )
        if not matches:
            raise KeyError(f"no child Wyckoff candidate for SG{child_sg} group {group}")
        child_row_offset, branch, child_branch = matches[0]
        occurrences: list[WyckoffSplitSourceOccurrence] = []
        orbit_flat_indices = tuple(int(index) for index in group["orbit_indices"])  # type: ignore[index]
        for flat_index in orbit_flat_indices:
            kernel_ordinal, parent_branch_ordinal = divmod(int(flat_index), branch_count)
            kernel_fraction = kernel_fractions[kernel_ordinal]
            parent_branch = parent_branches[parent_branch_ordinal]
            parent_affine_branch = (
                _fraction_add(parent_branch[0], kernel_fraction),
                parent_branch[1],
                parent_branch[2],
                parent_branch[3],
            )
            occurrences.append(
                WyckoffSplitSourceOccurrence(
                    flat_index=int(flat_index),
                    kernel_fraction_ordinal=int(kernel_ordinal),
                    parent_branch_ordinal=int(parent_branch_ordinal),
                    kernel_fraction=kernel_fraction,
                    parent_coset_record=tuple(int(value) for value in parent_cosets[parent_branch_ordinal]),
                    parent_affine_branch=parent_affine_branch,
                )
            )
        out.append(
            WyckoffSplitRowProvenance(
                source_group_ordinal=int(source_group_ordinal),
                child_row_offset=int(child_row_offset),
                formula15=formula15_for_child_row(
                    data,
                    int(child_sg),
                    child_row_offset,
                    branch,
                    child_branch,
                ),
                representative_flat_index=int(group["representative_index"]),
                orbit_flat_indices=orbit_flat_indices,
                occurrences=tuple(occurrences),
            )
        )
    return tuple(out)


def get_new_wyckoff_child_rows_formula15_provenance(
    data: SourceTables,
    parent_sg: int,
    child_sg: int,
    wyckoff_index: int,
    subgroup_basis: tuple[int, ...],
    subgroup_origin: FractionRecord,
) -> tuple[WyckoffSplitRowProvenance, ...]:
    """Return the frozen Source split ledger for one exact subgroup embedding.

    Formula15 topology depends only on the Source tables, parent/child rows,
    and the exact integer basis/origin records.  Site parameters are applied
    later, so sites sharing a parent Wyckoff row can safely reuse this result.
    """

    basis = tuple(int(value) for value in subgroup_basis)
    origin = tuple(int(value) for value in subgroup_origin)
    if len(basis) != 9 or len(origin) != 4:
        raise ValueError("Formula15 provenance requires basis9 and origin4")
    origin_record: FractionRecord = (origin[0], origin[1], origin[2], origin[3])
    return _get_new_wyckoff_child_rows_formula15_provenance_cached(
        data,
        int(parent_sg),
        int(child_sg),
        int(wyckoff_index),
        basis,
        origin_record,
    )


def _child_parameter_values(formula: Formula15, parent_params: dict[str, Any]) -> tuple[Fraction, Fraction, Fraction]:
    parent = {
        "x": Fraction(str(parent_params.get("x", 0))),
        "y": Fraction(str(parent_params.get("y", 0))),
        "z": Fraction(str(parent_params.get("z", 0))),
    }
    values: list[Fraction] = []
    for axis in range(3):
        start = 5 * axis
        constant, xcoef, ycoef, zcoef, den = (int(value) for value in formula[start:start + 5])
        if den == 0:
            values.append(Fraction(0))
            continue
        value = (
            Fraction(constant, den)
            + Fraction(xcoef, den) * parent["x"]
            + Fraction(ycoef, den) * parent["y"]
            + Fraction(zcoef, den) * parent["z"]
        )
        values.append(value % 1)
    return values[0], values[1], values[2]


def _evaluate_parent_affine_branch(
    branch: WyckoffBranch,
    parent_params: dict[str, Any],
) -> tuple[Fraction, Fraction, Fraction]:
    params = tuple(Fraction(str(parent_params.get(name, 0))) for name in ("x", "y", "z"))
    vectors = tuple(_fraction_values(record) for record in branch)
    return tuple(
        vectors[0][axis]
        + sum(params[index] * vectors[index + 1][axis] for index in range(3))
        for axis in range(3)
    )  # type: ignore[return-value]


def _expr_to_python(expr: str) -> str:
    expr = expr.strip().replace("−", "-")
    expr = re.sub(r"([+-])\s+(?=\d|[xyz])", r"\1", expr)
    return re.sub(
        r"(?<![A-Za-z0-9_])([+-]?(?:\d+(?:/\d+)?|\d*\.\d+))([xyz])\b",
        r"\1*\2",
        expr,
    )


def _eval_formula_expr(expr: str, params: dict[str, Fraction]) -> Fraction:
    return evaluate_fraction_expression(
        _expr_to_python(expr),
        params,
        f"unknown parameter {{name}} in {expr!r}",
        f"unsupported expression {expr!r}",
        decimal_names=False,
    ) % 1


def _representative_fraction(
    data: SourceTables,
    child_sg: int,
    row: WyckoffRow,
    child_params: tuple[Fraction, Fraction, Fraction],
) -> tuple[Fraction, Fraction, Fraction]:
    try:
        setting = int(data.default_inter_setting_record(int(child_sg))["id"])
        formula = data.inter_wyckoff_formula(int(child_sg), row, setting)
        parts = split_coordinate_expression3(str(formula["formula"]))
        if parts is None:
            raise ValueError("bad child Wyckoff formula")
        params = {"x": child_params[0], "y": child_params[1], "z": child_params[2]}
        folded = tuple(_eval_formula_expr(part, params) for part in parts)
    except Exception:
        vectors = data.wyckoff_fraction_vectors(row)
        xyz = list(vectors[0])
        for param_index in range(3):
            for axis in range(3):
                xyz[axis] += child_params[param_index] * vectors[param_index + 1][axis]
        folded = tuple(value % 1 for value in xyz)
    return folded  # type: ignore[return-value]


def undistorted_rows_from_wyckoff_split(
    data: SourceTables,
    *,
    parent_sg: int,
    child_sg: int | None,
    parent_wyckoff_row_id: Any,
    label_prefix: str,
    parent_params: dict[str, Any] | None,
    subgroup_basis: list[list[float]] | None,
    subgroup_origin: list[float] | tuple[float, float, float] | None,
) -> list[dict[str, Any]]:
    """Return split atom/site/xyz rows for one parent Wyckoff orbit.

    The returned coordinates are child-basis fractional coordinates.  The
    parent row identity is the Source row selected while parsing the parent;
    it is not reconstructed from a CIF display label.
    """

    if child_sg is None or subgroup_basis is None:
        raise ValueError("Wyckoff split requires a child group and basis")
    if isinstance(parent_wyckoff_row_id, bool) or not isinstance(
        parent_wyckoff_row_id, Integral
    ):
        raise TypeError("parent Wyckoff row id must be an exact integer")
    parent_params = parent_params or {}
    parent_rows = tuple(
        row
        for row in data.wyckoff_rows(int(parent_sg))
        if int(row.row_id) == int(parent_wyckoff_row_id)
    )
    if len(parent_rows) != 1:
        raise ValueError("parent Wyckoff row id does not identify one Source row")
    parent_row = parent_rows[0]
    basis = _basis_tuple(subgroup_basis)
    origin = _origin_record(subgroup_origin)
    splits = get_new_wyckoff_child_rows_formula15_provenance(
        data,
        int(parent_sg),
        int(child_sg),
        int(parent_row.offset0) + 1,
        basis,
        origin,
    )

    rows: list[dict[str, Any]] = []
    child_rows = data.wyckoff_rows(int(child_sg))
    for index, split in enumerate(splits, start=1):
        child_offset = int(split.child_row_offset)
        formula = split.formula15
        if child_offset < 1 or child_offset > len(child_rows):
            return []
        child_row = child_rows[child_offset - 1]
        child_params = _child_parameter_values(formula, parent_params)
        representative = _representative_fraction(data, int(child_sg), child_row, child_params)
        xyz = tuple(0.0 if value == 0 else float(value) for value in representative)
        rows.append(
            {
                "label": f"{label_prefix}_{index}",
                "site": f"{data.wyckoff_multiplicity(int(child_sg), child_row)}{child_row.label}",
                "xyz": list(xyz),
                "_wyckoff_formula15": formula,
                "_wyckoff_representative_fraction": _fraction_record(representative),
                "_wyckoff_source_group_ordinal": int(split.source_group_ordinal),
                "_wyckoff_source_representative_flat_index": int(split.representative_flat_index),
                "_wyckoff_source_occurrences": [
                    {
                        "flat_index": int(occurrence.flat_index),
                        "kernel_fraction_ordinal": int(occurrence.kernel_fraction_ordinal),
                        "parent_branch_ordinal": int(occurrence.parent_branch_ordinal),
                        "kernel_fraction": tuple(int(value) for value in occurrence.kernel_fraction),
                        "parent_coset_record": tuple(int(value) for value in occurrence.parent_coset_record),
                        "parent_point_fraction": _fraction_record(
                            _evaluate_parent_affine_branch(
                                occurrence.parent_affine_branch,
                                parent_params,
                            )
                        ),
                    }
                    for occurrence in split.occurrences
                ],
            }
        )
    return rows
