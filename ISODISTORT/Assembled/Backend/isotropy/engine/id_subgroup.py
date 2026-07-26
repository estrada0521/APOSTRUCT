"""Candidate mapping and origin-congruence logic for ``id_subgroup_``.

The module builds generator mappings, reduces candidate bases, constructs the
translation equations, and returns the first accepted subgroup and origin.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
import math
from typing import Sequence

from ISODISTORT.Assembled.Backend.exactmath import integer_determinant3
from ISODISTORT.Assembled.Backend.isotropy.engine.solve_eqs import solve_eqs_mod_int_first
from ISODISTORT.Assembled.Backend.isotropy.engine.source_data import SourceData


OperationRecord = tuple[int, int, int, int, int]


@dataclass(frozen=True)
class GeneratorMapping:
    candidate_subgroup: int
    candidate_generator_count: int
    per_generator_match_counts: tuple[int, ...]
    per_generator_match_slots: tuple[tuple[int, ...], ...]
    selected_input_generator_slots: tuple[int, ...]
    selected_candidate_generator_ops: tuple[int, ...]


@dataclass(frozen=True)
class OriginEquationSystem:
    matrix: tuple[tuple[int, int, int], ...]
    rhs: tuple[int, ...]
    denominator: int


@dataclass(frozen=True)
class IdSubgroupResult:
    subgroup: int
    basis: tuple[int, ...]
    origin: tuple[int, int, int, int]


def origin_equation_system(
    generator_matrices: Sequence[Sequence[Sequence[int]]],
    generator_diff_records: Sequence[Sequence[int]],
) -> OriginEquationSystem:
    """Build the origin congruences passed to ``solve_eqs_mod_int_``.

    For every accepted generator mapping, this adapter forms three equations

    ``(I - R_i^T) * origin = delta_i``

    and stores them in the column-stride-50 matrix handed to
    ``solve_eqs_mod_int_``.  The right-hand side is integerized with the
    product of the selected generator translation denominators. The transpose
    follows this adapter's row-vector buffer layout rather than introducing a
    separate mathematical convention.
    """

    if len(generator_matrices) != len(generator_diff_records):
        raise ValueError("matrix/diff generator count mismatch")
    denominator = 1
    for record in generator_diff_records:
        denominator *= int(record[3])
    rows: list[tuple[int, int, int]] = []
    rhs: list[int] = []
    identity = ((1, 0, 0), (0, 1, 0), (0, 0, 1))
    for matrix, diff in zip(generator_matrices, generator_diff_records):
        if len(matrix) != 3 or any(len(row) != 3 for row in matrix):
            raise ValueError("generator matrix must be 3x3")
        diff_den = int(diff[3])
        for row in range(3):
            rows.append(tuple(int(identity[row][col]) - int(matrix[col][row]) for col in range(3)))
            rhs.append((denominator * int(diff[row])) // diff_den)
    return OriginEquationSystem(tuple(rows), tuple(rhs), denominator)


def _reduce_fraction_record(values: Sequence[Fraction]) -> tuple[int, int, int, int]:
    denominator = 1
    for value in values:
        denominator = denominator * value.denominator // math.gcd(denominator, value.denominator)
    nums = [int(value * denominator) for value in values]
    gcd = abs(denominator)
    for value in nums:
        gcd = math.gcd(gcd, abs(value))
    if gcd > 1:
        nums = [value // gcd for value in nums]
        denominator //= gcd
    if denominator < 0:
        nums = [-value for value in nums]
        denominator = -denominator
    nums = [value % denominator for value in nums]
    return (nums[0], nums[1], nums[2], denominator)


def _integer_inverse_matrix(flat: Sequence[int]) -> tuple[tuple[int, int, int], ...]:
    inv = SourceData._inverse_fraction_matrix(tuple(int(value) for value in flat[:9]))
    out: list[tuple[int, int, int]] = []
    for row in inv:
        if any(value.denominator != 1 for value in row):
            raise ValueError(f"matrix inverse is not integral: {flat}")
        out.append(tuple(int(value) for value in row))  # type: ignore[arg-type]
    return tuple(out)


def _matinv_numerator_denominator(flat: Sequence[int]) -> tuple[tuple[int, ...], int]:
    inv = SourceData._inverse_fraction_matrix(tuple(int(value) for value in flat[:9]))
    denominator = 1
    for row in inv:
        for value in row:
            denominator = denominator * value.denominator // math.gcd(denominator, value.denominator)
    nums = [int(value * denominator) for row in inv for value in row]
    gcd = abs(denominator)
    for value in nums:
        gcd = math.gcd(gcd, abs(value))
    if gcd > 1:
        nums = [value // gcd for value in nums]
        denominator //= gcd
    if denominator < 0:
        nums = [-value for value in nums]
        denominator = -denominator
    return tuple(nums), denominator


def _vmlt_iso_vector_matrix(
    matrix: Sequence[Sequence[int]],
    record: Sequence[int],
) -> tuple[int, int, int, int]:
    x, y, z, denominator = (int(value) for value in record[:4])
    return (
        x * int(matrix[0][0]) + y * int(matrix[1][0]) + z * int(matrix[2][0]),
        x * int(matrix[0][1]) + y * int(matrix[1][1]) + z * int(matrix[2][1]),
        x * int(matrix[0][2]) + y * int(matrix[1][2]) + z * int(matrix[2][2]),
        denominator,
    )


def _matrix3_rows(flat: Sequence[int]) -> tuple[tuple[int, int, int], ...]:
    return tuple(
        tuple(int(flat[row * 3 + col]) for col in range(3))
        for row in range(3)
    )  # type: ignore[return-value]


def input_operation_records_from_inverse(
    inverse_matrix: Sequence[int],
    inverse_denominator: int,
    input_ops: Sequence[Sequence[int]],
) -> tuple[tuple[int, int, int, int, int], ...]:
    """Transform input operation translations through the inverse subgroup basis."""

    matrix = tuple(
        tuple(int(inverse_matrix[row * 3 + col]) for col in range(3))
        for row in range(3)
    )
    out: list[tuple[int, int, int, int, int]] = []
    for record in input_ops:
        transformed = _vmlt_iso_vector_matrix(matrix, record)
        denominator = int(inverse_denominator) * int(record[3])
        values = [
            Fraction(int(transformed[axis]), denominator)
            for axis in range(3)
        ]
        x, y, z, den = _reduce_fraction_record(values)
        out.append((x, y, z, den, int(record[4])))
    return tuple(out)


def id_subgroup_input_generator_matrices(
    data: SourceData,
    parent_sg: int,
    basis: Sequence[int],
    input_ops: Sequence[Sequence[int]],
) -> tuple[tuple[tuple[int, int, int], ...], ...]:
    """Return integer parent-operation matrices in the subgroup basis."""

    inverse, inverse_denominator = _matinv_numerator_denominator(basis)
    lattice = int(data.space["ispace_lattice"][int(parent_sg) - 1])
    setting = data.ml_lattice_setting_record(lattice)
    denom0, _denom1 = data.ml_lattice_setting_denominators(lattice)
    left = setting[:9]
    right = setting[9:18]
    out: list[tuple[tuple[int, int, int], ...]] = []
    for record in input_ops:
        raw = data._matmlt_iso_storage(
            data._matmlt_iso_storage(
                data._matmlt_iso_storage(inverse, left),
                data.point_operation_matrix(int(record[4])),
            ),
            right,
        )
        raw = data._matmlt_iso_storage(raw, tuple(int(value) for value in basis[:9]))
        reduced, denominator = data._reduce_integer_matrix_denominator(
            raw,
            int(denom0) * int(inverse_denominator),
        )
        if denominator != 1:
            raise ValueError(
                f"id_subgroup input matrix did not reduce to integer: "
                f"parent={parent_sg} op={record[4]} denominator={denominator}"
            )
        out.append(_matrix3_rows(reduced))
    return tuple(out)


_ID_SUBGROUP_SEARCH_VALUES = (0, 1, -1, 2, -2, 3, -3, 4, -4, 6, -6)


def _matrix_det3(matrix: Sequence[Sequence[int]]) -> int:
    return integer_determinant3(
        tuple(matrix[row][column] for row in range(3) for column in range(3))
    )


def _matmlt_flat(data: SourceData, left: Sequence[int], right: Sequence[Sequence[int]] | Sequence[int]) -> tuple[int, ...]:
    if right and isinstance(right[0], (tuple, list)):  # type: ignore[index]
        flat_right = tuple(int(value) for row in right for value in row)  # type: ignore[union-attr]
    else:
        flat_right = tuple(int(value) for value in right)  # type: ignore[arg-type]
    return data._matmlt_iso_storage(tuple(int(value) for value in left[:9]), flat_right)


def nice_lattice2_basis(data: SourceData, lattice: int, basis: Sequence[int]) -> tuple[int, ...]:
    """Normalize a candidate basis with the ``nice_lattice2_`` rules.

    The routine first maps the basis through the ML lattice's second 3x3 block,
    repeatedly shortens each basis vector by adding or subtracting another
    basis vector when that strictly lowers the squared length, then maps the
    result back through the first 3x3 block and applies ``reduc3_``.
    """

    setting = data.ml_lattice_setting_record(int(lattice))
    denom0, _denom1 = data.ml_lattice_setting_denominators(int(lattice))
    reduced_basis = list(_matmlt_flat(data, setting[9:18], basis))

    while True:
        previous = list(reduced_basis)
        for target in range(3):
            target_start = target * 3
            target_len = sum(reduced_basis[target_start + axis] ** 2 for axis in range(3))
            replaced = False
            for other in range(3):
                if other == target:
                    continue
                other_start = other * 3
                for sign in (-1, 1):
                    candidate = [
                        sign * reduced_basis[other_start + axis] + reduced_basis[target_start + axis]
                        for axis in range(3)
                    ]
                    if sum(value * value for value in candidate) < target_len:
                        reduced_basis[target_start : target_start + 3] = candidate
                        replaced = True
                        break
                if replaced:
                    break
        if reduced_basis == previous:
            break

    out = _matmlt_flat(data, setting[:9], reduced_basis)
    reduced, denominator = data._reduce_integer_matrix_denominator(out, denom0)
    if denominator != 1:
        raise ValueError(f"nice_lattice2_ did not reduce to integer basis: lattice={lattice} denominator={denominator}")
    return tuple(int(value) for value in reduced)


def _reduce_fraction_record_no_wrap(values: Sequence[Fraction]) -> tuple[int, int, int, int]:
    denominator = 1
    for value in values:
        denominator = denominator * value.denominator // math.gcd(denominator, value.denominator)
    nums = [int(value * denominator) for value in values]
    gcd = abs(denominator)
    for value in nums:
        gcd = math.gcd(gcd, abs(value))
    if gcd > 1:
        nums = [value // gcd for value in nums]
        denominator //= gcd
    if denominator < 0:
        nums = [-value for value in nums]
        denominator = -denominator
    return (nums[0], nums[1], nums[2], denominator)


def _origin_after_basis(basis: Sequence[int], solution: Sequence[int], denominator: int) -> tuple[int, int, int, int]:
    transformed = _vmlt_iso_vector_matrix(_matrix3_rows(basis), (int(solution[0]), int(solution[1]), int(solution[2]), int(denominator)))
    values = [Fraction(int(transformed[axis]), int(transformed[3])) for axis in range(3)]
    return _reduce_fraction_record_no_wrap(values)


def _nonpivot_count(rowop_matrix: Sequence[int], rows: int, cols: int) -> int:
    pivot_flags = [0] * (int(cols) + 1)
    for row in range(1, int(rows) + 1):
        for col in range(row, int(cols) + 1):
            index = (row - 1) + (col - 1) * 50
            if index < len(rowop_matrix) and int(rowop_matrix[index]) != 0:
                pivot_flags[col] = 1
                break
    return sum(1 for col in range(1, int(cols) + 1) if pivot_flags[col] == 0)


def id_subgroup_basis_candidate_from_rowop(
    rowop_matrix: Sequence[int],
    row_count: int,
    col_count: int,
    base_value_count: int,
    candidate_index: int,
) -> tuple[tuple[int, int, int], ...] | None:
    """Build one integer basis-transform candidate from a reduced row system.

    The reduced equation matrix uses a column-stride-50 layout. Non-pivot
    columns receive mixed-radix trial values from the fixed list
    `[0, 1, -1, 2, -2, 3, -3, 4, -4, 6, -6]`; pivot columns are solved by
    back-substitution. Only determinant `+1` candidates are returned.
    """

    rows = int(row_count)
    cols = int(col_count)
    base = int(base_value_count)
    if rows < 0 or cols < 0 or cols > 9:
        raise ValueError(f"unsupported row/column count: {rows}x{cols}")
    if not 0 < base <= len(_ID_SUBGROUP_SEARCH_VALUES):
        raise ValueError(f"unsupported base value count: {base}")

    def coeff(row: int, col: int) -> int:
        # `row` and `col` are 1-based; the workspace uses leading dimension 50.
        index = (row - 1) + (col - 1) * 50
        if index >= len(rowop_matrix):
            raise ValueError("rowop matrix dump is too narrow")
        return int(rowop_matrix[index])

    pivot_by_row: dict[int, int] = {}
    pivot_flags = [0] * (cols + 1)
    for row in range(1, rows + 1):
        for col in range(row, cols + 1):
            if coeff(row, col) != 0:
                pivot_by_row[row] = col
                pivot_flags[col] = 1
                break

    values = [0] * (cols + 1)
    residue = int(candidate_index) - 1
    for col in range(1, cols + 1):
        if pivot_flags[col] == 0:
            values[col] = _ID_SUBGROUP_SEARCH_VALUES[residue % base]
            residue //= base

    for row in range(1, rows + 1):
        pivot_col = pivot_by_row.get(row)
        if pivot_col is None:
            continue
        pivot = coeff(row, pivot_col)
        total = sum(coeff(row, col) * values[col] for col in range(pivot_col + 1, cols + 1))
        if pivot == 0 or total % pivot != 0:
            return None
        values[pivot_col] = -(total // pivot)

    if cols < 9:
        return None
    matrix = tuple(
        tuple(int(values[row + col * 3]) for col in range(3))
        for row in range(1, 4)
    )
    if _matrix_det3(matrix) != 1:
        return None
    return matrix


def id_subgroup_identify_with_generator_block(
    data: SourceData,
    parent_sg: int,
    basis: Sequence[int],
    input_ops: Sequence[Sequence[int]],
    flag: int,
    generators_conv: Sequence[int] | None = None,
) -> IdSubgroupResult | None:
    """Identify a subgroup in the candidate order used by ``id_subgroup_``.

    ``generators_conv`` optionally supplies an initialized generator block for
    compatibility fixtures. Normal callers read canonical generators from
    ``Source/data_space`` through ``SourceData.get_generators_records``.
    """

    parent = int(parent_sg)
    basis_in = tuple(int(value) for value in basis[:9])
    lattice = int(data.space["ispace_lattice"][parent - 1])
    basis_flat = nice_lattice2_basis(data, lattice, basis_in) if int(flag) >= 0 else basis_in
    active_ops = tuple(tuple(int(value) for value in record[:5]) for record in input_ops if len(record) >= 5 and int(record[4]) != 0)
    inverse, inverse_denominator = _matinv_numerator_denominator(basis_flat)
    local_input_ops = input_operation_records_from_inverse(inverse, inverse_denominator, active_ops)
    local_input_matrices = id_subgroup_input_generator_matrices(data, parent, basis_flat, local_input_ops)
    point_group_code = sum(int(data.space["ipoint_op_code"][int(record[4]) - 1]) for record in local_input_ops)
    try:
        point_group = [int(value) for value in data.space["ipoint_group_code"]].index(point_group_code) + 1
    except ValueError:
        return None

    if int(flag) == 0:
        candidate_range = range(1, 231)
    else:
        exact = abs(int(flag))
        candidate_range = range(exact, exact + 1)

    lattice_constraint: int | None = None
    for base_value_count in (3, 5, 7, 11):
        for candidate_sg in candidate_range:
            if int(data.space["ispace_point_group"][candidate_sg - 1]) != point_group:
                continue
            candidate_lattice = int(data.space["ispace_lattice"][candidate_sg - 1])
            if lattice_constraint is not None and candidate_lattice != lattice_constraint:
                continue
            if generators_conv is None:
                candidate_generators = data.get_generators_records(candidate_sg, "pml")
            else:
                candidate_generators = data.get_generators_records_from_initialized_block(
                    candidate_sg,
                    "pml",
                    tuple(int(value) for value in generators_conv),
                )
            mappings = candidate_generator_mappings(data, candidate_sg, local_input_ops, candidate_generators)
            if not mappings:
                continue
            candidate_matrices_by_op = {
                int(record[4]): data.id_subgroup_candidate_generator_matrix(candidate_sg, int(record[4]))
                for record in candidate_generators
            }
            for mapping in mappings:
                selected_input_matrices = [
                    local_input_matrices[int(slot) - 1]
                    for slot in mapping.selected_input_generator_slots
                ]
                candidate_matrices = [
                    candidate_matrices_by_op[int(op)]
                    for op in mapping.selected_candidate_generator_ops
                ]
                rowop_input = id_subgroup_rowop_setup_matrix(selected_input_matrices, candidate_matrices)
                rowop_matrix = rowop_reduce_matrix(rowop_input, len(candidate_matrices) * 9, 9)
                free_count = _nonpivot_count(rowop_matrix, len(candidate_matrices) * 9, 9)
                limit = int(base_value_count) ** int(free_count)
                for candidate_index in range(1, limit + 1):
                    basis_candidate = id_subgroup_basis_candidate_from_rowop(
                        rowop_matrix,
                        len(candidate_matrices) * 9,
                        9,
                        base_value_count,
                        candidate_index,
                    )
                    if basis_candidate is None:
                        continue
                    diffs = generator_diff_records(
                        data,
                        candidate_sg,
                        mapping.selected_input_generator_slots,
                        local_input_ops,
                        tuple(value for row in basis_candidate for value in row),
                        candidate_generators,
                    )
                    equations = origin_equation_system(candidate_matrices, diffs)
                    solve = solve_eqs_mod_int_first(equations.matrix, equations.rhs, equations.denominator)
                    lattice_constraint = candidate_lattice
                    if not solve.success:
                        continue
                    basis_out = _matmlt_flat(data, basis_flat, basis_candidate)
                    origin_out = _origin_after_basis(basis_out, solve.solution, solve.denominator)
                    return IdSubgroupResult(candidate_sg, basis_out, origin_out)
    return None


def _factor_row(values: list[int]) -> None:
    nonzero = [abs(value) for value in values if value != 0]
    if not nonzero:
        return
    divisor = min(nonzero)
    while divisor > 1:
        if all(value % divisor == 0 for value in values):
            for index, value in enumerate(values):
                values[index] = value // divisor
            return
        divisor -= 1


def rowop_reduce_matrix(
    matrix: Sequence[int],
    row_count: int,
    col_count: int,
    *,
    leading: int = 50,
) -> tuple[int, ...]:
    """Reduce an ISO column-stride integer matrix with ``rowop_`` rules.

    The routine performs integer row elimination over a `leading=50` matrix.
    It searches pivot columns from the diagonal with a skipped-column counter,
    swaps rows when needed, eliminates the pivot column from every other row,
    and reduces changed rows with `factor_`.
    """

    rows = int(row_count)
    cols = int(col_count)
    if rows < 0 or cols < 0:
        raise ValueError("negative matrix shape")
    out = [int(value) for value in matrix]
    min_size = max(0, (cols - 1) * leading + rows)
    if len(out) < min_size:
        out.extend([0] * (min_size - len(out)))

    def index(row: int, col: int) -> int:
        # `row` and `col` are 1-based.
        return (row - 1) + (col - 1) * leading

    skipped = 0
    pivot_row = 1
    while pivot_row <= rows:
        pivot_col = skipped + pivot_row
        if cols < pivot_col:
            break
        while out[index(pivot_row, pivot_col)] == 0:
            swap_row = None
            for row in range(pivot_row + 1, rows + 1):
                if out[index(row, pivot_col)] != 0:
                    swap_row = row
                    break
            if swap_row is not None:
                for col in range(1, cols + 1):
                    a = index(pivot_row, col)
                    b = index(swap_row, col)
                    out[a], out[b] = out[b], out[a]
                break
            skipped += 1
            pivot_col += 1
            if cols < pivot_col:
                return tuple(out)

        pivot = out[index(pivot_row, pivot_col)]
        for row in range(1, rows + 1):
            elim = out[index(row, pivot_col)]
            if row == pivot_row or elim == 0:
                continue
            for col in range(1, cols + 1):
                if col == pivot_col:
                    continue
                out[index(row, col)] = pivot * out[index(row, col)] - out[index(pivot_row, col)] * elim
            out[index(row, pivot_col)] = 0
            row_values = [out[index(row, col)] for col in range(1, cols + 1)]
            _factor_row(row_values)
            for col, value in enumerate(row_values, start=1):
                out[index(row, col)] = value
        pivot_row += 1
    return tuple(out)


def id_subgroup_rowop_setup_matrix(
    selected_input_matrices: Sequence[Sequence[Sequence[int]]],
    candidate_generator_matrices: Sequence[Sequence[Sequence[int]]],
    *,
    leading: int = 50,
) -> tuple[int, ...]:
    """Build the `id_subgroup_` matrix handed to ``rowop_``.

    For every accepted generator mapping, upstream constrains the trial basis
    transform `T` by

    ``C * T - T * R = 0``

    where `C` is the candidate subgroup generator matrix and `R` is the
    selected input generator matrix, both in the same integer setting used by
    the origin solve.  The nine unknowns of `T` are stored column-major, and
    the nine equations are emitted in column-major order.
    """

    if len(selected_input_matrices) != len(candidate_generator_matrices):
        raise ValueError("input/candidate generator matrix count mismatch")
    row_count = len(selected_input_matrices) * 9
    out = [0 for _ in range(max(row_count, 1) + (9 - 1) * leading)]

    def coeff_index(row: int, col: int) -> int:
        # `row` and `col` are 1-based; storage is column-stride-50.
        return (row - 1) + (col - 1) * leading

    def value(matrix: Sequence[Sequence[int]], row: int, col: int) -> int:
        return int(matrix[row][col])

    equation_row = 1
    for input_matrix, candidate_matrix in zip(selected_input_matrices, candidate_generator_matrices):
        if len(input_matrix) != 3 or len(candidate_matrix) != 3:
            raise ValueError("generator matrices must be 3x3")
        for col in range(3):
            for row in range(3):
                coefficients = [0 for _ in range(9)]
                for inner in range(3):
                    coefficients[inner + col * 3] += value(candidate_matrix, row, inner)
                    coefficients[row + inner * 3] -= value(input_matrix, inner, col)
                for variable, coefficient in enumerate(coefficients, start=1):
                    out[coeff_index(equation_row, variable)] = coefficient
                equation_row += 1
    return tuple(out)


def generator_diff_records(
    data: SourceData,
    candidate_subgroup: int,
    selected_input_generator_slots: Sequence[int],
    local_input_ops: Sequence[Sequence[int]],
    basis_transform: Sequence[int],
    candidate_generator_records: Sequence[Sequence[int]] | None = None,
) -> tuple[tuple[int, int, int, int], ...]:
    """Build the translation-difference records used by the origin equations.

    This boundary starts after the input operation records have already been
    transformed into the candidate lattice. This adapter applies the
    inverse of the accepted basis transform to each
    selected input generator translation, subtracts the corresponding
    candidate generator translation, reduces the fraction, and finally wraps
    numerators into `[0, denominator)`.
    """

    inverse = _integer_inverse_matrix(basis_transform)
    candidate_generators = (
        tuple(tuple(int(value) for value in record[:5]) for record in candidate_generator_records)
        if candidate_generator_records is not None
        else data.get_generators_records(int(candidate_subgroup), "pml")
    )
    out: list[tuple[int, int, int, int]] = []
    for index, slot in enumerate(selected_input_generator_slots):
        input_record = local_input_ops[int(slot) - 1]
        transformed = _vmlt_iso_vector_matrix(inverse, input_record)
        candidate = candidate_generators[index]
        values = [
            Fraction(int(transformed[axis]), int(transformed[3]))
            - Fraction(int(candidate[axis]), int(candidate[3]))
            for axis in range(3)
        ]
        out.append(_reduce_fraction_record(values))
    return tuple(out)


def _point_code(data: SourceData, point_op: int) -> int:
    """Return the point-operation code used by ``id_subgroup_``.

    This uses ``data_space:ipoint_op_code``, not the nearby
    ``ipoint_op_order`` table: point operations 27 and 28 compare equal here
    because both have code 81.
    """

    return int(data.space["ipoint_op_code"][int(point_op) - 1])


def _point_mul(data: SourceData, left: int, right: int) -> int:
    """Return the 1-based product from ``ipoint_op_mlt``.

    The table is row-major in 1-based point-operation labels. Products use
    ``(right, left)`` indexing.
    """

    return int(data.space["ipoint_op_mlt"][(int(right) - 1) * 72 + (int(left) - 1)])


def _closure_accepts(
    data: SourceData,
    input_ops: Sequence[OperationRecord],
    selected_input_slots: tuple[int, ...],
    selected_candidate_ops: tuple[int, ...],
) -> bool:
    input_slot_by_point_op = {int(record[4]): index + 1 for index, record in enumerate(input_ops)}
    mapped_index_by_input_slot: dict[int, int] = {}
    mapped_index_by_candidate_op: dict[int, int] = {}
    mapped_input_slots: list[int] = []
    mapped_candidate_ops: list[int] = []

    for slot, candidate_op in zip(selected_input_slots, selected_candidate_ops):
        if slot in mapped_index_by_input_slot:
            return False
        index = len(mapped_input_slots) + 1
        mapped_index_by_input_slot[int(slot)] = index
        mapped_index_by_candidate_op[int(candidate_op)] = index
        mapped_input_slots.append(int(slot))
        mapped_candidate_ops.append(int(candidate_op))

    changed = True
    while changed:
        changed = False
        mapped_count = len(mapped_input_slots)
        for left_index in range(mapped_count):
            for right_index in range(mapped_count):
                candidate_product = _point_mul(
                    data,
                    mapped_candidate_ops[left_index],
                    mapped_candidate_ops[right_index],
                )
                input_left_op = int(input_ops[mapped_input_slots[left_index] - 1][4])
                input_right_op = int(input_ops[mapped_input_slots[right_index] - 1][4])
                input_product_op = _point_mul(data, input_left_op, input_right_op)
                input_product_slot = input_slot_by_point_op.get(input_product_op)
                if input_product_slot is None:
                    return False

                existing_candidate_index = mapped_index_by_candidate_op.get(candidate_product)
                existing_input_index = mapped_index_by_input_slot.get(input_product_slot)
                if existing_candidate_index is None:
                    if existing_input_index is not None:
                        return False
                    index = len(mapped_input_slots) + 1
                    mapped_index_by_candidate_op[candidate_product] = index
                    mapped_index_by_input_slot[input_product_slot] = index
                    mapped_candidate_ops.append(candidate_product)
                    mapped_input_slots.append(input_product_slot)
                    changed = True
                    continue
                if mapped_input_slots[existing_candidate_index - 1] != input_product_slot:
                    return False
                if existing_input_index != existing_candidate_index:
                    return False
    return len(mapped_input_slots) >= len(input_ops)


def candidate_generator_mappings(
    data: SourceData,
    candidate_subgroup: int,
    input_ops: Sequence[OperationRecord],
    candidate_generators: Sequence[OperationRecord],
) -> tuple[GeneratorMapping, ...]:
    """Enumerate accepted generator mappings for one candidate subgroup."""

    match_slots: list[tuple[int, ...]] = []
    for generator in candidate_generators:
        code = _point_code(data, int(generator[4]))
        slots = tuple(
            index + 1
            for index, operation in enumerate(input_ops)
            if _point_code(data, int(operation[4])) == code
        )
        match_slots.append(slots)
    if any(not slots for slots in match_slots):
        return ()

    accepted: list[GeneratorMapping] = []
    candidate_ops = tuple(int(record[4]) for record in candidate_generators)
    total = 1
    for slots in match_slots:
        total *= len(slots)
    for ordinal in range(total):
        value = ordinal
        selected_list: list[int] = []
        for slots in match_slots:
            selected_list.append(int(slots[value % len(slots)]))
            value //= len(slots)
        selected = tuple(selected_list)
        if len(set(selected)) != len(selected):
            continue
        selected_slots = tuple(int(value) for value in selected)
        if not _closure_accepts(data, input_ops, selected_slots, candidate_ops):
            continue
        accepted.append(
            GeneratorMapping(
                candidate_subgroup=int(candidate_subgroup),
                candidate_generator_count=len(candidate_generators),
                per_generator_match_counts=tuple(len(slots) for slots in match_slots),
                per_generator_match_slots=tuple(match_slots),
                selected_input_generator_slots=selected_slots,
                selected_candidate_generator_ops=candidate_ops,
            )
        )
    return tuple(accepted)
