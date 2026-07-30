"""Modular equation solving for the ``id_subgroup_`` kernel.

The ``ifirst=1`` contract returns one representative solution of
``A x = rhs`` modulo an integer denominator, using the numerator/denominator
record convention shared by the Source operation tables.
"""

from __future__ import annotations

from dataclasses import dataclass
from math import gcd, lcm
from typing import Sequence


@dataclass(frozen=True)
class ModularSolveResult:
    success: bool
    denominator: int = 1
    solution: tuple[int, ...] = ()
    free_count: int = 0
    free_vectors: tuple[tuple[int, ...], ...] = ()


def _trunc_div(numerator: int, denominator: int) -> int:
    if denominator == 0:
        raise ZeroDivisionError("integer division by zero")
    sign = -1 if (numerator < 0) ^ (denominator < 0) else 1
    return sign * (abs(numerator) // abs(denominator))


def _factor(values: Sequence[int]) -> int:
    result = 0
    for value in values:
        result = gcd(result, abs(int(value)))
    return result or 1


def _mat_vec(matrix: Sequence[Sequence[int]], vector: Sequence[int]) -> tuple[int, ...]:
    return tuple(sum(int(row[col]) * int(vector[col]) for col in range(len(vector))) for row in matrix)


def _transpose(matrix: Sequence[Sequence[int]]) -> list[list[int]]:
    return [list(row) for row in zip(*matrix)] if matrix else []


def _has_offdiagonal(matrix: Sequence[Sequence[int]]) -> bool:
    row_count = len(matrix)
    col_count = len(matrix[0]) if matrix else 0
    for col in range(col_count):
        for row in range(row_count):
            if row != col and int(matrix[row][col]) != 0:
                return True
    return False


def _row_echelon_form_iso(
    matrix: Sequence[Sequence[int]],
    transform: Sequence[Sequence[int]],
) -> tuple[list[list[int]], list[list[int]]]:
    """Apply ``row_echelon_form_`` rules for the integer SNF solver.

    The ISO routine augments the active matrix with the row-transform matrix,
    repeatedly chooses the smallest nonzero pivot in the active column, applies
    C/Fortran truncating row elimination, and keeps the pivot positive.
    """

    row_count = len(matrix)
    col_count = len(matrix[0]) if matrix else 0
    transform_width = len(transform[0]) if transform else row_count
    augmented = [
        [int(value) for value in matrix[row]] + [int(value) for value in transform[row]]
        for row in range(row_count)
    ]
    for pivot_index in range(row_count):
        while pivot_index < row_count:
            nonzero_rows = [
                row
                for row in range(pivot_index, row_count)
                if pivot_index < col_count and augmented[row][pivot_index] != 0
            ]
            if not nonzero_rows:
                break
            best_row = min(nonzero_rows, key=lambda row: abs(augmented[row][pivot_index]))
            if best_row != pivot_index:
                augmented[pivot_index], augmented[best_row] = augmented[best_row], augmented[pivot_index]
            if augmented[pivot_index][pivot_index] < 0:
                augmented[pivot_index] = [-value for value in augmented[pivot_index]]
            if len(nonzero_rows) == 1:
                break
            pivot = augmented[pivot_index][pivot_index]
            for row in range(pivot_index + 1, row_count):
                value = augmented[row][pivot_index]
                if value == 0:
                    continue
                quotient = _trunc_div(value, pivot)
                for col in range(pivot_index, col_count + transform_width):
                    augmented[row][col] -= quotient * augmented[pivot_index][col]
    return (
        [row[:col_count] for row in augmented],
        [row[col_count:] for row in augmented],
    )


def _smith_normal_decomp_iso(
    rows: tuple[tuple[int, ...], ...],
) -> tuple[tuple[tuple[int, ...], ...], tuple[tuple[int, ...], ...], tuple[tuple[int, ...], ...]]:
    """Compute the ``smith_normal_form_`` decomposition needed for solution representatives."""

    matrix = [list(row) for row in rows]
    row_count = len(matrix)
    col_count = len(matrix[0]) if matrix else 0
    left = [[1 if row == col else 0 for col in range(row_count)] for row in range(row_count)]
    right = [[1 if row == col else 0 for col in range(col_count)] for row in range(col_count)]
    for _ in range(100):
        matrix, left = _row_echelon_form_iso(matrix, left)
        if not _has_offdiagonal(matrix):
            break
        transposed = _transpose(matrix)
        transposed, right = _row_echelon_form_iso(transposed, right)
        matrix = _transpose(transposed)
        if not _has_offdiagonal(matrix):
            break
    else:
            raise RuntimeError("Smith normal form did not converge")
    # Upstream transposes the column-transform buffer before returning.
    return (
        tuple(tuple(int(value) for value in row) for row in matrix),
        tuple(tuple(int(value) for value in row) for row in left),
        tuple(tuple(int(value) for value in row) for row in _transpose(right)),
    )


def solve_eqs_mod_int_first(
    matrix: Sequence[Sequence[int]],
    rhs: Sequence[int],
    modulus: int,
) -> ModularSolveResult:
    """Return the first ISO representative for ``matrix*x = rhs (mod modulus)``.

    The ISO routine converts the integer system to Smith normal form, picks the
    first solution branch for every nonzero diagonal entry, multiplies back by
    the right transform, reduces modulo the product denominator, then divides
    solution numerators and denominator by their common factor.
    """

    rows = tuple(tuple(int(value) for value in row) for row in matrix)
    if not rows:
        return ModularSolveResult(success=True, denominator=1, solution=())
    row_count = len(rows)
    col_count = len(rows[0]) if rows else 0
    if any(len(row) != col_count for row in rows):
        raise ValueError("ragged matrix")
    if len(rhs) != row_count:
        raise ValueError("rhs length does not match matrix row count")
    if modulus == 0:
        raise ValueError("modulus must be nonzero")

    diagonal, left, right = _smith_normal_decomp_iso(rows)
    transformed_rhs = tuple(
        0 if value != 0 and value % int(modulus) == 0 else value
        for value in _mat_vec(left, tuple(int(value) for value in rhs))
    )

    diag_values: list[int] = []
    for index in range(min(len(diagonal), col_count)):
        diag = int(diagonal[index][index])
        if diag != 0:
            diag_values.append(abs(diag))
    if not diag_values:
        for value in transformed_rhs:
            if int(value) % int(modulus) != 0:
                return ModularSolveResult(success=False)
        return ModularSolveResult(success=True, denominator=1, solution=tuple(0 for _ in range(col_count)))

    diag_lcm = 1
    for value in diag_values:
        diag_lcm = lcm(diag_lcm, abs(int(value)))
    product_denominator = diag_lcm * int(modulus)
    y = [0 for _ in range(col_count)]
    for index in range(min(len(diagonal), col_count)):
        diag = int(diagonal[index][index])
        value = int(transformed_rhs[index])
        if diag == 0:
            if value % int(modulus) != 0:
                return ModularSolveResult(success=False)
            continue
        # This is the p1=1 branch in solve_eqs_mod_int_: choose the first of
        # |diag| branches.  C/Fortran integer division truncates toward zero.
        y[index] = _trunc_div(_trunc_div(product_denominator * value, int(modulus)), diag)

    for index in range(min(len(diagonal), col_count), len(transformed_rhs)):
        if int(transformed_rhs[index]) % int(modulus) != 0:
            return ModularSolveResult(success=False)

    raw_solution = _mat_vec(right, tuple(y))
    solution = tuple(int(value) % product_denominator for value in raw_solution)
    common = _factor((product_denominator, *solution))
    return ModularSolveResult(
        success=True,
        denominator=product_denominator // common,
        solution=tuple(value // common for value in solution),
    )
