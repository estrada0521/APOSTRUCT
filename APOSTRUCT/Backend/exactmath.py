"""Shared exact three-dimensional linear algebra."""

from __future__ import annotations

from fractions import Fraction
import math
from typing import Sequence, TypeAlias


FractionVector3: TypeAlias = tuple[Fraction, Fraction, Fraction]
FractionMatrix3: TypeAlias = tuple[FractionVector3, FractionVector3, FractionVector3]


def integer_determinant3(matrix: Sequence[int]) -> int:
    values = tuple(int(value) for value in matrix)
    if len(values) != 9:
        raise ValueError(f"expected 9 matrix entries, got {len(values)}")
    return (
        values[0] * (values[4] * values[8] - values[5] * values[7])
        - values[1] * (values[3] * values[8] - values[5] * values[6])
        + values[2] * (values[3] * values[7] - values[4] * values[6])
    )


def integer_adjugate3(matrix: Sequence[int]) -> tuple[int, ...]:
    values = tuple(int(value) for value in matrix)
    if len(values) != 9:
        raise ValueError(f"expected 9 matrix entries, got {len(values)}")
    return (
        values[4] * values[8] - values[5] * values[7],
        values[2] * values[7] - values[1] * values[8],
        values[1] * values[5] - values[2] * values[4],
        values[5] * values[6] - values[3] * values[8],
        values[0] * values[8] - values[2] * values[6],
        values[2] * values[3] - values[0] * values[5],
        values[3] * values[7] - values[4] * values[6],
        values[1] * values[6] - values[0] * values[7],
        values[0] * values[4] - values[1] * values[3],
    )


def fraction_identity3() -> FractionMatrix3:
    return (
        (Fraction(1), Fraction(0), Fraction(0)),
        (Fraction(0), Fraction(1), Fraction(0)),
        (Fraction(0), Fraction(0), Fraction(1)),
    )


def fraction_record3(values: Sequence[Fraction]) -> tuple[int, int, int, int]:
    """Encode a three-dimensional rational vector as reduced integers."""

    denominator = math.lcm(*(Fraction(value).denominator for value in values), 1)
    numerators = [int(Fraction(value) * denominator) for value in values]
    divisor = math.gcd(
        denominator,
        math.gcd(abs(numerators[0]), math.gcd(abs(numerators[1]), abs(numerators[2]))),
    )
    if divisor > 1:
        denominator //= divisor
        numerators = [value // divisor for value in numerators]
    return numerators[0], numerators[1], numerators[2], denominator


def fraction_matrix_multiply3(
    left: Sequence[Sequence[Fraction]],
    right: Sequence[Sequence[Fraction]],
) -> FractionMatrix3:
    return tuple(
        tuple(
            sum(left[row][inner] * right[inner][column] for inner in range(3))
            for column in range(3)
        )
        for row in range(3)
    )  # type: ignore[return-value]


def fraction_row_multiply3(
    row: Sequence[Fraction],
    matrix: Sequence[Sequence[Fraction]],
) -> FractionVector3:
    return tuple(
        sum(row[inner] * matrix[inner][column] for inner in range(3))
        for column in range(3)
    )  # type: ignore[return-value]


def fraction_matrix_inverse3(
    matrix: Sequence[Sequence[Fraction]],
) -> FractionMatrix3:
    a, b, c = matrix[0]
    d, e, f = matrix[1]
    g, h, i = matrix[2]
    determinant = a * (e * i - f * h) - b * (d * i - f * g) + c * (d * h - e * g)
    if determinant == 0:
        raise ValueError(f"singular exact 3x3 matrix: {matrix}")
    return (
        ((e * i - f * h) / determinant, (c * h - b * i) / determinant, (b * f - c * e) / determinant),
        ((f * g - d * i) / determinant, (a * i - c * g) / determinant, (c * d - a * f) / determinant),
        ((d * h - e * g) / determinant, (b * g - a * h) / determinant, (a * e - b * d) / determinant),
    )


def fraction_rref(
    values: Sequence[Sequence[Fraction]],
    coefficient_columns: int | None = None,
) -> tuple[list[list[Fraction]], tuple[int, ...]]:
    """Return exact reduced row-echelon form with stable first-pivot choice."""

    matrix = [list(row) for row in values]
    columns = coefficient_columns if coefficient_columns is not None else (len(matrix[0]) if matrix else 0)
    pivot_columns: list[int] = []
    pivot_row = 0
    for column in range(columns):
        pivot = next(
            (row for row in range(pivot_row, len(matrix)) if matrix[row][column] != 0),
            None,
        )
        if pivot is None:
            continue
        matrix[pivot_row], matrix[pivot] = matrix[pivot], matrix[pivot_row]
        scale = matrix[pivot_row][column]
        matrix[pivot_row] = [value / scale for value in matrix[pivot_row]]
        for row in range(len(matrix)):
            if row == pivot_row or matrix[row][column] == 0:
                continue
            factor = matrix[row][column]
            matrix[row] = [
                value - factor * matrix[pivot_row][index]
                for index, value in enumerate(matrix[row])
            ]
        pivot_columns.append(column)
        pivot_row += 1
        if pivot_row == len(matrix):
            break
    return matrix, tuple(pivot_columns)


__all__ = [
    "FractionMatrix3",
    "FractionVector3",
    "fraction_identity3",
    "fraction_record3",
    "fraction_matrix_inverse3",
    "fraction_matrix_multiply3",
    "fraction_row_multiply3",
    "fraction_rref",
    "integer_adjugate3",
    "integer_determinant3",
]
