"""Shared exact three-dimensional linear algebra."""

from __future__ import annotations

from fractions import Fraction
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


__all__ = [
    "FractionMatrix3",
    "FractionVector3",
    "fraction_identity3",
    "fraction_matrix_inverse3",
    "fraction_matrix_multiply3",
    "fraction_row_multiply3",
    "integer_adjugate3",
    "integer_determinant3",
]
