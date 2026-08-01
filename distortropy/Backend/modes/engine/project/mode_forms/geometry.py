"""Unit-cell and vector-coordinate helpers used by final mode assembly."""

from __future__ import annotations

from fractions import Fraction
from typing import Iterable


def _apply_fractional_vector_transform(
    vector: Iterable[float],
    matrix: tuple[tuple[Fraction, Fraction, Fraction], ...],
) -> list[float]:
    values = tuple(float(value) for value in vector)
    return [
        sum(values[row] * float(matrix[row][col]) for row in range(3))
        for col in range(3)
    ]
