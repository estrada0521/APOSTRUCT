"""Unit-cell and vector-coordinate helpers used by final mode assembly."""

from __future__ import annotations

from fractions import Fraction
import math
from typing import Iterable

from ISODISTORT.Assembled.Backend.modes.engine.input import Case

def _unit_cell_rows(case: Case) -> tuple[tuple[float, float, float], ...]:
    a, b, c, alpha, beta, gamma = (float(value) for value in case.params)
    ar = math.radians(alpha)
    br = math.radians(beta)
    gr = math.radians(gamma)
    av = (a, 0.0, 0.0)
    bv = (b * math.cos(gr), b * math.sin(gr), 0.0)
    cx = c * math.cos(br)
    cy = c * (math.cos(ar) - math.cos(br) * math.cos(gr)) / math.sin(gr)
    cz = math.sqrt(max(c * c - cx * cx - cy * cy, 0.0))
    rows = (av, bv, (cx, cy, cz))
    unit_rows = []
    for row in rows:
        length = math.sqrt(sum(value * value for value in row))
        unit_rows.append(tuple(value / length if length else value for value in row))
    return tuple(unit_rows)  # type: ignore[return-value]
def _fractional_vector_to_cartesian(vector: Iterable[float], unit_rows: tuple[tuple[float, float, float], ...]) -> list[float]:
    values = tuple(float(value) for value in vector)
    return [
        sum(values[index] * unit_rows[index][axis] for index in range(3))
        for axis in range(3)
    ]
def _apply_fractional_vector_transform(
    vector: Iterable[float],
    matrix: tuple[tuple[Fraction, Fraction, Fraction], ...],
) -> list[float]:
    values = tuple(float(value) for value in vector)
    return [
        sum(values[row] * float(matrix[row][col]) for row in range(3))
        for col in range(3)
    ]
