"""Display-only geometry used by the interactive mode viewer."""

from __future__ import annotations

import math
from fractions import Fraction
from typing import Sequence

from distortropy.Backend.modes.presentation import (
    _inverse3 as _presentation_inverse3,
    _row_multiply as _presentation_row_multiply,
    centering_translations,
)


Vector3 = tuple[float, float, float]


def _nearest_integers(value: float) -> tuple[int, ...]:
    low = math.floor(value)
    high = math.ceil(value)
    return (low,) if low == high else (low, high)


def _origin3(origin: Sequence[float] | str) -> Vector3:
    if isinstance(origin, str):
        values = [part.strip() for part in origin.strip().strip("()").split(",")]
        if len(values) != 3:
            raise ValueError(f"invalid viewer origin: {origin!r}")
        return tuple(float(Fraction(value)) for value in values)  # type: ignore[return-value]
    values = list(origin)
    if len(values) >= 4 and float(values[3]) != 0:
        return tuple(float(values[index]) / float(values[3]) for index in range(3))  # type: ignore[return-value]
    if len(values) < 3:
        raise ValueError(f"invalid viewer origin: {origin!r}")
    return tuple(float(values[index]) for index in range(3))  # type: ignore[return-value]


def parent_cell_placement(
    basis: Sequence[Sequence[float]],
    origin: Sequence[float] | str,
    child_symbol: str,
) -> dict[str, list[list[float]] | list[float]]:
    """Return ISOVIZ-style parent basis/origin in child-cell coordinates.

    The exact algebraic origin is ``-origin @ inverse(basis)``. This adapter
    then chooses the lattice/centering-equivalent parent cell whose center is
    nearest to the displayed child-cell center.
    """

    try:
        parent_basis = _presentation_inverse3(basis)
    except ValueError as exc:
        if str(exc) != "mode presentation basis is singular":
            raise
        raise ValueError("viewer basis is singular") from exc
    base_origin = tuple(
        -value
        for value in _presentation_row_multiply(_origin3(origin), parent_basis)
    )
    parent_center_offset = tuple(0.5 * sum(parent_basis[row][axis] for row in range(3)) for axis in range(3))
    candidates: list[tuple[float, Vector3, int, tuple[int, int, int], Vector3]] = []
    for centering_index, centering in enumerate(centering_translations(child_symbol)):
        untranslated_center = tuple(base_origin[axis] + centering[axis] + parent_center_offset[axis] for axis in range(3))
        shifts = tuple(_nearest_integers(0.5 - untranslated_center[axis]) for axis in range(3))
        for x_shift in shifts[0]:
            for y_shift in shifts[1]:
                for z_shift in shifts[2]:
                    integer_shift = (x_shift, y_shift, z_shift)
                    candidate = tuple(base_origin[axis] + centering[axis] + integer_shift[axis] for axis in range(3))
                    center = tuple(candidate[axis] + parent_center_offset[axis] for axis in range(3))
                    score = sum((center[axis] - 0.5) ** 2 for axis in range(3))
                    candidates.append((score, center, centering_index, integer_shift, candidate))
    _, _, _, _, selected_origin = min(candidates, key=lambda item: (round(item[0], 12), item[1], item[2], item[3]))
    return {
        "parent_basis": [list(row) for row in parent_basis],
        "parent_origin": list(selected_origin),
    }
