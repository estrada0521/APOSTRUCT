"""Periodic fractional-cell point identity helpers."""

from __future__ import annotations

from fractions import Fraction
from typing import Any, Sequence


def periodic_numeric_close3(left: Sequence[Any], right: Sequence[Any], tolerance: float) -> bool:
    """Compare periodic coordinates without coercing operands before subtraction."""

    return all(abs(((left[i] - right[i] + 0.5) % 1.0) - 0.5) <= tolerance for i in range(3))


def periodic_float_close3(left: Sequence[Any], right: Sequence[Any], tolerance: float) -> bool:
    return all(abs(((float(left[i]) - float(right[i]) + 0.5) % 1.0) - 0.5) <= tolerance for i in range(3))


def periodic_fraction_float_close3(left: Sequence[Any], right: Sequence[Any], tolerance: float) -> bool:
    return all(
        abs(
            float(Fraction(left[axis]) - Fraction(right[axis]))
            - round(float(Fraction(left[axis]) - Fraction(right[axis])))
        )
        <= tolerance
        for axis in range(3)
    )
