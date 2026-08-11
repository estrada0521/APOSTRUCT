"""Periodic fractional-cell point identity helpers."""

from __future__ import annotations

from typing import Any, Sequence


def periodic_float_close3(left: Sequence[Any], right: Sequence[Any], tolerance: float) -> bool:
    return all(abs(((float(left[i]) - float(right[i]) + 0.5) % 1.0) - 0.5) <= tolerance for i in range(3))
