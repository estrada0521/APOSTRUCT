"""Normalization and small vector operations for printed mode vectors."""

from __future__ import annotations

def _normalize_mode_vectors(vectors: list[list[float]]) -> list[list[float]]:
    scale = max((abs(value) for vector in vectors for value in vector), default=0.0)
    if scale <= 1e-12:
        return vectors
    if abs(scale - 1.0) <= 1e-12:
        return vectors
    return [[value / scale for value in vector] for vector in vectors]


def _dominant_mode_component(mode: list[list[float]]) -> int:
    totals = [sum(abs(vector[axis]) for vector in mode) for axis in range(3)]
    return max(range(3), key=lambda axis: totals[axis])
def _same_mode(left: list[list[float]], right: list[list[float]], tol: float = 1e-9) -> bool:
    if len(left) != len(right):
        return False
    return all(
        abs(float(a[axis]) - float(b[axis])) <= tol
        for a, b in zip(left, right)
        for axis in range(3)
    )
def _active_mode_indices(mode: list[list[float]]) -> set[int]:
    return {
        index
        for index, vector in enumerate(mode)
        if any(abs(float(component)) > 1e-12 for component in vector)
    }
def _add_mode_vectors(left: list[list[float]], right: list[list[float]]) -> list[list[float]]:
    return [
        [float(a[axis]) + float(b[axis]) for axis in range(3)]
        for a, b in zip(left, right)
    ]
