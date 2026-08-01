"""Integer lattice arithmetic for Source operation records."""

from __future__ import annotations

from math import gcd
from functools import lru_cache, reduce

from distortropy.Backend.exactmath import integer_determinant3
from distortropy.Backend.lattice_quotient import (
    integer_inverse_denominator,
    integral_row_images_source_order,
)


def ndet(matrix: list[int] | tuple[int, ...]) -> int:
    """Return the 3x3 integer determinant used by ``ndet_`` callers."""

    values = tuple(int(value) for value in matrix[:9])
    if len(values) < 9:
        raise IndexError("list index out of range")
    return integer_determinant3(values)


def _gcd_many(values: list[int]) -> int:
    nonzero = [abs(value) for value in values if value != 0]
    if not nonzero:
        return 1
    return reduce(gcd, nonzero)


def _reduce_fraction_vector(vector: tuple[int, int, int, int]) -> tuple[int, int, int, int]:
    x, y, z, den = (int(item) for item in vector)
    if den == 0:
        raise ValueError("zero denominator")
    if den < 0:
        x, y, z, den = -x, -y, -z, -den
    factor = _gcd_many([x, y, z, den])
    return x // factor, y // factor, z // factor, den // factor


def matinv_denominator(matrix: list[int] | tuple[int, ...]) -> int:
    """Return the reduced denominator produced by ``matinv_``."""

    det = ndet(matrix)
    if det == 0:
        raise ValueError("singular matrix")
    return integer_inverse_denominator(tuple(int(value) for value in matrix[:9]))


@lru_cache(maxsize=2048)
def _get_new_fractionals_cached(
    matrix: tuple[int, ...],
) -> tuple[tuple[int, int, int, int], ...]:
    """Return Source-ordered quotient translations for an integer basis.

    The shared exact quotient kernel generates the same ``(x, y, z)``
    solutions and restores the Source lexicographic scan order without
    enumerating the full modular cube.
    """

    return integral_row_images_source_order(matrix, matinv_denominator(matrix))


def get_new_fractionals(matrix: list[int] | tuple[int, ...]) -> tuple[tuple[int, int, int, int], ...]:
    """Enumerate ``get_new_fractionals_`` representatives for a 3x3 integer lattice."""

    return _get_new_fractionals_cached(tuple(int(value) for value in matrix[:9]))
