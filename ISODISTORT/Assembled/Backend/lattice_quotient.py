"""Exact finite-quotient enumeration for three-dimensional integer lattices."""

from __future__ import annotations

from fractions import Fraction
from functools import lru_cache
from itertools import product
from math import gcd
from typing import Sequence

from ISODISTORT.Assembled.Backend.exactmath import (
    integer_adjugate3,
    integer_determinant3,
)


ModularVector3 = tuple[int, int, int]


def _cyclic_quotient_offsets(
    generator: ModularVector3,
    subgroup: set[ModularVector3],
    modulus: int,
) -> tuple[ModularVector3, ...]:
    """Return ``0,-g,-2g,...`` through the first subgroup recurrence."""

    negative = tuple((-value) % modulus for value in generator)
    offsets = [(0, 0, 0)]
    current = negative
    while current not in subgroup:
        offsets.append(current)
        current = tuple(
            (current[axis] + negative[axis]) % modulus for axis in range(3)
        )
    return tuple(offsets)


def reciprocal_mesh_source_order(
    basis: Sequence[int],
) -> tuple[tuple[Fraction, Fraction, Fraction], ...]:
    """Return the exact reciprocal quotient in Source first-occurrence order.

    Source scans each reciprocal coordinate as ``0,-1,-2,...`` modulo the
    determinant.  Factoring that scan through the cyclic subgroup chain
    ``<c> <= <b,c> <= <a,b,c>`` preserves its nested first/second/third order
    while visiting exactly ``abs(det)`` quotient points.
    """

    values = tuple(int(value) for value in basis)
    if len(values) < 9:
        raise IndexError("tuple index out of range")
    signed_det = integer_determinant3(values[:9])
    det = abs(signed_det)
    if det <= 0:
        return ()
    if det == 1:
        return ((Fraction(0), Fraction(0), Fraction(0)),)
    adjugate = integer_adjugate3(values[:9])
    sign = 1 if signed_det > 0 else -1
    generators = tuple(
        tuple(sign * adjugate[3 * row + axis] % det for row in range(3))
        for axis in range(3)
    )
    third = _cyclic_quotient_offsets(generators[2], {(0, 0, 0)}, det)
    third_group = set(third)
    second = _cyclic_quotient_offsets(generators[1], third_group, det)
    second_third_group = {
        tuple((second_row[axis] + third_row[axis]) % det for axis in range(3))
        for second_row in second
        for third_row in third
    }
    first = _cyclic_quotient_offsets(generators[0], second_third_group, det)
    return tuple(
        tuple(
            Fraction(
                (first_row[axis] + second_row[axis] + third_row[axis]) % det,
                det,
            )
            for axis in range(3)
        )
        for first_row in first
        for second_row in second
        for third_row in third
    )


def integer_inverse_denominator(matrix: Sequence[int]) -> int:
    """Return the reduced denominator of the exact 3x3 inverse numerator."""

    values = tuple(int(value) for value in matrix)
    if len(values) != 9:
        raise ValueError(f"expected 9 matrix entries, got {len(values)}")
    determinant = integer_determinant3(values)
    if determinant == 0:
        raise ValueError(f"singular basis matrix: {values}")
    adjugate = integer_adjugate3(values)
    factor = abs(determinant)
    for value in adjugate:
        factor = gcd(factor, abs(value))
    return abs(determinant) // factor


@lru_cache(maxsize=128)
def _integral_row_images_cached(
    matrix: tuple[int, ...],
    denominator: int,
) -> tuple[tuple[int, int, int, int], ...]:
    """Enumerate Source-ordered ``x M / denominator`` integral images.

    The ``get_new_fractionals_`` contract retains, in lexicographic order, the
    kernel of ``x M = 0 (mod denominator)``. Smith normal form generates that
    kernel directly; exact sorting then restores Source order. Work scales with
    the determinant plus sorting rather than an exhaustive denominator grid.
    """

    if len(matrix) != 9:
        raise ValueError(f"expected 9 matrix entries, got {len(matrix)}")
    modulus = int(denominator)
    if modulus <= 0:
        raise ValueError(f"denominator must be positive, got {denominator}")
    expected_count = abs(integer_determinant3(matrix))
    if expected_count == 0:
        raise ValueError(f"singular basis matrix: {matrix}")
    if expected_count == 1:
        return ((0, 0, 0, 1),)

    from sympy import Matrix, ZZ
    from sympy.matrices.normalforms import smith_normal_decomp

    transposed = Matrix(3, 3, lambda row, col: matrix[3 * col + row])
    diagonal, _left, right_transform = smith_normal_decomp(transposed, domain=ZZ)
    right = tuple(
        tuple(int(right_transform[row, col]) for col in range(3)) for row in range(3)
    )
    coordinate_values: list[tuple[int, ...]] = []
    for axis in range(3):
        solution_count = gcd(modulus, abs(int(diagonal[axis, axis])))
        step = modulus // solution_count
        coordinate_values.append(tuple(step * index for index in range(solution_count)))

    residues = sorted(
        {
            tuple(
                sum(right[row][col] * transformed[col] for col in range(3)) % modulus
                for row in range(3)
            )
            for transformed in product(*coordinate_values)
        }
    )
    if len(residues) != expected_count:
        raise ValueError(
            f"integer quotient count {len(residues)} != abs(det) "
            f"{expected_count} for {matrix}"
        )

    images: list[tuple[int, int, int, int]] = []
    for x, y, z in residues:
        numerators = tuple(
            x * matrix[col] + y * matrix[3 + col] + z * matrix[6 + col]
            for col in range(3)
        )
        quotients = tuple(divmod(value, modulus) for value in numerators)
        if any(remainder != 0 for _quotient, remainder in quotients):
            raise ValueError("Smith kernel produced a non-integral image")
        images.append(
            (
                quotients[0][0],
                quotients[1][0],
                quotients[2][0],
                1,
            )
        )
    return tuple(images)


def integral_row_images_source_order(
    matrix: Sequence[int],
    denominator: int,
) -> tuple[tuple[int, int, int, int], ...]:
    """Return exact quotient images in the established Source scan order."""

    return _integral_row_images_cached(
        tuple(int(value) for value in matrix),
        int(denominator),
    )
