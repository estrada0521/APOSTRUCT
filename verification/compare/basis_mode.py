"""Source-derived mode-row transport between equivalent child bases."""

from __future__ import annotations

from collections.abc import Sequence
from fractions import Fraction
import math

from verification.compare.basis_atom import exact_unimodular_matrix3
from verification.compare.basis_lattice import (
    child_coordinate_origin_shift,
)


IntMatrix3 = tuple[
    tuple[int, int, int],
    tuple[int, int, int],
    tuple[int, int, int],
]
FractionVector3 = tuple[Fraction, Fraction, Fraction]


def _fraction_vector3(value: object) -> FractionVector3 | None:
    if isinstance(value, (str, bytes, bytearray)):
        return None
    try:
        if not isinstance(value, Sequence) or len(value) != 3:
            return None
        out: list[Fraction] = []
        for item in value:
            if isinstance(item, bool):
                return None
            if isinstance(item, float) and not math.isfinite(item):
                return None
            out.append(item if isinstance(item, Fraction) else Fraction(str(item)))
    except (OverflowError, TypeError, ValueError, ZeroDivisionError):
        return None
    return tuple(out)  # type: ignore[return-value]


def _exact_integer(value: object) -> int | None:
    if isinstance(value, bool):
        return None
    try:
        fraction = value if isinstance(value, Fraction) else Fraction(str(value))
    except (OverflowError, TypeError, ValueError, ZeroDivisionError):
        return None
    return int(fraction) if fraction.denominator == 1 else None


def _determinant(matrix: IntMatrix3) -> int:
    return (
        matrix[0][0] * (matrix[1][1] * matrix[2][2] - matrix[1][2] * matrix[2][1])
        - matrix[0][1] * (matrix[1][0] * matrix[2][2] - matrix[1][2] * matrix[2][0])
        + matrix[0][2] * (matrix[1][0] * matrix[2][1] - matrix[1][1] * matrix[2][0])
    )


def _inverse_unimodular(matrix: IntMatrix3, determinant: int) -> IntMatrix3:
    return (
        (
            (matrix[1][1] * matrix[2][2] - matrix[1][2] * matrix[2][1]) // determinant,
            (matrix[0][2] * matrix[2][1] - matrix[0][1] * matrix[2][2]) // determinant,
            (matrix[0][1] * matrix[1][2] - matrix[0][2] * matrix[1][1]) // determinant,
        ),
        (
            (matrix[1][2] * matrix[2][0] - matrix[1][0] * matrix[2][2]) // determinant,
            (matrix[0][0] * matrix[2][2] - matrix[0][2] * matrix[2][0]) // determinant,
            (matrix[0][2] * matrix[1][0] - matrix[0][0] * matrix[1][2]) // determinant,
        ),
        (
            (matrix[1][0] * matrix[2][1] - matrix[1][1] * matrix[2][0]) // determinant,
            (matrix[0][1] * matrix[2][0] - matrix[0][0] * matrix[2][1]) // determinant,
            (matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0]) // determinant,
        ),
    )


def _condition_upper_bound(matrix: IntMatrix3, determinant: int) -> float | None:
    """Return a conservative Frobenius-norm condition bound."""

    inverse = _inverse_unimodular(matrix, determinant)
    try:
        left = math.sqrt(sum(value * value for row in matrix for value in row))
        right = math.sqrt(sum(value * value for row in inverse for value in row))
        condition = left * right
    except (OverflowError, ValueError):
        return None
    return condition if math.isfinite(condition) else None


def basis_equivalent_mode_frame(
    *,
    basis_change: object,
    request_origin: object,
    local_origin: object,
    request_child_space_group: object,
    local_child_space_group: object,
    condition_limit: float = 1e10,
    request_basis: object = None,
) -> tuple[IntMatrix3, FractionVector3] | None:
    """The frame that expresses local mode columns in the request child frame.

    Provenance: M.  ``B_local = U B_request`` gives both
    ``x_request = x_local U + shift`` and ``d_request = d_local U``; the
    columns themselves are carried and renormalized by the native comparison.
    """

    matrix = exact_unimodular_matrix3(basis_change)
    request_origin_value = _fraction_vector3(request_origin)
    local_origin_value = _fraction_vector3(local_origin)
    if request_origin_value is None or local_origin_value is None:
        return None
    if request_origin_value == local_origin_value:
        origin_shift: FractionVector3 = (
            Fraction(0),
            Fraction(0),
            Fraction(0),
        )
    elif isinstance(request_basis, Sequence) and not isinstance(
        request_basis, (str, bytes, bytearray)
    ):
        origin_shift = child_coordinate_origin_shift(
            request_basis,
            request_origin_value,
            local_origin_value,
        )
        if origin_shift is None:
            return None
    else:
        return None
    if matrix is None:
        return None
    determinant = _determinant(matrix)
    # The current direct binary/Web anchors all preserve orientation.  Keep
    # orientation-reversing magnetic and displacive print conventions closed
    # until an equally direct anchor fixes their presentation sign.
    if determinant != 1:
        return None
    condition = _condition_upper_bound(matrix, determinant)
    request_child = _exact_integer(request_child_space_group)
    local_child = _exact_integer(local_child_space_group)
    try:
        if isinstance(condition_limit, bool):
            return None
        limit = float(condition_limit)
    except (OverflowError, TypeError, ValueError):
        return None
    if (
        condition is None
        or not math.isfinite(limit)
        or limit < 1.0
        or condition > limit
        or request_child is None
        or local_child is None
        or request_child <= 0
        or request_child != local_child
    ):
        return None
    return matrix, origin_shift
