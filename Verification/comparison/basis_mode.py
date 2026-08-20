"""Source-derived mode-row transport between equivalent child bases."""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from fractions import Fraction
import math

from Verification.comparison.basis_atom import exact_unimodular_matrix3
from Verification.comparison.basis_lattice import (
    child_coordinate_origin_shift,
)


IntMatrix3 = tuple[
    tuple[int, int, int],
    tuple[int, int, int],
    tuple[int, int, int],
]
FloatVector3 = tuple[float, float, float]
FractionVector3 = tuple[Fraction, Fraction, Fraction]


@dataclass(frozen=True)
class BasisEquivalentModeRow:
    """One mode row transported without changing its Source ordinal."""

    source_row_index: int
    xyz: FloatVector3
    dxyz: FloatVector3


@dataclass(frozen=True)
class BasisEquivalentMode:
    """One complete mode column transported and normalized as a unit."""

    source_mode_index: int
    label: str
    normfactor: float | None
    scale: float
    rows: tuple[BasisEquivalentModeRow, ...]


@dataclass(frozen=True)
class BasisEquivalentModeFrameTransport:
    """A complete definition sequence expressed in the request child basis."""

    basis_change: IntMatrix3
    determinant: int
    condition_number: float
    request_origin: FractionVector3
    child_space_group: int
    vector_kind: str
    modes: tuple[BasisEquivalentMode, ...]
    provenance: str = "M:exact_gl3z_fractional_mode_transport"


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


def _float_vector3(value: object) -> FloatVector3 | None:
    if isinstance(value, (str, bytes, bytearray)):
        return None
    try:
        if not isinstance(value, Sequence) or len(value) != 3 or any(isinstance(item, bool) for item in value):
            return None
        out = tuple(float(item) for item in value)
    except (OverflowError, TypeError, ValueError):
        return None
    return out if all(math.isfinite(item) for item in out) else None  # type: ignore[return-value]


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


def _row_multiply(vector: FloatVector3, matrix: IntMatrix3) -> FloatVector3:
    return tuple(
        sum(vector[row] * matrix[row][column] for row in range(3))
        for column in range(3)
    )  # type: ignore[return-value]


def basis_equivalent_mode_frame_transport(
    local_modes: object,
    *,
    basis_change: object,
    request_origin: object,
    local_origin: object,
    request_child_space_group: object,
    local_child_space_group: object,
    vector_kind: str,
    condition_limit: float = 1e10,
    request_basis: object = None,
) -> BasisEquivalentModeFrameTransport | None:
    """Express complete local mode columns in the request child frame.

    Provenance: M.  ``B_local = U B_request`` gives both
    ``x_request = x_local U`` and ``d_request = d_local U``.  Unit-max
    normalization is applied once across every row in a mode column.  The
    Source order, labels, signs, and all nonnumeric metadata remain owned by
    the caller and are never inferred here.
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
    if (
        matrix is None
        or vector_kind not in {"dsp", "mag"}
        or isinstance(local_modes, (str, bytes, bytearray))
        or not isinstance(local_modes, Sequence)
        or not local_modes
    ):
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

    transformed_modes: list[BasisEquivalentMode] = []
    for mode_index, mode in enumerate(local_modes):
        if not isinstance(mode, Mapping):
            return None
        label = mode.get("label")
        rows = mode.get("rows")
        if (
            not isinstance(label, str)
            or not label
            or isinstance(rows, (str, bytes, bytearray))
            or not isinstance(rows, Sequence)
            or not rows
        ):
            return None
        raw_rows: list[tuple[FloatVector3, FloatVector3]] = []
        for row in rows:
            if not isinstance(row, Mapping):
                return None
            xyz = _float_vector3(row.get("xyz"))
            dxyz = _float_vector3(row.get("dxyz"))
            if xyz is None or dxyz is None:
                return None
            rotated_xyz = _row_multiply(xyz, matrix)
            request_xyz = tuple(
                (rotated_xyz[index] + float(origin_shift[index])) % 1.0
                for index in range(3)
            )
            request_dxyz = _row_multiply(dxyz, matrix)
            if not all(math.isfinite(value) for value in (*request_xyz, *request_dxyz)):
                return None
            raw_rows.append((request_xyz, request_dxyz))
        scale = max(abs(value) for _, dxyz in raw_rows for value in dxyz)
        if not math.isfinite(scale) or scale <= 0.0:
            return None
        normfactor_value = mode.get("normfactor")
        if normfactor_value is None:
            normfactor = None
        else:
            if isinstance(normfactor_value, bool):
                return None
            try:
                normfactor = float(normfactor_value) * scale
            except (OverflowError, TypeError, ValueError):
                return None
            if not math.isfinite(normfactor):
                return None
        transformed_modes.append(
            BasisEquivalentMode(
                source_mode_index=mode_index,
                label=label,
                normfactor=normfactor,
                scale=scale,
                rows=tuple(
                    BasisEquivalentModeRow(
                        source_row_index=row_index,
                        xyz=xyz,
                        dxyz=tuple(value / scale for value in dxyz),  # type: ignore[arg-type]
                    )
                    for row_index, (xyz, dxyz) in enumerate(raw_rows)
                ),
            )
        )
    return BasisEquivalentModeFrameTransport(
        basis_change=matrix,
        determinant=determinant,
        condition_number=condition,
        request_origin=request_origin_value,
        child_space_group=request_child,
        vector_kind=vector_kind,
        modes=tuple(transformed_modes),
    )
