"""Exact basis-equivalent lattice transport for validation frames.

The basis matrices in ISODISTORT are row-vector transforms from a parent
cell.  This module keeps the GL(3,Z) decision exact and uses floating point
only for the parent metric and displayed cell parameters.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
import math
from typing import Sequence


FractionMatrix3 = tuple[
    tuple[Fraction, Fraction, Fraction],
    tuple[Fraction, Fraction, Fraction],
    tuple[Fraction, Fraction, Fraction],
]
FloatMatrix3 = tuple[
    tuple[float, float, float],
    tuple[float, float, float],
    tuple[float, float, float],
]
IntMatrix3 = tuple[
    tuple[int, int, int],
    tuple[int, int, int],
    tuple[int, int, int],
]
CellParameters = tuple[float, float, float, float, float, float]
FractionVector3 = tuple[Fraction, Fraction, Fraction]


@dataclass(frozen=True)
class BasisEquivalentLatticeTransport:
    """Validated request/local metrics related by one unimodular basis map."""

    basis_change: IntMatrix3
    request_cell: CellParameters
    local_cell: CellParameters
    request_metric: FloatMatrix3
    local_metric: FloatMatrix3
    max_metric_residual: float
    condition_number: float


def _fraction_matrix(values: object) -> FractionMatrix3 | None:
    if isinstance(values, (str, bytes, bytearray)):
        return None
    try:
        if not isinstance(values, Sequence) or len(values) != 9:
            return None
        flat = tuple(Fraction(str(value)) for value in values)
    except (IndexError, OverflowError, TypeError, ValueError, ZeroDivisionError):
        return None
    return tuple(
        tuple(flat[3 * row + col] for col in range(3))
        for row in range(3)
    )  # type: ignore[return-value]


def _determinant(matrix: Sequence[Sequence[Fraction | float]]) -> Fraction | float:
    return (
        matrix[0][0]
        * (matrix[1][1] * matrix[2][2] - matrix[1][2] * matrix[2][1])
        - matrix[0][1]
        * (matrix[1][0] * matrix[2][2] - matrix[1][2] * matrix[2][0])
        + matrix[0][2]
        * (matrix[1][0] * matrix[2][1] - matrix[1][1] * matrix[2][0])
    )


def _inverse(matrix: FractionMatrix3) -> FractionMatrix3 | None:
    determinant = _determinant(matrix)
    if determinant == 0:
        return None
    inverse = (
        (
            (matrix[1][1] * matrix[2][2] - matrix[1][2] * matrix[2][1]) / determinant,
            (matrix[0][2] * matrix[2][1] - matrix[0][1] * matrix[2][2]) / determinant,
            (matrix[0][1] * matrix[1][2] - matrix[0][2] * matrix[1][1]) / determinant,
        ),
        (
            (matrix[1][2] * matrix[2][0] - matrix[1][0] * matrix[2][2]) / determinant,
            (matrix[0][0] * matrix[2][2] - matrix[0][2] * matrix[2][0]) / determinant,
            (matrix[0][2] * matrix[1][0] - matrix[0][0] * matrix[1][2]) / determinant,
        ),
        (
            (matrix[1][0] * matrix[2][1] - matrix[1][1] * matrix[2][0]) / determinant,
            (matrix[0][1] * matrix[2][0] - matrix[0][0] * matrix[2][1]) / determinant,
            (matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0]) / determinant,
        ),
    )
    return inverse  # type: ignore[return-value]


def _multiply(
    left: Sequence[Sequence[Fraction | float]],
    right: Sequence[Sequence[Fraction | float]],
) -> tuple[tuple[Fraction | float, ...], ...]:
    return tuple(
        tuple(sum(left[row][axis] * right[axis][col] for axis in range(3)) for col in range(3))
        for row in range(3)
    )


def _transpose(
    matrix: Sequence[Sequence[Fraction | float]],
) -> tuple[tuple[Fraction | float, ...], ...]:
    return tuple(tuple(matrix[col][row] for col in range(3)) for row in range(3))


def _condition_number(matrix: FractionMatrix3) -> float | None:
    inverse = _inverse(matrix)
    if inverse is None:
        return None
    try:
        norm = max(sum(abs(float(value)) for value in row) for row in matrix)
        inverse_norm = max(sum(abs(float(value)) for value in row) for row in inverse)
    except (OverflowError, ValueError):
        return None
    condition = norm * inverse_norm
    return condition if math.isfinite(condition) else None


def unimodular_basis_change(
    source_values: Sequence[object],
    target_values: Sequence[object],
) -> IntMatrix3 | None:
    """Return exact ``U`` for ``target = U source`` when ``U`` is unimodular."""

    source = _fraction_matrix(source_values)
    target = _fraction_matrix(target_values)
    if source is None or target is None:
        return None
    source_inverse = _inverse(source)
    if source_inverse is None:
        return None
    change = _multiply(target, source_inverse)
    if any(not isinstance(value, Fraction) or value.denominator != 1 for row in change for value in row):
        return None
    integral = tuple(tuple(int(value) for value in row) for row in change)
    if abs(int(_determinant(integral))) != 1:
        return None
    return integral  # type: ignore[return-value]


def child_coordinate_origin_shift(
    request_basis_values: Sequence[object],
    request_origin: Sequence[object],
    local_origin: Sequence[object],
) -> FractionVector3 | None:
    """Return the exact request-child shift caused by a parent-origin change."""

    request_basis = _fraction_matrix(request_basis_values)
    request_inverse = _inverse(request_basis) if request_basis is not None else None
    try:
        if (
            request_inverse is None
            or len(request_origin) != 3
            or len(local_origin) != 3
        ):
            return None
        request = tuple(Fraction(str(value)) for value in request_origin)
        local = tuple(Fraction(str(value)) for value in local_origin)
    except (OverflowError, TypeError, ValueError, ZeroDivisionError):
        return None
    delta = tuple(local[index] - request[index] for index in range(3))
    return tuple(
        sum(delta[row] * request_inverse[row][column] for row in range(3))
        for column in range(3)
    )  # type: ignore[return-value]


def _metric_from_cell(cell: CellParameters) -> FloatMatrix3 | None:
    a, b, c, alpha, beta, gamma = cell
    if (
        not all(math.isfinite(value) for value in cell)
        or min(a, b, c) <= 0.0
        or not all(0.0 < angle < 180.0 for angle in (alpha, beta, gamma))
    ):
        return None
    ca = math.cos(math.radians(alpha))
    cb = math.cos(math.radians(beta))
    cg = math.cos(math.radians(gamma))
    metric = (
        (a * a, a * b * cg, a * c * cb),
        (a * b * cg, b * b, b * c * ca),
        (a * c * cb, b * c * ca, c * c),
    )
    leading_two = metric[0][0] * metric[1][1] - metric[0][1] * metric[1][0]
    determinant = float(_determinant(metric))
    if metric[0][0] <= 0.0 or leading_two <= 0.0 or determinant <= 0.0:
        return None
    return metric


def _basis_metric(
    basis: FractionMatrix3,
    parent_metric: FloatMatrix3,
) -> FloatMatrix3 | None:
    try:
        basis_float = tuple(tuple(float(value) for value in row) for row in basis)
    except (OverflowError, ValueError):
        return None
    product = _multiply(_multiply(basis_float, parent_metric), _transpose(basis_float))
    try:
        result = tuple(tuple(float(value) for value in row) for row in product)
    except (OverflowError, ValueError):
        return None
    if any(not math.isfinite(value) for row in result for value in row):
        return None
    return result  # type: ignore[return-value]


def _cell_from_metric(metric: FloatMatrix3) -> CellParameters | None:
    diagonal = tuple(metric[index][index] for index in range(3))
    if any(not math.isfinite(value) or value <= 0.0 for value in diagonal):
        return None
    lengths = tuple(math.sqrt(value) for value in diagonal)

    def angle(left: int, right: int) -> float:
        denominator = lengths[left] * lengths[right]
        cosine = metric[left][right] / denominator
        if not math.isfinite(cosine):
            return math.nan
        return math.degrees(math.acos(max(-1.0, min(1.0, cosine))))

    cell = (
        lengths[0],
        lengths[1],
        lengths[2],
        angle(1, 2),
        angle(0, 2),
        angle(0, 1),
    )
    return cell if all(math.isfinite(value) for value in cell) else None


def basis_equivalent_lattice_transport(
    parent_cell: Sequence[float],
    request_basis_values: Sequence[object],
    local_basis_values: Sequence[object],
    *,
    residual_tolerance: float = 1e-10,
    condition_limit: float = 1e10,
) -> BasisEquivalentLatticeTransport | None:
    """Transport a child metric from a local basis to an equivalent request basis.

    This function proves only the lattice relation.  Callers remain responsible
    for exact child-group, origin, slot, and nonbasis subgroup identity gates.
    """

    if (
        len(parent_cell) != 6
        or not math.isfinite(residual_tolerance)
        or residual_tolerance <= 0.0
        or not math.isfinite(condition_limit)
        or condition_limit < 1.0
    ):
        return None
    try:
        cell = tuple(float(value) for value in parent_cell)
    except (TypeError, ValueError, OverflowError):
        return None
    parent_metric = _metric_from_cell(cell)  # type: ignore[arg-type]
    request_basis = _fraction_matrix(request_basis_values)
    local_basis = _fraction_matrix(local_basis_values)
    if parent_metric is None or request_basis is None or local_basis is None:
        return None

    basis_change = unimodular_basis_change(request_basis_values, local_basis_values)
    if basis_change is None:
        return None
    change_fraction = tuple(tuple(Fraction(value) for value in row) for row in basis_change)
    change_inverse = _inverse(change_fraction)  # type: ignore[arg-type]
    if change_inverse is None or _multiply(change_fraction, request_basis) != local_basis:
        return None

    conditions = tuple(
        _condition_number(matrix)
        for matrix in (request_basis, local_basis, change_fraction)
    )
    if any(value is None for value in conditions):
        return None
    condition = max(float(value) for value in conditions if value is not None)
    if condition > condition_limit:
        return None

    request_metric = _basis_metric(request_basis, parent_metric)
    local_metric = _basis_metric(local_basis, parent_metric)
    if request_metric is None or local_metric is None:
        return None
    try:
        inverse_float = tuple(tuple(float(value) for value in row) for row in change_inverse)
    except (OverflowError, ValueError):
        return None
    transported = _multiply(_multiply(inverse_float, local_metric), _transpose(inverse_float))
    residual = max(
        abs(float(transported[row][col]) - request_metric[row][col])
        for row in range(3)
        for col in range(3)
    )
    scale = max(1.0, *(abs(value) for row in request_metric for value in row))
    if not math.isfinite(residual) or residual > residual_tolerance * scale:
        return None

    request_cell = _cell_from_metric(request_metric)
    local_cell = _cell_from_metric(local_metric)
    if request_cell is None or local_cell is None:
        return None
    return BasisEquivalentLatticeTransport(
        basis_change=basis_change,
        request_cell=request_cell,
        local_cell=local_cell,
        request_metric=request_metric,
        local_metric=local_metric,
        max_metric_residual=residual,
        condition_number=condition,
    )
