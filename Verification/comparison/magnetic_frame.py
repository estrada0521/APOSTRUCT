"""Exact coordinate bridge for magnetic atom comparison."""

from __future__ import annotations

from fractions import Fraction
from typing import Iterable, Sequence

from sympy import Matrix

from APOSTRUCT.Backend.source.tables import source_tables


FractionPoint = tuple[Fraction, Fraction, Fraction]
FractionMatrix = tuple[
    tuple[Fraction, Fraction, Fraction],
    tuple[Fraction, Fraction, Fraction],
    tuple[Fraction, Fraction, Fraction],
]


def _matrix(values: Iterable[object]) -> FractionMatrix:
    flat = tuple(Fraction(str(value)) for value in values)
    if len(flat) != 9:
        raise ValueError(f"expected 3x3 matrix, got {len(flat)} values")
    return tuple(
        tuple(flat[row * 3 + col] for col in range(3)) for row in range(3)
    )  # type: ignore[return-value]


def _inverse(matrix: FractionMatrix) -> FractionMatrix:
    inverse = Matrix(matrix).inv()
    return tuple(
        tuple(Fraction(inverse[row, col]) for col in range(3))
        for row in range(3)
    )  # type: ignore[return-value]


def _row(point: Sequence[Fraction], matrix: FractionMatrix) -> FractionPoint:
    return tuple(
        sum(Fraction(point[row]) * matrix[row][col] for row in range(3))
        for col in range(3)
    )  # type: ignore[return-value]


def _add(left: Sequence[Fraction], right: Sequence[Fraction]) -> FractionPoint:
    return tuple(
        Fraction(left[index]) + Fraction(right[index]) for index in range(3)
    )  # type: ignore[return-value]


def _sub(left: Sequence[Fraction], right: Sequence[Fraction]) -> FractionPoint:
    return tuple(
        Fraction(left[index]) - Fraction(right[index]) for index in range(3)
    )  # type: ignore[return-value]


def _origin(value: object) -> FractionPoint:
    if isinstance(value, str):
        parts = [part.strip() for part in value.strip().strip("()").split(",")]
        if len(parts) != 3:
            raise ValueError(f"bad origin: {value!r}")
        return tuple(Fraction(part) for part in parts)  # type: ignore[return-value]
    if not isinstance(value, (list, tuple)):
        raise TypeError(f"bad origin: {value!r}")
    if len(value) == 4:
        denominator = int(value[3])
        if denominator == 0:
            raise ZeroDivisionError("zero origin denominator")
        return tuple(
            Fraction(int(value[index]), denominator) for index in range(3)
        )  # type: ignore[return-value]
    if len(value) == 3:
        return tuple(Fraction(str(item)) for item in value)  # type: ignore[return-value]
    raise ValueError(f"bad origin length: {len(value)}")


def presentation_to_standard_child_cinter(
    xyz: Sequence[float | Fraction],
    *,
    parent_sg: int,
    child_sg: int,
    presentation_basis: Iterable[object],
    presentation_origin: object,
    source_basis: Iterable[object],
    source_origin: object,
    parent_setting_id: int | None = None,
) -> tuple[float, float, float]:
    """Map a displayed child point to the selected subgroup standard frame."""

    data = source_tables()
    displayed = tuple(Fraction(str(value)) for value in xyz)
    display_matrix = _matrix(presentation_basis)
    source_matrix = _matrix(source_basis)
    parent_cinter = _add(_row(displayed, display_matrix), _origin(presentation_origin))

    if parent_setting_id is not None:
        default_id = int(data.default_inter_setting_record(int(parent_sg))["id"])
        if int(parent_setting_id) != default_id:
            default_matrix = data.cml_to_cinter_matrix(int(parent_sg), default_id)
            selected_matrix = data.cml_to_cinter_matrix(
                int(parent_sg), int(parent_setting_id)
            )
            transform = tuple(
                _row(row, selected_matrix) for row in _inverse(default_matrix)
            )
            default_origin = data.cml_to_cinter_origin(int(parent_sg), default_id)
            selected_origin = data.cml_to_cinter_origin(
                int(parent_sg), int(parent_setting_id)
            )
            shift = tuple(
                selected_origin[col]
                - sum(
                    default_origin[row] * transform[row][col]
                    for row in range(3)
                )
                for col in range(3)
            )
            parent_cinter = _row(_sub(parent_cinter, shift), _inverse(transform))

    parent_pml_to_cinter = data.pml_to_cinter_matrix(int(parent_sg))
    parent_setting_origin = data.cml_to_cinter_origin(int(parent_sg))
    parent_pml = _row(
        _sub(parent_cinter, parent_setting_origin),
        _inverse(parent_pml_to_cinter),
    )
    child_pml = _row(
        _sub(parent_pml, _origin(source_origin)),
        _inverse(source_matrix),
    )
    child_pml_to_cinter = data.pml_to_cinter_matrix(int(child_sg))
    child_setting_origin = data.cml_to_cinter_origin(int(child_sg))
    child_cinter = _add(_row(child_pml, child_pml_to_cinter), child_setting_origin)
    return tuple(float(value % 1) for value in child_cinter)  # type: ignore[return-value]


def selected_magnetic_frame(isotropy: dict[str, object]) -> dict[str, object] | None:
    """Return a validated numeric frame payload or fail closed."""

    subgroup = isotropy.get("magnetic_subgroup") or isotropy.get("subgroup")
    if not isinstance(subgroup, dict) or subgroup.get("ordinary_number") is None:
        return None
    basis = (
        isotropy.get("presentation_basis")
        or isotropy.get("display_basis")
        or isotropy.get("basis")
    )
    source_basis = isotropy.get("source_basis_values")
    source_origin = isotropy.get("source_origin_values")
    if not (
        isinstance(basis, list)
        and len(basis) == 3
        and all(isinstance(row, list) and len(row) == 3 for row in basis)
        and isinstance(source_basis, list)
        and len(source_basis) == 9
        and isinstance(source_origin, list)
        and len(source_origin) == 4
    ):
        return None
    try:
        return {
            "magnetic_group": int(subgroup["number"]),
            "child_sg": int(subgroup["ordinary_number"]),
            "presentation_basis": tuple(value for row in basis for value in row),
            "presentation_origin": isotropy.get("origin")
            or isotropy.get("origin_values")
            or (0, 0, 0),
            "source_basis": tuple(source_basis),
            "source_origin": tuple(source_origin),
            "parent_setting_id": (
                int(isotropy["parent_inter_setting_id"])
                if isotropy.get("parent_inter_setting_id") is not None
                else None
            ),
        }
    except (KeyError, TypeError, ValueError):
        return None


__all__ = ["presentation_to_standard_child_cinter", "selected_magnetic_frame"]
