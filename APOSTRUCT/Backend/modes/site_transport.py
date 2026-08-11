"""Mode-site transport helpers."""

from __future__ import annotations

from fractions import Fraction
from numbers import Integral
from typing import Any
import gemmi
from APOSTRUCT.Backend.cif_numbers import float_cif_number
from APOSTRUCT.Backend.exactmath import (
    fraction_matrix_inverse3,
    fraction_matrix_multiply3,
)
from APOSTRUCT.Backend.source.tables import source_tables
from APOSTRUCT.Backend.modes.common import (
    _fraction_row_multiply,
)



def _exact_integer(value: Any, *, field: str) -> int:
    if isinstance(value, bool) or not isinstance(value, Integral):
        raise TypeError(f"{field} must be an exact integer")
    return int(value)


def _parent_setting_bridge(
    sg: int,
    selected_setting_id: int | None,
) -> tuple[
    tuple[tuple[Fraction, Fraction, Fraction], ...],
    tuple[Fraction, Fraction, Fraction],
] | None:
    """Return the Source-default -> selected-parent affine setting bridge."""

    if selected_setting_id is None:
        return None
    sg = _exact_integer(sg, field="parent space group")
    selected_setting_id = _exact_integer(selected_setting_id, field="parent setting")
    data = source_tables()
    if int(data.inter_setting_record(selected_setting_id)["space_group"]) != sg:
        raise ValueError("parent setting does not belong to the parent space group")
    default_id = int(data.default_inter_setting_record(sg)["id"])
    if selected_setting_id == default_id:
        return None
    default_matrix = data.cml_to_cinter_matrix(sg, default_id)
    selected_matrix = data.cml_to_cinter_matrix(sg, selected_setting_id)
    transform = fraction_matrix_multiply3(
        fraction_matrix_inverse3(default_matrix),
        selected_matrix,
    )
    default_origin = data.cml_to_cinter_origin(sg, default_id)
    selected_origin = data.cml_to_cinter_origin(sg, selected_setting_id)
    shift = tuple(
        selected_origin[col]
        - sum(default_origin[row] * transform[row][col] for row in range(3))
        for col in range(3)
    )
    return transform, shift  # type: ignore[return-value]


def _parent_frame_point(point: Any) -> tuple[Fraction, Fraction, Fraction]:
    values = tuple(Fraction(str(value)) for value in point)
    if len(values) != 3:
        raise ValueError("parent point must have three coordinates")
    return values  # type: ignore[return-value]


def _parent_point_to_default(
    point: Any,
    bridge: Any,
) -> tuple[Fraction, Fraction, Fraction]:
    values = _parent_frame_point(point)
    if bridge is None:
        return values
    transform, shift = bridge
    return _fraction_row_multiply(
        tuple(values[axis] - shift[axis] for axis in range(3)),
        fraction_matrix_inverse3(transform),
    )


def _source_default_wyckoff_params(
    sg: int,
    site: dict[str, Any],
    selected_setting_id: int | None,
    symmetry_operations: list[str],
) -> dict[str, float] | None:
    """Express selected-setting Wyckoff parameters in Source's default setting."""

    raw_params = site.get("wyckoff_params")
    fallback = dict(raw_params) if isinstance(raw_params, dict) and raw_params else None
    bridge = _parent_setting_bridge(int(sg), selected_setting_id)
    if bridge is None:
        return fallback
    try:
        row_id = int(site["wyckoff_row_id"])
    except (KeyError, TypeError, ValueError):
        return fallback
    source_values = tuple(float_cif_number(value) for value in site.get("fract") or ())
    if len(source_values) != 3 or any(value is None for value in source_values):
        return fallback
    source_point = tuple(float(value) for value in source_values)
    data = source_tables()
    default_id = int(data.default_inter_setting_record(int(sg))["id"])
    try:
        row = next(item for item in data.wyckoff_rows(int(sg)) if int(item.row_id) == row_id)
    except StopIteration:
        return fallback
    vectors = data._inter_wyckoff_fraction_vectors(int(sg), row, default_id)
    candidates = [source_point]
    for triplet in symmetry_operations:
        try:
            candidates.append(tuple(float(value) for value in gemmi.Op(str(triplet)).apply_to_xyz(source_point)))
        except (RuntimeError, TypeError, ValueError):
            continue
    for candidate in candidates:
        default_point = tuple(
            float(value % 1) for value in _parent_point_to_default(candidate, bridge)
        )
        solved = data._solve_wyckoff_params_from_vectors(
            vectors, default_point, tol=1e-5
        )
        if solved is not None:
            return {str(key): float(value) for key, value in solved.items()}
    return fallback
