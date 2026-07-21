"""Mode-site transport helpers.

Extracted mechanically from the former monolithic runtime.
"""

from __future__ import annotations

from fractions import Fraction
from itertools import permutations
from numbers import Integral
from typing import Any
import gemmi
from ISODISTORT.Assembled.Backend.exactmath import (
    fraction_matrix_inverse3,
    fraction_matrix_multiply3,
)
from ISODISTORT.Assembled.Backend.source.tables import source_tables
from ISODISTORT.Assembled.Backend.modes.engine.decoder import ModeDataDecoder
from ISODISTORT.Assembled.Backend.modes.engine.subgroup_structure.presentation_transport import (
    representative_operation_record_from_cinter,
)
from ISODISTORT.Assembled.Backend.modes.common import (
    _fraction_matmul,
    _fraction_matrix_inverse_3,
    _fraction_row_multiply,
    _site_params,
)

def _site_representative_operation_record(
    decoder: ModeDataDecoder,
    sg: int,
    site: dict[str, Any],
) -> tuple[int, int, int, int, int] | None:
    """Return Stage1's exact source-to-canonical operation in parent PML."""

    text = site.get("representative_operation")
    lattice = site.get("representative_lattice_translation")
    if not isinstance(text, str) or not isinstance(lattice, (list, tuple)) or len(lattice) != 3:
        return None
    try:
        operation = gemmi.Op(text)
        denominator = int(operation.DEN)
        rotation = [
            [Fraction(int(operation.rot[col][row]), denominator) for col in range(3)]
            for row in range(3)
        ]
        translation = [
            Fraction(int(operation.tran[axis]), denominator) - int(lattice[axis])
            for axis in range(3)
        ]
        return representative_operation_record_from_cinter(
            decoder,
            int(sg),
            rotation,
            translation,
        )
    except (KeyError, TypeError, ValueError, ZeroDivisionError):
        return None



def _operation_row_affine(
    text: str,
    lattice_translation: list[int] | tuple[int, ...],
) -> tuple[
    tuple[tuple[Fraction, Fraction, Fraction], ...],
    tuple[Fraction, Fraction, Fraction],
] | None:
    try:
        operation = gemmi.Op(text)
        denominator = int(operation.DEN)
        rotation = tuple(
            tuple(Fraction(int(operation.rot[col][row]), denominator) for col in range(3))
            for row in range(3)
        )
        translation = tuple(
            Fraction(int(operation.tran[axis]), denominator) - int(lattice_translation[axis])
            for axis in range(3)
        )
        return rotation, translation  # type: ignore[return-value]
    except (IndexError, RuntimeError, TypeError, ValueError, ZeroDivisionError):
        return None



def _presentation_site_relative_operation_record(
    decoder: ModeDataDecoder,
    sg: int,
    site: dict[str, Any],
    symmetry_operations: list[str],
    selected_setting_id: int | None,
) -> tuple[int, int, int, int, int] | None:
    """Bridge the CIF-visible Wyckoff representative to Source's canonical one."""

    try:
        row_id = int(site["wyckoff_row_id"])
        source_point = tuple(float(value) for value in site["fract"])
        source_affine = _operation_row_affine(
            str(site["representative_operation"]),
            tuple(int(value) for value in site["representative_lattice_translation"]),
        )
    except (KeyError, TypeError, ValueError):
        return None
    if source_affine is None:
        return None
    data = source_tables()
    default_id = int(data.default_inter_setting_record(int(sg))["id"])
    setting_id = default_id if selected_setting_id is None else int(selected_setting_id)
    try:
        row = next(item for item in data.wyckoff_rows(int(sg)) if int(item.row_id) == row_id)
    except StopIteration:
        return None
    vectors = data._inter_wyckoff_fraction_vectors(int(sg), row, setting_id)
    setting_bridge = _parent_setting_bridge(int(sg), selected_setting_id)
    parameter_count = sum(1 for vector in vectors[1:] if any(value != 0 for value in vector))
    if setting_bridge is None or parameter_count == 0:
        presentation_affine = source_affine
    else:
        candidates: list[tuple[float, int, str, tuple[int, int, int]]] = []
        for operation_index, text in enumerate(symmetry_operations):
            try:
                applied = tuple(float(value) for value in gemmi.Op(str(text)).apply_to_xyz(source_point))
            except (RuntimeError, TypeError, ValueError):
                continue
            solved = data._solve_wyckoff_params_from_vectors(vectors, applied, tol=2e-5)
            params = _site_params({"wyckoff_params": solved}) if solved is not None else None
            if params is None or len(params) != parameter_count:
                continue
            score = min(
                sum(
                    abs(((float(params[index]) - float(source_point[axis]) + 0.5) % 1.0) - 0.5)
                    for index, axis in enumerate(axes)
                )
                for axes in permutations(range(3), parameter_count)
            )
            target = tuple(
                float(vectors[0][axis])
                + sum(float(params[index]) * float(vectors[index + 1][axis]) for index in range(parameter_count))
                for axis in range(3)
            )
            lattice = tuple(int(round(applied[axis] - target[axis])) for axis in range(3))
            if not all(abs(applied[axis] - target[axis] - lattice[axis]) <= 2e-5 for axis in range(3)):
                continue
            candidates.append((float(score), int(operation_index), str(text), lattice))
        if not candidates:
            return None
        _score, _index, text, lattice = min(candidates, key=lambda item: (item[0], item[1]))
        presentation_affine = _operation_row_affine(text, lattice)
    if presentation_affine is None:
        return None
    presentation_rotation, presentation_translation = presentation_affine
    source_rotation, source_translation = source_affine
    presentation_inverse = _fraction_matrix_inverse_3(
        [[float(value) for value in row] for row in presentation_rotation]
    )
    if presentation_inverse is None:
        return None
    inverse_translation = tuple(
        -value
        for value in _fraction_row_multiply(
            presentation_translation,
            presentation_inverse,
        )
    )
    relative_rotation = _fraction_matmul(presentation_inverse, source_rotation)
    relative_translation = tuple(
        _fraction_row_multiply(inverse_translation, source_rotation)[axis]
        + source_translation[axis]
        for axis in range(3)
    )
    if setting_bridge is not None:
        bridge_matrix, bridge_shift = setting_bridge
        bridge_inverse = fraction_matrix_inverse3(bridge_matrix)
        selected_rotation = relative_rotation
        relative_rotation = fraction_matrix_multiply3(
            fraction_matrix_multiply3(bridge_matrix, selected_rotation),
            bridge_inverse,
        )
        relative_translation = _fraction_row_multiply(
            tuple(
                _fraction_row_multiply(bridge_shift, selected_rotation)[axis]
                + relative_translation[axis]
                - bridge_shift[axis]
                for axis in range(3)
            ),
            bridge_inverse,
        )
    try:
        return representative_operation_record_from_cinter(
            decoder,
            int(sg),
            relative_rotation,
            relative_translation,
        )
    except (KeyError, TypeError, ValueError):
        return None



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


def _parent_point_from_default(
    point: Any,
    bridge: Any,
) -> tuple[Fraction, Fraction, Fraction]:
    values = _parent_frame_point(point)
    if bridge is None:
        return values
    transform, shift = bridge
    transformed = _fraction_row_multiply(values, transform)
    return tuple(transformed[axis] + shift[axis] for axis in range(3))  # type: ignore[return-value]


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
        source_point = tuple(float(value) for value in site["fract"])
    except (KeyError, TypeError, ValueError):
        return fallback
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


def _source_default_site_params(
    sg: int,
    site: dict[str, Any],
    selected_setting_id: int | None,
    symmetry_operations: list[str],
) -> tuple[float, ...] | None:
    """Return Source-default Wyckoff parameters in positional kernel order."""

    params = _source_default_wyckoff_params(
        int(sg),
        site,
        selected_setting_id,
        symmetry_operations,
    )
    return _site_params({"wyckoff_params": params}) if params else None
