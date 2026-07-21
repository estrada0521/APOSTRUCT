"""Source-only magnetic orbit bridge for selected presentation grids."""

from __future__ import annotations

from fractions import Fraction
import math
from typing import Any, Callable, Sequence

from ISODISTORT.Assembled.Backend.exactmath import (
    fraction_matrix_inverse3 as _matrix_inverse,
    fraction_row_multiply3 as _row_multiply,
)
from ISODISTORT.Assembled.Backend.source.magnetic import data as magnetic_data
from ISODISTORT.Assembled.Backend.source.tables import source_tables
from ISODISTORT.Assembled.Backend.modes.presentation import (
    _bucket_candidates,
    _bucket_width,
    _periodic_bucket_key,
)
from ISODISTORT.Assembled.Backend.modes.site_transport import (
    _parent_point_from_default,
    _parent_point_to_default,
    _parent_setting_bridge,
)
from ISODISTORT.Assembled.Backend.modes.structure.magnetic_wyckoff import (
    _input_fraction,
    _pml_point_operation_matrix,
    identify_magnetic_wyckoff_branch,
    magnetic_orbit_points,
)


FractionPoint = tuple[Fraction, Fraction, Fraction]
FractionMatrix = tuple[FractionPoint, FractionPoint, FractionPoint]
PeriodicPointIndex = tuple[
    float,
    int,
    dict[tuple[int, int, int], tuple[int, ...]],
]


def _periodic_point_index(
    points: Sequence[Sequence[Any]],
    tolerance: float,
) -> PeriodicPointIndex | None:
    """Index periodic points without becoming a matching authority."""

    width = _bucket_width(float(tolerance))
    if width is None:
        return None
    bucket_count = max(1, int(math.floor(1.0 / width)))
    mutable: dict[tuple[int, int, int], list[int]] = {}
    for index, point in enumerate(points):
        try:
            key = _periodic_bucket_key(point, width, bucket_count)
        except (OverflowError, TypeError, ValueError):
            return None
        if key is None:
            return None
        mutable.setdefault(key, []).append(index)
    return (
        width,
        bucket_count,
        {key: tuple(indices) for key, indices in mutable.items()},
    )


def _periodic_match_indices(
    point: Sequence[Any],
    candidates: Sequence[Sequence[Any]],
    tolerance: float,
    predicate: Callable[[Any, Any], bool],
    point_index: PeriodicPointIndex | None,
) -> list[int]:
    """Return the same ascending matches as a full periodic scan."""

    narrowed: Sequence[int] = range(len(candidates))
    if point_index is not None:
        width, bucket_count, buckets = point_index
        try:
            key = _periodic_bucket_key(point, width, bucket_count)
        except (OverflowError, TypeError, ValueError):
            key = None
        if key is not None:
            narrowed = _bucket_candidates(
                buckets,
                key,
                periodic_count=bucket_count,
            )
    return [
        index
        for index in narrowed
        if predicate(point, candidates[index])
    ]


def _periodic_match_owner_indices(
    point: Sequence[Any],
    candidates: Sequence[Sequence[Any]],
    owners: Sequence[int],
    tolerance: float,
    predicate: Callable[[Any, Any], bool],
    point_index: PeriodicPointIndex | None,
) -> tuple[int, ...]:
    """Return candidate owners in Source order without becoming authority."""

    if len(candidates) != len(owners):
        return ()
    return tuple(
        sorted(
            {
                owners[index]
                for index in _periodic_match_indices(
                    point,
                    candidates,
                    tolerance,
                    predicate,
                    point_index,
                )
            }
        )
    )


def _fraction_matrix(values: Sequence[Any]) -> FractionMatrix:
    rows = tuple(values)
    if len(rows) == 3 and all(isinstance(row, Sequence) and len(row) == 3 for row in rows):
        return tuple(
            tuple(_input_fraction(value) for value in row)
            for row in rows
        )  # type: ignore[return-value]
    if len(rows) != 9:
        raise ValueError(f"expected 3x3 basis, got {values!r}")
    return tuple(
        tuple(_input_fraction(rows[3 * row + col]) for col in range(3))
        for row in range(3)
    )  # type: ignore[return-value]


def _fraction_origin(values: Sequence[Any] | str) -> FractionPoint:
    if isinstance(values, str):
        text = values.strip()
        if text.startswith("(") and text.endswith(")"):
            text = text[1:-1]
        parts = tuple(part.strip() for part in text.split(","))
        if len(parts) != 3:
            raise ValueError(f"expected origin triplet, got {values!r}")
        return tuple(Fraction(part) for part in parts)  # type: ignore[return-value]
    raw = tuple(values)
    if len(raw) == 4:
        denominator = int(raw[3])
        if denominator == 0:
            raise ValueError("zero origin denominator")
        return tuple(Fraction(int(raw[axis]), denominator) for axis in range(3))  # type: ignore[return-value]
    if len(raw) != 3:
        raise ValueError(f"expected origin triplet or record, got {values!r}")
    return tuple(_input_fraction(value) for value in raw)  # type: ignore[return-value]


def _fold(point: Sequence[Fraction]) -> FractionPoint:
    return tuple(Fraction(value) % 1 for value in point)  # type: ignore[return-value]


def _periodic_close(left: Sequence[Fraction], right: Sequence[Fraction], tolerance: float) -> bool:
    return all(
        abs(float(Fraction(left[axis]) - Fraction(right[axis]))
            - round(float(Fraction(left[axis]) - Fraction(right[axis])))) <= tolerance
        for axis in range(3)
    )


def _correspondence_value(item: Any, name: str) -> Any:
    if isinstance(item, dict):
        return item.get(name)
    return getattr(item, name, None)


def _operation_record(item: Any, name: str) -> tuple[int, int, int, int, int]:
    raw = _correspondence_value(item, name)
    if not isinstance(raw, Sequence) or len(raw) != 5:
        raise ValueError(f"missing {name} in magnetic correspondence")
    record = tuple(int(value) for value in raw)
    if record[3] == 0:
        raise ValueError(f"zero denominator in {name}")
    return record  # type: ignore[return-value]


def _apply_pml_operation(
    sg: int,
    record: tuple[int, int, int, int, int],
    point: Sequence[Fraction],
) -> FractionPoint:
    table = magnetic_data().table
    ordinary_point_op = int(table["mag_point_op_mag2nonmag"][int(record[4]) - 1])
    rotation = _pml_point_operation_matrix(int(sg), ordinary_point_op)
    rotated = _row_multiply(tuple(Fraction(value) for value in point), rotation)
    translation = tuple(Fraction(int(record[axis]), int(record[3])) for axis in range(3))
    return tuple(rotated[axis] + translation[axis] for axis in range(3))  # type: ignore[return-value]


def selected_magnetic_correspondence_standard_points(
    *,
    parent_sg: int,
    child_sg: int,
    presentation_points: Sequence[Sequence[Fraction | int | float]],
    presentation_basis: Sequence[Any],
    presentation_origin: Sequence[Any] | str,
    selected_basis: Sequence[Any],
    selected_origin: Sequence[Any] | str,
    full_operation_correspondence: Sequence[Any],
    parent_setting_id: int | None = None,
    tolerance: float = 1e-8,
) -> tuple[FractionPoint, ...]:
    """Map a presentation grid through one selected magnetic BNS embedding.

    ``presentation_basis/origin`` map presentation child-cinter coordinates to
    the selected parent cinter setting. ``selected_basis/origin`` are the
    accepted subgroup result and map canonical child-PML coordinates to the
    Source-default parent PML. ``parent_setting_id`` supplies the intervening
    default-to-selected parent setting bridge. The retained full operation
    correspondence then validates the affine map on every grid row. Search-only
    candidate basis/origin values are deliberately not accepted.
    """

    if not presentation_points or not full_operation_correspondence:
        return ()
    try:
        direct_points = tuple(
            tuple(_input_fraction(value) for value in point)
            for point in presentation_points
        )
        if any(len(point) != 3 for point in direct_points):
            return ()
        display_basis = _fraction_matrix(presentation_basis)
        display_inverse = _matrix_inverse(display_basis)
        display_origin = _fraction_origin(presentation_origin)
        result_basis = _fraction_matrix(selected_basis)
        result_inverse = _matrix_inverse(result_basis)
        result_origin = _fraction_origin(selected_origin)
        data = source_tables()
        parent_pml_to_cinter = tuple(
            tuple(Fraction(value) for value in row)
            for row in data.pml_to_cinter_matrix(int(parent_sg))
        )
        parent_cinter_to_pml = _matrix_inverse(parent_pml_to_cinter)
        parent_setting_origin = tuple(
            Fraction(value) for value in data.cml_to_cinter_origin(int(parent_sg))
        )
        setting_bridge = _parent_setting_bridge(int(parent_sg), parent_setting_id)
        child_pml_to_cinter = tuple(
            tuple(Fraction(value) for value in row)
            for row in data.pml_to_cinter_matrix(int(child_sg))
        )
        child_setting_origin = tuple(
            Fraction(value) for value in data.cml_to_cinter_origin(int(child_sg))
        )
        correspondences = tuple(
            (
                int(_correspondence_value(item, "input_slot") or index),
                _operation_record(item, "input_record"),
                _operation_record(item, "canonical_record"),
            )
            for index, item in enumerate(full_operation_correspondence, start=1)
        )
        if len({slot for slot, _raw, _canonical in correspondences}) != len(correspondences):
            return ()
    except (IndexError, KeyError, TypeError, ValueError, ZeroDivisionError):
        return ()

    def presentation_to_parent_pml(point: FractionPoint) -> FractionPoint:
        parent_cinter_offset = _row_multiply(point, display_basis)
        parent_cinter = tuple(
            parent_cinter_offset[axis] + display_origin[axis]
            for axis in range(3)
        )
        parent_default_cinter = _parent_point_to_default(
            parent_cinter, setting_bridge
        )
        return _row_multiply(
            tuple(
                parent_default_cinter[axis] - parent_setting_origin[axis]
                for axis in range(3)
            ),
            parent_cinter_to_pml,
        )

    def parent_pml_to_presentation(point: FractionPoint) -> FractionPoint:
        cinter_offset = _row_multiply(point, parent_pml_to_cinter)
        parent_cinter = tuple(
            cinter_offset[axis] + parent_setting_origin[axis]
            for axis in range(3)
        )
        parent_cinter = _parent_point_from_default(parent_cinter, setting_bridge)
        return _fold(
            _row_multiply(
                tuple(parent_cinter[axis] - display_origin[axis] for axis in range(3)),
                display_inverse,
            )
        )

    def parent_to_canonical_pml(point: FractionPoint) -> FractionPoint:
        return _row_multiply(
            tuple(point[axis] - result_origin[axis] for axis in range(3)),
            result_inverse,
        )

    def canonical_to_parent_pml(point: FractionPoint) -> FractionPoint:
        offset = _row_multiply(point, result_basis)
        return tuple(offset[axis] + result_origin[axis] for axis in range(3))  # type: ignore[return-value]

    def canonical_pml_to_cinter(point: FractionPoint) -> FractionPoint:
        offset = _row_multiply(point, child_pml_to_cinter)
        return _fold(
            tuple(offset[axis] + child_setting_origin[axis] for axis in range(3))
        )

    anchors: list[FractionPoint] = []
    canonical_points: list[FractionPoint] = []
    parent_points: list[FractionPoint] = []
    for point in direct_points:
        parent_point = presentation_to_parent_pml(point)  # type: ignore[arg-type]
        canonical_point = parent_to_canonical_pml(parent_point)
        if not _periodic_close(
            canonical_to_parent_pml(canonical_point), parent_point, tolerance
        ):
            return ()
        parent_points.append(parent_point)
        canonical_points.append(canonical_point)
        anchors.append(canonical_pml_to_cinter(canonical_point))

    assigned: list[FractionPoint | None] = [None] * len(direct_points)
    covered: set[int] = set()
    for parent_point in parent_points:
        canonical_seed = parent_to_canonical_pml(parent_point)
        for _slot, input_record, canonical_record in correspondences:
            raw_image = _apply_pml_operation(int(parent_sg), input_record, parent_point)
            canonical_image = _apply_pml_operation(int(child_sg), canonical_record, canonical_seed)
            if not _periodic_close(
                parent_to_canonical_pml(raw_image), canonical_image, tolerance
            ):
                return ()
            direct_image = parent_pml_to_presentation(raw_image)
            matches = [
                index
                for index, candidate in enumerate(direct_points)
                if _periodic_close(direct_image, candidate, tolerance)
            ]
            if len(matches) != 1:
                return ()
            target = matches[0]
            # PML integer lifts become centering translations in a conventional
            # child setting.  Compare the paired action in primitive PML, then
            # retain the target row's direct cinter anchor rather than choosing
            # a path-dependent conventional lift.
            if not _periodic_close(canonical_image, canonical_points[target], tolerance):
                return ()
            standard_image = anchors[target]
            if assigned[target] is not None and not _periodic_close(
                assigned[target], standard_image, tolerance
            ):
                return ()
            assigned[target] = standard_image
            covered.add(target)
    if covered != set(range(len(direct_points))) or any(point is None for point in assigned):
        return ()
    return tuple(point for point in assigned if point is not None)


def presentation_branch_labels_from_correspondence(
    *,
    magnetic_group: int,
    branches: Sequence[tuple[str, Sequence[Fraction | int | float], int]],
) -> tuple[str, ...]:
    """Map canonical branch labels to one presentation child frame.

    The selected operation correspondence establishes the orbit components;
    this helper only transports their printed Wyckoff-row identities.  A
    mapping is accepted when direct Source-table identification preserves each
    component multiplicity and induces a bijection between canonical and
    presentation labels.
    """

    canonical_to_presentation: dict[str, str] = {}
    presentation_to_canonical: dict[str, str] = {}
    labels: list[str] = []
    for canonical_label, point, multiplicity in branches:
        try:
            identification = identify_magnetic_wyckoff_branch(
                int(magnetic_group), point, setting="cinter"
            )
            orbit = magnetic_orbit_points(
                int(magnetic_group), point, record_setting="cinter"
            )
        except (IndexError, KeyError, TypeError, ValueError, ZeroDivisionError):
            return ()
        presentation_label = str(identification.row.label)
        canonical_label = str(canonical_label)
        if (
            len(orbit) != int(multiplicity)
            or canonical_to_presentation.get(canonical_label, presentation_label)
            != presentation_label
            or presentation_to_canonical.get(presentation_label, canonical_label)
            != canonical_label
        ):
            return ()
        canonical_to_presentation[canonical_label] = presentation_label
        presentation_to_canonical[presentation_label] = canonical_label
        labels.append(presentation_label)
    return tuple(labels)


def presentation_grid_rows_for_magnetic_groups(
    *,
    magnetic_group: int,
    label_prefix: str,
    presentation_points: Sequence[Sequence[float]],
    groups: Sequence[Any],
    to_standard: Callable[[Any], tuple[Fraction, Fraction, Fraction]],
    tolerance: float = 1e-7,
) -> list[dict[str, Any]]:
    """Return direct presentation branches under a unique full-orbit bijection."""

    if not presentation_points:
        return []
    standard_points = [to_standard(point) for point in presentation_points]

    def periodic_match(left: Any, right: Any) -> bool:
        return all(
            abs(
                float(Fraction(left[axis]) - Fraction(right[axis]))
                - round(float(Fraction(left[axis]) - Fraction(right[axis])))
            )
            <= tolerance
            for axis in range(3)
        )

    def same_point_set(
        left: Any,
        right: Any,
        right_index: PeriodicPointIndex | None,
    ) -> bool:
        if len(left) != len(right):
            return False
        unmatched = set(range(len(right)))
        for point in left:
            matches = [
                index
                for index in _periodic_match_indices(
                    point,
                    right,
                    tolerance,
                    periodic_match,
                    right_index,
                )
                if index in unmatched
            ]
            if len(matches) != 1:
                return False
            unmatched.remove(matches[0])
        return not unmatched

    standard_point_index = _periodic_point_index(standard_points, tolerance)
    components: list[tuple[tuple[int, ...], tuple[Any, ...]]] = []
    unused = set(range(len(standard_points)))
    while unused:
        seed = min(unused)
        try:
            orbit = magnetic_orbit_points(
                int(magnetic_group), standard_points[seed], record_setting="cinter"
            )
        except (IndexError, KeyError, TypeError, ValueError, ZeroDivisionError):
            return []
        component: set[int] = set()
        for point in orbit:
            matches = _periodic_match_indices(
                point,
                standard_points,
                tolerance,
                periodic_match,
                standard_point_index,
            )
            if len(matches) != 1:
                return []
            component.add(matches[0])
        if len(component) != len(orbit) or seed not in component or not component <= unused:
            return []
        components.append((tuple(sorted(component)), tuple(orbit)))
        unused.difference_update(component)

    group_orbits: list[tuple[Any, ...]] = []
    for group in groups:
        try:
            orbit = magnetic_orbit_points(
                int(magnetic_group), group.standard_representative, record_setting="cinter"
            )
        except (IndexError, KeyError, TypeError, ValueError, ZeroDivisionError):
            return []
        if len(orbit) != int(group.multiplicity):
            return []
        group_orbits.append(tuple(orbit))
    group_orbit_indexes = [
        _periodic_point_index(orbit, tolerance)
        for orbit in group_orbits
    ]
    group_points = tuple(point for orbit in group_orbits for point in orbit)
    group_point_owners = tuple(
        group_index
        for group_index, orbit in enumerate(group_orbits)
        for _point in orbit
    )
    group_point_index = _periodic_point_index(group_points, tolerance)

    candidates: list[list[int]] = []
    for _indices, component_orbit in components:
        candidate_group_indices = (
            _periodic_match_owner_indices(
                component_orbit[0],
                group_points,
                group_point_owners,
                tolerance,
                periodic_match,
                group_point_index,
            )
            if component_orbit
            else tuple(range(len(groups)))
        )
        matching_groups = []
        for group_index in candidate_group_indices:
            group = groups[group_index]
            group_orbit = group_orbits[group_index]
            if len(component_orbit) != int(group.multiplicity):
                continue
            if same_point_set(
                component_orbit,
                group_orbit,
                group_orbit_indexes[group_index],
            ):
                matching_groups.append(group_index)
        candidates.append(matching_groups)
    if any(not row for row in candidates):
        return []
    assignments: list[tuple[int, ...]] = []

    def assign(component_index: int, used_groups: set[int], selected: list[int]) -> None:
        if len(assignments) > 1:
            return
        if component_index == len(candidates):
            assignments.append(tuple(selected))
            return
        for group_index in candidates[component_index]:
            if group_index in used_groups:
                continue
            used_groups.add(group_index)
            selected.append(group_index)
            assign(component_index + 1, used_groups, selected)
            selected.pop()
            used_groups.remove(group_index)

    assign(0, set(), [])
    if (
        len(assignments) != 1
        or len(assignments[0]) != len(groups)
        or set(assignments[0]) != set(range(len(groups)))
    ):
        return []

    rows: list[dict[str, Any]] = []
    covered_direct: set[int] = set()
    presentation_point_index = _periodic_point_index(
        presentation_points,
        tolerance,
    )
    for indices, group_index in zip(
        (component[0] for component in components), assignments[0], strict=True
    ):
        group = groups[group_index]
        seed = min(indices)
        direct_point = presentation_points[seed]
        try:
            direct_orbit = magnetic_orbit_points(
                int(magnetic_group), direct_point, record_setting="cinter"
            )
            identification = identify_magnetic_wyckoff_branch(
                int(magnetic_group), direct_point, setting="cinter"
            )
        except (IndexError, KeyError, TypeError, ValueError, ZeroDivisionError):
            return []
        direct_indices: set[int] = set()
        for point in direct_orbit:
            matches = _periodic_match_indices(
                point,
                presentation_points,
                tolerance,
                periodic_match,
                presentation_point_index,
            )
            if len(matches) != 1:
                return []
            direct_indices.add(matches[0])
        mapped_direct = tuple(to_standard(point) for point in direct_orbit)
        if (
            len(direct_orbit) != int(group.multiplicity)
            or direct_indices != set(indices)
            or direct_indices & covered_direct
            or not same_point_set(
                mapped_direct,
                group_orbits[group_index],
                group_orbit_indexes[group_index],
            )
        ):
            return []
        covered_direct.update(direct_indices)
        rows.append(
            {
                "label": f"{label_prefix}_{len(rows) + 1}",
                "site": f"{int(group.multiplicity)}{identification.row.label}",
                "xyz": list(direct_point),
                "_mode_row_orbit_points": [
                    list(presentation_points[index])
                    for index in sorted(direct_indices)
                ],
            }
        )
    return rows if covered_direct == set(range(len(presentation_points))) else []
