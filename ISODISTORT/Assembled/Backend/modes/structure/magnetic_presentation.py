"""Source-only magnetic orbit bridge for selected presentation grids."""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from typing import Any, Sequence

from ISODISTORT.Assembled.Backend.exactmath import (
    fraction_matrix_inverse3 as _matrix_inverse,
    fraction_row_multiply3 as _row_multiply,
)
from ISODISTORT.Assembled.Backend.source.magnetic import data as magnetic_data
from ISODISTORT.Assembled.Backend.source.tables import source_tables
from ISODISTORT.Assembled.Backend.source.magnetic_operations import (
    _pml_point_operation_matrix,
    generate_magnetic_space_group_records,
)
from ISODISTORT.Assembled.Backend.modes.structure.magnetic_wyckoff import (
    _input_fraction,
    identify_magnetic_wyckoff_branch,
    magnetic_orbit_points,
)


FractionPoint = tuple[Fraction, Fraction, Fraction]
FractionMatrix = tuple[FractionPoint, FractionPoint, FractionPoint]


@dataclass(frozen=True)
class MagneticAtomAction:
    """Exact magnetic action on one canonical child atom population."""

    components: tuple[tuple[int, ...], ...]
    canonical_cinter_points: tuple[FractionPoint, ...]


def _fraction_matrix(values: Sequence[Any]) -> FractionMatrix:
    rows = tuple(values)
    if len(rows) == 3 and all(
        isinstance(row, Sequence) and len(row) == 3 for row in rows
    ):
        return tuple(tuple(_input_fraction(value) for value in row) for row in rows)  # type: ignore[return-value]
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


def _apply_pml_operation(
    sg: int,
    record: tuple[int, int, int, int, int],
    point: Sequence[Fraction],
) -> FractionPoint:
    table = magnetic_data().table
    ordinary_point_op = int(table["mag_point_op_mag2nonmag"][int(record[4]) - 1])
    rotation = _pml_point_operation_matrix(int(sg), ordinary_point_op)
    rotated = _row_multiply(tuple(Fraction(value) for value in point), rotation)
    translation = tuple(
        Fraction(int(record[axis]), int(record[3])) for axis in range(3)
    )
    return tuple(rotated[axis] + translation[axis] for axis in range(3))  # type: ignore[return-value]


def selected_magnetic_atom_action(
    *,
    magnetic_group: int,
    child_sg: int,
    parent_points: Sequence[Sequence[Fraction | int]],
    ordinary_orbits: Sequence[Sequence[int]],
    selected_basis: Sequence[Any],
    selected_origin: Sequence[Any] | str,
) -> MagneticAtomAction:
    """Build the exact atom permutation induced by the Source BNS group."""

    if not parent_points:
        raise ValueError("magnetic atom action requires canonical atoms")
    points = tuple(tuple(Fraction(value) for value in point) for point in parent_points)
    if any(len(point) != 3 for point in points):
        raise ValueError("magnetic atom action requires point triplets")

    result_basis = _fraction_matrix(selected_basis)
    result_inverse = _matrix_inverse(result_basis)
    result_origin = _fraction_origin(selected_origin)
    data = source_tables()
    child_pml_to_cinter = tuple(
        tuple(Fraction(value) for value in row)
        for row in data.pml_to_cinter_matrix(int(child_sg))
    )
    child_setting_origin = tuple(
        Fraction(value) for value in data.cml_to_cinter_origin(int(child_sg))
    )
    records = generate_magnetic_space_group_records(
        int(magnetic_group), setting="binary"
    )
    if not records:
        raise ValueError(f"magnetic group {magnetic_group} has no Source operations")

    canonical_pml = tuple(
        _row_multiply(
            tuple(point[axis] - result_origin[axis] for axis in range(3)),
            result_inverse,
        )
        for point in points
    )

    def canonical_pml_to_cinter(point: FractionPoint) -> FractionPoint:
        offset = _row_multiply(point, child_pml_to_cinter)
        return _fold(
            tuple(offset[axis] + child_setting_origin[axis] for axis in range(3))
        )

    canonical_cinter = tuple(canonical_pml_to_cinter(point) for point in canonical_pml)
    cinter_owner: dict[FractionPoint, int] = {}
    for index, point in enumerate(canonical_cinter):
        previous = cinter_owner.setdefault(point, index)
        if previous != index:
            raise ValueError(
                f"canonical atoms share one child-cinter point: {previous}, {index}"
            )
    orbit_indices = tuple(
        tuple(int(index) for index in orbit) for orbit in ordinary_orbits
    )
    if (
        not orbit_indices
        or sorted(index for orbit in orbit_indices for index in orbit)
        != list(range(len(points)))
        or any(not orbit for orbit in orbit_indices)
    ):
        raise ValueError("ordinary child sites do not partition magnetic atoms")
    adjacency: list[set[int]] = [{index} for index in range(len(points))]
    for orbit in orbit_indices:
        for index in orbit:
            adjacency[index].update(orbit)
    for operation_index, canonical_record in enumerate(records, start=1):
        targets: set[int] = set()
        for source, canonical_point in enumerate(canonical_pml):
            canonical_image = canonical_pml_to_cinter(
                _apply_pml_operation(int(child_sg), canonical_record, canonical_point)
            )
            target = cinter_owner.get(canonical_image)
            if target is None:
                raise KeyError(
                    "magnetic Source operation leaves the canonical atom population: "
                    f"operation={operation_index}, source={source}, "
                    f"image={canonical_image}"
                )
            targets.add(target)
            adjacency[source].add(target)
            adjacency[target].add(source)
        if targets != set(range(len(points))):
            raise ValueError(
                f"magnetic operation {operation_index} is not an atom permutation"
            )

    visited: set[int] = set()
    components: list[tuple[int, ...]] = []
    for start in range(len(points)):
        if start in visited:
            continue
        component: set[int] = set()
        pending = [start]
        while pending:
            current = pending.pop()
            if current in component:
                continue
            component.add(current)
            pending.extend(adjacency[current] - component)
        visited.update(component)
        components.append(tuple(sorted(component)))
    return MagneticAtomAction(
        components=tuple(components),
        canonical_cinter_points=canonical_cinter,
    )


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
