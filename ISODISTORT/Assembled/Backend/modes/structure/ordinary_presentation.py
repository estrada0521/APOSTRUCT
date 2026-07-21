"""Source-only matching of Formula15 rows to presentation-grid orbits."""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
import re
from typing import Any, Mapping, Sequence

import gemmi

from ISODISTORT.Assembled.Backend.exactmath import (
    fraction_matrix_inverse3,
    fraction_row_multiply3 as _row_multiply,
)


ExactPoint = tuple[Fraction, Fraction, Fraction]
ExactMatrix3 = tuple[ExactPoint, ExactPoint, ExactPoint]


@dataclass(frozen=True)
class FormulaGridOrbitAssignment:
    """One unambiguous Formula15 row to complete-grid orbit assignment."""

    formula_index: int
    formula_label: str
    formula_site: str
    formula15: tuple[int, ...]
    formula_representative: ExactPoint
    grid_component_index: int
    grid_indices: tuple[int, ...]
    representative_grid_index: int


def match_formula_assignments_to_geometric_grid_orbits(
    *,
    child_sg: int,
    formula_assignments: Sequence[FormulaGridOrbitAssignment],
    presentation_grid_points: Sequence[Sequence[Any]],
    geometric_rows: Sequence[Mapping[str, Any]],
) -> tuple[int, ...] | None:
    """Return geometric-row indices in Formula15 component order.

    Formula15 owns the split topology and Source row order.  The geometric
    rows own presentation-setting representatives and Wyckoff letters.  The
    two surfaces are joined only by a complete, exact child-orbit partition;
    site letters are not compared because a setting change can relabel the
    same orbit.  Multiplicity and component membership remain exact.
    """

    if (
        not formula_assignments
        or len(formula_assignments) != len(geometric_rows)
        or not presentation_grid_points
    ):
        return None
    grid = tuple(_point(point) for point in presentation_grid_points)
    if any(point is None for point in grid):
        return None
    exact_grid = tuple(point for point in grid if point is not None)
    if len(set(exact_grid)) != len(exact_grid):
        return None
    try:
        space_group = gemmi.find_spacegroup_by_number(int(child_sg))
        if space_group is None:
            return None
        operations = tuple(space_group.operations())
    except (RuntimeError, TypeError, ValueError):
        return None
    if not operations:
        return None

    grid_index = {point: index for index, point in enumerate(exact_grid)}
    geometric_components: list[tuple[int, ...]] = []
    geometric_coverage: set[int] = set()
    for row in geometric_rows:
        representative = _point(row.get("xyz"))
        multiplicity = _site_multiplicity(str(row.get("site") or ""))
        if representative is None or multiplicity is None:
            return None
        orbit = _orbit(operations, representative)
        if not orbit or not orbit <= grid_index.keys():
            return None
        indices = tuple(sorted(grid_index[point] for point in orbit))
        if len(indices) != multiplicity or set(indices) & geometric_coverage:
            return None
        source_index = row.get("_source_index")
        if source_index is not None:
            try:
                if int(source_index) not in indices:
                    return None
            except (TypeError, ValueError):
                return None
        geometric_components.append(indices)
        geometric_coverage.update(indices)
    if geometric_coverage != set(range(len(exact_grid))):
        return None

    formula_coverage: set[int] = set()
    valid_indices = set(range(len(exact_grid)))
    for assignment in formula_assignments:
        multiplicity = _site_multiplicity(assignment.formula_site)
        indices = set(assignment.grid_indices)
        if (
            multiplicity is None
            or len(indices) != multiplicity
            or indices & formula_coverage
            or not indices <= valid_indices
        ):
            return None
        formula_coverage.update(indices)
    if formula_coverage != valid_indices:
        return None
    ordered_rows: list[int] = []
    used_rows: set[int] = set()
    for assignment in formula_assignments:
        candidates = [
            index
            for index, component in enumerate(geometric_components)
            if index not in used_rows and component == assignment.grid_indices
        ]
        if len(candidates) != 1:
            return None
        used_rows.add(candidates[0])
        ordered_rows.append(candidates[0])
    return tuple(ordered_rows) if len(used_rows) == len(geometric_rows) else None


def _mod1(value: Fraction) -> Fraction:
    return value % 1


def _input_fraction(value: Any) -> Fraction:
    """Recover the exact decimal carried through presentation float math."""

    if isinstance(value, float):
        decimal = Fraction(str(round(value, 9)))
        crystallographic = decimal.limit_denominator(384)
        if abs(crystallographic - decimal) <= Fraction(1, 200_000_000):
            return crystallographic
        return decimal
    return Fraction(str(value))


def _point(values: Any) -> ExactPoint | None:
    if not isinstance(values, (list, tuple)):
        return None
    if len(values) == 4:
        try:
            denominator = int(values[3])
            if denominator == 0:
                return None
            return tuple(
                _mod1(Fraction(int(values[axis]), denominator))
                for axis in range(3)
            )  # type: ignore[return-value]
        except (TypeError, ValueError, ZeroDivisionError):
            return None
    if len(values) != 3:
        return None
    try:
        return tuple(_mod1(_input_fraction(values[axis])) for axis in range(3))  # type: ignore[return-value]
    except (TypeError, ValueError, ZeroDivisionError):
        return None


def _unfolded_point(values: Any) -> ExactPoint | None:
    if isinstance(values, str):
        text = values.strip()
        if text.startswith("(") and text.endswith(")"):
            text = text[1:-1]
        values = tuple(part.strip() for part in text.split(","))
    if not isinstance(values, (list, tuple)):
        return None
    if len(values) == 4:
        try:
            denominator = int(values[3])
            if denominator == 0:
                return None
            return tuple(
                Fraction(int(values[axis]), denominator)
                for axis in range(3)
            )  # type: ignore[return-value]
        except (TypeError, ValueError, ZeroDivisionError):
            return None
    if len(values) != 3:
        return None
    try:
        return tuple(_input_fraction(values[axis]) for axis in range(3))  # type: ignore[return-value]
    except (TypeError, ValueError, ZeroDivisionError):
        return None


def _matrix3(values: Any) -> ExactMatrix3 | None:
    if not isinstance(values, (list, tuple)) or len(values) != 3:
        return None
    rows = tuple(_unfolded_point(row) for row in values)
    if any(row is None for row in rows):
        return None
    return tuple(row for row in rows if row is not None)  # type: ignore[return-value]


def _inverse3(matrix: ExactMatrix3) -> ExactMatrix3 | None:
    try:
        return fraction_matrix_inverse3(matrix)
    except ValueError:
        return None


def _site_multiplicity(site: str) -> int | None:
    match = re.fullmatch(r"\s*(\d+)\s*[^\s]+\s*", str(site or ""))
    return int(match.group(1)) if match else None


def _operation_image(operation: gemmi.Op, point: ExactPoint) -> ExactPoint:
    denominator = int(gemmi.Op.DEN)
    return tuple(
        _mod1(
            sum(
                Fraction(int(operation.rot[axis][column]), denominator) * point[column]
                for column in range(3)
            )
            + Fraction(int(operation.tran[axis]), denominator)
        )
        for axis in range(3)
    )  # type: ignore[return-value]


def _orbit(operations: Sequence[gemmi.Op], point: ExactPoint) -> frozenset[ExactPoint]:
    return frozenset(_operation_image(operation, point) for operation in operations)


def _centering_translations(child_sg: int) -> tuple[ExactPoint, ...] | None:
    try:
        space_group = gemmi.find_spacegroup_by_number(int(child_sg))
        if space_group is None:
            return None
        operations = tuple(space_group.operations())
    except (RuntimeError, TypeError, ValueError):
        return None
    denominator = int(gemmi.Op.DEN)
    identity = (
        (denominator, 0, 0),
        (0, denominator, 0),
        (0, 0, denominator),
    )
    translations: list[ExactPoint] = []
    for operation in operations:
        rotation = tuple(tuple(int(value) for value in row) for row in operation.rot)
        if rotation != identity:
            continue
        translation = tuple(
            _mod1(Fraction(int(operation.tran[axis]), denominator))
            for axis in range(3)
        )
        if translation not in translations:
            translations.append(translation)  # type: ignore[arg-type]
    zero = (Fraction(0), Fraction(0), Fraction(0))
    if zero not in translations:
        return None
    return (zero, *(translation for translation in translations if translation != zero))


def _formula_representative(row: Mapping[str, Any]) -> ExactPoint | None:
    exact = row.get("_wyckoff_representative_fraction")
    return _point(exact if exact is not None else row.get("xyz"))


def assign_formula_source_occurrences_to_presentation_grid(
    *,
    child_sg: int | None = None,
    formula_rows: Sequence[Mapping[str, Any]],
    presentation_grid_points: Sequence[Sequence[Any]],
    parent_pml_to_cinter_matrix: Sequence[Sequence[Any]],
    parent_pml_to_cinter_origin: Sequence[Any],
    presentation_basis: Sequence[Sequence[Any]],
    presentation_origin: Sequence[Any],
) -> tuple[FormulaGridOrbitAssignment, ...] | None:
    """Transport Source split occurrences directly into a presentation grid.

    Formula15 provenance supplies complete child-orbit occurrence sets as
    exact parent-PML points.  Each point is moved through parent cinter into
    the selected child basis.  A result is returned only when those sets form
    a unique, multiplicity-preserving partition of the complete grid.
    """

    if not formula_rows or not presentation_grid_points:
        return None
    pml_to_cinter = _matrix3(parent_pml_to_cinter_matrix)
    cinter_origin = _unfolded_point(parent_pml_to_cinter_origin)
    selected_basis = _matrix3(presentation_basis)
    selected_origin = _unfolded_point(presentation_origin)
    selected_inverse = None if selected_basis is None else _inverse3(selected_basis)
    if (
        pml_to_cinter is None
        or cinter_origin is None
        or selected_origin is None
        or selected_inverse is None
    ):
        return None
    grid = tuple(_point(point) for point in presentation_grid_points)
    if any(point is None for point in grid):
        return None
    exact_grid = tuple(point for point in grid if point is not None)
    if len(set(exact_grid)) != len(exact_grid):
        return None
    grid_index = {point: index for index, point in enumerate(exact_grid)}

    assignments: list[FormulaGridOrbitAssignment] = []
    used_indices: set[int] = set()
    component_order: dict[int, int] = {}
    for formula_index, row in enumerate(formula_rows):
        site = str(row.get("site") or "")
        multiplicity = _site_multiplicity(site)
        representative = _formula_representative(row)
        occurrences = row.get("_wyckoff_source_occurrences")
        if (
            multiplicity is None
            or representative is None
            or not isinstance(occurrences, (list, tuple))
            or not occurrences
        ):
            return None
        occurrence_count = len(occurrences)
        if occurrence_count <= 0 or multiplicity % occurrence_count:
            return None
        expansion = multiplicity // occurrence_count
        if expansion == 1:
            translations = ((Fraction(0), Fraction(0), Fraction(0)),)
        else:
            if child_sg is None:
                return None
            all_translations = _centering_translations(int(child_sg))
            if all_translations is None or len(all_translations) != expansion:
                return None
            translations = all_translations
        transported: dict[int, tuple[int, int]] = {}
        for occurrence in occurrences:
            if not isinstance(occurrence, Mapping):
                return None
            parent_pml = _unfolded_point(occurrence.get("parent_point_fraction"))
            try:
                flat_index = int(occurrence["flat_index"])
            except (KeyError, TypeError, ValueError):
                return None
            if parent_pml is None:
                return None
            parent_cinter_linear = _row_multiply(parent_pml, pml_to_cinter)
            parent_cinter = tuple(
                parent_cinter_linear[axis] + cinter_origin[axis]
                for axis in range(3)
            )
            selected = _row_multiply(
                tuple(parent_cinter[axis] - selected_origin[axis] for axis in range(3)),
                selected_inverse,
            )
            for centering_ordinal, translation in enumerate(translations):
                target = tuple(
                    _mod1(_input_fraction(float(selected[axis] + translation[axis])))
                    for axis in range(3)
                )
                target_index = grid_index.get(target)
                if target_index is None:
                    return None
                identity = (flat_index, centering_ordinal)
                previous = transported.get(target_index)
                if previous is not None and previous != identity:
                    return None
                transported[target_index] = identity
        indices = tuple(sorted(transported))
        if len(indices) != multiplicity or set(indices) & used_indices:
            return None
        try:
            representative_flat_index = int(row["_wyckoff_source_representative_flat_index"])
            representative_index = next(
                index for index, identity in transported.items()
                if identity == (representative_flat_index, 0)
            )
            formula15 = tuple(int(value) for value in row.get("_wyckoff_formula15") or ())
        except (KeyError, StopIteration, TypeError, ValueError):
            return None
        used_indices.update(indices)
        component_order[min(indices)] = formula_index
        assignments.append(
            FormulaGridOrbitAssignment(
                formula_index=formula_index,
                formula_label=str(row.get("label") or ""),
                formula_site=site,
                formula15=formula15,
                formula_representative=representative,
                grid_component_index=0,
                grid_indices=indices,
                representative_grid_index=representative_index,
            )
        )

    if used_indices != set(range(len(exact_grid))):
        return None
    ordered_components = {
        formula_index: component_index
        for component_index, (_minimum, formula_index) in enumerate(sorted(component_order.items()))
    }
    return tuple(
        FormulaGridOrbitAssignment(
            formula_index=item.formula_index,
            formula_label=item.formula_label,
            formula_site=item.formula_site,
            formula15=item.formula15,
            formula_representative=item.formula_representative,
            grid_component_index=ordered_components[item.formula_index],
            grid_indices=item.grid_indices,
            representative_grid_index=item.representative_grid_index,
        )
        for item in assignments
    )
