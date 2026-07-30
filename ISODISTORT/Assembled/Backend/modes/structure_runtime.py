"""Undistorted-structure construction helpers for modes."""

from __future__ import annotations

from dataclasses import replace
from fractions import Fraction
from itertools import product
import math
from numbers import Integral
import re
from typing import Any
from ISODISTORT.Assembled.Backend.exactmath import (
    fraction_matrix_inverse3,
    fraction_matrix_multiply3,
    fraction_row_multiply3,
    integer_determinant3,
)
from ISODISTORT.Assembled.Backend.source.tables import SourceTables
from ISODISTORT.Assembled.Backend.modes.structure.wyckoff_split import (
    undistorted_rows_from_wyckoff_split,
)
from ISODISTORT.Assembled.Backend.modes.structure.ordinary_presentation import (
    _centering_translations,
    _input_fraction,
    formula_child_sites_in_presentation,
)
from ISODISTORT.Assembled.Backend.modes.presentation import present_mode_rows
from ISODISTORT.Assembled.Backend.modes.periodic import periodic_float_close3
from ISODISTORT.Assembled.Backend.modes.structure.child_atom_layout import (
    ChildAtomLayout,
    ChildAtomPresentationRow,
    child_atom_layout_in_presentation_order,
    child_atom_layout_from_formula_sites,
    exact_operation_record,
    operation_record_row_assignments,
    regroup_child_atom_layout,
)
from ISODISTORT.Assembled.Backend.modes.structure.magnetic_wyckoff import (
    identify_magnetic_wyckoff,
    magnetic_group_setting,
    magnetic_orbit_points,
)
from ISODISTORT.Assembled.Backend.modes.structure.magnetic_presentation import (
    presentation_branch_labels_from_correspondence,
    selected_magnetic_atom_action,
)
from ISODISTORT.Assembled.Backend.modes.engine.decoder import ModeDataDecoder

from ISODISTORT.Assembled.Backend.modes.common import (
    _assembled_data,
    _float_matrix_inverse_3,
    _fraction_matmul,
    _fraction_matrix_inverse_3,
    _fraction_row_multiply,
    _fraction_vecadd,
    _fraction_vecsub,
    _isotropy_from_opd_row,
    _matrix_from_basis_tuple,
    _origin_record_vector,
    _origin_vector,
    _row_multiply,
)
from ISODISTORT.Assembled.Backend.modes.site_transport import _parent_setting_bridge


def _web_magnetic_occurrence_gauge(
    decoder: ModeDataDecoder,
    sg: int,
    rows: list[dict[str, Any]],
) -> list[dict[str, Any]]:
    """Bridge faithful axial rows to Web's occurrence-frame convention.

    Each operation record is the Source-canonical-site to occurrence witness.
    Its determinant is invariant under the parent-to-child setting conjugation,
    so no point-operation ordinal is part of this presentation rule.
    """

    out: list[dict[str, Any]] = []
    for row in rows:
        record = row.get("_operation_record")
        if not isinstance(record, (list, tuple)) or len(record) != 5:
            out.append(row)
            continue
        lattice = int(decoder.iso.space["ispace_lattice"][int(sg) - 1])
        code = int(
            decoder.iso.space["ipoint_op_psettings"][
                (lattice - 1) * 72 + int(record[4]) - 1
            ]
        )
        matrix = []
        for _ in range(9):
            matrix.append(code % 3 - 1)
            code //= 3
        determinant = integer_determinant3(matrix)
        if determinant >= 0:
            out.append(row)
        else:
            out.append(
                {
                    **row,
                    "dxyz": [
                        -float(value) for value in row.get("dxyz") or (0.0, 0.0, 0.0)
                    ],
                }
            )
    return out


def _exact_space_group_number(value: Any, name: str) -> int:
    if isinstance(value, bool) or not isinstance(value, Integral):
        raise TypeError(f"{name} must be an exact integer")
    number = int(value)
    if not 1 <= number <= 230:
        raise ValueError(f"{name} must be between 1 and 230")
    return number


def _source_child_atom_layout_for_site(
    *,
    sg: int,
    child_sg: int | None,
    site: dict[str, Any],
    label_prefix: str,
    split_basis: list[list[float]] | None,
    split_origin: Any,
    parent_params: dict[str, float] | None,
) -> ChildAtomLayout:
    """Construct the complete child atom layout from Source occurrences."""

    if child_sg is None or split_basis is None:
        raise ValueError("child atom layout requires a subgroup basis")
    parent_space_group = _exact_space_group_number(sg, "parent space group")
    child_space_group = _exact_space_group_number(child_sg, "child space group")
    data = _assembled_data()
    formula_rows = undistorted_rows_from_wyckoff_split(
        data,
        parent_sg=parent_space_group,
        child_sg=child_space_group,
        parent_wyckoff_row_id=site.get("wyckoff_row_id"),
        label_prefix=label_prefix,
        parent_params=parent_params,
        subgroup_basis=split_basis,
        subgroup_origin=split_origin,
    )
    if not formula_rows:
        raise ValueError("Source Wyckoff split produced no Formula15 rows")
    formula_sites = formula_child_sites_in_presentation(
        child_sg=child_space_group,
        formula_rows=formula_rows,
        subgroup_basis=split_basis,
        subgroup_origin=split_origin,
        child_pml_to_cinter_matrix=data.pml_to_cinter_matrix(child_space_group),
        child_pml_to_cinter_origin=data.cml_to_cinter_origin(child_space_group),
    )
    if formula_sites is None:
        raise ValueError("Formula15 did not partition the presentation atoms")
    return child_atom_layout_from_formula_sites(
        formula_sites,
        label_prefix=label_prefix,
    )


ExactMatrix3 = tuple[
    tuple[Fraction, Fraction, Fraction],
    tuple[Fraction, Fraction, Fraction],
    tuple[Fraction, Fraction, Fraction],
]
ExactVector3 = tuple[Fraction, Fraction, Fraction]
AffineOperation = tuple[ExactMatrix3, ExactVector3]


def _exact_integer_tuple(value: Any, size: int, name: str) -> tuple[int, ...]:
    if not isinstance(value, (list, tuple)) or len(value) != size:
        raise ValueError(f"{name} requires {size} entries")
    if any(isinstance(item, bool) or not isinstance(item, Integral) for item in value):
        raise TypeError(f"{name} must contain exact integers")
    return tuple(int(item) for item in value)


def _selected_source_subgroup_transform(
    selected_opd: dict[str, Any] | None,
) -> tuple[tuple[int, ...], tuple[int, int, int, int]]:
    iso = _isotropy_from_opd_row(selected_opd)
    if iso is None:
        raise ValueError("selected OPD has no Source subgroup transform")
    basis = _exact_integer_tuple(
        iso.get("source_basis_values"), 9, "Source subgroup basis"
    )
    raw_origin = _exact_integer_tuple(
        iso.get("source_origin_values"), 4, "Source subgroup origin"
    )
    if raw_origin[3] <= 0:
        raise ValueError("Source subgroup origin denominator must be positive")
    return basis, raw_origin  # type: ignore[return-value]


def _exact_presentation_origin(origin: Any) -> ExactVector3:
    if origin is None:
        return Fraction(0), Fraction(0), Fraction(0)
    if isinstance(origin, str):
        text = origin.strip()
        if not (text.startswith("(") and text.endswith(")")):
            raise ValueError(f"invalid presentation origin: {origin!r}")
        parts = tuple(part.strip() for part in text[1:-1].split(","))
        if len(parts) != 3:
            raise ValueError(f"invalid presentation origin: {origin!r}")
        return tuple(Fraction(part) for part in parts)  # type: ignore[return-value]
    record = _exact_integer_tuple(origin, 4, "presentation origin")
    if record[3] <= 0:
        raise ValueError("presentation origin denominator must be positive")
    return tuple(Fraction(record[axis], record[3]) for axis in range(3))  # type: ignore[return-value]


def _fraction_matrix_determinant3(matrix: ExactMatrix3) -> Fraction:
    return (
        matrix[0][0] * (matrix[1][1] * matrix[2][2] - matrix[1][2] * matrix[2][1])
        - matrix[0][1] * (matrix[1][0] * matrix[2][2] - matrix[1][2] * matrix[2][0])
        + matrix[0][2] * (matrix[1][0] * matrix[2][1] - matrix[1][1] * matrix[2][0])
    )


def _primitive_lattice_offsets(matrix: ExactMatrix3) -> tuple[ExactVector3, ...]:
    denominator = math.lcm(*(value.denominator for row in matrix for value in row), 1)
    offsets = {
        tuple(
            sum(Fraction(values[row]) * matrix[row][column] for row in range(3)) % 1
            for column in range(3)
        )
        for values in product(range(denominator), repeat=3)
    }
    determinant = abs(_fraction_matrix_determinant3(matrix))
    if determinant == 0:
        raise ValueError("child presentation transform is singular")
    expected = Fraction(1, 1) / determinant
    if expected.denominator != 1 or len(offsets) != int(expected):
        raise ValueError(
            "child presentation transform does not define a conventional lattice"
        )
    return tuple(sorted(offsets))  # type: ignore[return-value]


def _space_group_action_in_affine_setting(
    data: SourceTables,
    *,
    sg: int,
    matrix: ExactMatrix3,
    origin: ExactVector3,
) -> frozenset[AffineOperation]:
    inverse = fraction_matrix_inverse3(matrix)
    offsets = _primitive_lattice_offsets(matrix)
    operations: set[AffineOperation] = set()
    for x, y, z, denominator, point_op in data.generate_space_group_records(int(sg)):
        rotation = _point_rotation_matrix(data, int(sg), int(point_op))
        transformed_rotation = fraction_matrix_multiply3(
            inverse, fraction_matrix_multiply3(rotation, matrix)
        )
        translated = fraction_row_multiply3(
            (
                Fraction(int(x), int(denominator)),
                Fraction(int(y), int(denominator)),
                Fraction(int(z), int(denominator)),
            ),
            matrix,
        )
        rotated_origin = fraction_row_multiply3(origin, transformed_rotation)
        base_translation = tuple(
            (translated[axis] + origin[axis] - rotated_origin[axis]) % 1
            for axis in range(3)
        )
        for offset in offsets:
            operations.add(
                (
                    transformed_rotation,
                    tuple(
                        (base_translation[axis] + offset[axis]) % 1 for axis in range(3)
                    ),
                )
            )
    expected = len(data.generate_space_group_records(int(sg))) * len(offsets)
    if len(operations) != expected:
        raise ValueError(
            "child space-group action contains duplicate affine operations"
        )
    return frozenset(operations)


def _transform_affine_action(
    operations: frozenset[AffineOperation],
    *,
    matrix: ExactMatrix3,
    origin: ExactVector3,
) -> frozenset[AffineOperation]:
    inverse = fraction_matrix_inverse3(matrix)
    transformed: set[AffineOperation] = set()
    for rotation, translation in operations:
        output_rotation = fraction_matrix_multiply3(
            inverse, fraction_matrix_multiply3(rotation, matrix)
        )
        output_translation = fraction_row_multiply3(translation, matrix)
        rotated_origin = fraction_row_multiply3(origin, output_rotation)
        transformed.add(
            (
                output_rotation,
                tuple(
                    (output_translation[axis] + origin[axis] - rotated_origin[axis]) % 1
                    for axis in range(3)
                ),
            )
        )
    if len(transformed) != len(operations):
        raise ValueError("child setting transform collapses affine operations")
    return frozenset(transformed)


def _child_presentation_setting(
    *,
    data: SourceTables,
    parent_sg: int,
    child_sg: int,
    parent_setting_id: int | None,
    source_basis: tuple[int, ...],
    source_origin: tuple[int, int, int, int],
    presentation_basis_pml: tuple[int, ...],
    presentation_origin: Any,
) -> tuple[
    tuple[int, ...],
    ExactMatrix3,
    ExactVector3,
    tuple[tuple[int, frozenset[AffineOperation]], ...],
]:
    """Identify the displayed child setting and its exact coordinate map."""

    default_id = int(data.default_inter_setting_record(int(child_sg))["id"])
    changed = data.subgroup_change_setting_cinter(
        int(parent_sg),
        int(child_sg),
        source_basis,
        source_origin,
        parent_setting_id=parent_setting_id,
        subgroup_setting_id=default_id,
    )
    basis_denominator = int(changed["basis_denominator"])
    if basis_denominator == 0:
        raise ValueError("child setting transform has a zero basis denominator")
    child_to_parent = tuple(
        tuple(Fraction(int(value), basis_denominator) for value in row)
        for row in changed["basis"]
    )
    changed_origin = _exact_integer_tuple(changed["origin"], 4, "child setting origin")
    if changed_origin[3] <= 0:
        raise ValueError("child setting origin denominator must be positive")
    child_to_parent_origin = tuple(
        Fraction(changed_origin[axis], changed_origin[3]) for axis in range(3)
    )
    presentation_pml_values = _exact_integer_tuple(
        presentation_basis_pml, 9, "presentation PML basis"
    )
    presentation_pml = tuple(
        tuple(
            Fraction(presentation_pml_values[3 * row + column]) for column in range(3)
        )
        for row in range(3)
    )
    parent_pml_to_selected = data.pml_to_cinter_matrix(
        int(parent_sg), parent_setting_id
    )
    display_basis = fraction_matrix_multiply3(presentation_pml, parent_pml_to_selected)
    display_basis_inverse = fraction_matrix_inverse3(display_basis)
    child_to_display = fraction_matrix_multiply3(child_to_parent, display_basis_inverse)
    display_origin = _exact_presentation_origin(presentation_origin)
    child_to_display_origin = fraction_row_multiply3(
        tuple(child_to_parent_origin[axis] - display_origin[axis] for axis in range(3)),
        display_basis_inverse,
    )
    identity_basis = (1, 0, 0, 0, 1, 0, 0, 0, 1)
    zero_origin = (0, 0, 0, 1)
    coordinate_matches: list[int] = []
    for setting_id in data.inter_setting_ids_for_space_group(int(child_sg)):
        setting_change = data.subgroup_change_setting_cinter(
            int(child_sg),
            int(child_sg),
            identity_basis,
            zero_origin,
            parent_setting_id=default_id,
            subgroup_setting_id=int(setting_id),
        )
        inverse_basis_values = _exact_integer_tuple(
            tuple(value for row in setting_change["inverse_basis"] for value in row),
            9,
            "child inter-setting inverse basis",
        )
        inverse_denominator = setting_change["inverse_basis_denominator"]
        if isinstance(inverse_denominator, bool) or not isinstance(
            inverse_denominator, Integral
        ):
            raise TypeError("child inter-setting basis denominator must be an integer")
        if inverse_denominator <= 0:
            raise ValueError("child inter-setting basis denominator must be positive")
        inverse_matrix = tuple(
            tuple(
                Fraction(
                    inverse_basis_values[3 * row + column],
                    int(inverse_denominator),
                )
                for column in range(3)
            )
            for row in range(3)
        )
        inverse_origin_record = _exact_integer_tuple(
            setting_change["inverse_origin"],
            4,
            "child inter-setting inverse origin",
        )
        if inverse_origin_record[3] <= 0:
            raise ValueError("child inter-setting origin denominator must be positive")
        inverse_origin = tuple(
            Fraction(inverse_origin_record[axis], inverse_origin_record[3])
            for axis in range(3)
        )
        if inverse_matrix == child_to_display and all(
            (inverse_origin[axis] - child_to_display_origin[axis]) % 1 == 0
            for axis in range(3)
        ):
            coordinate_matches.append(int(setting_id))
    candidate_setting_ids = (
        tuple(coordinate_matches)
        if coordinate_matches
        else tuple(data.inter_setting_ids_for_space_group(int(child_sg)))
    )
    target = _transform_affine_action(
        _space_group_action_in_affine_setting(
            data,
            sg=int(child_sg),
            matrix=data.pml_to_cinter_matrix(int(child_sg), default_id),
            origin=data.cml_to_cinter_origin(int(child_sg), default_id),
        ),
        matrix=child_to_display,
        origin=child_to_display_origin,
    )
    matching_actions = tuple(
        (
            setting_id,
            _space_group_action_in_affine_setting(
                data,
                sg=int(child_sg),
                matrix=data.pml_to_cinter_matrix(int(child_sg), int(setting_id)),
                origin=data.cml_to_cinter_origin(int(child_sg), int(setting_id)),
            ),
        )
        for setting_id in candidate_setting_ids
    )
    matches = tuple(
        setting_id for setting_id, action in matching_actions if action == target
    )
    setting_ids = (default_id,) if default_id in matches else matches
    if not setting_ids:
        raise ValueError("displayed child action matches no Source inter setting")
    selected_actions = tuple(
        (setting_id, action)
        for setting_id, action in matching_actions
        if setting_id in setting_ids
    )
    return (
        setting_ids,
        child_to_display,
        child_to_display_origin,
        selected_actions,
    )


def _fraction_floor(value: Fraction) -> int:
    return value.numerator // value.denominator


def _fraction_ceil(value: Fraction) -> int:
    return -((-value.numerator) // value.denominator)


def _solve_exact_wyckoff_parameters(
    vectors: tuple[ExactVector3, ...],
    target: ExactVector3,
) -> tuple[Fraction, ...] | None:
    """Solve one Source Wyckoff representative exactly modulo the lattice."""

    base, *raw_parameters = vectors
    parameters = tuple(vector for vector in raw_parameters if any(vector))
    if not parameters:
        return (
            ()
            if all((target[axis] - base[axis]) % 1 == 0 for axis in range(3))
            else None
        )

    column_count = len(parameters)
    ranges: list[range] = []
    for axis in range(3):
        lower = sum(
            (
                coefficient
                for coefficient in (vector[axis] for vector in parameters)
                if coefficient < 0
            ),
            Fraction(0),
        )
        upper = sum(
            (
                coefficient
                for coefficient in (vector[axis] for vector in parameters)
                if coefficient > 0
            ),
            Fraction(0),
        )
        delta = target[axis] - base[axis]
        ranges.append(
            range(
                _fraction_ceil(lower - delta),
                _fraction_floor(upper - delta) + 1,
            )
        )

    for lattice_shift in product(*ranges):
        augmented = [
            [parameters[column][row] for column in range(column_count)]
            + [target[row] + lattice_shift[row] - base[row]]
            for row in range(3)
        ]
        pivot_row = 0
        pivot_columns: list[int] = []
        for column in range(column_count):
            pivot = next(
                (row for row in range(pivot_row, 3) if augmented[row][column] != 0),
                None,
            )
            if pivot is None:
                continue
            augmented[pivot_row], augmented[pivot] = (
                augmented[pivot],
                augmented[pivot_row],
            )
            divisor = augmented[pivot_row][column]
            augmented[pivot_row] = [value / divisor for value in augmented[pivot_row]]
            for row in range(3):
                if row == pivot_row or augmented[row][column] == 0:
                    continue
                factor = augmented[row][column]
                augmented[row] = [
                    augmented[row][index] - factor * augmented[pivot_row][index]
                    for index in range(column_count + 1)
                ]
            pivot_columns.append(column)
            pivot_row += 1
        if len(pivot_columns) != column_count:
            raise ValueError("Source Wyckoff parameters are linearly dependent")
        if any(
            all(row[column] == 0 for column in range(column_count))
            and row[column_count] != 0
            for row in augmented
        ):
            continue
        solution = [Fraction(0)] * column_count
        for row, column in enumerate(pivot_columns):
            solution[column] = augmented[row][column_count]
        if not all(Fraction(0) <= value < 1 for value in solution):
            continue
        if all(
            sum(
                parameters[column][axis] * solution[column]
                for column in range(column_count)
            )
            == target[axis] + lattice_shift[axis] - base[axis]
            for axis in range(3)
        ):
            return tuple(solution)
    return None


def _source_wyckoff_orbit_matches(
    *,
    vectors: tuple[ExactVector3, ...],
    operations: frozenset[AffineOperation],
    points: frozenset[ExactVector3],
) -> bool:
    active_vectors = tuple(vector for vector in vectors[1:] if any(vector))
    for target in points:
        parameters = _solve_exact_wyckoff_parameters(vectors, target)
        if parameters is None:
            continue
        representative = tuple(
            (
                vectors[0][axis]
                + sum(
                    active_vectors[index][axis] * parameters[index]
                    for index in range(len(parameters))
                )
            )
            % 1
            for axis in range(3)
        )
        orbit = frozenset(
            tuple(
                (
                    fraction_row_multiply3(representative, rotation)[axis]
                    + translation[axis]
                )
                % 1
                for axis in range(3)
            )
            for rotation, translation in operations
        )
        if orbit == points:
            return True
    return False


def _ordinary_layout_in_public_setting(
    *,
    child_sg: int,
    setting_ids: tuple[int, ...],
    setting_actions: tuple[tuple[int, frozenset[AffineOperation]], ...],
    coordinate_matrix: ExactMatrix3,
    coordinate_origin: ExactVector3,
    layout: ChildAtomLayout,
    label_correspondence: dict[str, str],
) -> ChildAtomLayout:
    """Name each Formula15 partition from its complete exact Source orbit."""

    data = _assembled_data()
    if not setting_ids:
        raise ValueError("child presentation has no Source inter setting")
    atom_by_id = {atom.atom_id: atom for atom in layout.atoms}
    if len(atom_by_id) != len(layout.atoms):
        raise ValueError("child atom layout contains duplicate atom identities")
    operations_by_setting = dict(setting_actions)
    if len(operations_by_setting) != len(setting_actions) or set(
        operations_by_setting
    ) != set(setting_ids):
        raise ValueError("child presentation setting actions are incomplete")
    labels: list[str] = []
    for site in layout.sites:
        source_site = site.wyckoff_site
        known_label = label_correspondence.get(source_site)
        if known_label is not None:
            labels.append(known_label)
            continue
        source_label = re.fullmatch(r"\s*\d+\s*(\S+)\s*", source_site)
        if source_label is None:
            raise ValueError(f"invalid Formula15 Wyckoff site: {source_site!r}")
        source_site_pg = int(
            data.wyckoff_row_by_label(int(child_sg), source_label.group(1)).site_pg
        )
        points = frozenset(
            tuple(
                (
                    fraction_row_multiply3(
                        atom_by_id[atom_id].structure_xyz,
                        coordinate_matrix,
                    )[axis]
                    + coordinate_origin[axis]
                )
                % 1
                for axis in range(3)
            )
            for atom_id in site.atom_ids
        )
        if len(points) != len(site.atom_ids):
            raise ValueError("Formula15 child site contains duplicate points")
        setting_labels: list[str] = []
        for setting_id in setting_ids:
            matches = tuple(
                row
                for row in data.wyckoff_rows(int(child_sg))
                if data.wyckoff_multiplicity(int(child_sg), row) == len(points)
                and int(row.site_pg) == source_site_pg
                and _source_wyckoff_orbit_matches(
                    vectors=data.inter_wyckoff_fraction_vectors(
                        int(child_sg), row, int(setting_id)
                    ),
                    operations=operations_by_setting[setting_id],
                    points=points,
                )
            )
            if len(matches) != 1:
                raise ValueError(
                    "Source child-setting Wyckoff identity is not unique for "
                    f"{site.child_site_id} {site.wyckoff_site}: "
                    f"setting={setting_id}, rows={[row.row_id for row in matches]}"
                )
            label = str(matches[0].label)
            if not label:
                raise ValueError("Source child-setting Wyckoff identity is incomplete")
            setting_labels.append(f"{len(points)}{label}")
        if len(set(setting_labels)) != 1:
            raise ValueError(
                "equivalent Source child settings disagree on one Wyckoff identity: "
                f"settings={setting_ids}, source={source_site}"
            )
        label_correspondence[source_site] = setting_labels[0]
        labels.append(setting_labels[0])
    return ChildAtomLayout(
        sites=tuple(
            replace(site, wyckoff_site=label)
            for site, label in zip(layout.sites, labels, strict=True)
        ),
        atoms=layout.atoms,
    )


def _magnetic_child_atom_layout_for_site(
    *,
    magnetic_group: int,
    child_sg: int,
    label_prefix: str,
    layout: ChildAtomLayout,
    presentation_positions: dict[str, tuple[float, float, float]],
    source_basis: list[list[float]],
    source_origin: Any,
) -> ChildAtomLayout:
    """Regroup canonical atoms by the exact selected magnetic group action."""

    atom_index_by_id = {atom.atom_id: index for index, atom in enumerate(layout.atoms)}
    action = selected_magnetic_atom_action(
        magnetic_group=int(magnetic_group),
        child_sg=int(child_sg),
        parent_points=tuple(atom.parent_xyz for atom in layout.atoms),
        ordinary_orbits=tuple(
            tuple(atom_index_by_id[atom_id] for atom_id in site.atom_ids)
            for site in layout.sites
        ),
        selected_basis=source_basis,
        selected_origin=source_origin,
    )
    if set(presentation_positions) != set(atom_index_by_id):
        raise ValueError("magnetic presentation positions do not cover canonical atoms")
    displayed_points = tuple(
        tuple(
            _input_fraction(value) % 1 for value in presentation_positions[atom.atom_id]
        )
        for atom in layout.atoms
    )
    if len(set(displayed_points)) != len(displayed_points):
        raise ValueError("magnetic presentation atoms share one position")
    components = action.components
    site_index_by_atom_id = {
        atom_id: site_index
        for site_index, site in enumerate(layout.sites)
        for atom_id in site.atom_ids
    }
    regrouping = []
    presentation_branches = []
    for component in components:
        component_atom_ids = tuple(layout.atoms[index].atom_id for index in component)
        site_indices = tuple(
            dict.fromkeys(
                site_index_by_atom_id[atom_id] for atom_id in component_atom_ids
            )
        )
        site_atom_ids = tuple(
            atom_id
            for site_index in site_indices
            for atom_id in layout.sites[site_index].atom_ids
        )
        if set(site_atom_ids) != set(component_atom_ids):
            raise ValueError("magnetic atom component splits an ordinary child site")
        seed_index = min(
            component,
            key=lambda index: (
                layout.atoms[index].source_raw_index,
                layout.atoms[index].centering_ordinal,
            ),
        )
        seed_point = action.canonical_cinter_points[seed_index]
        row, standard = identify_magnetic_wyckoff(int(magnetic_group), seed_point)
        orbit = magnetic_orbit_points(int(magnetic_group), seed_point)
        component_points = {
            action.canonical_cinter_points[index] for index in component
        }
        if len(orbit) != len(component) or set(orbit) != component_points:
            raise ValueError("canonical magnetic orbit does not equal atom membership")
        target = tuple(Fraction(value) % 1 for value in standard)
        matches = [
            index
            for index in component
            if action.canonical_cinter_points[index] == target
        ]
        if len(matches) != 1:
            raise ValueError(
                "magnetic representative does not identify one canonical atom: "
                f"matches={matches}"
            )
        displayed_representative = displayed_points[matches[0]]
        regrouping.append((site_indices, len(component), displayed_representative))
        presentation_branches.append(
            (str(row.label), displayed_representative, len(component))
        )
    presentation_labels = presentation_branch_labels_from_correspondence(
        magnetic_group=int(magnetic_group),
        branches=presentation_branches,
    )
    if presentation_labels and len(presentation_labels) != len(regrouping):
        raise ValueError("magnetic presentation branch map is incomplete")
    # Display coordinates need not themselves use the standard BNS setting.
    # Canonical labels remain authoritative unless a complete relabeling exists.
    output_labels = (
        presentation_labels
        if presentation_labels
        else tuple(branch[0] for branch in presentation_branches)
    )
    if len(output_labels) != len(regrouping):
        raise ValueError("magnetic child sites have no Wyckoff branch labels")
    return regroup_child_atom_layout(
        layout,
        label_prefix=label_prefix,
        groups=tuple(
            (
                site_indices,
                f"{multiplicity}{output_label}",
                displayed_representative,
            )
            for (
                site_indices,
                multiplicity,
                displayed_representative,
            ), output_label in zip(regrouping, output_labels, strict=True)
        ),
    )


def _child_centering_parent_translation_record(
    *,
    child_sg: int,
    presentation_basis_pml: tuple[int, ...],
    centering_ordinal: int,
) -> tuple[int, int, int, int]:
    child_space_group = _exact_space_group_number(child_sg, "child space group")
    if isinstance(centering_ordinal, bool) or not isinstance(
        centering_ordinal, Integral
    ):
        raise TypeError("child centering ordinal must be an exact integer")
    exact_centering_ordinal = int(centering_ordinal)
    translations = _centering_translations(child_space_group)
    if translations is None or not 0 <= exact_centering_ordinal < len(translations):
        raise ValueError(f"child centering ordinal out of range: {centering_ordinal}")
    basis_values = _exact_integer_tuple(
        presentation_basis_pml, 9, "presentation PML basis"
    )
    translation = translations[exact_centering_ordinal]
    parent_pml = tuple(
        sum(translation[row] * basis_values[3 * row + column] for row in range(3))
        for column in range(3)
    )
    denominator = math.lcm(*(value.denominator for value in parent_pml), 1)
    return tuple(int(value * denominator) for value in parent_pml) + (denominator,)  # type: ignore[return-value]


def _child_atom_layout_operation_records(
    decoder: ModeDataDecoder,
    *,
    child_sg: int,
    presentation_basis_pml: tuple[int, ...],
    layout: ChildAtomLayout,
) -> tuple[tuple[int, int, int, int, int], ...]:
    child_space_group = _exact_space_group_number(child_sg, "child space group")
    point_operation_count = len(decoder.iso.space["ipoint_op_inverse"])
    records: list[tuple[int, int, int, int, int]] = []
    for atom in layout.atoms:
        source_record = decoder.add_translation_to_operation_record(
            atom.source_kernel_fraction,
            atom.source_parent_coset_record,
        )
        centering_record = _child_centering_parent_translation_record(
            child_sg=child_space_group,
            presentation_basis_pml=presentation_basis_pml,
            centering_ordinal=atom.centering_ordinal,
        )
        records.append(
            exact_operation_record(
                decoder.add_translation_to_operation_record(
                    centering_record, source_record
                ),
                point_operation_count=point_operation_count,
            )
        )
    return tuple(records)


def _mode_row_assignments(
    decoder: ModeDataDecoder,
    *,
    child_sg: int,
    presentation_basis_pml: tuple[int, ...],
    layout: ChildAtomLayout,
    rows: list[dict[str, Any]],
) -> tuple[tuple[int, ...], ...]:
    """Map presented atom rows to canonical atoms by affine identity."""

    canonical_records = _child_atom_layout_operation_records(
        decoder,
        child_sg=child_sg,
        presentation_basis_pml=presentation_basis_pml,
        layout=layout,
    )
    point_operation_count = len(decoder.iso.space["ipoint_op_inverse"])
    row_records: list[tuple[int, int, int, int, int]] = []
    for row in rows:
        raw_record = row.get("_operation_record")
        if not isinstance(raw_record, (list, tuple)):
            raise ValueError("mode atom row is missing its Source operation record")
        record = exact_operation_record(
            raw_record,
            point_operation_count=point_operation_count,
        )
        centering_ordinal = row.get("_presentation_centering_ordinal")
        if centering_ordinal is not None:
            if isinstance(centering_ordinal, bool) or not isinstance(
                centering_ordinal, int
            ):
                raise TypeError("presentation centering ordinal must be an integer")
            centering_record = _child_centering_parent_translation_record(
                child_sg=child_sg,
                presentation_basis_pml=presentation_basis_pml,
                centering_ordinal=centering_ordinal,
            )
            record = exact_operation_record(
                decoder.add_translation_to_operation_record(centering_record, record),
                point_operation_count=point_operation_count,
            )
        row_records.append(record)  # type: ignore[arg-type]
    return operation_record_row_assignments(
        canonical_records,
        row_records,
        presentation_basis_pml=presentation_basis_pml,
        point_operation_count=point_operation_count,
    )


def _compile_child_atom_mode_topology(
    decoder: ModeDataDecoder,
    *,
    child_sg: int,
    presentation_basis_pml: tuple[int, ...],
    layout: ChildAtomLayout,
    atom_operation_records: list[Any],
) -> dict[tuple[int, int, tuple[int, int, int, int, int]], str]:
    """Compile every Source row/centering representation to one child atom."""

    child_space_group = _exact_space_group_number(child_sg, "child space group")
    centering_count = len(_centering_translations(child_space_group) or ())
    if centering_count <= 0 or not atom_operation_records:
        raise ValueError("mode topology requires Source rows and centerings")
    point_operation_count = len(decoder.iso.space["ipoint_op_inverse"])
    rows: list[dict[str, Any]] = []
    identities: list[tuple[int, int, tuple[int, int, int, int, int]]] = []
    for source_raw_index, raw_record in enumerate(atom_operation_records):
        if not isinstance(raw_record, (list, tuple)):
            raise TypeError("mode topology operation record must be a sequence")
        operation_record = exact_operation_record(
            raw_record,
            point_operation_count=point_operation_count,
        )
        for centering_ordinal in range(centering_count):
            rows.append(
                {
                    "_operation_record": operation_record,
                    "_presentation_centering_ordinal": centering_ordinal,
                }
            )
            identities.append((source_raw_index, centering_ordinal, operation_record))
    assignments = _mode_row_assignments(
        decoder,
        child_sg=child_space_group,
        presentation_basis_pml=presentation_basis_pml,
        layout=layout,
        rows=rows,
    )
    topology: dict[tuple[int, int, tuple[int, int, int, int, int]], str] = {}
    for atom, candidate_indices in zip(layout.atoms, assignments, strict=True):
        for candidate_index in candidate_indices:
            identity = identities[candidate_index]
            previous = topology.setdefault(identity, atom.atom_id)
            if previous != atom.atom_id:
                raise ValueError("mode topology identity maps to multiple child atoms")
    if len(topology) != len(identities):
        raise ValueError("mode topology does not resolve every Source row")
    if set(topology.values()) != {atom.atom_id for atom in layout.atoms}:
        raise ValueError("mode topology does not cover every canonical child atom")
    return topology


def _present_child_atom_layout(
    decoder: ModeDataDecoder,
    *,
    parent_sg: int,
    child_sg: int,
    parent_wyckoff: str,
    site_params: tuple[float, ...] | None,
    presentation_basis: list[list[float]],
    presentation_basis_pml: tuple[int, ...],
    presentation_origin: Any,
    child_symbol: str,
    label_prefix: str,
    layout: ChildAtomLayout,
    parent_setting_bridge: tuple[Any, Any] | None,
) -> tuple[
    ChildAtomLayout,
    dict[str, tuple[float, float, float]],
]:
    """Bind Formula15 identities to the Source display atom table once."""

    source_rows = [
        row
        for row in decoder.wyckoff_rows(int(parent_sg))
        if row.label == parent_wyckoff
    ]
    if len(source_rows) != 1:
        raise ValueError(
            f"expected one Source Wyckoff row for SG{parent_sg} "
            f"{parent_wyckoff}, got {len(source_rows)}"
        )
    source_row = source_rows[0]
    raw_fractionals = decoder.display_distortion_atom_fractionals(
        int(parent_sg),
        source_row,
        site_params,
        presentation_basis_pml,
    )
    raw_records = decoder.supercell_atom_operation_records(
        int(parent_sg),
        source_row,
        presentation_basis_pml,
    )
    if len(raw_fractionals) != len(raw_records) or not raw_fractionals:
        raise ValueError("Source atom positions and operation records differ")
    if parent_setting_bridge is None:
        displayed_fractionals = raw_fractionals
    else:
        setting_matrix, setting_origin = parent_setting_bridge
        displayed_fractionals = tuple(
            tuple(
                sum(
                    position[row] * Fraction(setting_matrix[row][column])
                    for row in range(3)
                )
                + Fraction(setting_origin[column])
                for column in range(3)
            )
            for position in raw_fractionals
        )
    presented = present_mode_rows(
        (
            {
                "xyz": position,
                "dxyz": (0.0, 0.0, 0.0),
                "_source_raw_index": index,
                "_operation_record": raw_records[index],
            }
            for index, position in enumerate(displayed_fractionals)
        ),
        basis=presentation_basis,
        origin=_origin_vector(presentation_origin),
        centering_symbol=child_symbol,
        include_centering_ordinal=True,
    )
    rows = [dict(row) for row in presented["rows"]]
    assignments = _mode_row_assignments(
        decoder,
        child_sg=child_sg,
        presentation_basis_pml=presentation_basis_pml,
        layout=layout,
        rows=rows,
    )
    if any(len(candidate_indices) != 1 for candidate_indices in assignments):
        raise ValueError(
            "Source atom table does not bijectively cover canonical child atoms"
        )
    atom_ids_by_row: list[str | None] = [None] * len(rows)
    presentation_rows: list[ChildAtomPresentationRow | None] = [None] * len(rows)
    positions: dict[str, tuple[float, float, float]] = {}
    for atom, candidate_indices in zip(layout.atoms, assignments, strict=True):
        row_index = candidate_indices[0]
        row = rows[row_index]
        raw_xyz = row.get("xyz")
        if not isinstance(raw_xyz, (list, tuple)) or len(raw_xyz) != 3:
            raise ValueError("Source atom table contains no presented position")
        xyz = tuple(float(value) % 1.0 for value in raw_xyz)
        if not all(math.isfinite(value) for value in xyz):
            raise ValueError("Source atom table contains a nonfinite position")
        atom_ids_by_row[row_index] = atom.atom_id
        positions[atom.atom_id] = xyz  # type: ignore[assignment]
        source_raw_index = row.get("_source_raw_index")
        centering_ordinal = row.get("_presentation_centering_ordinal")
        if (
            isinstance(source_raw_index, bool)
            or not isinstance(source_raw_index, Integral)
            or int(source_raw_index) < 0
        ):
            raise TypeError("Source display row has no exact raw index")
        if (
            isinstance(centering_ordinal, bool)
            or not isinstance(centering_ordinal, Integral)
            or int(centering_ordinal) < 0
        ):
            raise TypeError("Source display row has no exact centering ordinal")
        raw_record = row.get("_operation_record")
        if not isinstance(raw_record, (list, tuple)):
            raise ValueError("Source display row has no operation record")
        presentation_rows[row_index] = ChildAtomPresentationRow(
            atom_id=atom.atom_id,
            source_raw_index=int(source_raw_index),
            centering_ordinal=int(centering_ordinal),
            operation_record=exact_operation_record(
                raw_record,
                point_operation_count=len(decoder.iso.space["ipoint_op_inverse"]),
            ),
        )
    if len(positions) != len(layout.atoms):
        raise ValueError("Source atom table contains duplicate atom identities")
    if any(atom_id is None for atom_id in atom_ids_by_row):
        raise ValueError("Source atom table contains an unassigned presentation row")
    if any(row is None for row in presentation_rows):
        raise ValueError("Source atom table contains an unbound presentation row")
    ordered_layout = child_atom_layout_in_presentation_order(
        layout,
        label_prefix=label_prefix,
        atom_ids=tuple(atom_id for atom_id in atom_ids_by_row if atom_id is not None),
        atom_positions=positions,
    )
    ordered_layout = replace(
        ordered_layout,
        presentation_rows=tuple(row for row in presentation_rows if row is not None),
        point_operation_count=len(decoder.iso.space["ipoint_op_inverse"]),
    )
    return (
        ordered_layout,
        positions,
    )


def _mode_rows_on_child_atom_layout(
    *,
    layout: ChildAtomLayout,
    mode_topology: dict[tuple[int, int, tuple[int, int, int, int, int]], str],
    rows: list[dict[str, Any]],
    mode_atom_positions: dict[str, tuple[float, float, float]],
) -> list[dict[str, Any]]:
    """Attach vectors to the Source display-row identities bound by the layout."""

    if (
        layout.point_operation_count <= 0
        or len(layout.presentation_rows) != len(layout.atoms)
        or len(rows) != len(layout.presentation_rows)
    ):
        raise ValueError("mode rows do not cover the canonical presentation layout")
    atom_by_id = {atom.atom_id: atom for atom in layout.atoms}
    if len(atom_by_id) != len(layout.atoms):
        raise ValueError("child atom layout contains duplicate atom identities")
    source_by_atom_id: dict[str, dict[str, Any]] = {}
    for source in rows:
        source_raw_index = source.get("_source_raw_index")
        raw_centering_ordinal = source.get("_presentation_centering_ordinal")
        centering_ordinal = (
            0 if raw_centering_ordinal is None else raw_centering_ordinal
        )
        raw_record = source.get("_operation_record")
        if (
            isinstance(source_raw_index, bool)
            or not isinstance(source_raw_index, Integral)
            or isinstance(centering_ordinal, bool)
            or not isinstance(centering_ordinal, Integral)
            or not isinstance(raw_record, (list, tuple))
        ):
            raise ValueError("mode row has no exact Source display-row identity")
        identity = (
            int(source_raw_index),
            int(centering_ordinal),
            exact_operation_record(
                raw_record,
                point_operation_count=layout.point_operation_count,
            ),
        )
        atom_id = mode_topology.get(identity)
        if atom_id is None or atom_id in source_by_atom_id:
            raise ValueError("mode row differs from its Source display-row identity")
        source_by_atom_id[atom_id] = source
    if set(source_by_atom_id) != set(atom_by_id):
        raise ValueError("mode rows do not bijectively cover canonical child atoms")
    site_by_id = {site.child_site_id: site for site in layout.sites}
    aligned: list[dict[str, Any]] = []
    for atom_index, atom in enumerate(layout.atoms):
        source = source_by_atom_id[atom.atom_id]
        source_fields = {
            key: value
            for key, value in source.items()
            if key != "_presentation_centering_ordinal"
        }
        raw_xyz = source.get("xyz")
        raw_dxyz = source.get("dxyz")
        if not (isinstance(raw_dxyz, (list, tuple)) and len(raw_dxyz) == 3):
            raise ValueError("mode atom row requires one vector")
        dxyz = tuple(float(value) for value in raw_dxyz)
        expected_xyz = mode_atom_positions.get(atom.atom_id)
        if expected_xyz is None:
            raise ValueError(f"canonical atom has no mode coordinate: {atom.atom_id}")
        if not all(math.isfinite(value) for value in (*expected_xyz, *dxyz)):
            raise ValueError("mode atom row contains a nonfinite position or vector")
        if raw_xyz is not None:
            if not isinstance(raw_xyz, (list, tuple)) or len(raw_xyz) != 3:
                raise ValueError("mode atom row contains an invalid position")
            mode_xyz = tuple(float(value) % 1.0 for value in raw_xyz)
            if not all(math.isfinite(value) for value in mode_xyz):
                raise ValueError("mode atom row contains a nonfinite position")
        else:
            mode_xyz = expected_xyz
        if not periodic_float_close3(expected_xyz, mode_xyz, 1e-12):
            raise ValueError(
                "mode row differs from its canonical presentation coordinate: "
                f"{atom.atom_id}: {expected_xyz!r} != {mode_xyz!r}"
            )
        aligned.append(
            {
                **source_fields,
                "atom": atom.child_site_id if atom.member_order == 0 else None,
                "atom_id": atom.atom_id,
                "child_site": atom.child_site_id,
                "wyckoff_site": site_by_id[atom.child_site_id].wyckoff_site,
                "xyz": list(expected_xyz),
                "dxyz": list(dxyz),
                "_presentation_grid_index": atom_index,
            }
        )
    return aligned


def _basis_from_opd_row(
    selected_opd: dict[str, Any] | None,
) -> list[list[float]] | None:
    if not isinstance(selected_opd, dict):
        return None
    iso = selected_opd.get("isotropy") or selected_opd
    basis = iso.get("basis") if isinstance(iso, dict) else None
    if not (
        isinstance(basis, list)
        and len(basis) == 3
        and all(isinstance(row, list) and len(row) == 3 for row in basis)
    ):
        return None
    try:
        return [[float(value) for value in row] for row in basis]
    except (TypeError, ValueError):
        return None


def _source_split_basis_from_opd_row(
    selected_opd: dict[str, Any] | None,
) -> list[list[float]] | None:
    if not isinstance(selected_opd, dict):
        return None
    iso = selected_opd.get("isotropy") or selected_opd
    basis = iso.get("source_basis_values") if isinstance(iso, dict) else None
    if not isinstance(basis, list) or len(basis) != 9:
        return None
    try:
        return [[float(basis[row * 3 + col]) for col in range(3)] for row in range(3)]
    except (TypeError, ValueError):
        return None


def _source_split_origin_from_opd_row(
    selected_opd: dict[str, Any] | None,
) -> tuple[int, int, int, int] | None:
    if not isinstance(selected_opd, dict):
        return None
    iso = selected_opd.get("isotropy") or selected_opd
    origin = iso.get("source_origin_values") if isinstance(iso, dict) else None
    if not isinstance(origin, list) or len(origin) != 4:
        return None
    try:
        return tuple(int(value) for value in origin)  # type: ignore[return-value]
    except (TypeError, ValueError):
        return None


def _point_rotation_matrix(
    data: SourceTables,
    child_sg: int,
    point_op: int,
) -> tuple[tuple[Fraction, Fraction, Fraction], ...]:
    units = (
        (Fraction(1), Fraction(0), Fraction(0)),
        (Fraction(0), Fraction(1), Fraction(0)),
        (Fraction(0), Fraction(0), Fraction(1)),
    )
    return tuple(
        tuple(data.vrot_fraction(int(child_sg), int(point_op), unit)) for unit in units
    )  # type: ignore[return-value]


def _fraction_vector_record(
    values: tuple[Fraction, Fraction, Fraction],
    point_op: int,
) -> tuple[int, int, int, int, int]:
    denominator = math.lcm(*(value.denominator for value in values), 1)
    return tuple(int(value * denominator) for value in values) + (
        denominator,
        int(point_op),
    )  # type: ignore[return-value]


def _subgroup_parent_operation_records(
    data: Any,
    parent_sg: int,
    child_sg: int,
    basis: tuple[int, ...],
    origin: tuple[int, int, int, int],
) -> tuple[tuple[int, int, int, int, int], ...]:
    """Map child operations into the selected raw parent embedding."""

    basis_matrix = _matrix_from_basis_tuple(basis)
    basis_inverse = _fraction_matrix_inverse_3(
        [[float(value) for value in row] for row in basis_matrix]
    )
    if basis_inverse is None:
        return ()
    origin_vector = _origin_record_vector(origin)
    parent_point_ops: dict[tuple[tuple[Fraction, Fraction, Fraction], ...], int] = {}
    for record in data.generate_space_group_records(int(parent_sg)):
        parent_point_ops[
            _point_rotation_matrix(data, int(parent_sg), int(record[4]))
        ] = int(record[4])

    out: list[tuple[int, int, int, int, int]] = []
    for x, y, z, den, child_point_op in data.generate_space_group_records(
        int(child_sg)
    ):
        rotation = _point_rotation_matrix(data, int(child_sg), int(child_point_op))
        translation = (
            Fraction(int(x), int(den)),
            Fraction(int(y), int(den)),
            Fraction(int(z), int(den)),
        )
        parent_rotation = _fraction_matmul(
            basis_inverse, _fraction_matmul(rotation, basis_matrix)
        )
        point_op = parent_point_ops.get(parent_rotation)
        if point_op is None:
            raise ValueError(
                f"child SG{child_sg} operation {child_point_op} does not map to a point operation in SG{parent_sg}"
            )
        parent_translation = _fraction_vecadd(
            _fraction_vecsub(
                origin_vector, _fraction_row_multiply(origin_vector, parent_rotation)
            ),
            _fraction_row_multiply(translation, basis_matrix),
        )
        out.append(_fraction_vector_record(parent_translation, point_op))
    return tuple(out)


def _split_basis_origin_for_wyckoff(
    *,
    selected_opd: dict[str, Any] | None,
) -> tuple[list[list[float]], tuple[int, int, int, int]]:
    """Return the transform consumed by Source-only ``get_new_wyckoff_`` splitting.

    Complete-mode Wyckoff splitting follows ISO's raw ``data_isotropy``
    subgroup basis/origin, not the public/VALUE CELL presentation basis.  The
    two transforms can describe the same subgroup physically while only the raw
    transform closes the child Wyckoff orbit buffer used by ``get_new_wyckoff_``.
    """

    basis, origin = _selected_source_subgroup_transform(selected_opd)
    return (
        [[float(basis[row * 3 + column]) for column in range(3)] for row in range(3)],
        origin,
    )


def _selected_subgroup_number(selected_opd: dict[str, Any] | None) -> int | None:
    if not isinstance(selected_opd, dict):
        return None
    iso = selected_opd.get("isotropy") or selected_opd
    if not isinstance(iso, dict):
        return None
    subgroup = iso.get("subgroup")
    if isinstance(subgroup, dict):
        subgroup = subgroup.get("ordinary_number") or subgroup.get("number")
    try:
        return int(subgroup)
    except (TypeError, ValueError):
        return None


def _selected_magnetic_group_number(
    selected_opd: dict[str, Any] | None,
    child_sg: int | None,
) -> int | None:
    iso = _isotropy_from_opd_row(selected_opd)
    if not iso:
        return None
    subgroup = iso.get("subgroup")
    if not isinstance(subgroup, dict) or subgroup.get("ordinary_number") is None:
        return None
    try:
        group = int(subgroup["number"])
        setting = magnetic_group_setting(group)
    except (IndexError, KeyError, TypeError, ValueError):
        return None
    if child_sg is None or int(setting.ordinary_space_group) != int(child_sg):
        return None
    return group


def _presentation_basis_candidate(
    basis: list[list[float]] | None,
) -> tuple[list[list[float]] | None, str]:
    if basis is None:
        return None, "none"
    identity = ((1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, 1.0))
    is_identity = all(
        abs(float(basis[row][col]) - identity[row][col]) <= 1e-12
        for row in range(3)
        for col in range(3)
    )
    return basis, "identity" if is_identity else "child_basis"


def _basis_cinter_to_pml(
    decoder: ModeDataDecoder,
    sg: int,
    basis_cinter: list[list[float]] | None,
    selected_setting_id: int | None = None,
) -> tuple[int, ...] | None:
    if basis_cinter is None:
        return None
    rows_list: list[tuple[Fraction, Fraction, Fraction]] = []
    for row in basis_cinter:
        values: list[Fraction] = []
        for value in row:
            raw = Fraction(str(value))
            snapped = raw.limit_denominator(384)
            values.append(snapped if abs(float(raw - snapped)) <= 5e-9 else raw)
        rows_list.append(tuple(values))
    rows = tuple(rows_list)
    setting_bridge = _parent_setting_bridge(int(sg), selected_setting_id)
    if setting_bridge is not None:
        selected_to_default = fraction_matrix_inverse3(setting_bridge[0])
        rows = fraction_matrix_multiply3(rows, selected_to_default)
    inverse = fraction_matrix_inverse3(decoder.pml_to_cinter_matrix(sg))
    basis_pml = fraction_matrix_multiply3(rows, inverse)
    if any(value.denominator != 1 for row in basis_pml for value in row):
        return None
    return tuple(int(value) for row in basis_pml for value in row)


def _presentation_mode_vectors(
    modes: list[Any],
    basis: list[list[float]] | None,
) -> list[list[list[float]]]:
    inverse = _float_matrix_inverse_3(basis) if basis is not None else None
    return [
        [
            (
                [float(value) for value in vector[:3]]
                if inverse is None
                else _row_multiply([float(value) for value in vector[:3]], inverse)
            )
            for vector in mode
        ]
        for mode in modes
        if isinstance(mode, list)
    ]
