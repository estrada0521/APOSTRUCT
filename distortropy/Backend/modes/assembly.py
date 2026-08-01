"""Assemble mode definitions and their undistorted structures."""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from fractions import Fraction
from typing import Any
import gemmi
import numpy as np
from distortropy.Backend.modes.engine.input import Case
from distortropy.Backend.modes.engine.dynamic_subduction import (
    _strict_integral_values,
)
from distortropy.Backend.modes.engine.trace.pipeline_trace import (
    project_mode_block,
)
from distortropy.Backend.modes.engine.project.display_distortion import (
    occurrence_sublattice,
)
from distortropy.Backend.modes.presentation import (
    child_lattice_cartesian,
    present_mode_rows,
)
from distortropy.Backend.modes.viewer import parent_cell_placement
from distortropy.Backend.parameter_names import parameter_name

from distortropy.Backend.modes.definition_presentation import (
    _apply_dynamic_source_family_presentation,
    _complete_mode_normfactor,
    _mode_label,
    _orthogonalize_definition_modes,
    _strain_mode_definitions,
)
from distortropy.Backend.modes.common import (
    _source_tables,
    _cell_from_basis,
    _cell_params,
    _freeparam_from_opd_row,
    _integer_basis_tuple,
    _is_full_parameter_opd,
    _isotropy_from_opd_row,
    _mapped_sites,
    _mode_decoder,
    _origin_from_opd_row,
    _origin_vector,
    _site_params,
)
from distortropy.Backend.modes.structure_runtime import (
    _basis_cinter_to_pml,
    _basis_from_opd_row,
    _child_presentation_setting,
    _compile_child_atom_mode_topology,
    _magnetic_child_atom_layout_for_site,
    _mode_rows_on_child_atom_layout,
    _ordinary_layout_in_public_setting,
    _present_child_atom_layout,
    _presentation_basis_candidate,
    _presentation_mode_vectors,
    _selected_magnetic_group_number,
    _selected_subgroup_number,
    _selected_source_subgroup_transform,
    _source_split_basis_from_opd_row,
    _source_split_origin_from_opd_row,
    _source_child_atom_layout_for_site,
    _split_basis_origin_for_wyckoff,
    _web_magnetic_occurrence_gauge,
)
from distortropy.Backend.modes.structure.child_atom_layout import (
    ChildAtomLayout,
    child_atom_layout_in_presentation_order,
    exact_operation_record,
)
from distortropy.Backend.modes.request_context import (
    _case_k_params,
    _evaluated_kvector_text,
    _little_star_arm_width,
    _opd_direction_text,
    _selected_dynamic_gid,
    _should_emit_all_source_opd_rows,
    _spec_opd_direction_text,
    _spec_display_kvector,
)
from distortropy.Backend.modes.subduction_specs import (
    _coupled_render_specs,
    _magnetic_subduced_mode_specs,
    _subduced_mode_specs,
)
from distortropy.Backend.modes.print_layout import (
    _mode_site_irrep_labels,
    _mode_source_layout,
    _mode_source_metadata_layout,
)
from distortropy.Backend.modes.print_intertwiner import (
    _rank1_source_print_intertwiner,
)
from distortropy.Backend.modes.occurrence_aliases import (
    OccurrenceSiteEmission,
    admitted_occurrence_alias_spec_orders,
    admitted_self_complete_occurrence_alias_spec_orders,
    direction_distinct_occurrence_alias_spec_orders,
)
from distortropy.Backend.modes.site_transport import (
    _parent_setting_bridge,
    _source_default_wyckoff_params,
)


def _secondary_invariant_factors(
    *,
    render_specs: list[tuple[dict[str, Any], str, int, str]],
    emitted_spec_orders: set[int],
    primary_slot_count: int,
    coupled: bool,
    magnetic_primary: bool,
) -> list[dict[str, Any]]:
    """Expose rendered static secondary subspaces for Source-domain resolution."""

    factors: list[dict[str, Any]] = []
    for spec_order, (spec, mode_kind, _setting, _display_label) in enumerate(
        render_specs
    ):
        if (
            spec_order not in emitted_spec_orders
            or spec.get("primary") is not False
            or int(spec.get("old_id") or 0) <= 0
        ):
            continue
        gid = spec.get("gid")
        if type(gid) is not int or gid <= 0:
            raise ValueError("rendered static secondary lost its Source irrep identity")
        factor = {
            "slot": int(primary_slot_count) + len(factors) + 1,
            "gid": gid,
            "label": str(spec.get("display_label") or spec.get("label") or ""),
            "magnetic": mode_kind == "mag",
            "k_parameters": [],
            "role": "secondary",
        }
        source_free_count = spec.get("source_free_count")
        if (
            not isinstance(source_free_count, bool)
            and isinstance(source_free_count, int)
            and source_free_count > 0
        ):
            factor["parameter_count"] = source_free_count
        row_id = spec.get("row_id")
        domain = spec.get("domain")
        opd = spec.get("opd")
        if (
            type(row_id) is int
            and row_id > 0
            and type(domain) is int
            and domain > 0
            and isinstance(opd, str)
            and opd
        ):
            factor.update({"opd": opd, "domain": domain})
        else:
            direction_matrix = spec.get("direction_matrix")
            if (
                not coupled
                and not magnetic_primary
                and mode_kind == "dsp"
            ):
                raise ValueError(
                    "rendered static secondary lost its Source OPD/domain identity"
                )
            if not isinstance(direction_matrix, list) or not direction_matrix:
                raise ValueError("rendered static secondary lost its direction subspace")
            factor["direction_matrix"] = [list(row) for row in direction_matrix]
        factors.append(factor)
    return factors


def _attach_definition_factor_metadata(
    definitions: Sequence[dict[str, Any]],
    *,
    render_specs: Sequence[tuple[dict[str, Any], str, int, str]],
    selected_factors: Sequence[Mapping[str, Any]],
) -> None:
    """Attach stable primary-factor handles without parsing display labels."""

    primary_factors: list[dict[str, Any]] = []
    parameter_offset = 0
    for raw_factor in selected_factors:
        parameter_count = raw_factor["parameter_count"]
        factor = {key: raw_factor[key] for key in ("slot", "gid")}
        factor["parameters"] = [
            parameter_name(parameter_offset + index)
            for index in range(parameter_count)
        ]
        parameter_offset += parameter_count
        primary_factors.append(factor)

    for definition in definitions:
        spec_order = definition.get("_spec_order")
        if type(spec_order) is not int or not 0 <= spec_order < len(render_specs):
            continue
        spec = render_specs[spec_order][0]
        if spec.get("primary") is not True:
            definition["role"] = "secondary"
            continue
        definition["role"] = "primary"
        primary_order = spec.get("_primary_slot_order")
        if primary_order is None and len(primary_factors) == 1:
            primary_order = 0
        if (
            type(primary_order) is int
            and 0 <= primary_order < len(primary_factors)
        ):
            definition["factor"] = dict(primary_factors[primary_order])


def _exact_source_kparam(value: Any) -> tuple[int, int, int, int] | None:
    try:
        if isinstance(value, (str, bytes, bytearray)):
            return None
        values = tuple(Fraction(str(item)) for item in value)
    except (TypeError, ValueError, ZeroDivisionError):
        return None
    if len(values) != 4 or any(item.denominator != 1 for item in values):
        return None
    values = tuple(int(item) for item in values)
    if values[3] == 0:
        return None
    return values  # type: ignore[return-value]


def _freeze_occurrence_source_value(value: Any) -> Any:
    if isinstance(value, Mapping):
        return tuple(
            (str(key), _freeze_occurrence_source_value(value[key]))
            for key in sorted(value, key=lambda item: str(item))
        )
    if isinstance(value, Sequence) and not isinstance(value, (str, bytes, bytearray)):
        return tuple(_freeze_occurrence_source_value(item) for item in value)
    if isinstance(value, float) and not np.isfinite(value):
        raise ValueError("non-finite Source identity")
    hash(value)
    return value


def _freeze_integral_record(record: Any) -> tuple[int, ...]:
    try:
        if isinstance(record, (str, bytes, bytearray)):
            raise ValueError
        values = tuple(Fraction(str(item)) for item in record)
    except (TypeError, ValueError, ZeroDivisionError):
        raise ValueError("malformed Source operation record") from None
    if (
        len(values) != 5
        or any(item.denominator != 1 for item in values)
        or values[3] <= 0
        or values[4] <= 0
    ):
        raise ValueError("malformed Source operation record")
    return tuple(int(item) for item in values)


def _freeze_fractional_xyz(value: Any) -> tuple[tuple[int, int], ...]:
    try:
        if isinstance(value, (str, bytes, bytearray)):
            raise ValueError
        values = tuple(Fraction(str(item)) for item in value)
    except (TypeError, ValueError, ZeroDivisionError):
        raise ValueError("malformed Source fractional coordinate") from None
    if len(values) != 3:
        raise ValueError("malformed Source fractional coordinate")
    return tuple(
        (item.numerator % item.denominator, item.denominator) for item in values
    )


def _occurrence_source_signature(
    source_metadata_layout: Sequence[Mapping[str, Any]] | None,
    mode_count: int,
) -> tuple[Any, ...] | None:
    """Freeze every emitted Source column identity without reordering it."""

    if (
        source_metadata_layout is None
        or int(mode_count) <= 0
        or len(source_metadata_layout) != int(mode_count)
    ):
        return None
    signature: list[Any] = []
    try:
        for identity in source_metadata_layout:
            opd_groups = identity.get("opd_groups")
            if not isinstance(opd_groups, list) or not opd_groups:
                return None
            required = _strict_integral_values(
                [
                    identity[field]
                    for field in (
                        "gid",
                        "site_pg",
                        "pg_irrep",
                        "family",
                        "component",
                        "print_component",
                        "family_width",
                        "source_row_count",
                        "little_type",
                    )
                ],
                size=9,
            )
            if required is None:
                return None
            (
                gid,
                site_pg,
                pg_irrep,
                family,
                component,
                print_component,
                family_width,
                source_row_count,
                little_type,
            ) = required
            if (
                gid <= 0
                or site_pg <= 0
                or pg_irrep <= 0
                or family_width <= 0
                or source_row_count <= 0
                or little_type <= 0
                or family < 0
                or component < 0
                or print_component < 0
            ):
                return None
            signature.append(
                (
                    *required,
                    _freeze_occurrence_source_value(opd_groups),
                )
            )
    except (KeyError, TypeError, ValueError):
        return None
    return tuple(signature)


def _occurrence_modes_are_nonzero(mode_vectors: Sequence[Any]) -> bool:
    if not mode_vectors:
        return False
    for vectors in mode_vectors:
        try:
            values = np.asarray(vectors, dtype=float)
        except (TypeError, ValueError):
            return False
        if (
            values.size == 0
            or not np.all(np.isfinite(values))
            or float(np.linalg.norm(values)) <= 1e-12
        ):
            return False
    return True


def _occurrence_alias_spec_orders(
    render_specs: Sequence[tuple[dict[str, Any], str, int, str]],
) -> tuple[frozenset[int], dict[int, int]]:
    """Join speculative aliases to one exact non-alias representative spec."""

    candidate_orders = frozenset(
        index
        for index, (spec, mode_kind, _setting, _label) in enumerate(render_specs)
        if spec.get("_occurrence_alias_anchor") is not None
    )
    anchors: dict[tuple[str, int, tuple[int, int, int, int]], list[int]] = {}
    for index, (spec, mode_kind, _setting, _label) in enumerate(render_specs):
        if index in candidate_orders:
            continue
        source_kparam = _exact_source_kparam(spec.get("source_kparam"))
        try:
            gid = int(spec.get("gid") or 0)
        except (TypeError, ValueError):
            continue
        if gid <= 0 or source_kparam is None:
            continue
        anchors.setdefault((mode_kind, gid, source_kparam), []).append(index)
    candidate_to_anchor: dict[int, int] = {}
    for candidate_order in candidate_orders:
        spec = render_specs[candidate_order][0]
        raw_anchor = spec.get("_occurrence_alias_anchor")
        try:
            gid = int(raw_anchor[0])
            source_kparam = _exact_source_kparam(raw_anchor[1])
        except (IndexError, TypeError, ValueError):
            continue
        mode_kind = render_specs[candidate_order][1]
        matches = (
            anchors.get((mode_kind, gid, source_kparam), ())
            if source_kparam is not None
            else ()
        )
        if len(matches) == 1:
            candidate_to_anchor[candidate_order] = matches[0]
    return candidate_orders, candidate_to_anchor


def _build_child_structure(
    *,
    cif_info: dict[str, Any],
    selected_opd: dict[str, Any] | None = None,
    parent_inter_setting_id: int | None = None,
    decoder: Any,
) -> tuple[
    dict[str, Any],
    tuple[ChildAtomLayout, ...],
    dict[str, tuple[float, float, float]],
    tuple[tuple[float, ...] | None, ...],
]:
    """Build the one child-site/atom identity shared by every output surface."""

    sites = _mapped_sites(cif_info)
    if not sites:
        raise ValueError("mode kernel requires at least one mapped displacive site")
    sg = int(cif_info["parent"]["number"])
    child_sg = _selected_subgroup_number(selected_opd)
    if child_sg is None:
        raise ValueError("selected subgroup has no ordinary space-group identity")
    parent_cell = _cell_params(cif_info)
    display_basis = _basis_from_opd_row(selected_opd)
    presentation_basis, presentation_rule = _presentation_basis_candidate(display_basis)
    if presentation_basis is None:
        raise ValueError("selected subgroup has no presentation basis")
    child_origin = _origin_from_opd_row(selected_opd)
    presentation_basis_pml = _basis_cinter_to_pml(
        decoder,
        int(sg),
        presentation_basis,
        parent_inter_setting_id,
    )
    if presentation_basis_pml is None:
        raise ValueError("selected subgroup has no exact presentation basis")
    split_basis, split_origin = _split_basis_origin_for_wyckoff(
        selected_opd=selected_opd,
    )
    magnetic_group = _selected_magnetic_group_number(selected_opd, child_sg)
    child_symbol = gemmi.find_spacegroup_by_number(int(child_sg)).hm
    parent_setting_bridge = _parent_setting_bridge(int(sg), parent_inter_setting_id)
    source_basis_record, source_origin_record = _selected_source_subgroup_transform(
        selected_opd
    )
    (
        child_setting_ids,
        child_coordinate_matrix,
        child_coordinate_origin,
        child_setting_actions,
    ) = _child_presentation_setting(
        data=_source_tables(),
        parent_sg=int(sg),
        child_sg=int(child_sg),
        parent_setting_id=parent_inter_setting_id,
        source_basis=source_basis_record,
        source_origin=source_origin_record,
        presentation_basis_pml=presentation_basis_pml,
        presentation_origin=child_origin,
    )
    undistorted_atoms: list[dict[str, Any]] = []
    layouts: list[ChildAtomLayout] = []
    mode_atom_positions: dict[str, tuple[float, float, float]] = {}
    site_params_by_site: list[tuple[float, ...] | None] = []
    public_wyckoff_labels: dict[str, str] = {}
    symmetry_operations = [
        str(value) for value in (cif_info.get("symmetry_operations") or [])
    ]
    for site_index, site in enumerate(sites):
        atom_label = str(site.get("type") or site.get("label") or "X")
        label_prefix = str(site.get("label") or atom_label)
        parent_params = _source_default_wyckoff_params(
            int(sg),
            site,
            parent_inter_setting_id,
            symmetry_operations,
        )
        site_params = (
            _site_params({"wyckoff_params": parent_params}) if parent_params else None
        )
        layout = _source_child_atom_layout_for_site(
            sg=int(sg),
            child_sg=int(child_sg),
            site=site,
            label_prefix=label_prefix,
            split_basis=split_basis,
            split_origin=split_origin,
            parent_params=parent_params,
        )
        layout = _ordinary_layout_in_public_setting(
            child_sg=int(child_sg),
            setting_ids=child_setting_ids,
            setting_actions=child_setting_actions,
            coordinate_matrix=child_coordinate_matrix,
            coordinate_origin=child_coordinate_origin,
            layout=layout,
            label_correspondence=public_wyckoff_labels,
        )
        layout, site_mode_positions = _present_child_atom_layout(
            decoder,
            parent_sg=int(sg),
            child_sg=int(child_sg),
            parent_wyckoff=str(site["wyckoff"]),
            site_params=site_params,
            presentation_basis=presentation_basis,
            presentation_basis_pml=presentation_basis_pml,
            presentation_origin=child_origin,
            child_symbol=child_symbol,
            label_prefix=label_prefix,
            layout=layout,
            parent_setting_bridge=parent_setting_bridge,
        )
        if magnetic_group is not None:
            if split_basis is None:
                raise ValueError("magnetic child layout has no Source subgroup basis")
            layout = _magnetic_child_atom_layout_for_site(
                magnetic_group=int(magnetic_group),
                child_sg=int(child_sg),
                label_prefix=label_prefix,
                layout=layout,
                presentation_positions=site_mode_positions,
                source_basis=split_basis,
                source_origin=split_origin,
            )
            layout = child_atom_layout_in_presentation_order(
                layout,
                label_prefix=label_prefix,
                atom_ids=tuple(row.atom_id for row in layout.presentation_rows),
                atom_positions=site_mode_positions,
            )
        overlap = set(mode_atom_positions).intersection(site_mode_positions)
        if overlap:
            raise ValueError(
                f"canonical atom identities cross parent sites: {sorted(overlap)[:8]!r}"
            )
        mode_atom_positions.update(site_mode_positions)
        layouts.append(layout)
        site_params_by_site.append(site_params)
        undistorted_atoms.extend(
            {
                "label": child_site.child_site_id,
                "child_site": child_site.child_site_id,
                "site": child_site.wyckoff_site,
                "multiplicity": len(child_site.atom_ids),
                "atom_ids": list(child_site.atom_ids),
                "xyz": [float(value) for value in child_site.representative_xyz],
                "_presentation_orbit_points": [
                    list(mode_atom_positions[atom_id]) for atom_id in child_site.atom_ids
                ],
            }
            for child_site in layout.sites
        )
    return (
        {
            "status": "ok",
            "lattice": _cell_from_basis(parent_cell, presentation_basis),
            "subgroup_details": {
                **(_isotropy_from_opd_row(selected_opd) or {}),
                "parent_inter_setting_id": parent_inter_setting_id,
            },
            "undistorted_atoms": undistorted_atoms,
            "_presentation_rule": presentation_rule,
        },
        tuple(layouts),
        mode_atom_positions,
        tuple(site_params_by_site),
    )


def build_mode_details(
    *,
    cif_info: dict[str, Any],
    selected_k: dict[str, Any],
    selected_irrep: dict[str, Any],
    selected_opd: dict[str, Any] | None = None,
    selected_slots: list[dict[str, Any]] | None = None,
    selected_factors: list[dict[str, Any]] | None = None,
    k_params: dict[str, str] | None = None,
    include_strain: bool = False,
    include_magnetic: bool = False,
    displacive_row_ids: list[int] | None = None,
    magnetic_row_ids: list[int] | None = None,
    displacive_site_indexes: list[int] | None = None,
    magnetic_site_indexes: list[int] | None = None,
    parent_inter_setting_id: int | None = None,
    include_structure: bool = True,
) -> dict[str, Any]:
    """Return complete mode details for the selected OPD."""

    sites = _mapped_sites(cif_info)
    if not sites:
        return {
            "status": "unsupported",
            "reason": "mode kernel requires at least one mapped displacive site",
        }
    selected_slots = list(selected_slots or [])
    selected_factors = list(selected_factors or [])
    coupled = bool(
        (_isotropy_from_opd_row(selected_opd) or {}).get("coupled") and selected_slots
    )
    if (
        not coupled
        and int(selected_irrep.get("old_id") or 0) <= 0
        and _selected_dynamic_gid(selected_irrep, selected_opd) is None
    ):
        return {
            "status": "unsupported",
            "reason": "mode kernel requires either a static old_id or a Source dynamic_gid",
            "irrep": selected_irrep,
        }
    sg = int(cif_info["parent"]["number"])
    decoder = _mode_decoder()
    parent_cell = _cell_params(cif_info)
    child_sg = _selected_subgroup_number(selected_opd)
    target_label = str(
        selected_irrep.get("symbol") or selected_irrep.get("label") or ""
    )
    if not target_label:
        return {
            "status": "unsupported",
            "reason": "mode kernel requires a selected irrep label",
            "irrep": selected_irrep,
        }
    render_specs = (
        _coupled_render_specs(
            decoder,
            sg,
            selected_slots,
            selected_opd,
            include_displacive=displacive_row_ids is None or bool(displacive_row_ids),
            include_magnetic=include_magnetic
            and (magnetic_row_ids is None or bool(magnetic_row_ids)),
        )
        if coupled
        else []
    )
    if not render_specs and not coupled:
        magnetic_selected = bool(include_magnetic and selected_irrep.get("magnetic"))
        mode_specs = _subduced_mode_specs(decoder, sg, selected_irrep, selected_opd)
        render_specs = [(spec, "dsp", 1, str(spec["label"])) for spec in mode_specs]
        if magnetic_selected:
            magnetic_specs = _magnetic_subduced_mode_specs(
                decoder, sg, selected_irrep, selected_opd
            )
            render_specs.extend(
                (
                    spec,
                    "mag",
                    2,
                    str(spec.get("display_label") or "m" + str(spec["label"])),
                )
                for spec in magnetic_specs
            )
    (
        occurrence_alias_candidate_orders,
        occurrence_alias_candidate_to_anchor,
    ) = _occurrence_alias_spec_orders(render_specs)
    occurrence_direction_distinct_orders = (
        direction_distinct_occurrence_alias_spec_orders(
            render_specs,
            occurrence_alias_candidate_to_anchor,
        )
    )
    occurrence_observation_orders = frozenset(
        index
        for index, (spec, _mode_kind, _setting, _label) in enumerate(render_specs)
        if bool(spec.get("_occurrence_alias_observation_only"))
    )
    occurrence_observed_orders = frozenset(
        set(occurrence_alias_candidate_to_anchor)
        | set(occurrence_alias_candidate_to_anchor.values())
    )
    occurrence_eligible_sites: dict[int, list[int]] = {
        order: [] for order in occurrence_observed_orders
    }
    occurrence_emissions: dict[tuple[int, int], OccurrenceSiteEmission] = {}
    pending_occurrence_definitions: dict[int, list[dict[str, Any]]] = {}
    pending_occurrence_counts: dict[tuple[int, int], int] = {}
    strain_mode_specs = [
        {**spec, "_spec_order": spec_order}
        for spec_order, (spec, mode_kind, _setting, _label) in enumerate(render_specs)
        if mode_kind == "dsp"
    ]
    selected_source_rows = (
        (_isotropy_from_opd_row(selected_opd) or {}).get("source_numeric_rows")
        or (_isotropy_from_opd_row(selected_opd) or {}).get("source_matrix")
        or []
    )
    if selected_source_rows:
        for strain_spec in strain_mode_specs:
            if not strain_spec.get("primary") or strain_spec.get("direction_matrix"):
                continue
            try:
                rows = [[float(value) for value in row] for row in selected_source_rows]
                width = len(rows[0])
                if width == 0 or any(len(row) != width for row in rows):
                    continue
                strain_spec["source_numeric_rows"] = rows
                strain_spec["direction_matrix"] = [
                    [rows[free][coordinate] for free in range(len(rows))]
                    for coordinate in range(width)
                ]
            except (TypeError, ValueError):
                continue
    for strain_spec in strain_mode_specs:
        if strain_spec.get("gid") is not None:
            continue
        try:
            strain_spec["gid"] = int(
                decoder.little_record(int(sg), str(strain_spec.get("label") or "")).gid
            )
        except (KeyError, TypeError, ValueError):
            continue
    if include_strain and coupled and not strain_mode_specs:
        strain_mode_specs = [
            spec
            for spec, mode_kind, _setting, _label in _coupled_render_specs(
                decoder,
                sg,
                selected_slots,
                selected_opd,
                include_displacive=True,
                include_magnetic=False,
            )
            if mode_kind == "dsp"
        ]
    parent_symbol = _source_tables().default_setting_space_symbol(int(sg))
    display_kvector = _evaluated_kvector_text(selected_k, k_params)
    selected_gid = _selected_dynamic_gid(selected_irrep, selected_opd)
    opd_direction = _opd_direction_text(
        selected_opd,
        selected_irrep,
        group_width=(
            _little_star_arm_width(decoder, int(selected_gid))
            if selected_gid is not None
            else None
        ),
    )
    display_basis = _basis_from_opd_row(selected_opd)
    presentation_basis, presentation_rule = _presentation_basis_candidate(display_basis)
    child_lattice = _cell_from_basis(parent_cell, presentation_basis)
    child_cartesian = child_lattice_cartesian(
        tuple(
            float(child_lattice[key])
            for key in ("a", "b", "c", "alpha", "beta", "gamma")
        ),
        ((1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, 1.0)),
    )
    split_basis, split_origin = _split_basis_origin_for_wyckoff(
        selected_opd=selected_opd,
    )
    selected_freeparam = _freeparam_from_opd_row(decoder, selected_opd)
    selected_full_parameter_opd = _is_full_parameter_opd(
        selected_opd, selected_freeparam
    )
    selected_origin_override = _origin_from_opd_row(selected_opd)
    parent_setting_bridge = _parent_setting_bridge(sg, parent_inter_setting_id)
    child_symbol = (
        gemmi.find_spacegroup_by_number(int(child_sg)).hm
        if child_sg is not None
        else "P1"
    )
    viewer_placement = parent_cell_placement(
        presentation_basis or [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
        _origin_vector(selected_origin_override),
        child_symbol,
    )
    source_basis_override = _integer_basis_tuple(
        _source_split_basis_from_opd_row(selected_opd)
    )
    selected_basis_pml_override = (
        _basis_cinter_to_pml(
            decoder,
            sg,
            presentation_basis,
            parent_inter_setting_id,
        )
        if presentation_basis is not None
        else None
    )
    if (
        child_sg is None
        or presentation_basis is None
        or selected_basis_pml_override is None
    ):
        raise ValueError(
            "canonical child atom layout requires an exact presentation basis"
        )
    definitions: list[dict[str, Any]] = []
    magnetic_definitions: list[dict[str, Any]] = []
    (
        structure_payload,
        child_atom_layouts,
        mode_atom_positions,
        parent_site_params,
    ) = _build_child_structure(
        cif_info=cif_info,
        selected_opd=selected_opd,
        parent_inter_setting_id=parent_inter_setting_id,
        decoder=decoder,
    )
    undistorted_atoms: list[dict[str, Any]] = [
        dict(atom)
        for atom in (
            structure_payload.get("undistorted_atoms") or []
            if include_structure
            else []
        )
        if isinstance(atom, dict)
    ]
    selected_opd_isotropy: dict[str, Any] | None = _isotropy_from_opd_row(selected_opd)
    selected_isotropy: dict[str, Any] | None = selected_opd_isotropy
    per_site_status: list[dict[str, Any]] = []
    for site_index, site in enumerate(sites):
        atom_label = str(site.get("type") or site.get("label") or "X")
        child_atom_layout = child_atom_layouts[site_index]
        site_params = parent_site_params[site_index]
        site_definition_count = 0
        site_atom_count = 0
        mode_topologies_by_records: dict[
            tuple[tuple[int, int, int, int, int], ...],
            dict[tuple[int, int, tuple[int, int, int, int, int]], str],
        ] = {}
        for spec_index, (
            spec,
            mode_kind,
            vector_setting,
            display_irrep_label,
        ) in enumerate(render_specs):
            allowed_site_indexes = (
                magnetic_site_indexes if mode_kind == "mag" else displacive_site_indexes
            )
            if allowed_site_indexes is not None and site_index + 1 not in {
                int(value) for value in allowed_site_indexes
            }:
                continue
            site_row_id = int(site.get("wyckoff_row_id") or 0)
            allowed_row_ids = (
                magnetic_row_ids if mode_kind == "mag" else displacive_row_ids
            )
            if allowed_row_ids is not None and site_row_id not in {
                int(value) for value in allowed_row_ids
            }:
                continue
            occurrence_alias_candidate = spec_index in occurrence_alias_candidate_orders
            occurrence_observation_only = spec_index in occurrence_observation_orders
            if spec_index in occurrence_observed_orders:
                occurrence_eligible_sites[spec_index].append(site_index)
            spec_k_label = str(spec["k_label"])
            spec_row_id = spec.get("row_id")
            case = Case(
                sg=sg,
                wyckoff=str(site["wyckoff"]),
                k_label=spec_k_label,
                params=parent_cell,
                atom_label=atom_label,
                site_params=site_params,
                k_params=_case_k_params(
                    decoder,
                    selected_k=selected_k,
                    selected_irrep=selected_irrep,
                    selected_opd=selected_opd,
                    spec=spec,
                    k_params=k_params,
                ),
            )
            spec_display_kvector = _spec_display_kvector(
                decoder, case, spec, display_kvector
            )
            spec_opd_direction = _spec_opd_direction_text(
                decoder,
                spec=spec,
                selected_direction=opd_direction,
            )
            trace_uses_raw_basis = bool(
                source_basis_override is not None
                and (mode_kind == "mag" or spec.get("source_numeric_rows"))
            )
            emit_full_spec_rows = _should_emit_all_source_opd_rows(
                decoder,
                coupled=coupled,
                spec=spec,
            )
            source_orderparam_rows = (
                spec.get("source_numeric_rows")
                or (
                    (_isotropy_from_opd_row(selected_opd) or {}).get(
                        "source_numeric_rows"
                    )
                    if mode_kind == "mag"
                    else []
                )
                or []
            )
            trace_gid = int(spec.get("gid") or selected_irrep.get("gid") or 0)
            if trace_gid <= 0:
                raise ValueError("mode spec has no Source gid")
            orderparam_rows_override = (
                {
                    trace_gid: [
                        [float(value) for value in row]
                        for row in source_orderparam_rows
                    ]
                }
                if trace_gid > 0 and source_orderparam_rows
                else None
            )
            occurrence_basis_proof = occurrence_sublattice(
                source_basis_override,
                selected_basis_pml_override,
            )
            occurrence_basis_override = (
                selected_basis_pml_override
                if (
                    trace_uses_raw_basis
                    and mode_kind == "dsp"
                    and case.k_params
                    and orderparam_rows_override is not None
                    and occurrence_basis_proof is not None
                )
                else None
            )
            trace_kwargs: dict[str, Any] = {
                "selected_gid": trace_gid,
                "selected_isotropy_row_id": (
                    int(spec_row_id) if spec_row_id is not None else None
                ),
                "selected_basis_pml_override": (
                    source_basis_override
                    if trace_uses_raw_basis
                    else selected_basis_pml_override
                ),
                "emit_all_opd_rows": selected_full_parameter_opd or emit_full_spec_rows,
                "emit_freeparam_opd_groups": (
                    not selected_full_parameter_opd and not emit_full_spec_rows
                ),
                "vector_setting": vector_setting,
                "orderparam_rows_by_gid": orderparam_rows_override,
                "occurrence_basis_pml_override": occurrence_basis_override,
                "presentation_setting_matrix": (
                    None if parent_setting_bridge is None else parent_setting_bridge[0]
                ),
                "presentation_setting_origin": (
                    None if parent_setting_bridge is None else parent_setting_bridge[1]
                ),
            }
            trace = project_mode_block(decoder, case, **trace_kwargs)
            block = trace["little_irreps"][0]

            site_isotropy = (
                block.get("selected_isotropy_row") if isinstance(block, dict) else None
            )
            if spec.get("primary") and isinstance(site_isotropy, dict):
                selected_isotropy = site_isotropy
            atom_fractionals = block.get("atom_fractionals") or []
            atom_operation_records = block.get("atom_operation_records") or []
            mode_vectors = block.get("mode_vectors") or []
            source_display_projection = bool(block.get("source_display_projection"))
            child_origin = (
                selected_origin_override
                if selected_origin_override is not None
                else (
                    selected_isotropy.get("origin")
                    if isinstance(selected_isotropy, dict)
                    else None
                )
            )
            source_metadata_layout = _mode_source_metadata_layout(
                block,
                len(mode_vectors),
                site_pg=int((trace.get("wyckoff") or {}).get("site_pg") or 0),
            )
            if spec_index in occurrence_observed_orders:
                source_signature = _occurrence_source_signature(
                    source_metadata_layout,
                    len(mode_vectors),
                )
                if source_signature is not None:
                    occurrence_emissions[(spec_index, site_index)] = (
                        OccurrenceSiteEmission(
                            site_index=site_index,
                            row_id=site_row_id,
                            source_signature=source_signature,
                            mode_count=len(mode_vectors),
                            all_modes_nonzero=_occurrence_modes_are_nonzero(
                                mode_vectors
                            ),
                            atom_operation_records=tuple(
                                _freeze_integral_record(record)
                                for record in atom_operation_records
                            ),
                            atom_fractionals=tuple(
                                _freeze_fractional_xyz(xyz) for xyz in atom_fractionals
                            ),
                        )
                    )
            if occurrence_observation_only:
                continue
            printable_mode_vectors = (
                _presentation_mode_vectors(mode_vectors, presentation_basis)
                if not trace_uses_raw_basis or source_display_projection
                else mode_vectors
            )
            site_irrep_labels = _mode_site_irrep_labels(
                decoder,
                trace,
                block,
                len(printable_mode_vectors),
            )
            display_site_irrep_labels = (
                _mode_site_irrep_labels(
                    decoder,
                    trace,
                    block,
                    len(printable_mode_vectors),
                    opd_rows_as_components=True,
                )
                if any(
                    int(summary.get("little_type") or 0) == 3
                    or len(summary.get("opd_groups") or []) > 1
                    for summary in block.get("mode_block_summaries") or []
                    if isinstance(summary, dict)
                )
                else None
            )
            source_layout = _mode_source_layout(block, len(printable_mode_vectors))
            spec_gid = int(spec.get("gid") or 0)
            child_symbol = gemmi.find_spacegroup_by_number(int(child_sg)).hm
            operation_record_key = tuple(
                exact_operation_record(
                    record,
                    point_operation_count=len(decoder.iso.space["ipoint_op_inverse"]),
                )
                for record in atom_operation_records
            )
            mode_topology = mode_topologies_by_records.get(operation_record_key)
            if mode_topology is None:
                mode_topology = _compile_child_atom_mode_topology(
                    decoder,
                    child_sg=child_sg,
                    presentation_basis_pml=selected_basis_pml_override,
                    layout=child_atom_layout,
                    atom_operation_records=list(operation_record_key),
                )
                mode_topologies_by_records[operation_record_key] = mode_topology
            for mode_index, vectors in enumerate(printable_mode_vectors, start=1):
                if not isinstance(vectors, list):
                    continue
                site_irrep_label = (
                    None
                    if site_irrep_labels is None
                    else site_irrep_labels[mode_index - 1]
                )
                display_site_irrep_label = (
                    None
                    if display_site_irrep_labels is None
                    else display_site_irrep_labels[mode_index - 1]
                )
                public_site_irrep = (
                    display_site_irrep_label
                    or site_irrep_label
                    or ("mode" if mode_index == 1 else f"mode_{mode_index}")
                )
                rows = [
                    {
                        "dxyz": [float(value) for value in vector],
                        "_source_raw_index": atom_index,
                        "_operation_record": (
                            atom_operation_records[atom_index]
                            if atom_index < len(atom_operation_records)
                            else None
                        ),
                    }
                    for atom_index, vector in enumerate(vectors)
                    if isinstance(vector, list)
                ]
                if (
                    trace_uses_raw_basis
                    and not source_display_projection
                    and presentation_basis is not None
                ):
                    raw_rows = [
                        {
                            "xyz": [
                                float(Fraction(str(value)))
                                for value in atom_fractionals[index]
                            ],
                            "dxyz": [float(value) for value in vectors[index]],
                            "_source_raw_index": index,
                            **(
                                {"_operation_record": atom_operation_records[index]}
                                if index < len(atom_operation_records)
                                else {}
                            ),
                        }
                        for index in range(min(len(atom_fractionals), len(vectors)))
                        if isinstance(vectors[index], list)
                    ]
                    presented = present_mode_rows(
                        raw_rows,
                        basis=presentation_basis,
                        origin=_origin_vector(child_origin),
                        centering_symbol=child_symbol,
                        include_centering_ordinal=True,
                    )
                    rows = [
                        {
                            "atom": f"{site.get('label') or atom_label}_1"
                            if index == 0
                            else None,
                            "xyz": list(row["xyz"]),
                            "dxyz": list(row["dxyz"]),
                            **(
                                {"_operation_record": row["_operation_record"]}
                                if row.get("_operation_record") is not None
                                else {}
                            ),
                            **(
                                {"_source_raw_index": row["_source_raw_index"]}
                                if row.get("_source_raw_index") is not None
                                else {}
                            ),
                            **(
                                {
                                    "_presentation_centering_ordinal": row[
                                        "_presentation_centering_ordinal"
                                    ]
                                }
                                if row.get("_presentation_centering_ordinal")
                                is not None
                                else {}
                            ),
                        }
                        for index, row in enumerate(presented["rows"])
                    ]
                rows = _mode_rows_on_child_atom_layout(
                    layout=child_atom_layout,
                    mode_topology=mode_topology,
                    rows=rows,
                    mode_atom_positions=mode_atom_positions,
                )
                source_print_identity = (
                    None
                    if source_metadata_layout is None
                    else source_metadata_layout[mode_index - 1]
                )
                source_print_scalar = (
                    _rank1_source_print_intertwiner(
                        decoder,
                        sg=sg,
                        child_sg=child_sg,
                        case=case,
                        spec=spec,
                        trace=trace,
                        block=block,
                        site=site,
                        rows=rows,
                        identity=source_print_identity,
                        presentation_basis=presentation_basis,
                        vector_setting=vector_setting,
                        subgroup_basis=source_basis_override,
                        subgroup_origin=_source_split_origin_from_opd_row(selected_opd),
                    )
                    if mode_kind == "mag"
                    and spec_gid > 0
                    and int(spec.get("old_id") or 0) <= 0
                    else None
                )
                if mode_kind == "mag":
                    rows = _web_magnetic_occurrence_gauge(decoder, sg, rows)
                target_definitions = (
                    pending_occurrence_definitions.setdefault(spec_index, [])
                    if occurrence_alias_candidate
                    else magnetic_definitions
                    if mode_kind == "mag"
                    else definitions
                )
                target_definitions.append(
                    {
                        "label": _mode_label(
                            parent_symbol,
                            spec_display_kvector,
                            display_irrep_label,
                            spec_opd_direction,
                            site,
                            site_irrep_label,
                            mode_index,
                            mode_kind,
                        ),
                        "mode": {
                            "kind": (
                                "magnetic" if mode_kind == "mag" else "displacive"
                            ),
                            "k": spec_k_label,
                            "kvector": f"({spec_display_kvector})",
                            "irrep": display_irrep_label,
                            "gid": spec_gid,
                            "direction": spec_opd_direction,
                            "site_type": str(site.get("type") or atom_label),
                            "site_label": str(site.get("label") or atom_label),
                            "wyckoff": str(site.get("wyckoff") or ""),
                            "site_irrep": public_site_irrep,
                        },
                        "normfactor": _complete_mode_normfactor(
                            rows, child_cartesian, child_sg
                        ),
                        "rows": rows,
                        "_spec_order": spec_index,
                        "_site_order": site_index,
                        "_mode_kind": mode_kind,
                        "_source_family_dynamic": bool(
                            spec_gid > 0 and int(spec.get("old_id") or 0) <= 0
                        ),
                        "_source_family": (
                            None
                            if source_layout is None
                            else source_layout[mode_index - 1][0]
                        ),
                        "_source_family_component": (
                            None
                            if source_layout is None
                            else source_layout[mode_index - 1][1]
                        ),
                        "_source_family_width": (
                            None
                            if source_layout is None
                            else source_layout[mode_index - 1][2]
                        ),
                        "_source_family_phase": bool(
                            spec_gid > 0
                            and int(spec.get("old_id") or 0) <= 0
                            and int(decoder.little_record_by_gid(spec_gid).irrep_type)
                            == 1
                            and (
                                mode_kind == "mag"
                                or decoder.little_transform_block_count(spec_gid) > 1
                            )
                        ),
                        "_source_family_positive_lead": bool(
                            mode_kind == "mag"
                            and spec_gid > 0
                            and int(spec.get("old_id") or 0) <= 0
                            and int(decoder.little_record_by_gid(spec_gid).irrep_type)
                            == 1
                        ),
                        "_source_print_identity": (
                            None
                            if source_print_identity is None
                            else {
                                **source_print_identity,
                                "source_kparam": tuple(
                                    int(value)
                                    for value in spec.get("carrier_source_kparam")
                                    or spec.get("source_kparam")
                                    or ()
                                ),
                                "direction_matrix": [
                                    list(row)
                                    for row in spec.get("direction_matrix") or []
                                ],
                            }
                        ),
                        "_display_site_irrep_label": (
                            display_site_irrep_label
                        ),
                        "_source_print_scalar": source_print_scalar,
                    }
                )
            if occurrence_alias_candidate:
                pending_occurrence_counts[(spec_index, site_index)] = len(mode_vectors)
            else:
                site_definition_count += len(mode_vectors)
                site_atom_count = max(site_atom_count, len(atom_fractionals))
        per_site_status.append(
            {
                "site": site.get("label"),
                "status": "ok",
                "definition_count": site_definition_count,
                "atom_count": site_atom_count,
                "missing": [],
            }
        )
    admitted_occurrence_orders = admitted_occurrence_alias_spec_orders(
        occurrence_alias_candidate_to_anchor,
        eligible_sites_by_spec=occurrence_eligible_sites,
        emissions=occurrence_emissions,
    ) | admitted_self_complete_occurrence_alias_spec_orders(
        occurrence_direction_distinct_orders,
        eligible_sites_by_spec=occurrence_eligible_sites,
        emissions=occurrence_emissions,
    )
    for spec_order in sorted(admitted_occurrence_orders):
        target_definitions = (
            magnetic_definitions
            if render_specs[spec_order][1] == "mag"
            else definitions
        )
        target_definitions.extend(pending_occurrence_definitions.get(spec_order, ()))
        for site_index in occurrence_eligible_sites.get(spec_order, ()):
            count = int(pending_occurrence_counts.get((spec_order, site_index), 0))
            if count > 0 and 0 <= site_index < len(per_site_status):
                per_site_status[site_index]["definition_count"] = (
                    int(per_site_status[site_index].get("definition_count") or 0)
                    + count
                )
    definitions.sort(
        key=lambda item: (
            int(item.get("_spec_order", 0)),
            int(item.get("_site_order", 0)),
        )
    )
    magnetic_definitions.sort(
        key=lambda item: (
            int(item.get("_spec_order", 0)),
            int(item.get("_site_order", 0)),
        )
    )
    emitted_spec_orders = {
        int(item["_spec_order"])
        for item in (*definitions, *magnetic_definitions)
        if type(item.get("_spec_order")) is int
    }
    secondary_opd_factors = _secondary_invariant_factors(
        render_specs=render_specs,
        emitted_spec_orders=emitted_spec_orders,
        primary_slot_count=max(1, len(selected_slots)),
        coupled=coupled,
        magnetic_primary=bool(selected_irrep.get("magnetic")),
    )
    definitions = _orthogonalize_definition_modes(
        definitions, child_cartesian, child_sg
    )
    magnetic_definitions = _orthogonalize_definition_modes(
        magnetic_definitions, child_cartesian, child_sg
    )
    definitions = _apply_dynamic_source_family_presentation(definitions)
    magnetic_definitions = _apply_dynamic_source_family_presentation(
        magnetic_definitions
    )
    strain_definitions = (
        _strain_mode_definitions(
            sg,
            parent_symbol,
            decoder,
            strain_mode_specs,
        )
        if include_strain
        else []
    )
    _attach_definition_factor_metadata(
        (*definitions, *magnetic_definitions, *strain_definitions),
        render_specs=render_specs,
        selected_factors=selected_factors,
    )
    viewer_atoms: list[dict[str, Any]] = []
    if include_structure:
        expected_atom_ids = {
            atom.atom_id for layout in child_atom_layouts for atom in layout.atoms
        }
        if set(mode_atom_positions) != expected_atom_ids:
            missing = sorted(expected_atom_ids - set(mode_atom_positions))
            extra = sorted(set(mode_atom_positions) - expected_atom_ids)
            raise ValueError(
                "mode coordinates do not cover the canonical child layout: "
                f"missing={missing[:8]!r}, extra={extra[:8]!r}"
            )
        for site_index, (site, layout) in enumerate(
            zip(sites, child_atom_layouts, strict=True)
        ):
            element = str(site.get("type") or site.get("label") or "X")
            for atom_index, atom in enumerate(layout.atoms):
                viewer_atoms.append(
                    {
                        "label": atom.child_site_id,
                        "atom_id": atom.atom_id,
                        "child_site": atom.child_site_id,
                        "element": element,
                        "site_order": site_index,
                        "atom_index": atom_index,
                        "xyz": list(mode_atom_positions[atom.atom_id]),
                    }
                )
    if not definitions and not magnetic_definitions and not strain_definitions:
        return {
            "status": "missing",
            "reason": f"selected irrep {target_label!r} produced no local site definitions",
            "sites": per_site_status,
        }
    return {
        "status": "ok",
        "sites": per_site_status,
        "lattice": child_lattice,
        "viewer_parent_basis": viewer_placement["parent_basis"],
        "viewer_parent_origin": viewer_placement["parent_origin"],
        "subgroup_details": (
            None
            if not isinstance(selected_isotropy, dict)
            else {
                "row_id": selected_isotropy.get("row_id"),
                "subgroup": selected_isotropy.get("subgroup"),
                "arms": selected_isotropy.get("arms"),
                "direction": selected_isotropy.get("direction"),
                "basis": selected_isotropy.get("basis"),
                "basis_pml_to_cinter": selected_isotropy.get("basis_pml_to_cinter"),
                "display_basis": display_basis,
                "presentation_basis": presentation_basis,
                "presentation_rule": presentation_rule,
                "origin": selected_origin_override,
                "magnetic_subgroup": (
                    selected_opd_isotropy.get("subgroup")
                    if isinstance(selected_opd_isotropy, dict)
                    else None
                ),
                "source_basis_values": (
                    selected_opd_isotropy.get("source_basis_values")
                    if isinstance(selected_opd_isotropy, dict)
                    else None
                ),
                "source_origin_values": (
                    selected_opd_isotropy.get("source_origin_values")
                    if isinstance(selected_opd_isotropy, dict)
                    else None
                ),
                "parent_inter_setting_id": parent_inter_setting_id,
            }
        ),
        "undistorted_atoms": undistorted_atoms,
        "viewer_atoms": viewer_atoms,
        "displacive_definitions": definitions,
        "strain_definitions": strain_definitions,
        "magnetic_definitions": magnetic_definitions,
        "secondary_opd_factors": secondary_opd_factors,
        "displacive_amplitudes": [],
    }
