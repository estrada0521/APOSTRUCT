"""Expose the mode-kernel atom-displacement projection trace.

The trace follows Source projection rows through coefficient preparation and
the final vector expansion while preserving the runtime call boundaries.
"""

from __future__ import annotations

from fractions import Fraction
from typing import Any, Sequence

from distortropy.Backend.modes.engine.decoder import ModeDataDecoder
from distortropy.Backend.modes.engine.input import Case
from distortropy.Backend.modes.engine.project.mode_forms import (
    _add_mode_vectors,
    _apply_fractional_vector_transform,
    _normalize_mode_vectors,
    _regroup_repeated_component_modes,
    _type1_parametric_scalar_plane_print_order,
    _type3_parametric_component_print_order,
    _type3_parametric_kdim2_print_basis,
    _type3_real_print_modes,
)
from distortropy.Backend.modes.engine.project.main_bridge import (
    dense_orderparam_source,
    dense_source_e8,
    dense_type3_parametric_source,
    direct_bridge_basis_function,
    direct_bridge_coefficients,
    effective_orderparam_rows,
)
from distortropy.Backend.modes.engine.project.display_distortion import (
    DisplayProjectionBlock,
    resolve_display_projection,
    source_display_projection_plan,
)
from distortropy.Backend.modes.engine.project.mode_counts import little_records_for_k, mode_totals
from distortropy.Backend.modes.engine.project.entry import (
    project_entry_trace,
    project_local408_trace,
    site_get_irreps_trace,
)
from distortropy.Backend.modes.engine.trace.project_sources import (
    parent_get_irrep4_trace,
    project_return_basis_type1_trace,
    project_vector_bridge_source_trace,
    project_vector_orderparam_source_trace,
    project_vector_prep_trace,
    site_ssgn_matrix_trace,
)

def pipeline_trace(
    decoder: ModeDataDecoder,
    case: Case,
    *,
    selected_gid: int | None = None,
    selected_isotropy_row_id: int | None = None,
    selected_basis_pml_override: tuple[int, ...] | None = None,
    emit_all_opd_rows: bool = False,
    emit_freeparam_opd_groups: bool = False,
    emit_all_project_families: bool = False,
    vector_setting: int = 1,
    orderparam_rows_by_gid: dict[int, list[list[float]]] | None = None,
    occurrence_basis_pml_override: tuple[int, ...] | None = None,
    occurrence_translations_cinter: Sequence[Sequence[float]] | None = None,
    presentation_setting_matrix: Sequence[Sequence[Fraction | float | int]] | None = None,
    presentation_setting_origin: Sequence[Fraction | float | int] | None = None,
    _diagnostics: bool = True,
) -> dict[str, object]:
    """Return the mode-projection intermediate state for a parsed input case."""

    if occurrence_basis_pml_override is not None and (
        len(occurrence_basis_pml_override) != 9
        or any(
            isinstance(value, bool) or not isinstance(value, int)
            for value in occurrence_basis_pml_override
        )
    ):
        raise ValueError("occurrence basis override must contain nine exact integers")

    rows = [row for row in decoder.wyckoff_rows(case.sg) if row.label == case.wyckoff]
    if len(rows) != 1:
        raise KeyError(f"expected one Wyckoff row for SG{case.sg} {case.wyckoff}, got {len(rows)}")
    row = rows[0]
    all_littles = little_records_for_k(decoder, case.sg, case.k_label)
    if selected_gid is None:
        littles = all_littles
    else:
        if isinstance(selected_gid, bool) or not isinstance(selected_gid, int) or selected_gid <= 0:
            raise ValueError("selected gid must be a positive exact integer")
        littles = [little for little in all_littles if little.gid == selected_gid]
        if len(littles) != 1:
            raise ValueError(
                f"gid {selected_gid} does not identify one SG{case.sg} {case.k_label} irrep"
            )
    totals = mode_totals(decoder, case) if _diagnostics else {}
    if int(vector_setting) not in (1, 2):
        raise ValueError(f"unsupported vector setting: {vector_setting}")
    project_entries = project_entry_trace(
        decoder,
        case,
        row,
        vector_setting=int(vector_setting),
        selected_gid=selected_gid,
    )
    project_local408 = project_local408_trace(decoder, case, row, project_entries)
    if _diagnostics:
        site_ssgn = site_ssgn_matrix_trace(
            decoder, row, project_entries, project_local408
        )
        parent_get_irrep4 = parent_get_irrep4_trace(
            decoder, case, row, project_entries
        )
    project_basis_type1 = project_return_basis_type1_trace(decoder, case, row, project_entries, project_local408)
    project_vector_prep = project_vector_prep_trace(decoder, case, row, project_entries, project_local408)
    project_vector_bridge_source = project_vector_bridge_source_trace(decoder, project_basis_type1)
    if _diagnostics:
        k_vector = decoder.k_vector_from_case(case)
        display_k_vector = decoder.display_k_vector_from_case(case)
        newlat_input = decoder.newlat_input_records_for_case(case) if case.k_params else ()
        find_isotropy_orderparam = (
            decoder.find_isotropy_orderparam_records_for_case(case)
            if case.k_params
            else ()
        )
    find_isotropy_basis = decoder.find_isotropy_initial_basis_for_case(case) if case.k_params else ()
    if _diagnostics:
        k_supercell_index = decoder.k_supercell_index(case)
        wyckoff_multiplicity = decoder.wyckoff_multiplicity(case.sg, row)
        supercell_atom_count = decoder.supercell_atom_count(case, row)
    pml_to_cinter = decoder.pml_to_cinter_matrix(case.sg)
    isotropy_rows_by_gid = {
        little.gid: decoder.isotropy_rows_for_little_gid(little.gid)
        for little in littles
    }
    selected_isotropy_row_by_gid = {}
    selected_basis_by_gid = {}
    project_entry_by_gid = {
        int(item["gid"]): item
        for item in project_entries
        if isinstance(item, dict) and "gid" in item
    }
    canonical_atom_basis = None
    for little in all_littles:
        rows_for_little = isotropy_rows_by_gid.get(little.gid)
        if rows_for_little is None:
            rows_for_little = decoder.isotropy_rows_for_little_gid(little.gid)
        if rows_for_little:
            canonical_atom_basis = tuple(rows_for_little[0].basis)
            break
    if canonical_atom_basis is None and find_isotropy_basis:
        canonical_atom_basis = tuple(int(value) for value in find_isotropy_basis)
    for little in littles:
        rows_for_little = isotropy_rows_by_gid[little.gid]
        project_entry = project_entry_by_gid.get(little.gid, {})
        target_direction = int(project_entry["param5_second"]) if "param5_second" in project_entry else None
        explicit_selected_row = (
            next((item for item in rows_for_little if int(item.row_id) == int(selected_isotropy_row_id)), None)
            if selected_isotropy_row_id is not None
            else None
        )
        selected_isotropy_row_by_gid[little.gid] = explicit_selected_row or (
            decoder.canonical_isotropy_row(rows_for_little)
            if not case.k_params
            else (
                next((item for item in rows_for_little if int(item.direction) == target_direction), None)
                if target_direction is not None
                else None
            )
        )
        if selected_isotropy_row_by_gid[little.gid] is None:
            selected_isotropy_row_by_gid[little.gid] = next(
                (
                    item
                    for item in rows_for_little
                    if tuple(int(value) for value in item.origin) == (0, 0, 0, 1)
                ),
                rows_for_little[0] if rows_for_little else None,
            )
        if selected_basis_pml_override is not None and (
            explicit_selected_row is not None
            or int(little.gid) in (orderparam_rows_by_gid or {})
        ):
            selected_basis_by_gid[little.gid] = tuple(int(value) for value in selected_basis_pml_override)
        elif selected_isotropy_row_by_gid[little.gid] is not None:
            selected_basis_by_gid[little.gid] = tuple(selected_isotropy_row_by_gid[little.gid].basis)
        elif find_isotropy_basis:
            selected_basis_by_gid[little.gid] = tuple(int(value) for value in find_isotropy_basis)
        else:
            selected_basis_by_gid[little.gid] = None
    project_vector_orderparam_source = project_vector_orderparam_source_trace(
        decoder,
        littles,
        selected_isotropy_row_by_gid,
        case,
        orderparam_rows_by_gid,
    )
    baseline_atom_basis_by_gid = {}
    occurrence_atom_basis_by_gid = {}
    canonical_rows = None
    if canonical_atom_basis is not None:
        canonical_rows = {
            tuple(canonical_atom_basis[3 * row_index:3 * row_index + 3])
            for row_index in range(3)
        }
    for little in littles:
        selected_basis = selected_basis_by_gid[little.gid]
        atom_basis = selected_basis
        if selected_basis is not None and canonical_atom_basis is not None and canonical_rows is not None:
            selected_rows = {
                tuple(selected_basis[3 * row_index:3 * row_index + 3])
                for row_index in range(3)
            }
            if selected_rows == canonical_rows:
                atom_basis = canonical_atom_basis
        baseline_atom_basis_by_gid[little.gid] = atom_basis
        occurrence_atom_basis_by_gid[little.gid] = (
            occurrence_basis_pml_override
            if occurrence_basis_pml_override is not None
            and int(little.gid) in (orderparam_rows_by_gid or {})
            else None
        )
    atom_fractionals_by_gid = {}
    atom_operation_records_by_gid = {}
    atom_order_by_gid = {}
    mode_vectors_by_gid = {}
    mode_block_summaries_by_gid: dict[int, list[dict[str, object]]] = {}
    source_display_projection_by_gid: dict[int, dict[str, int]] = {}
    project_vector_prep_by_gid: dict[int, list[dict[str, object]]] = {}
    for item in project_vector_prep:
        if isinstance(item, dict) and "gid" in item:
            project_vector_prep_by_gid.setdefault(int(item["gid"]), []).append(item)
    project_basis_by_branch = {
        (int(item["gid"]), int(item["pg_irrep"])): item
        for item in project_basis_type1
        if isinstance(item, dict) and "gid" in item and "pg_irrep" in item
    }
    bridge_sources_by_branch: dict[tuple[int, int], list[dict[str, object]]] = {}
    for item in project_vector_bridge_source:
        if isinstance(item, dict) and "gid" in item and "pg_irrep" in item:
            bridge_sources_by_branch.setdefault((int(item["gid"]), int(item["pg_irrep"])), []).append(item)
    for values in bridge_sources_by_branch.values():
        values.sort(key=lambda item: int(item["family"]))
    orderparam_source_by_gid = {
        int(item["gid"]): item
        for item in project_vector_orderparam_source
        if isinstance(item, dict) and "gid" in item
    }
    project_weight_buffers: dict[tuple[object, ...], Sequence[float]] = {}

    def project_weight_buffer(
        gid: int,
        record: tuple[int, int, int, int, int],
    ) -> Sequence[float]:
        key = (
            int(gid),
            tuple(int(value) for value in record),
        )
        if key not in project_weight_buffers:
            project_weight_buffers[key] = (
                decoder.project_vector_bridge_weight_view_for_record(
                    int(gid),
                    key[1],
                    case,
                )
            )
        return project_weight_buffers[key]

    final_vector_transform = decoder.cml_to_cinter_matrix(case.sg)
    setting_matrix = (
        tuple(tuple(Fraction(value) for value in row) for row in presentation_setting_matrix)
        if presentation_setting_matrix is not None
        else None
    )
    setting_origin = (
        tuple(Fraction(value) for value in presentation_setting_origin)
        if presentation_setting_origin is not None
        else (Fraction(0), Fraction(0), Fraction(0))
    )
    setting_handedness = 1
    if setting_matrix is not None:
        determinant = (
            setting_matrix[0][0]
            * (setting_matrix[1][1] * setting_matrix[2][2] - setting_matrix[1][2] * setting_matrix[2][1])
            - setting_matrix[0][1]
            * (setting_matrix[1][0] * setting_matrix[2][2] - setting_matrix[1][2] * setting_matrix[2][0])
            + setting_matrix[0][2]
            * (setting_matrix[1][0] * setting_matrix[2][1] - setting_matrix[1][1] * setting_matrix[2][0])
        )
        setting_handedness = -1 if determinant < 0 else 1

    def transform_setting_point(
        point: Sequence[Fraction],
    ) -> tuple[Fraction, Fraction, Fraction]:
        if setting_matrix is None:
            return tuple(point)  # type: ignore[return-value]
        return tuple(
            sum(point[row] * setting_matrix[row][col] for row in range(3)) + setting_origin[col]
            for col in range(3)
        )  # type: ignore[return-value]

    def finalize_vector(vector: list[float]) -> list[float]:
        if final_vector_transform is not None:
            vector = _apply_fractional_vector_transform(vector, final_vector_transform)
        if setting_matrix is not None:
            vector = _apply_fractional_vector_transform(vector, setting_matrix)
            if int(vector_setting) == 2 and setting_handedness < 0:
                vector = [-value for value in vector]
        return vector

    def materialize_atom_layout(
        gid: int,
        atom_basis: Sequence[int],
    ) -> tuple[
        Sequence[Sequence[Fraction]],
        Sequence[Sequence[int]],
        Sequence[Sequence[Fraction]],
        Sequence[int],
    ]:
        raw_atom_fractionals = decoder.display_distortion_atom_fractionals(
            case.sg,
            row,
            case.site_params,
            atom_basis,
        )
        raw_atom_records = decoder.supercell_atom_operation_records(
            case.sg,
            row,
            atom_basis,
        )
        if occurrence_translations_cinter and gid not in (orderparam_rows_by_gid or {}):
            expanded_fractionals = []
            expanded_records = []
            for translation_values in occurrence_translations_cinter:
                translation = tuple(Fraction(str(value)) for value in translation_values)
                if len(translation) != 3:
                    raise ValueError(
                        "expected cinter occurrence translation triplet, "
                        f"got {translation_values!r}"
                    )
                translation_record = decoder.cinter_translation_to_pml_record(
                    case.sg,
                    translation,
                )
                for position, operation in zip(
                    raw_atom_fractionals,
                    raw_atom_records,
                    strict=True,
                ):
                    expanded_fractionals.append(
                        tuple(position[axis] + translation[axis] for axis in range(3))
                    )
                    expanded_records.append(
                        decoder.add_translation_to_operation_record(
                            translation_record[:4],
                            operation,
                        )
                    )
            raw_atom_fractionals = tuple(expanded_fractionals)
            raw_atom_records = tuple(expanded_records)
        displayed_atom_fractionals = tuple(
            transform_setting_point(position) for position in raw_atom_fractionals
        )
        atom_order = decoder.final_atom_order_indices(raw_atom_records)
        return (
            raw_atom_fractionals,
            raw_atom_records,
            displayed_atom_fractionals,
            atom_order,
        )

    for little in littles:
        rows_for_little = isotropy_rows_by_gid[little.gid]
        baseline_atom_basis = baseline_atom_basis_by_gid[little.gid]
        occurrence_atom_basis = occurrence_atom_basis_by_gid[little.gid]
        if occurrence_atom_basis is not None or baseline_atom_basis is not None:
            candidate_layout = None
            if occurrence_atom_basis is not None:
                try:
                    candidate_layout = materialize_atom_layout(
                        int(little.gid),
                        occurrence_atom_basis,
                    )
                except (
                    IndexError,
                    KeyError,
                    OverflowError,
                    TypeError,
                    ValueError,
                    ZeroDivisionError,
                ):
                    candidate_layout = None
            raw_preps = project_vector_prep_by_gid.get(little.gid, [])
            direct_direction_rows = (orderparam_rows_by_gid or {}).get(int(little.gid))

            def build_candidate_projection(
                layout: tuple[
                    Sequence[Sequence[Fraction]],
                    Sequence[Sequence[int]],
                    Sequence[Sequence[Fraction]],
                    Sequence[int],
                ],
            ) -> tuple[DisplayProjectionBlock, ...] | None:
                if not (
                    little.irrep_type == 1
                    and bool(case.k_params)
                    and int(vector_setting) == 1
                    and direct_direction_rows
                ):
                    return None
                _fractionals, candidate_records, _displayed, candidate_order = layout
                return source_display_projection_plan(
                    decoder,
                    case,
                    gid=int(little.gid),
                    full_dim=int(little.full_dim),
                    vector_setting=int(vector_setting),
                    direction_rows=direct_direction_rows,
                    occurrence_records=[
                        candidate_records[index] for index in candidate_order
                    ],
                    raw_preps=raw_preps,
                    project_basis_by_branch=project_basis_by_branch,
                    bridge_sources_by_branch=bridge_sources_by_branch,
                )

            def build_baseline_layout() -> tuple[
                Sequence[Sequence[Fraction]],
                Sequence[Sequence[int]],
                Sequence[Sequence[Fraction]],
                Sequence[int],
            ] | None:
                if baseline_atom_basis is None:
                    return None
                return materialize_atom_layout(int(little.gid), baseline_atom_basis)

            layout, direct_projection, projection_adopted = resolve_display_projection(
                candidate_layout,
                build_candidate_projection,
                build_baseline_layout,
            )
            if layout is None:
                atom_fractionals_by_gid[little.gid] = []
                mode_vectors_by_gid[little.gid] = []
                continue
            (
                raw_atom_fractionals,
                raw_atom_records,
                displayed_atom_fractionals,
                atom_order,
            ) = layout
            direct_projection_by_prep = {
                block.prep_index: block for block in direct_projection
            }
            if projection_adopted:
                source_display_projection_by_gid[int(little.gid)] = {
                    "occurrence_count": len(atom_order),
                    "prep_count": len(direct_projection),
                    "family_count": sum(
                        len(block.families) for block in direct_projection
                    ),
                }

            atom_fractionals_by_gid[little.gid] = [
                [str(value) for value in position]
                for position in (displayed_atom_fractionals[index] for index in atom_order)
            ]
            atom_operation_records_by_gid[little.gid] = [
                [int(value) for value in raw_atom_records[index]]
                for index in atom_order
            ]
            atom_order_by_gid[little.gid] = [int(index) for index in atom_order]
            mode_vectors_by_gid[little.gid] = []
            mode_block_summaries_by_gid[little.gid] = []
            for prep_index, prep in enumerate(raw_preps):
                if not (
                    isinstance(prep, dict)
                    and little.irrep_type in {1, 2, 3}
                    and int(prep.get("vector_dim", 0)) in {1, 2, 3}
                ):
                    continue
                pg_irrep = int(prep["pg_irrep"])
                project_count = int(prep.get("project_count", 0))
                vector_rep_values = decoder.site_vector_reps(row.site_pg)
                vector_start = (int(vector_setting) - 1) * 3
                vector_slots = tuple(int(slot) for slot in vector_rep_values[vector_start:vector_start + 3])
                matching_components = [
                    component_start
                    for component_start, slot in enumerate(vector_slots, start=1)
                    if slot == pg_irrep
                ]

                direct_block = direct_projection_by_prep.get(prep_index)
                if direct_block is not None:
                    direct_groups = [
                        {
                            "index": int(index),
                            "row_count": 1,
                            "rows": [[float(value) for value in row_values]],
                        }
                        for index, row_values in enumerate(direct_direction_rows or [])
                    ]
                    for family in direct_block.families:
                        modes = [
                            [finalize_vector(list(vector)) for vector in mode]
                            for mode in family.modes
                        ]
                        mode_vectors_by_gid[little.gid].extend(modes)
                        mode_block_summaries_by_gid[little.gid].append(
                            {
                                "gid": int(little.gid),
                                "label": little.label,
                                "pg_irrep": int(direct_block.pg_irrep),
                                "source_family": int(family.family),
                                "source_row_count": int(len(direct_groups)),
                                "opd_group_count": int(len(direct_groups)),
                                "opd_groups": direct_groups,
                                "block_count": int(len(direct_groups)),
                                "block_mode_counts": [1] * len(direct_groups),
                                "flat_mode_count": int(len(modes)),
                                "emitted_mode_count": int(direct_block.project_count),
                                "project_count": int(direct_block.project_count),
                                "bridge_project_count": int(
                                    direct_block.bridge_project_counts[family.family]
                                ),
                                "vector_dim": int(direct_block.vector_dim),
                                "little_type": int(little.irrep_type),
                                "little_full_dim": int(little.full_dim),
                            }
                        )
                    continue

                def recordwise_project_vector_modes() -> list[list[list[float]]]:
                    sources = bridge_sources_by_branch.get((little.gid, pg_irrep), [])
                    orderparam_source = orderparam_source_by_gid.get(little.gid)
                    if not sources or not isinstance(orderparam_source, dict):
                        return []
                    source_rows = orderparam_source.get("rows", [])
                    if not isinstance(source_rows, list) or not source_rows:
                        return []
                    source_row_count = effective_orderparam_rows(source_rows)  # type: ignore[arg-type]
                    if source_row_count <= 0:
                        return []
                    displayed_records = [raw_atom_records[index] for index in atom_order]
                    type1_dim4_identity_pair = (
                        little.irrep_type == 1
                        and bool(case.k_params)
                        and int(little.full_dim) == 4
                        and int(prep.get("vector_dim", 0)) in (1, 2)
                        and source_row_count == 2
                        and len(source_rows[0]) >= 2  # type: ignore[index]
                        and len(source_rows[1]) >= 2  # type: ignore[index]
                        and abs(float(source_rows[0][0]) - 1.0) <= 1e-12  # type: ignore[index]
                        and abs(float(source_rows[1][1]) - 1.0) <= 1e-12  # type: ignore[index]
                        and all(abs(float(value)) <= 1e-12 for index, value in enumerate(source_rows[0]) if index != 0)  # type: ignore[index]
                        and all(abs(float(value)) <= 1e-12 for index, value in enumerate(source_rows[1]) if index != 1)  # type: ignore[index]
                    )
                    source_iter = list(reversed(sources)) if type1_dim4_identity_pair else sources
                    modes: list[list[list[float]]] = []
                    for source in source_iter:
                        term_count = int(source["bridge_project_count"]) + 1
                        two_row_slots = (
                            little.irrep_type == 3
                            and source_row_count >= 2
                            and (
                                source_row_count == 2
                                or
                                len(displayed_records) <= 4
                                or (
                                    source_row_count >= 4
                                    and project_count in (1, 2)
                                    and int(prep.get("vector_dim", 0)) in (1, 2)
                                )
                            )
                        )
                        source_row_id = int(orderparam_source.get("row_id", 0))
                        freeparam_group_count = 0
                        if emit_freeparam_opd_groups and source_row_id > 0:
                            try:
                                freeparam = int(decoder.isotropy_orderparam_freeparam(source_row_id))
                                freeparam_group_count = max(0, freeparam // int(decoder.TYPE_DIVISOR[little.irrep_type]))
                            except Exception:
                                freeparam_group_count = 0
                        if emit_freeparam_opd_groups and source_row_id > 0 and little.irrep_type == 2 and source_row_count > 1:
                            opd_groups = [([source_rows[i]], 1) for i in range(source_row_count)]  # type: ignore[index]
                        elif freeparam_group_count > 0 and little.irrep_type == 1:
                            opd_groups = [([source_rows[i]], 1) for i in range(min(source_row_count, freeparam_group_count))]  # type: ignore[index]
                        elif freeparam_group_count > 0 and little.irrep_type == 3:
                            opd_groups = [
                                (list(source_rows[index:index + 2]), 2)  # type: ignore[index]
                                for index in range(0, min(source_row_count, freeparam_group_count * 2), 2)
                                if index + 1 < source_row_count
                            ]
                        elif freeparam_group_count > 0:
                            opd_groups = [([source_rows[i]], 1) for i in range(min(source_row_count, freeparam_group_count))]  # type: ignore[index]
                        elif emit_all_opd_rows and little.irrep_type == 1:
                            opd_groups = [([source_rows[i]], 1) for i in range(source_row_count)]  # type: ignore[index]
                        elif (
                            emit_all_opd_rows
                            and little.irrep_type == 3
                            and source_row_count >= 2
                            and source_row_count % 2 == 0
                        ):
                            opd_groups = [
                                (list(source_rows[index:index + 2]), 2)  # type: ignore[index]
                                for index in range(0, source_row_count, 2)
                            ]
                        elif emit_all_opd_rows:
                            opd_groups = [([source_rows[i]], 1) for i in range(source_row_count)]  # type: ignore[index]
                        elif (
                            source_row_id <= 0
                            and int(little.old_id) > 0
                            and little.irrep_type == 1
                            and source_row_count > 1
                        ):
                            opd_groups = [([source_rows[i]], 1) for i in range(source_row_count)]  # type: ignore[index]
                        elif little.irrep_type == 1:
                            opd_groups = [([source_rows[0]], 1)]  # type: ignore[index]
                        elif (
                            two_row_slots
                            and source_row_count >= 4
                            and case.k_params
                            and little.irrep_type == 3
                            and int(little.full_dim) == 8
                        ):
                            pair_count = max(1, source_row_count // 2)
                            row_start = 2 * (int(source.get("family", 0)) % pair_count)
                            opd_groups = [
                                (list(source_rows[row_start:row_start + 2]), 2)  # type: ignore[index]
                            ]
                        elif two_row_slots and source_row_count >= 4 and case.k_params:
                            opd_groups = [
                                (list(source_rows[index:index + 2]), 2)  # type: ignore[index]
                                for index in range(0, source_row_count, 2)
                            ]
                        elif two_row_slots:
                            opd_groups = [(list(source_rows[:2]), 2)]  # type: ignore[index]
                        elif little.irrep_type == 3 and case.k_params and source_row_count >= 4:
                            opd_groups = [
                                (list(source_rows[index:index + 2]), 2)  # type: ignore[index]
                                for index in range(0, source_row_count, 2)
                            ]
                        elif little.irrep_type == 3 and source_row_count >= 4:
                            opd_groups = [(list(source_rows[0:2]), 2)]  # type: ignore[index]
                        else:
                            opd_groups = [([source_rows[i]], 1) for i in range(source_row_count)]  # type: ignore[index]
                        opd_mode_blocks: list[list[list[list[float]]]] = []
                        for opd_group_index, (opd_rows, opd_row_count) in enumerate(opd_groups):
                            source_178 = (
                                dense_type3_parametric_source(term_count)
                                if (
                                    little.irrep_type == 3
                                    and case.k_params
                                    and int(opd_row_count) == 2
                                    and int(prep.get("site_pg", 0)) == 4
                                    and int(prep.get("vector_basis_id", 0)) == 26
                                )
                                else None
                            )
                            if source_178 is None:
                                source_178 = dense_orderparam_source(
                                    opd_rows,  # type: ignore[arg-type]
                                    term_count,
                                    opd_row_count,
                                )
                            source_modes: list[list[list[float]]] = []
                            for atom_index, record in enumerate(displayed_records):
                                weight_buffer = project_weight_buffer(
                                    little.gid,
                                    record,
                                )
                                direct_coefficients = direct_bridge_coefficients(
                                    project_basis_by_branch,
                                    little.gid,
                                    pg_irrep,
                                    source,
                                    weight_buffer,
                                    int(source["vector_loop_count"]),
                                    int(source["bridge_project_count"]),
                                    181,
                                )
                                coefficients = (
                                    direct_coefficients
                                    if direct_coefficients is not None
                                    else decoder.project_vector_bridge_first_stage(
                                        vector_loop_count=int(source["vector_loop_count"]),
                                        supercell_atom_count=int(source["bridge_project_count"]) + 1,
                                        project_count=int(source["bridge_project_count"]),
                                        source_e8_minus512=dense_source_e8(source),
                                        stack_0490=weight_buffer,
                                    )
                                )
                                basis_function = (
                                    direct_bridge_basis_function(
                                        opd_rows,  # type: ignore[arg-type]
                                        int(source["vector_loop_count"]),
                                        opd_row_count,
                                        int(source["bridge_project_count"]),
                                        coefficients,
                                        360,
                                    )
                                    if direct_coefficients is not None
                                    else decoder.project_vector_bridge_second_stage(
                                        vector_loop_count=int(source["vector_loop_count"]),
                                        atom_count=opd_row_count,
                                        project_count=int(source["bridge_project_count"]),
                                        source_minus256=source_178,
                                        coefficients=coefficients,
                                        output_length=360,
                                    )
                                )
                                project_vector_count = (
                                    int(project_count)
                                    if emit_all_project_families
                                    else 1
                                )
                                emitted, output = decoder.project_vector_from_boundary(
                                    site_pg=row.site_pg,
                                    vector_basis_id=int(prep["vector_basis_id"]),
                                    target_vector_rep=pg_irrep,
                                    atom_count=opd_row_count,
                                    vector_dim=int(prep["vector_dim"]),
                                    vector_setting=int(vector_setting),
                                    point_op=int(record[4]),
                                    project_count=project_vector_count,
                                    basis_function=basis_function,
                                    output_length=720,
                                )
                                split_row_slots = little.irrep_type == 3 and opd_row_count > 1
                                if split_row_slots:
                                    while len(source_modes) < max(0, emitted) * opd_row_count:
                                        source_modes.append([[0.0, 0.0, 0.0] for _ in displayed_records])
                                    for emitted_index in range(max(0, emitted)):
                                        for row_index in range(opd_row_count):
                                            mode_index = emitted_index * opd_row_count + row_index
                                            offset = emitted_index * 144 + 3 * row_index
                                            vector = [float(output[offset + axis]) for axis in range(3)]
                                            source_modes[mode_index][atom_index] = finalize_vector(vector)
                                else:
                                    while len(source_modes) < emitted:
                                        source_modes.append([[0.0, 0.0, 0.0] for _ in displayed_records])
                                    for emitted_index in range(max(0, emitted)):
                                        offset = emitted_index * 144
                                        vector = [float(output[offset + axis]) for axis in range(3)]
                                        source_modes[emitted_index][atom_index] = finalize_vector(vector)
                            opd_mode_blocks.append([_normalize_mode_vectors(mode) for mode in source_modes])
                        mode_block_summaries_by_gid[little.gid].append({
                            "gid": int(little.gid),
                            "label": little.label,
                            "pg_irrep": int(pg_irrep),
                            "source_family": int(source["family"]),
                            "source_row_count": int(source_row_count),
                            "opd_group_count": int(len(opd_groups)),
                            "opd_groups": [
                                {
                                    "index": int(group_index),
                                    "row_count": int(group_row_count),
                                    "rows": [[float(value) for value in row_values] for row_values in group_rows],
                                }
                                for group_index, (group_rows, group_row_count) in enumerate(opd_groups)
                            ],
                            "block_count": int(len(opd_mode_blocks)),
                            "block_mode_counts": [int(len(block)) for block in opd_mode_blocks],
                            "flat_mode_count": int(sum(len(block) for block in opd_mode_blocks)),
                            "emitted_mode_count": None,
                            "project_count": int(project_count),
                            "bridge_project_count": int(source["bridge_project_count"]),
                            "vector_dim": int(prep.get("vector_dim", 0)),
                            "little_type": int(little.irrep_type),
                            "little_full_dim": int(little.full_dim),
                        })
                        if (
                            little.irrep_type == 2
                            and source_row_count > 2
                            and source_row_count % 2 == 0
                            and len(opd_mode_blocks) == source_row_count
                            and opd_mode_blocks
                        ):
                            print_group_count = source_row_count if case.k_params else min(4, source_row_count)
                            emitted_count = min(len(block) for block in opd_mode_blocks[:print_group_count])
                            for mode_index in range(emitted_count):
                                for block in opd_mode_blocks[:print_group_count]:
                                    modes.append(block[mode_index])
                        elif (
                            little.irrep_type == 3
                            and case.k_params
                            and source_row_count > 2
                            and source_row_count % 2 == 0
                            and len(opd_mode_blocks) == source_row_count
                            and opd_mode_blocks
                        ):
                            for block in opd_mode_blocks:
                                modes.extend(block)
                        elif (
                            little.irrep_type != 1
                            and source_row_count > 2
                            and source_row_count % 2 == 0
                            and len(opd_mode_blocks) == source_row_count
                        ):
                            for block_index in range(0, len(opd_mode_blocks), 2):
                                left = opd_mode_blocks[block_index]
                                right = opd_mode_blocks[block_index + 1]
                                for left_mode, right_mode in zip(left, right):
                                    modes.append(_normalize_mode_vectors(_add_mode_vectors(left_mode, right_mode)))
                        elif (
                            little.irrep_type == 3
                            and source_row_count == 2
                            and len(opd_mode_blocks) == 2
                            and len(opd_mode_blocks[0]) == len(opd_mode_blocks[1])
                            and len(opd_mode_blocks[0]) > 1
                        ):
                            for mode_index in range(len(opd_mode_blocks[0])):
                                modes.append(opd_mode_blocks[0][mode_index])
                                modes.append(opd_mode_blocks[1][mode_index])
                        else:
                            for block in opd_mode_blocks:
                                modes.extend(block)
                    if not modes:
                        return []
                    if little.irrep_type == 3:
                        modes = _type3_real_print_modes(modes)
                    return _regroup_repeated_component_modes(modes)

                def local_vector_for(component_start: int, vector_setting: int) -> list[float]:
                    vector_basis_id = int(prep["vector_basis_id"])
                    code_base = (
                        (vector_basis_id - 1) * 18
                        + (vector_setting - 1) * 9
                        + (component_start - 1) * 3
                    )
                    return [
                        float(decoder.iso.const[int(decoder.iso.wyckoff["iwyckoff_pg_vector_basis"][code_base + axis])])
                        for axis in range(3)
                    ]

                def vector_for_terms(
                    terms: tuple[tuple[int, int, int, float], ...],
                ) -> list[list[float]]:
                    atom_records = raw_atom_records
                    uses_column_fallback = decoder.atom_expansion_uses_column_fallback(
                        case.sg,
                        row,
                        case.site_params,
                        atom_basis,
                    )
                    vectors = []
                    for record in atom_records:
                        weight_buffer = project_weight_buffer(little.gid, record)
                        point = decoder.iso.point_ops[int(record[4]) - 1]
                        vector = [0.0, 0.0, 0.0]
                        for component_start, vector_setting, weight_index, coefficient in terms:
                            local_vector = local_vector_for(component_start, vector_setting)
                            weight = float(weight_buffer[weight_index]) if weight_index < len(weight_buffer) else 0.0
                            weight *= coefficient
                            if uses_column_fallback and pg_irrep in {2, 4} and int(record[4]) != 1:
                                weight *= -1.0
                            for output_axis in range(3):
                                vector[output_axis] += sum(
                                    float(point[inner, output_axis])
                                    * weight
                                    * local_vector[inner]
                                    for inner in range(3)
                                )
                        vectors.append(finalize_vector(vector))
                    return _normalize_mode_vectors(vectors)

                recordwise_modes = recordwise_project_vector_modes()
                if recordwise_modes:
                    mode_vectors_by_gid[little.gid].extend(recordwise_modes)
                    continue

                if project_count == 1:
                    basis_item = project_basis_by_branch.get((little.gid, pg_irrep), {})
                    selected_indices = basis_item.get("selected_indices", []) if isinstance(basis_item, dict) else []
                    weight_index = int(selected_indices[0]) if selected_indices else 0
                    if isinstance(basis_item, dict):
                        active_rows = {
                            int(active[0]) // 3
                            for active in basis_item.get("active_values", [])
                            if (
                                isinstance(active, list)
                                and len(active) >= 2
                                and abs(float(active[1])) > 1e-12
                            )
                        }
                        if len(active_rows) > 1 and weight_index == 0:
                            partner_row = max(active_rows)
                            weight_index = partner_row * 48 + partner_row
                    for component_start in matching_components:
                        coefficient = -1.0 if little.full_dim > 1 and weight_index % little.full_dim else 1.0
                        mode_vectors_by_gid[little.gid].append(
                            vector_for_terms(((component_start, 1, weight_index, coefficient),))
                        )
                elif project_count == 2 and int(prep.get("vector_dim", 0)) == 1 and little.full_dim == 2:
                    basis_item = project_basis_by_branch.get((little.gid, pg_irrep), {})
                    family_stride = 144
                    active_component_by_family: dict[int, int] = {}
                    if isinstance(basis_item, dict):
                        for active in basis_item.get("active_values", []):
                            if (
                                isinstance(active, list)
                                and len(active) >= 2
                                and abs(float(active[1])) > 1e-12
                            ):
                                active_index = int(active[0])
                                active_component_by_family[active_index // family_stride] = (
                                    active_index % family_stride
                                ) // 3
                    for family in range(project_count):
                        base_index = family * 48
                        active_component = active_component_by_family.get(family, family % little.full_dim)
                        if active_component == 0:
                            weight_indices = (base_index, base_index + 1)
                        elif active_component == 1:
                            weight_indices = (base_index + 1, base_index)
                        else:
                            weight_indices = (base_index + active_component,)
                        for component_start in matching_components:
                            mode_vectors_by_gid[little.gid].append(
                                vector_for_terms(
                                    tuple((component_start, 1, weight_index, 1.0) for weight_index in weight_indices)
                                )
                            )
                elif project_count == 2 and int(prep.get("vector_dim", 0)) == 2 and little.full_dim == 2:
                    if len(matching_components) >= 2:
                        first_weight = 0
                        second_weight = 49
                        mode_vectors_by_gid[little.gid].append(
                            vector_for_terms((
                                (matching_components[0], 1, first_weight, 1.0),
                                (matching_components[1], 2, second_weight, 1.0),
                            ))
                        )
                        mode_vectors_by_gid[little.gid].append(
                            vector_for_terms((
                                (matching_components[1], 1, first_weight, 1.0),
                                (matching_components[0], 2, second_weight, 1.0),
                            ))
                        )
            if little.irrep_type == 3 and case.k_params:
                mode_vectors_by_gid[little.gid] = _type3_parametric_component_print_order(
                    mode_vectors_by_gid[little.gid]
                )
            mode_vectors_by_gid[little.gid] = _regroup_repeated_component_modes(mode_vectors_by_gid[little.gid])
            if little.irrep_type == 1 and case.k_params:
                mode_vectors_by_gid[little.gid] = _type1_parametric_scalar_plane_print_order(
                    mode_vectors_by_gid[little.gid],
                    allow_support_pairing=int(vector_setting) == 1,
                )
            if little.irrep_type == 3 and case.k_params:
                if 143 <= int(case.sg) <= 194 and int(little.full_dim) == 12:
                    mode_vectors_by_gid[little.gid] = _type3_parametric_kdim2_print_basis(
                        mode_vectors_by_gid[little.gid]
                    )
        else:
            atom_fractionals_by_gid[little.gid] = []
            mode_vectors_by_gid[little.gid] = []
    little_irreps = []
    for little in littles:
        selected_row = selected_isotropy_row_by_gid[little.gid]
        selected_basis = selected_basis_by_gid[little.gid]
        block = {
            "gid": little.gid,
            "old_id": little.old_id,
            "label": little.label,
            "full_dim": little.full_dim,
            "irrep_type": little.irrep_type,
            "lif": little.lif,
        }
        if _diagnostics:
            block["mode_total"] = totals.get(little.label, 0)
            block["isotropy_rows"] = [
                {
                    "row_id": item.row_id,
                    "subgroup": item.subgroup,
                    "arms": item.arms,
                    "direction": item.direction,
                    "basis": list(item.basis),
                    "basis_pml_to_cinter": decoder.format_fraction_matrix(
                        decoder.transform_basis_rows(item.basis, pml_to_cinter)
                    ),
                    "origin": list(item.origin),
                    "new_fractionals": [
                        list(record)
                        for record in decoder.get_new_fractionals_from_basis(item.basis)
                    ],
                }
                for item in isotropy_rows_by_gid[little.gid]
            ]
        block.update(
            {
                "selected_isotropy_row": None
                if selected_basis is None
                else {
                    "row_id": 0 if selected_row is None else selected_row.row_id,
                    "subgroup": 0 if selected_row is None else selected_row.subgroup,
                    "arms": 0 if selected_row is None else selected_row.arms,
                    "direction": 0 if selected_row is None else selected_row.direction,
                    "basis": list(selected_basis),
                    "basis_pml_to_cinter": decoder.format_fraction_matrix(
                        decoder.transform_basis_rows(selected_basis, pml_to_cinter)
                    ),
                    "origin": [] if selected_row is None else list(selected_row.origin),
                },
                "atom_fractionals": atom_fractionals_by_gid[little.gid],
                "atom_operation_records": atom_operation_records_by_gid.get(little.gid, []),
                "atom_order": atom_order_by_gid.get(little.gid, []),
                "mode_vectors": mode_vectors_by_gid[little.gid],
                "mode_block_summaries": mode_block_summaries_by_gid.get(little.gid, []),
                "source_display_projection": source_display_projection_by_gid.get(
                    int(little.gid)
                ),
            }
        )
        little_irreps.append(block)

    if not _diagnostics:
        if len(little_irreps) != 1:
            raise ValueError("production projection requires exactly one selected irrep")
        return {
            "wyckoff": {"site_pg": row.site_pg},
            "little_irreps": little_irreps,
            "project_return_basis_type1": project_basis_type1,
            "project_vector_bridge_source": project_vector_bridge_source,
        }

    result: dict[str, object] = {
        "input": {
            "title": case.title,
            "sg": case.sg,
            "atom_label": case.atom_label,
            "wyckoff": case.wyckoff,
            "k_label": case.k_label,
            "k_direction": case.k_direction,
            "lattice_parameters": list(case.params),
            "site_params": None if case.site_params is None else list(case.site_params),
            "k_params": [str(value) for value in case.k_params],
        },
        "k_vector": [str(value) for value in k_vector],
        "display_k_vector": [str(value) for value in display_k_vector],
        "find_isotropy_orderparam_records": [list(record) for record in find_isotropy_orderparam],
        "newlat_input_records": [list(record) for record in newlat_input],
        "find_isotropy_basis_pml": list(find_isotropy_basis),
        "k_supercell_index": int(k_supercell_index),
        "wyckoff_multiplicity": int(wyckoff_multiplicity),
        "supercell_atom_count": int(supercell_atom_count),
        "setting_transforms": {
            "pml_to_cinter": decoder.format_fraction_matrix(pml_to_cinter),
        },
        "wyckoff": {
            "row_id": row.row_id,
            "offset0": row.offset0,
            "label": row.label,
            "site_pg": row.site_pg,
            "vector_reps": list(
                decoder.site_vector_reps(row.site_pg)[
                    (int(vector_setting) - 1) * 3:(int(vector_setting) - 1) * 3 + 3
                ]
            ),
            "vector_setting": int(vector_setting),
        },
        "little_irreps": little_irreps,
        "mode_totals": totals,
        "project_entry": project_entries,
        "project_local408": project_local408,
        "site_get_irreps": site_get_irreps_trace(project_local408),
        "site_get_ssgn_matrix": site_ssgn,
        "parent_get_irrep4": parent_get_irrep4,
        "project_return_basis_type1": project_basis_type1,
        "project_vector_prep": project_vector_prep,
        "project_vector_bridge_source": project_vector_bridge_source,
        "project_vector_orderparam_source": project_vector_orderparam_source,
    }
    return result


def project_mode_block(
    decoder: ModeDataDecoder,
    case: Case,
    **kwargs: Any,
) -> dict[str, object]:
    """Return the selected projection without constructing diagnostic fields."""

    return pipeline_trace(decoder, case, _diagnostics=False, **kwargs)
