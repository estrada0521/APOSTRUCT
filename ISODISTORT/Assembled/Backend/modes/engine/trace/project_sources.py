"""Build decoded Source intermediates consumed by the mode projection trace."""

from __future__ import annotations

from typing import Any

import numpy as np

from ISODISTORT.Assembled.Backend.modes.engine.input import Case
from ISODISTORT.Assembled.Backend.modes.engine.project.mode_forms import _parametric_type1_real_orderparam_rows
from ISODISTORT.Assembled.Backend.modes.engine.records import WyckoffRow


def _matrix_json(matrix: np.ndarray) -> list[list[float]]:
    return [[float(value) for value in row] for row in np.asarray(matrix).tolist()]


def site_ssgn_matrix_trace(
    decoder: Any,
    row: WyckoffRow,
    project_entries: list[dict[str, int]],
    project_local408: list[dict[str, object]],
) -> list[dict[str, object]]:
    out: list[dict[str, object]] = []
    for entry, local in zip(project_entries, project_local408):
        pg_irrep = int(entry["pg_irrep"])
        try:
            matrices = decoder.site_little_sparse_matrices(row.site_pg, pg_irrep)
            matrix_items = {
                str(op): _matrix_json(matrices[int(op)])
                for op in local["local_ops"]  # type: ignore[index]
            }
            out.append({
                "gid": int(entry["gid"]),
                "pg_irrep": pg_irrep,
                "param5_second": int(entry["param5_second"]),
                "site_old_id": int(local["site_old_id"]),
                "local_ops": [int(op) for op in local["local_ops"]],  # type: ignore[index]
                "matrices": matrix_items,
            })
        except Exception as exc:
            out.append({
                "gid": int(entry["gid"]),
                "pg_irrep": pg_irrep,
                "param5_second": int(entry["param5_second"]),
                "site_old_id": int(local["site_old_id"]),
                "local_ops": [int(op) for op in local.get("local_ops", [])],  # type: ignore[union-attr]
                "error": str(exc),
            })
    return out


def parent_get_irrep4_trace(
    decoder: Any,
    case: Case,
    row: WyckoffRow,
    project_entries: list[dict[str, int]],
) -> list[dict[str, object]]:
    out: list[dict[str, object]] = []
    seen: set[tuple[int, tuple[int, int, int, int, int]]] = set()
    for entry in project_entries:
        gid = int(entry["gid"])
        little = decoder.little_record_by_gid(gid)
        try:
            records = decoder.project_parent_operation_records_for_little_row(case.sg, little, row, case.site_params)
            for record in records:
                record_tuple = tuple(int(x) for x in record)
                key = (gid, record_tuple)
                if key in seen:
                    continue
                seen.add(key)
                if case.k_params:
                    phases = decoder.operation_record_phases_for_case(gid, record_tuple, case)
                    matrix = decoder.little_phase_matrix_by_gid_record_for_case(gid, record_tuple, case)
                else:
                    phases = decoder.operation_record_phases(gid, record_tuple)
                    matrix = decoder.little_phase_matrix_by_gid_record(gid, record_tuple, phases)
                out.append({
                    "gid": gid,
                    "op_record": list(record_tuple),
                    "phases": [str(phase) for phase in phases],
                    "matrix": _matrix_json(matrix),
                })
        except Exception as exc:
            out.append({
                "gid": gid,
                "pg_irrep": int(entry["pg_irrep"]),
                "param5_second": int(entry["param5_second"]),
                "error": str(exc),
            })
    return out


def project_return_basis_type1_trace(
    decoder: Any,
    case: Case,
    row: WyckoffRow,
    project_entries: list[dict[str, int]],
    project_local408: list[dict[str, object]],
    length: int = 720,
) -> list[dict[str, object]]:
    family_stride = 0x480 // 8
    out: list[dict[str, object]] = []
    for entry, local in zip(project_entries, project_local408):
        gid = int(entry["gid"])
        pg_irrep = int(entry["pg_irrep"])
        little = decoder.little_record_by_gid(gid)
        item: dict[str, object] = {
            "gid": gid,
            "pg_irrep": pg_irrep,
            "param5_second": int(entry["param5_second"]),
            "count": int(local["project_count"]),
            "irrep_type": int(little.irrep_type),
            "family_stride": family_stride,
        }
        if little.irrep_type not in {1, 2, 3}:
            item["skipped"] = "irrep_type_not_1_2_or_3"
            out.append(item)
            continue
        try:
            records = decoder.project_parent_operation_records_for_little_row(case.sg, little, row, case.site_params)
            basis_length = max(int(length), int(local["project_count"]) * family_stride)
            if little.irrep_type == 1:
                selection = decoder.project_select_real_blocks_type1(
                    gid,
                    row.site_pg,
                    pg_irrep,
                    records,
                    int(local["project_count"]),
                    case,
                )
            elif little.irrep_type == 2:
                selection = decoder.project_select_real_blocks_type2_from_records(
                    gid,
                    row.site_pg,
                    pg_irrep,
                    records,
                    int(local["project_count"]),
                    case,
                )
            else:
                selection = decoder.project_select_real_blocks_type3(
                    gid,
                    row.site_pg,
                    pg_irrep,
                    records,
                    int(local["project_count"]),
                    case,
                )
            active_values: list[list[object]] = []
            for family, block in enumerate(selection.real_blocks):
                start = family * family_stride
                if start >= basis_length:
                    break
                real_block = block.real
                parent_dim, site_dim = real_block.shape
                for parent_row in range(parent_dim):
                    for site_col in range(site_dim):
                        index = start + parent_row * 3 + site_col
                        if index < basis_length:
                            active_values.append([int(index), float(real_block[parent_row, site_col])])
            item["active_values"] = active_values
            item["selected_indices"] = [int(index) for index in selection.selected_indices]
        except Exception as exc:
            item["error"] = str(exc)
        out.append(item)
    return out


def project_vector_prep_trace(
    decoder: Any,
    case: Case,
    row: WyckoffRow,
    project_entries: list[dict[str, int]],
    project_local408: list[dict[str, object]],
) -> list[dict[str, object]]:
    """Trace MAIN__ values prepared around `project_` for vector expansion."""

    out: list[dict[str, object]] = []
    for entry, local in zip(project_entries, project_local408):
        gid = int(entry["gid"])
        pg_irrep = int(entry["pg_irrep"])
        little = decoder.little_record_by_gid(gid)
        site_old_id = decoder.site_pg_irrep_old_id(row.site_pg, pg_irrep)
        site_dim = decoder.image_record(site_old_id).dimension
        try:
            records = decoder.project_parent_operation_records_for_little_row(
                case.sg,
                little,
                row,
                case.site_params,
            )
            vector_basis_id = decoder.vector_basis_id_for_site_records(row.site_pg, records)
            parent_ops = [int(record[4]) for record in records]
        except Exception as exc:
            vector_basis_id = None
            parent_ops = []
            error = str(exc)
        else:
            error = None
        item: dict[str, object] = {
            "gid": gid,
            "old_id": int(little.old_id),
            "site_pg": int(row.site_pg),
            "pg_irrep": pg_irrep,
            "param5_second": int(entry["param5_second"]),
            "project_count": int(local["project_count"]),
            "site_old_id": int(site_old_id),
            "vector_dim": int(site_dim),
            "vector_basis_id": vector_basis_id,
            "project_param_r8": [pg_irrep, int(entry["param5_second"])],
            "parent_ops": parent_ops,
        }
        if error is not None:
            item["error"] = error
        out.append(item)
    return out


def project_vector_bridge_source_trace(
    decoder: Any,
    project_basis: list[dict[str, object]],
) -> list[dict[str, object]]:
    """Expose MAIN__ bridge source buffers synthesized from `project_` basis."""

    out: list[dict[str, object]] = []
    for item in project_basis:
        if not isinstance(item, dict) or item.get("error") or "active_values" not in item:
            continue
        for source in decoder.project_vector_bridge_source_e8_from_active_values(
            item["active_values"],  # type: ignore[arg-type]
            family_stride=int(item.get("family_stride", 144)),
        ):
            source_values = source["source_e8_minus512"]
            active_values = [
                [index, float(value)]
                for index, value in enumerate(source_values)  # type: ignore[arg-type]
                if abs(float(value)) > 1e-12
            ]
            out.append({
                "gid": int(item["gid"]),
                "pg_irrep": int(item["pg_irrep"]),
                "param5_second": int(item["param5_second"]),
                "family": int(source["family"]),
                "bridge_project_count": int(source["bridge_project_count"]),
                "vector_loop_count": int(source["vector_loop_count"]),
                "source_active_values": active_values,
            })
    return out


def project_vector_orderparam_source_trace(
    decoder: Any,
    littles: list,
    selected_isotropy_row_by_gid: dict[int, object],
    case: Case | None = None,
    orderparam_rows_by_gid: dict[int, list[list[float]]] | None = None,
) -> list[dict[str, object]]:
    """Expose the order-parameter source matrix used before `project_vector_`.

    Fixed-k builds it from the selected isotropy row's order-parameter matrix.
    Parametric-k has no isotropy row_id; the source is instead built from the
    `subgroup_to_orderparam_` direction, whose operations are the local
    `newlat_` supercell basis rows fed as pure translations `{I|row}`.
    """

    out: list[dict[str, object]] = []
    for little in littles:
        row = selected_isotropy_row_by_gid.get(int(little.gid))
        row_id = int(getattr(row, "row_id", 0) or 0)
        override = (orderparam_rows_by_gid or {}).get(int(little.gid))
        if override is not None:
            matrix = [[float(value) for value in values] for values in override]
            row_id = 0
        elif row_id > 0:
            matrix = decoder.isotropy_orderparam_matrix(row_id, int(little.full_dim))
        elif case is not None and case.k_params:
            try:
                basis = decoder.newlat_basis(1, decoder.newlat_input_records_for_case(case))
                identity_op = int(decoder.generate_space_group_records(case.sg)[0][4])
                ops = [
                    (basis[3 * i], basis[3 * i + 1], basis[3 * i + 2], 1, identity_op)
                    for i in range(3)
                ]
                matrices = [decoder.get_irreps_matrix_for_case(int(little.gid), op, case) for op in ops]
                reduced = decoder.xrowop2_like(decoder.subgroup_constraint_matrix_from_irrep_matrices(matrices))
                count, direction = decoder.subgroup_direction_from_reduced_matrix(reduced)
                use_find_isotropy_candidates = bool(case.k_params)
                skip_trial_direction = False
                if (
                    int(little.irrep_type) == 3
                    and int(little.full_dim) == 4
                    and count == 4
                    and not use_find_isotropy_candidates
                ):
                    direction = direction[[0, 2], :]
                    count = 2
                    skip_trial_direction = True
                if count > 1 and not skip_trial_direction:
                    if use_find_isotropy_candidates:
                        candidate_records = decoder.find_isotropy_candidate_operation_records_for_case(case)
                    else:
                        seen_ops = {identity_op}
                        candidate_list: list[tuple[int, int, int, int, int]] = []
                        for record in decoder.generate_space_group_records(case.sg):
                            point_op = int(record[4])
                            if point_op in seen_ops:
                                continue
                            seen_ops.add(point_op)
                            candidate_list.append(record)
                        candidate_records = tuple(candidate_list)
                    accepted_records: list[tuple[int, int, int, int, int]] = []
                    for record in candidate_records:
                        try:
                            if use_find_isotropy_candidates:
                                trial_records = accepted_records + [record]
                                trial_matrices = [
                                    decoder.get_irreps_matrix_for_case(int(little.gid), trial_record, case)
                                    for trial_record in (ops + trial_records)
                                ]
                            else:
                                trial_matrices = matrices + [
                                    decoder.get_irreps_matrix_for_case(int(little.gid), record, case)
                                ]
                        except Exception:
                            continue
                        trial_reduced = decoder.xrowop2_like(
                            decoder.subgroup_constraint_matrix_from_irrep_matrices(trial_matrices)
                        )
                        trial_count, trial_direction = decoder.subgroup_direction_from_reduced_matrix(trial_reduced)
                        if 0 < trial_count < count:
                            if use_find_isotropy_candidates:
                                accepted_records.append(record)
                            count = trial_count
                            direction = trial_direction
                            if count == 1:
                                break
                matrix = [[float(value) for value in dir_row] for dir_row in direction]
                if int(little.irrep_type) == 1:
                    matrix = (
                        _parametric_type1_real_orderparam_rows(matrix, case, int(little.full_dim))
                        or matrix
                    )
            except Exception:
                continue
        else:
            continue
        out.append({
            "gid": int(little.gid),
            "label": little.label,
            "row_id": row_id,
            "full_dim": int(little.full_dim),
            "rows": [[float(value) for value in row_values] for row_values in matrix],
        })
    return out
