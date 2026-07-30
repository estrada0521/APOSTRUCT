"""Mode presentation intertwiner helpers."""

from __future__ import annotations

from collections import OrderedDict
from dataclasses import replace
from fractions import Fraction
import re
from typing import Any
import gemmi
import numpy as np
from ISODISTORT.Assembled.Backend.modes.engine.decoder import ModeDataDecoder
from ISODISTORT.Assembled.Backend.modes.engine.input import Case
from ISODISTORT.Assembled.Backend.modes.engine.dynamic_subduction import (
    kvec_standard_type1_project_surface,
)
from ISODISTORT.Assembled.Backend.modes.engine.subgroup_structure.presentation_transport import (
    child_print_aligned_subduction_basis,
    child_print_subduction_embedding,
    child_site_irrep_projectors,
    canonical_child_print_factors,
    expand_occurrence_print_columns,
    factor_presentation_occurrences,
    factor_carrier_coefficients_on_target_rows,
    factor_resolved_parent_site_vectors,
    factor_target_row_indices,
    printed_project_block_sources,
    printed_project_blocks_from_trace,
    presentation_carrier_coefficients,
    compose_seitz_records,
    site_vector_print_columns,
    site_print_basis_intertwiner,
    solve_canonical_rank1_print_scalar,
    solve_source_print_intertwiner,
)
from ISODISTORT.Assembled.Backend.modes.presentation import (
    present_mode_rows,
)

from ISODISTORT.Assembled.Backend.modes.common import (
    _float_matrix_inverse_3,
    _origin_vector,
    _row_multiply,
)
from ISODISTORT.Assembled.Backend.modes.print_layout import (
    _parallel_source_child_vector,
    _source_family_vector_matrix,
)
from ISODISTORT.Assembled.Backend.modes.request_context import (
    _pml_vector_to_case_k_params,
)


def _factor_rows_and_target_indices(
    decoder: ModeDataDecoder,
    *,
    sg: int,
    child_sg: int,
    parent_row: Any,
    basis: tuple[int, ...],
    origin: tuple[int, int, int, int],
    rows: list[dict[str, Any]],
) -> tuple[
    tuple[Any, ...],
    tuple[int, ...],
    tuple[tuple[int, int, int, int, int], ...],
    tuple[tuple[float, ...], ...],
]:
    """Share exact occurrence factorizations in a bounded decoder-local LRU."""

    records = tuple(
        tuple(int(value) for value in row["_operation_record"]) for row in rows
    )
    points = tuple(
        tuple(float(value) for value in row.get("xyz") or ()) for row in rows
    )
    key = (
        int(sg),
        int(child_sg),
        parent_row,
        tuple(int(value) for value in basis),
        tuple(int(value) for value in origin),
        records,
        points,
    )
    cache = getattr(decoder, "_factor_rows_and_target_indices_cache", None)
    if cache is None:
        cache = OrderedDict()
        decoder._factor_rows_and_target_indices_cache = cache
    cached = cache.get(key)
    if cached is not None:
        cache.move_to_end(key)
        factors, target_rows = cached
        return factors, target_rows, records, points
    factors = factor_presentation_occurrences(
        decoder,
        parent_sg=int(sg),
        parent_wyckoff_row=parent_row,
        child_sg=int(child_sg),
        subgroup_basis=basis,
        subgroup_origin=origin,
        occurrence_records=records,
    )
    target_rows = factor_target_row_indices(
        decoder,
        child_sg=int(child_sg),
        factors=factors,
        raw_child_points=points,
        occurrence_records=records,
    )
    cache[key] = (factors, target_rows)
    if len(cache) > 64:
        cache.popitem(last=False)
    return factors, target_rows, records, points


def _paired_type1_presentation_project_matrices(
    decoder: ModeDataDecoder,
    *,
    case: Case,
    gid: int,
    pg_irrep: int,
    source_kparam: tuple[int, int, int, int],
    request_k_params: tuple[Fraction, ...] | None,
    project_item: dict[str, Any],
    sources: tuple[Any, ...],
) -> tuple[Case, tuple[np.ndarray, ...]] | None:
    """Join trace blocks to the paired faithful type-1 standard surface.

    The trace remains the authority for the source-side printed blocks.  A
    presentation matrix is returned only when every family retains the exact
    selected candidate and numerical source block exposed by the paired
    projection surface.  Family ordinal then joins the independently evaluated
    source and standard cases in Source print order.
    """

    if not sources or len(source_kparam) != 4 or int(source_kparam[3]) == 0:
        return None
    family_order = tuple(
        dict.fromkeys(int(source.project_family_index) for source in sources)
    )
    if family_order != tuple(range(len(family_order))):
        return None
    project_count = int(project_item.get("count") or len(family_order))
    if project_count != len(family_order):
        return None
    try:
        surface = kvec_standard_type1_project_surface(
            decoder,
            case,
            int(gid),
            int(pg_irrep),
            project_count,
            tuple(int(value) for value in source_kparam),
        )
    except (KeyError, RuntimeError, TypeError, ValueError, np.linalg.LinAlgError):
        return None
    if (
        tuple(surface.provenance.source_kparam) != tuple(source_kparam)
        or tuple(surface.provenance.standard_kparam) == tuple(source_kparam)
        or tuple(case.k_params or ()) != tuple(surface.source_case.k_parameters)
        or len(surface.source_case.selected_indices) != project_count
        or len(surface.standard_case.selected_indices) != project_count
    ):
        return None
    arm_permutation = tuple(int(value) for value in surface.provenance.arm_permutation)
    arm_identity = arm_permutation == tuple(range(len(arm_permutation)))
    arm_conjugated = any(bool(value) for value in surface.provenance.arm_conjugated)
    request_differs = request_k_params is not None and tuple(
        surface.source_case.k_parameters
    ) != tuple(request_k_params)
    if arm_identity and not arm_conjugated and not request_differs:
        return None
    source_blocks = np.asarray(surface.source_case.real_blocks, dtype=complex)
    standard_blocks = np.asarray(surface.standard_case.real_blocks, dtype=complex)
    if (
        source_blocks.shape != standard_blocks.shape
        or source_blocks.shape[0] != project_count
    ):
        return None
    presentation: list[np.ndarray] = []
    for source in sources:
        family = int(source.project_family_index)
        candidate = source.selected_candidate_index
        if (
            candidate is None
            or int(candidate) != int(surface.source_case.selected_indices[family])
            or source_blocks[family].shape != np.asarray(source.project_matrix).shape
            or not np.allclose(
                source_blocks[family],
                np.asarray(source.project_matrix, dtype=complex),
                atol=1e-10,
                rtol=0.0,
            )
        ):
            return None
        presentation.append(np.asarray(standard_blocks[family], dtype=complex))
    presentation_case = replace(
        case, k_params=tuple(surface.standard_case.k_parameters)
    )
    return presentation_case, tuple(presentation)


def _request_presentation_case(
    decoder: ModeDataDecoder,
    *,
    gid: int,
    case: Case,
    spec: dict[str, Any],
    carrier_representative_operation_record: tuple[int, int, int, int, int],
) -> Case | None:
    """Evaluate non-central carrier phases at the requested input parameters."""

    raw = spec.get("request_k_params", spec.get("case_k_params"))
    if raw is None or isinstance(raw, (str, bytes, bytearray)):
        return None
    try:
        requested = tuple(Fraction(str(value)) for value in raw)
    except (TypeError, ValueError, ZeroDivisionError):
        return None
    if len(requested) != len(tuple(case.k_params or ())):
        return None
    presentation_case = replace(case, k_params=requested)
    try:
        source_phases = tuple(
            Fraction(value) % 1
            for value in decoder.operation_record_phases_for_case(
                gid, carrier_representative_operation_record, case
            )
        )
        presentation_phases = tuple(
            Fraction(value) % 1
            for value in decoder.operation_record_phases_for_case(
                gid, carrier_representative_operation_record, presentation_case
            )
        )
    except (IndexError, KeyError, TypeError, ValueError, ZeroDivisionError):
        return None
    if not source_phases or len(source_phases) != len(presentation_phases):
        return None
    phase_quotient = tuple(
        (target - source) % 1
        for source, target in zip(source_phases, presentation_phases, strict=True)
    )
    if len(phase_quotient) > 1 and len(set(phase_quotient)) == 1:
        return None
    return presentation_case


def _rank1_factor_target(
    decoder: ModeDataDecoder,
    *,
    sg: int,
    child_sg: int,
    parent_row: Any,
    basis: tuple[int, ...],
    origin: tuple[int, int, int, int],
    rows: list[dict[str, Any]],
    parent_column: Any,
    pg_irrep: int,
    vector_setting: int,
    presentation_basis: list[list[float]],
    gid: int,
    presentation_case: Case,
    direction_matrix: list[list[float]],
    source_vectors: np.ndarray,
) -> tuple[np.ndarray, tuple[int, ...], tuple[int, ...]] | None:
    """Build a Source-aligned rank-one target from every exact child factor."""
    try:
        factors, target_rows, records, points = _factor_rows_and_target_indices(
            decoder,
            sg=int(sg),
            child_sg=int(child_sg),
            parent_row=parent_row,
            basis=basis,
            origin=origin,
            rows=rows,
        )
        inverse = _float_matrix_inverse_3(presentation_basis)
        final_setting = decoder.cml_to_cinter_matrix(int(sg))
        if inverse is None or final_setting is None:
            return None
        transform = np.asarray(final_setting, dtype=complex) @ np.asarray(
            inverse, dtype=complex
        )
        columns = np.zeros((len(factors), 3, 1), dtype=complex)
        parent_vectors = np.asarray(parent_column.vector, dtype=complex).reshape(3, 1)
        orbit_representatives = sorted(
            {int(factor.orbit_representative_index) for factor in factors}
        )
        for representative_index in orbit_representatives:
            factor_indices = [
                index
                for index, factor in enumerate(factors)
                if int(factor.orbit_representative_index) == representative_index
            ]
            orbit_factors = [factors[index] for index in factor_indices]
            match = re.fullmatch(
                r"\s*\d+\s*([A-Za-z]+)\s*",
                str(rows[representative_index].get("wyckoff_site") or ""),
            )
            if match is None:
                return None
            child_row = next(
                item
                for item in decoder.wyckoff_rows(int(child_sg))
                if item.label == match.group(1)
            )
            child_records = decoder.wyc_pg_elements_records(int(child_sg), child_row)
            projectors = child_site_irrep_projectors(
                decoder,
                parent_site_pg=int(parent_row.site_pg),
                parent_pg_irrep=int(pg_irrep),
                parent_site_records=decoder.wyc_pg_elements_records(
                    int(sg), parent_row
                ),
                child_site_pg=int(child_row.site_pg),
                child_stabilizer_point_ops=tuple(
                    int(record[4]) for record in child_records
                ),
                factors=orbit_factors,
            )
            aligned_candidates = []
            for projector in projectors:
                child_columns = site_vector_print_columns(
                    decoder,
                    site_pg=int(child_row.site_pg),
                    pg_irrep=int(projector.child_pg_irrep),
                    site_operation_records=child_records,
                    vector_setting=int(vector_setting),
                )
                blocks: dict[int, list[Any]] = {}
                for column in child_columns:
                    blocks.setdefault(int(column.component_start), []).append(column)
                for block in blocks.values():
                    if len(block) != int(projector.rank):
                        continue
                    aligned = child_print_aligned_subduction_basis(
                        decoder,
                        parent_sg=int(sg),
                        representative_record=records[representative_index],
                        child_basis_rows_in_parent_cinter=presentation_basis,
                        parent_site_vectors=parent_vectors,
                        projector_basis=projector.basis,
                        child_print_vectors=np.asarray(
                            [column.vector for column in block], dtype=complex
                        ).T,
                        vector_setting=int(vector_setting),
                    )
                    if aligned is not None:
                        aligned_candidates.append(aligned)
            aligned_basis = (
                aligned_candidates[0].basis if len(aligned_candidates) == 1 else None
            )
            if aligned_basis is None:
                return None
            site_vectors = factor_resolved_parent_site_vectors(
                decoder,
                parent_site_pg=int(parent_row.site_pg),
                parent_pg_irrep=int(pg_irrep),
                parent_site_records=decoder.wyc_pg_elements_records(
                    int(sg), parent_row
                ),
                occurrence_records=records,
                factors=orbit_factors,
                parent_site_vectors=parent_vectors,
                child_subduction_basis=aligned_basis,
                vector_setting=int(vector_setting),
            )
            site_vectors = np.einsum("ab,fbc->fac", transform.T, site_vectors)
            carrier = presentation_carrier_coefficients(
                decoder,
                gid=int(gid),
                presentation_records=[
                    factor.presentation_record for factor in orbit_factors
                ],
                presentation_case=presentation_case,
                direction_matrix=direction_matrix,
                source_vectors=source_vectors,
            )
            orbit_columns = expand_occurrence_print_columns(
                carrier,
                site_vectors.transpose(0, 2, 1),
            ).reshape(len(orbit_factors), 3, 1)
            for local_index, factor_index in enumerate(factor_indices):
                columns[factor_index] = orbit_columns[local_index]
    except (IndexError, KeyError, StopIteration, ValueError, np.linalg.LinAlgError):
        return None
    try:
        canonical = canonical_child_print_factors(
            decoder,
            child_sg=int(child_sg),
            factors=factors,
        )
    except (IndexError, KeyError, ValueError):
        return None
    indices = tuple(int(value) for value in canonical.source_factor_indices)
    target = np.asarray([columns[index] for index in indices], dtype=complex).reshape(
        len(indices) * 3,
        1,
    )
    row_indices = tuple(int(target_rows[index]) for index in indices)
    return (
        target,
        row_indices,
        tuple(int(value) for value in canonical.orbit_representative_indices),
    )


def _selected_site_print_component_is_identity(
    decoder: ModeDataDecoder,
    *,
    sg: int,
    parent_row: Any,
    pg_irrep: int,
    vector_setting: int,
    representative_operation_record: tuple[int, int, int, int, int],
    component_start: int,
    expected_canonical_columns: np.ndarray,
) -> bool:
    """Accept a dependent-site fallback only for the literal identity operation."""

    identity_record = tuple(
        int(value) for value in decoder.generate_space_group_records(int(sg))[0]
    )
    if (
        tuple(int(value) for value in representative_operation_record)
        != identity_record
    ):
        return False
    selected_change = site_print_basis_intertwiner(
        decoder,
        parent_sg=int(sg),
        parent_wyckoff_row=parent_row,
        pg_irrep=int(pg_irrep),
        vector_setting=int(vector_setting),
        representative_operation_record=representative_operation_record,
        component_start=int(component_start),
        expected_canonical_columns=expected_canonical_columns,
    )
    if selected_change is None:
        return False
    selected_matrix = np.real_if_close(selected_change.matrix, tol=1000)
    expected_columns = np.asarray(expected_canonical_columns)
    if expected_columns.ndim != 2:
        return False
    expected_width = int(expected_columns.shape[1])
    return bool(
        not np.iscomplexobj(selected_matrix)
        and selected_matrix.shape == (expected_width, expected_width)
        and np.allclose(
            selected_matrix,
            np.eye(expected_width),
            atol=1e-10,
            rtol=0.0,
        )
    )


def _factor_print_block_target(
    decoder: ModeDataDecoder,
    *,
    sg: int,
    child_sg: int,
    parent_row: Any,
    basis: tuple[int, ...],
    origin: tuple[int, int, int, int],
    rows: list[dict[str, Any]],
    pg_irrep: int,
    vector_setting: int,
    presentation_basis: list[list[float]],
    gid: int,
    source_case: Case,
    presentation_case: Case,
    direction_matrix: list[list[float]],
    project_matrix: np.ndarray,
    presentation_project_matrix: np.ndarray | None,
    parent_site_vectors: np.ndarray,
    site_vector_component: int,
    representative_operation_record: tuple[int, int, int, int, int] | None = None,
    carrier_representative_operation_record: tuple[int, int, int, int, int]
    | None = None,
    defer_matrix_site_change: bool = False,
) -> tuple[np.ndarray, tuple[int, ...], np.ndarray] | None:
    """Transport one faithful printed block into the selected child layout."""

    target_project_matrix = (
        np.asarray(project_matrix, dtype=complex)
        if presentation_project_matrix is None
        else np.asarray(presentation_project_matrix, dtype=complex)
    )
    if target_project_matrix.shape != np.asarray(project_matrix).shape:
        return None
    try:
        factors, target_rows, records, points = _factor_rows_and_target_indices(
            decoder,
            sg=int(sg),
            child_sg=int(child_sg),
            parent_row=parent_row,
            basis=basis,
            origin=origin,
            rows=rows,
        )
        target_records = records
        target_factors = factors
        if carrier_representative_operation_record is not None:
            target_factors = tuple(
                replace(
                    factor,
                    presentation_record=compose_seitz_records(
                        decoder,
                        int(sg),
                        first=carrier_representative_operation_record,
                        then=factor.presentation_record,
                    ),
                )
                for factor in factors
            )
        inverse = _float_matrix_inverse_3(presentation_basis)
        final_setting = decoder.cml_to_cinter_matrix(int(sg))
        if inverse is None or final_setting is None:
            return None
        transform = np.asarray(final_setting, dtype=complex) @ np.asarray(
            inverse, dtype=complex
        )
        print_row_count = len(direction_matrix)
        if print_row_count == 0:
            return None
        factor_columns = np.zeros((len(factors), 3, print_row_count), dtype=complex)
        source_factor_columns = np.zeros(
            (len(factors), 3, print_row_count), dtype=complex
        )
        orbit_representatives = sorted(
            {int(factor.orbit_representative_index) for factor in factors}
        )
        for representative_index in orbit_representatives:
            factor_indices = [
                index
                for index, factor in enumerate(factors)
                if int(factor.orbit_representative_index) == representative_index
            ]
            orbit_factors = [factors[index] for index in factor_indices]
            target_orbit_factors = [target_factors[index] for index in factor_indices]
            match = re.fullmatch(
                r"\s*\d+\s*([A-Za-z]+)\s*",
                str(rows[representative_index].get("wyckoff_site") or ""),
            )
            if match is None:
                return None
            child_row = next(
                item
                for item in decoder.wyckoff_rows(int(child_sg))
                if item.label == match.group(1)
            )
            child_records = decoder.wyc_pg_elements_records(int(child_sg), child_row)
            projectors = child_site_irrep_projectors(
                decoder,
                parent_site_pg=int(parent_row.site_pg),
                parent_pg_irrep=int(pg_irrep),
                parent_site_records=decoder.wyc_pg_elements_records(
                    int(sg), parent_row
                ),
                child_site_pg=int(child_row.site_pg),
                child_stabilizer_point_ops=tuple(
                    int(record[4]) for record in child_records
                ),
                factors=orbit_factors,
            )
            embedding_candidates = []
            for projector in projectors:
                child_columns = site_vector_print_columns(
                    decoder,
                    site_pg=int(child_row.site_pg),
                    pg_irrep=int(projector.child_pg_irrep),
                    site_operation_records=child_records,
                    vector_setting=int(vector_setting),
                )
                embedding = child_print_subduction_embedding(
                    decoder,
                    parent_sg=int(sg),
                    representative_record=records[representative_index],
                    child_basis_rows_in_parent_cinter=presentation_basis,
                    parent_site_vectors=parent_site_vectors,
                    projector_basis=projector.basis,
                    child_print_vectors=np.asarray(
                        [column.vector for column in child_columns], dtype=complex
                    ).T,
                    vector_setting=int(vector_setting),
                )
                if embedding is not None:
                    embedding_candidates.append(embedding)
            if len(embedding_candidates) != 1:
                return None
            aligned_basis = np.asarray(embedding_candidates[0].basis, dtype=complex)
            if aligned_basis.shape[1] != project_matrix.shape[1]:
                return None
            site_vectors = factor_resolved_parent_site_vectors(
                decoder,
                parent_site_pg=int(parent_row.site_pg),
                parent_pg_irrep=int(pg_irrep),
                parent_site_records=decoder.wyc_pg_elements_records(
                    int(sg), parent_row
                ),
                occurrence_records=records,
                factors=orbit_factors,
                parent_site_vectors=parent_site_vectors,
                child_subduction_basis=aligned_basis,
                vector_setting=int(vector_setting),
            )
            site_vectors = np.einsum("ab,fbc->fac", transform.T, site_vectors)
            carrier = factor_carrier_coefficients_on_target_rows(
                decoder,
                gid=int(gid),
                factors=target_orbit_factors,
                source_factors=orbit_factors,
                target_row_indices=[
                    target_rows[factor_index] for factor_index in factor_indices
                ],
                target_occurrence_records=target_records,
                presentation_case=presentation_case,
                source_case=source_case,
                direction_matrix=direction_matrix,
                source_vectors=target_project_matrix,
            ).coefficients
            source_carrier = factor_carrier_coefficients_on_target_rows(
                decoder,
                gid=int(gid),
                factors=orbit_factors,
                target_row_indices=[
                    target_rows[factor_index] for factor_index in factor_indices
                ],
                target_occurrence_records=records,
                presentation_case=source_case,
                direction_matrix=direction_matrix,
                source_vectors=project_matrix,
            ).coefficients
            expanded = expand_occurrence_print_columns(
                carrier,
                site_vectors.transpose(0, 2, 1),
            ).reshape(len(orbit_factors), 3, print_row_count, -1)
            if expanded.shape[3] != target_project_matrix.shape[1]:
                return None
            orbit_columns = np.sum(expanded, axis=3)
            source_expanded = expand_occurrence_print_columns(
                source_carrier,
                site_vectors.transpose(0, 2, 1),
            ).reshape(len(orbit_factors), 3, print_row_count, -1)
            source_orbit_columns = np.sum(source_expanded, axis=3)
            for local_index, factor_index in enumerate(factor_indices):
                factor_columns[factor_index] = orbit_columns[local_index]
                source_factor_columns[factor_index] = source_orbit_columns[local_index]
    except (IndexError, KeyError, StopIteration, ValueError, np.linalg.LinAlgError):
        return None
    if representative_operation_record is not None:
        site_change = site_print_basis_intertwiner(
            decoder,
            parent_sg=int(sg),
            parent_wyckoff_row=parent_row,
            pg_irrep=int(pg_irrep),
            vector_setting=int(vector_setting),
            representative_operation_record=representative_operation_record,
        )
        if site_change is None:
            if not _selected_site_print_component_is_identity(
                decoder,
                sg=int(sg),
                parent_row=parent_row,
                pg_irrep=int(pg_irrep),
                vector_setting=int(vector_setting),
                representative_operation_record=representative_operation_record,
                component_start=int(site_vector_component) + 1,
                expected_canonical_columns=parent_site_vectors,
            ):
                return None
            site_change = None
        if site_change is not None and site_change.matrix.shape != (1, 1):
            if defer_matrix_site_change:
                site_change = None
            else:
                return None
        if site_change is None:
            scalar = complex(1.0)
        else:
            scalar = complex(site_change.matrix[0, 0])
        if abs(scalar.imag) > 1e-10:
            return None
        factor_columns *= float(scalar.real)
    try:
        canonical = canonical_child_print_factors(
            decoder,
            child_sg=int(child_sg),
            factors=factors,
        )
    except (IndexError, KeyError, ValueError):
        return None
    indices = tuple(int(value) for value in canonical.source_factor_indices)
    target = np.asarray(
        [factor_columns[index] for index in indices], dtype=complex
    ).reshape(
        len(indices) * 3,
        print_row_count,
    )
    source_target = np.asarray(
        [source_factor_columns[index] for index in indices], dtype=complex
    ).reshape(len(indices) * 3, print_row_count)
    row_indices = tuple(int(target_rows[index]) for index in indices)
    return target, row_indices, source_target


def _rank1_source_print_intertwiner(
    decoder: ModeDataDecoder,
    *,
    sg: int,
    child_sg: int | None,
    case: Case,
    spec: dict[str, Any],
    trace: dict[str, Any],
    block: dict[str, Any],
    site: dict[str, Any],
    rows: list[dict[str, Any]],
    identity: dict[str, Any] | None,
    presentation_basis: list[list[float]] | None,
    vector_setting: int,
    subgroup_basis: tuple[int, ...] | None,
    subgroup_origin: tuple[int, int, int, int] | None,
) -> float | None:
    """Solve a nondegenerate rank-one Source printed-column gauge."""

    if (
        identity is None
        or child_sg is None
        or presentation_basis is None
        or int(identity.get("family_width") or 0) != 1
        or not spec.get("source_numeric_rows")
    ):
        return None
    records = [row.get("_operation_record") for row in rows]
    if not records or any(
        not isinstance(record, (list, tuple)) or len(record) != 5 for record in records
    ):
        return None
    gid = int(identity.get("gid") or 0)
    pg_irrep = int(identity.get("pg_irrep") or 0)
    full_dim = int(block.get("full_dim") or 0)
    source_vectors = _source_family_vector_matrix(
        trace,
        gid=gid,
        pg_irrep=pg_irrep,
        family=int(identity.get("family") or 0),
        full_dim=full_dim,
    )
    if source_vectors is None or source_vectors.shape[1] != 1:
        return None
    try:
        parent_row = next(
            item
            for item in decoder.wyckoff_rows(int(sg))
            if item.label == str(site["wyckoff"])
        )
    except (KeyError, StopIteration):
        return None
    parent_columns = site_vector_print_columns(
        decoder,
        site_pg=int(parent_row.site_pg),
        pg_irrep=pg_irrep,
        site_operation_records=decoder.wyc_pg_elements_records(int(sg), parent_row),
        vector_setting=int(vector_setting),
    )
    if len(parent_columns) != 1:
        return None
    inverse = _float_matrix_inverse_3(presentation_basis)
    if inverse is None:
        return None
    final_setting = decoder.cml_to_cinter_matrix(int(sg))
    if final_setting is None:
        return None
    reference = np.asarray(
        _row_multiply(
            _row_multiply(list(parent_columns[0].vector), final_setting),
            inverse,
        ),
        dtype=float,
    )
    target_site_vectors: list[np.ndarray] | None = []
    for row in rows:
        target = _parallel_source_child_vector(
            decoder,
            child_sg=int(child_sg),
            child_wyckoff_site=str(row.get("wyckoff_site") or ""),
            reference=reference,
            vector_setting=int(vector_setting),
        )
        if target is None:
            target_site_vectors = None
            break
        target_site_vectors.append(target.reshape(1, 3))
    reciprocal = spec.get("reciprocal_vector_pml")
    if not isinstance(reciprocal, (list, tuple)) or len(reciprocal) != 3:
        return None
    presentation_params = _pml_vector_to_case_k_params(
        decoder,
        gid=gid,
        reciprocal_vector_pml=tuple(Fraction(value) for value in reciprocal),
    )
    if presentation_params is None:
        return None
    presentation_case = Case(
        sg=case.sg,
        wyckoff=case.wyckoff,
        k_label=case.k_label,
        params=case.params,
        title=case.title,
        atom_label=case.atom_label,
        k_direction=case.k_direction,
        site_params=case.site_params,
        k_params=tuple(presentation_params),
    )
    factor_target = (
        _rank1_factor_target(
            decoder,
            sg=int(sg),
            child_sg=int(child_sg),
            parent_row=parent_row,
            basis=subgroup_basis,
            origin=subgroup_origin,
            rows=rows,
            parent_column=parent_columns[0],
            pg_irrep=pg_irrep,
            vector_setting=int(vector_setting),
            presentation_basis=presentation_basis,
            gid=gid,
            presentation_case=presentation_case,
            direction_matrix=spec.get("source_numeric_rows") or [],
            source_vectors=source_vectors,
        )
        if subgroup_basis is not None and subgroup_origin is not None
        else None
    )
    if factor_target is not None:
        target, row_indices, orbit_indices = factor_target
        scale = float(np.max(np.abs(target))) if target.size else 0.0
        if scale <= 1e-12 or target.shape != (len(orbit_indices) * 3, 1):
            return None
        try:
            current_columns = [
                np.asarray(rows[index].get("dxyz") or (), dtype=complex)
                for index in row_indices
            ]
        except IndexError:
            return None
        target_columns = [
            column.reshape(3)
            for column in (target / scale).reshape(len(orbit_indices), 3, 1)
        ]
        solved = solve_canonical_rank1_print_scalar(
            current_columns,
            target_columns,
            orbit_representative_indices=orbit_indices,
        )
        return None if solved is None else float(solved.scalar)
    target: np.ndarray
    if target_site_vectors is None:
        return None
    try:
        carrier = presentation_carrier_coefficients(
            decoder,
            gid=gid,
            presentation_records=records,  # type: ignore[arg-type]
            presentation_case=presentation_case,
            direction_matrix=spec.get("source_numeric_rows") or [],
            source_vectors=source_vectors,
        )
        target = expand_occurrence_print_columns(
            carrier, np.asarray(target_site_vectors)
        )
    except (ValueError, np.linalg.LinAlgError):
        return None
    scale = float(np.max(np.abs(target))) if target.size else 0.0
    if scale <= 1e-12:
        return None
    target = target / scale
    source = np.asarray(
        [float(value) for row in rows for value in row.get("dxyz") or ()],
        dtype=complex,
    ).reshape(-1, 1)
    solved = solve_source_print_intertwiner(source, target)
    if solved is None or solved.matrix.shape != (1, 1):
        return None
    value = complex(solved.matrix[0, 0])
    return float(value.real) if abs(value.imag) <= 1e-10 else None


def _remix_site_print_columns(
    target_columns: list[np.ndarray],
    sources: list[Any],
    print_rows: list[int],
    site_matrix: np.ndarray,
) -> list[np.ndarray] | None:
    """Apply a Source site-basis intertwiner within each family/OPD row."""

    if len(target_columns) != len(sources) or len(sources) != len(print_rows):
        return None
    remixed_columns = list(target_columns)
    for family in dict.fromkeys(int(source.project_family_index) for source in sources):
        for print_row in dict.fromkeys(print_rows):
            positions = [
                index
                for index, (source, row_index) in enumerate(
                    zip(sources, print_rows, strict=True)
                )
                if int(source.project_family_index) == family
                and int(row_index) == int(print_row)
            ]
            positions.sort(
                key=lambda index: (
                    int(sources[index].vector_component),
                    int(sources[index].printed_block_index),
                )
            )
            components = [int(sources[index].vector_component) for index in positions]
            if len(positions) != site_matrix.shape[0] or len(set(components)) != len(
                components
            ):
                return None
            block_matrix = np.asarray(
                [remixed_columns[index] for index in positions],
                dtype=complex,
            ).T
            remixed = block_matrix @ site_matrix
            for column, position in enumerate(positions):
                remixed_columns[position] = remixed[:, column]
    return remixed_columns


def _multi_source_print_modes(
    decoder: ModeDataDecoder,
    *,
    sg: int,
    child_sg: int | None,
    case: Case,
    spec: dict[str, Any],
    trace: dict[str, Any],
    block: dict[str, Any],
    site: dict[str, Any],
    mode_vectors: list[Any],
    identities: list[dict[str, Any]] | None,
    atom_fractionals: list[Any],
    atom_operation_records: list[Any],
    presentation_basis: list[list[float]] | None,
    child_origin: Any,
    vector_setting: int,
    subgroup_basis: tuple[int, ...] | None,
    subgroup_origin: tuple[int, int, int, int] | None,
    representative_operation_record: tuple[int, int, int, int, int] | None = None,
    carrier_representative_operation_record: tuple[int, int, int, int, int]
    | None = None,
    allow_single_columns: bool = False,
    use_type1_project_surface: bool = False,
) -> tuple[list[Any], set[int], np.ndarray | None]:
    """Apply Source-only multi-column print intertwiners before layout joins."""

    unchanged = (mode_vectors, set(), None)
    if (
        child_sg is None
        or presentation_basis is None
        or subgroup_basis is None
        or subgroup_origin is None
        or identities is None
        or len(identities) != len(mode_vectors)
        or not spec.get("source_numeric_rows")
    ):
        return unchanged
    gid = int(spec.get("gid") or 0)
    full_dim = int(block.get("full_dim") or 0)
    site_pg = int((trace.get("wyckoff") or {}).get("site_pg") or 0)
    if gid <= 0 or full_dim <= 0 or site_pg <= 0:
        return unchanged
    try:
        parent_row = next(
            item
            for item in decoder.wyckoff_rows(int(sg))
            if item.label == str(site["wyckoff"])
        )
        trace_blocks = printed_project_blocks_from_trace(
            decoder,
            gid=gid,
            site_pg=site_pg,
            vector_setting=int(vector_setting),
            project_return_basis_type1=trace.get("project_return_basis_type1") or [],
        )
        project_item_by_irrep = {
            int(item.get("pg_irrep") or 0): item
            for item in trace.get("project_return_basis_type1") or []
            if isinstance(item, dict) and int(item.get("gid") or 0) == gid
        }
        presentation_case = (
            _request_presentation_case(
                decoder,
                gid=gid,
                case=case,
                spec=spec,
                carrier_representative_operation_record=(
                    carrier_representative_operation_record
                ),
            )
            if carrier_representative_operation_record is not None
            else case
        )
        if presentation_case is None:
            return unchanged
        child_symbol = gemmi.find_spacegroup_by_number(int(child_sg)).hm
        presented_modes: list[list[dict[str, Any]]] = []
        raw_mode_rows: list[list[dict[str, Any]]] = []
        for vectors in mode_vectors:
            if not isinstance(vectors, list):
                return unchanged
            raw_rows = [
                {
                    "xyz": [
                        float(Fraction(str(value))) for value in atom_fractionals[index]
                    ],
                    "dxyz": [float(value) for value in vectors[index]],
                    "_operation_record": atom_operation_records[index],
                    "_source_raw_index": index,
                }
                for index in range(
                    min(
                        len(atom_fractionals), len(vectors), len(atom_operation_records)
                    )
                )
                if isinstance(vectors[index], list)
            ]
            raw_mode_rows.append(raw_rows)
            presented_modes.append(
                list(
                    present_mode_rows(
                        raw_rows,
                        basis=presentation_basis,
                        origin=_origin_vector(child_origin),
                        centering_symbol=child_symbol,
                    )["rows"]
                )
            )
        if not presented_modes or any(
            len(rows) != len(presented_modes[0]) for rows in presented_modes
        ):
            return unchanged
        calibration_rows = raw_mode_rows[0]
        if any(row.get("_operation_record") is None for row in calibration_rows):
            return unchanged
    except (IndexError, KeyError, StopIteration, TypeError, ValueError):
        return unchanged

    transformed = [
        [[float(value) for value in vector] for vector in vectors]
        for vectors in mode_vectors
    ]
    applied: set[int] = set()
    overall_matrix = np.eye(len(mode_vectors), dtype=float)
    parent_identity_point_op = int(decoder.generate_space_group_records(int(sg))[0][4])
    for pg_irrep in dict.fromkeys(
        int(identity.get("pg_irrep") or 0) for identity in identities
    ):
        mode_indices = [
            index
            for index, identity in enumerate(identities)
            if int(identity.get("pg_irrep") or 0) == pg_irrep
        ]
        if not mode_indices or (
            len(mode_indices) == 1
            and parent_identity_point_op == 1
            and not allow_single_columns
        ):
            continue
        item = project_item_by_irrep.get(pg_irrep)
        blocks_for_irrep = [
            block_item
            for block_item in trace_blocks
            if block_item.target_vector_rep == pg_irrep
        ]
        if item is None or not blocks_for_irrep:
            continue
        try:
            sources = printed_project_block_sources(
                decoder,
                site_pg=site_pg,
                pg_irrep=pg_irrep,
                site_operation_records=decoder.wyc_pg_elements_records(
                    int(sg), parent_row
                ),
                vector_setting=int(vector_setting),
                full_dim=full_dim,
                project_basis_active_values=item.get("active_values") or [],
                printed_blocks=blocks_for_irrep,
                family_stride=int(item.get("family_stride") or 144),
            )
            paired_type1 = (
                _paired_type1_presentation_project_matrices(
                    decoder,
                    case=case,
                    gid=gid,
                    pg_irrep=pg_irrep,
                    source_kparam=tuple(
                        int(value) for value in spec.get("source_kparam") or ()
                    ),
                    request_k_params=(
                        None
                        if spec.get("request_k_params", spec.get("case_k_params"))
                        is None
                        else tuple(
                            Fraction(str(value))
                            for value in spec.get(
                                "request_k_params", spec.get("case_k_params")
                            )
                            or ()
                        )
                    ),
                    project_item=item,
                    sources=sources,
                )
                if use_type1_project_surface
                else None
            )
            if paired_type1 is None:
                irrep_presentation_case = presentation_case
                presentation_project_matrices = tuple(
                    np.asarray(source.project_matrix, dtype=complex)
                    for source in sources
                )
            else:
                irrep_presentation_case, presentation_project_matrices = paired_type1
            target_columns = []
            source_block_columns = []
            expanded_sources = []
            expanded_print_rows = []
            site_change_matrix: np.ndarray | None = None
            if representative_operation_record is not None and int(
                decoder.little_record_by_gid(gid).irrep_type
            ) in {1, 3}:
                site_change = site_print_basis_intertwiner(
                    decoder,
                    parent_sg=int(sg),
                    parent_wyckoff_row=parent_row,
                    pg_irrep=int(pg_irrep),
                    vector_setting=int(vector_setting),
                    representative_operation_record=representative_operation_record,
                )
                if site_change is not None:
                    matrix = np.real_if_close(site_change.matrix, tol=1000)
                    if (
                        not np.iscomplexobj(matrix)
                        and matrix.ndim == 2
                        and matrix.shape[0] == matrix.shape[1]
                        and matrix.shape[0] > 1
                    ):
                        site_change_matrix = np.asarray(matrix, dtype=complex)
            canonical_rows: tuple[int, ...] | None = None
            for source_index, source in enumerate(sources):
                target_result = _factor_print_block_target(
                    decoder,
                    sg=int(sg),
                    child_sg=int(child_sg),
                    parent_row=parent_row,
                    basis=subgroup_basis,
                    origin=subgroup_origin,
                    rows=calibration_rows,
                    pg_irrep=pg_irrep,
                    vector_setting=int(vector_setting),
                    presentation_basis=presentation_basis,
                    gid=gid,
                    source_case=case,
                    presentation_case=irrep_presentation_case,
                    direction_matrix=spec.get("source_numeric_rows") or [],
                    project_matrix=np.asarray(source.project_matrix, dtype=complex),
                    presentation_project_matrix=presentation_project_matrices[
                        source_index
                    ],
                    parent_site_vectors=np.asarray(source.site_vectors, dtype=complex),
                    site_vector_component=int(source.vector_component),
                    representative_operation_record=representative_operation_record,
                    carrier_representative_operation_record=(
                        None
                        if paired_type1 is not None
                        else carrier_representative_operation_record
                    ),
                    defer_matrix_site_change=site_change_matrix is not None,
                )
                if target_result is None:
                    target_columns = []
                    break
                target, target_rows, source_block = target_result
                if canonical_rows is None:
                    canonical_rows = target_rows
                elif canonical_rows != target_rows:
                    target_columns = []
                    break
                for print_row in range(target.shape[1]):
                    target_columns.append(target[:, print_row])
                    source_block_columns.append(source_block[:, print_row])
                    expanded_sources.append(source)
                    expanded_print_rows.append(print_row)
            if not target_columns or canonical_rows is None:
                continue
            if site_change_matrix is not None:
                remixed = _remix_site_print_columns(
                    target_columns,
                    expanded_sources,
                    expanded_print_rows,
                    site_change_matrix,
                )
                if remixed is None:
                    continue
                target_columns = remixed
            sources = tuple(expanded_sources)
            desired_by_family: dict[int, int] = {}
            for index in mode_indices:
                family = int(identities[index].get("family") or 0)
                desired_by_family[family] = desired_by_family.get(family, 0) + 1
            contracted_targets: list[np.ndarray] = []
            contracted_sources: list[Any] = []
            contraction_valid = True
            for family in dict.fromkeys(
                source.project_family_index for source in sources
            ):
                positions = [
                    index
                    for index, source in enumerate(sources)
                    if int(source.project_family_index) == int(family)
                ]
                desired = desired_by_family.get(int(family), 0)
                if desired == len(positions):
                    for position in positions:
                        contracted_targets.append(target_columns[position])
                        contracted_sources.append(sources[position])
                elif desired == 1 and positions:
                    source_blocks = np.asarray(
                        [source_block_columns[position] for position in positions],
                        dtype=complex,
                    ).T
                    family_mode_indices = [
                        index
                        for index in mode_indices
                        if int(identities[index].get("family") or 0) == int(family)
                    ]
                    local_family = np.asarray(
                        [
                            [
                                complex(float(value))
                                for row_index in canonical_rows
                                for value in next(
                                    row
                                    for row in presented_modes[index]
                                    if int(row.get("_source_raw_index", -1))
                                    == int(row_index)
                                ).get("dxyz")
                                or ()
                            ]
                            for index in family_mode_indices
                        ],
                        dtype=complex,
                    ).T
                    contraction, *_ = np.linalg.lstsq(
                        source_blocks, local_family, rcond=None
                    )
                    if (
                        contraction.shape != (len(positions), 1)
                        or np.linalg.matrix_rank(source_blocks, tol=1e-10) < 1
                        or float(
                            np.max(np.abs(source_blocks @ contraction - local_family))
                        )
                        > 1e-10
                    ):
                        contraction_valid = False
                        break
                    target_blocks = np.asarray(
                        [target_columns[position] for position in positions],
                        dtype=complex,
                    ).T
                    contracted_targets.append((target_blocks @ contraction).reshape(-1))
                    contracted_sources.append(sources[positions[0]])
                else:
                    contraction_valid = False
                    break
            if not contraction_valid:
                continue
            target_columns = contracted_targets
            sources = tuple(contracted_sources)
            target_all = np.asarray(target_columns, dtype=complex).T
            if target_all.shape[1] != len(mode_indices):
                continue
            target_basis = target_all
            source_matrix = np.asarray(
                [
                    [
                        complex(float(value))
                        for row_index in canonical_rows
                        for value in next(
                            row
                            for row in presented_modes[index]
                            if int(row.get("_source_raw_index", -1)) == int(row_index)
                        ).get("dxyz")
                        or ()
                    ]
                    for index in mode_indices
                ],
                dtype=complex,
            ).T
            solved = solve_source_print_intertwiner(source_matrix, target_basis)
            if solved is None or solved.matrix.shape != (
                len(mode_indices),
                len(mode_indices),
            ):
                continue
            matrix = np.real_if_close(solved.matrix, tol=1000)
            if np.iscomplexobj(matrix):
                continue
            original = [
                np.asarray(mode_vectors[index], dtype=float) for index in mode_indices
            ]
            for target_local, target_index in enumerate(mode_indices):
                values = sum(
                    original[source_local] * float(matrix[source_local, target_local])
                    for source_local in range(len(mode_indices))
                )
                transformed[target_index] = values.tolist()
            overall_matrix[np.ix_(mode_indices, mode_indices)] = np.asarray(
                matrix, dtype=float
            )
            applied.update(mode_indices)
        except (IndexError, KeyError, TypeError, ValueError, np.linalg.LinAlgError):
            continue
    return transformed, applied, (overall_matrix if applied else None)
