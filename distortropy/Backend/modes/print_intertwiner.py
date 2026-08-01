"""Mode presentation intertwiner helpers."""

from __future__ import annotations

from collections import OrderedDict
from fractions import Fraction
import re
from typing import Any
import numpy as np
from distortropy.Backend.modes.engine.decoder import ModeDataDecoder
from distortropy.Backend.modes.engine.input import Case
from distortropy.Backend.modes.engine.subgroup_structure.presentation_transport import (
    child_print_aligned_subduction_basis,
    child_site_irrep_projectors,
    canonical_child_print_factors,
    expand_occurrence_print_columns,
    factor_presentation_occurrences,
    factor_resolved_parent_site_vectors,
    factor_target_row_indices,
    presentation_carrier_coefficients,
    site_vector_print_columns,
    solve_canonical_rank1_print_scalar,
    solve_source_print_intertwiner,
)

from distortropy.Backend.modes.common import (
    _float_matrix_inverse_3,
    _row_multiply,
)
from distortropy.Backend.modes.print_layout import (
    _parallel_source_child_vector,
    _source_family_vector_matrix,
)
from distortropy.Backend.modes.request_context import (
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
