"""Generate and select site-irrep projection blocks."""

from __future__ import annotations

from typing import Iterable

import numpy as np

from APOSTRUCT.Backend.modes.engine.input import Case
from APOSTRUCT.Backend.modes.engine.records import ProjectCandidateRows, ProjectSelection


class ProjectSelectionMixin:
    def project_real_blocks_from_records(
        self,
        gid: int,
        site_pg: int,
        pg_irrep: int,
        parent_records: Iterable[tuple[int, int, int, int, int]],
        case: Case | None = None,
    ) -> ProjectCandidateRows:
        """Build normalized real candidate blocks using mode-kernel full records.

        This is the common numerical core before the type-specific conversion
        into the complex `zrowop2_` work matrix.  Parent matrices are read via
        `get_irrep4_` semantics: full operation record, generated-operation
        sparse matrix, and k-star phase.  Site matrices use the local setting
        selected by the site operation records.  Their matrix orientation is
        the transpose consumed by the type-2 reduction below.
        """

        little = self.little_record_by_gid(gid)
        records = tuple(tuple(int(x) for x in record) for record in parent_records)
        site_records = self.map_parent_ops_to_site_project_records(site_pg, [record[4] for record in records])
        site_ops = tuple(int(record[4]) for record in site_records)
        site_mats = self.site_little_sparse_matrices(site_pg, pg_irrep)
        parent_dim = little.full_dim
        site_old_id = self.site_pg_irrep_old_id(site_pg, pg_irrep)
        site_dim = self.image_record(site_old_id).dimension

        parent_mats = []
        for record in records:
            if case is not None and case.k_params:
                parent_mats.append(self.little_phase_matrix_by_gid_record_for_case(gid, record, case))
            else:
                phases = self.operation_record_phases(gid, record)
                try:
                    parent_mats.append(self.little_phase_matrix_by_gid_record(gid, record, phases))
                except KeyError:
                    parent_mats.append(self.little_phase_matrix_by_gid_op(gid, record[4], phases))
        blocks = []
        for candidate_row in range(parent_dim):
            accum = np.zeros((parent_dim, site_dim), dtype=complex)
            for parent_matrix, site_op in zip(parent_mats, site_ops):
                accum += np.outer(
                    parent_matrix[candidate_row, :],
                    site_mats[site_op][0, :],
                )
            scale = np.max(np.abs(accum))
            if scale > 1e-12:
                accum = accum / scale
            blocks.append(accum)
        return ProjectCandidateRows(
            parent_ops=tuple(int(record[4]) for record in records),
            site_ops=site_ops,
            pg_irrep=pg_irrep,
            rows=np.array(blocks, dtype=complex),
        )

    def project_select_real_blocks_type1(
        self,
        gid: int,
        site_pg: int,
        pg_irrep: int,
        parent_records: Iterable[tuple[int, int, int, int, int]],
        project_count: int,
        case: Case | None = None,
        tol: float = 1e-10,
    ) -> ProjectSelection:
        """Select independent real blocks for a type-1 little irrep."""

        candidates = self.project_real_blocks_from_records(gid, site_pg, pg_irrep, parent_records, case)
        kept_rows: np.ndarray | None = None
        selected_indices: list[int] = []
        selected_blocks: list[np.ndarray] = []
        for index, block in enumerate(candidates.rows):
            if np.max(np.abs(block)) <= tol:
                continue
            trial_row = block.reshape(1, -1)
            trial = trial_row if kept_rows is None else np.vstack([kept_rows, trial_row])
            reduced = self.zrowop2_like(trial, tol=tol)
            independent = bool(np.max(np.abs(reduced[-1, :].real) + np.abs(reduced[-1, :].imag)) > tol)
            if not independent:
                continue
            kept_rows = reduced
            selected_indices.append(index)
            selected_blocks.append(block)
            if len(selected_blocks) == project_count:
                break
        if len(selected_blocks) != project_count:
            raise RuntimeError(
                f"project_ type-1 selection found {len(selected_blocks)} blocks, expected {project_count}"
            )
        return ProjectSelection(
            parent_ops=candidates.parent_ops,
            site_ops=candidates.site_ops,
            pg_irrep=pg_irrep,
            selected_indices=tuple(selected_indices),
            real_blocks=np.array(selected_blocks, dtype=complex),
            reduced_rows=np.array(kept_rows, dtype=complex),
        )

    def project_select_real_blocks_type3(
        self,
        gid: int,
        site_pg: int,
        pg_irrep: int,
        parent_records: Iterable[tuple[int, int, int, int, int]],
        project_count: int,
        case: Case | None = None,
        tol: float = 1e-10,
    ) -> ProjectSelection:
        """Select independent real blocks for a type-3 little irrep.

        Return the selected normalized real blocks while using the transformed
        row-reduction state only to test independence.  The support gates below
        preserve distinct Source block representatives before that reduction.
        """

        candidates = self.project_real_blocks_from_records(gid, site_pg, pg_irrep, parent_records, case)
        _transform, inverse = self.transform_irrep2complex_matrices(gid)
        first_block_count = self.little_transform_block_count(gid)
        support_keys: list[tuple[int, ...]] = []
        use_block_support = False
        if len(candidates.rows) and first_block_count > 1:
            support_keys = []
            support_ok = True
            for candidate_block in candidates.rows:
                half_dim0 = candidate_block.shape[0] // 2
                candidate_first = candidate_block.T @ inverse[: candidate_block.shape[0], :half_dim0]
                flat = candidate_first.reshape(-1)
                if len(flat) % first_block_count:
                    support_ok = False
                    break
                block_width0 = len(flat) // first_block_count
                support = tuple(
                    block_index
                    for block_index in range(first_block_count)
                    if np.max(
                        np.abs(flat[block_index * block_width0:(block_index + 1) * block_width0].real)
                        + np.abs(flat[block_index * block_width0:(block_index + 1) * block_width0].imag)
                    )
                    > tol
                )
                support_keys.append(support)
            nonempty_supports = {support for support in support_keys if support}
            use_block_support = (
                support_ok
                and len(nonempty_supports) >= int(project_count)
            )
        singleton_support_keys: list[int] = []
        use_singleton_block_support = False
        if len(candidates.rows) and candidates.rows[0].shape[1] == 1 and first_block_count > 1:
            half_dim0 = candidates.rows[0].shape[0] // 2
            if half_dim0 % first_block_count == 0:
                block_width0 = half_dim0 // first_block_count
                singleton_support_keys = []
                singleton_supports = True
                for candidate_block in candidates.rows:
                    candidate_first = candidate_block.T @ inverse[: candidate_block.shape[0], :half_dim0]
                    flat = candidate_first.reshape(-1)
                    support = [
                        block_index
                        for block_index in range(first_block_count)
                        if np.max(
                            np.abs(flat[block_index * block_width0:(block_index + 1) * block_width0].real)
                            + np.abs(flat[block_index * block_width0:(block_index + 1) * block_width0].imag)
                        )
                        > tol
                    ]
                    if len(support) != 1:
                        singleton_supports = False
                        break
                    singleton_support_keys.append(support[0])
                nonempty_singleton_supports = set(singleton_support_keys)
                use_singleton_block_support = (
                    singleton_supports
                    and len(nonempty_singleton_supports) >= int(project_count)
                )
        if support_keys:
            exact_support_indices: dict[tuple[int, ...], int] = {}
            for index, candidate_block in enumerate(candidates.rows):
                if np.max(np.abs(candidate_block)) <= tol:
                    continue
                half_dim0 = candidate_block.shape[0] // 2
                candidate_first = candidate_block.T @ inverse[: candidate_block.shape[0], :half_dim0]
                exact_support = tuple(
                    support_index
                    for support_index, value in enumerate(candidate_first.reshape(-1))
                    if abs(value.real) + abs(value.imag) > tol
                )
                if not exact_support:
                    continue
                exact_support_indices.setdefault(exact_support, index)
            if len(exact_support_indices) >= int(project_count):
                selected_indices = []
                selected_blocks = []
                reduced_rows: np.ndarray | None = None
                for index in exact_support_indices.values():
                    block = candidates.rows[index]
                    transformed_row = (
                        block.T
                        @ inverse[: block.shape[0], : block.shape[0] // 2]
                    ).reshape(1, -1)
                    trial = (
                        transformed_row
                        if reduced_rows is None
                        else np.vstack([reduced_rows, transformed_row])
                    )
                    reduced = self.zrowop2_like(trial, tol=tol)
                    rank = int(
                        np.count_nonzero(
                            np.max(
                                np.abs(reduced.real) + np.abs(reduced.imag), axis=1
                            )
                            > tol
                        )
                    )
                    if rank <= len(selected_indices):
                        continue
                    selected_indices.append(index)
                    selected_blocks.append(block)
                    reduced_rows = reduced
                    if len(selected_indices) == int(project_count):
                        break
                if len(selected_indices) == int(project_count):
                    return ProjectSelection(
                        parent_ops=candidates.parent_ops,
                        site_ops=candidates.site_ops,
                        pg_irrep=pg_irrep,
                        selected_indices=tuple(selected_indices),
                        real_blocks=np.array(selected_blocks, dtype=complex),
                        reduced_rows=np.asarray(reduced_rows, dtype=complex),
                    )
            grouped_supports: list[tuple[tuple[int, ...], list[int]]] = []
            for index, support in enumerate(support_keys):
                if not support:
                    continue
                for existing_support, indices in grouped_supports:
                    if existing_support == support:
                        indices.append(index)
                        break
                else:
                    grouped_supports.append((support, [index]))
            if len(grouped_supports) > 1:
                site_dim = int(candidates.rows[0].shape[1])
                distributed_indices: list[int] = []
                for support, indices in grouped_supports:
                    take = 1 if site_dim > 1 else len(support)
                    distributed_indices.extend(indices[:take])
                if len(distributed_indices) == int(project_count):
                    selected_blocks = [candidates.rows[index] for index in distributed_indices]
                    return ProjectSelection(
                        parent_ops=candidates.parent_ops,
                        site_ops=candidates.site_ops,
                        pg_irrep=pg_irrep,
                        selected_indices=tuple(distributed_indices),
                        real_blocks=np.array(selected_blocks, dtype=complex),
                        reduced_rows=np.empty((0, 0), dtype=complex),
                    )
        kept_rows: np.ndarray | None = None
        selected_indices: list[int] = []
        selected_blocks: list[np.ndarray] = []
        for index, block in enumerate(candidates.rows):
            if np.max(np.abs(block)) <= tol:
                continue
            half_dim = block.shape[0] // 2
            first = block.T @ inverse[: block.shape[0], :half_dim]
            block_count = self.little_transform_block_count(gid)
            if use_block_support:
                candidate_row = np.zeros((1, first_block_count), dtype=complex)
                for block_index in support_keys[index]:
                    candidate_row[0, block_index] = 1.0
            elif use_singleton_block_support and block.shape[1] == 1 and block_count > 1 and half_dim % block_count == 0:
                block_width = half_dim // block_count
                support = np.zeros((1, block_count), dtype=complex)
                flat = first.reshape(-1)
                for block_index in range(block_count):
                    start = block_index * block_width
                    stop = start + block_width
                    if np.max(np.abs(flat[start:stop].real) + np.abs(flat[start:stop].imag)) > tol:
                        support[0, block_index] = 1.0
                candidate_row = support
            elif block.shape[1] > 1 and block_count == 2:
                candidate_row = first.reshape(1, -1)[:, ::2]
            else:
                candidate_row = first.reshape(1, -1)
            if np.max(np.abs(candidate_row.real) + np.abs(candidate_row.imag)) <= tol:
                continue
            rank_before = 0 if kept_rows is None else len(kept_rows)
            trial = candidate_row if kept_rows is None else np.vstack([kept_rows, candidate_row])
            reduced = self.zrowop2_like(trial, tol=tol)
            rank_after = int(
                np.count_nonzero(
                    np.max(np.abs(reduced.real) + np.abs(reduced.imag), axis=1) > tol
                )
            )
            independent = rank_after > rank_before
            if not independent:
                continue
            kept_rows = reduced
            selected_indices.append(index)
            selected_blocks.append(block)
            if len(selected_blocks) == project_count:
                break
        if len(selected_blocks) != project_count:
            raise RuntimeError(
                f"project_ type-3 selection found {len(selected_blocks)} blocks, expected {project_count}"
            )
        return ProjectSelection(
            parent_ops=candidates.parent_ops,
            site_ops=candidates.site_ops,
            pg_irrep=pg_irrep,
            selected_indices=tuple(selected_indices),
            real_blocks=np.array(selected_blocks, dtype=complex),
            reduced_rows=np.array(kept_rows, dtype=complex),
        )

    def project_select_real_blocks_type2_from_records(
        self,
        gid: int,
        site_pg: int,
        pg_irrep: int,
        parent_records: Iterable[tuple[int, int, int, int, int]],
        project_count: int,
        case: Case | None = None,
        tol: float = 1e-10,
    ) -> ProjectSelection:
        """Mirror the type-2 ``project_`` returned-block selection.

        Type-2 irreps return the original normalized real block, but test
        independence after mapping the block through the inverse matrix emitted
        by ``transform_irrep2complex_``.  Each candidate contributes two
        complex rows: one from the first half of the inverse columns and one
        from the conjugate partner half.
        """

        candidates = self.project_real_blocks_from_records(gid, site_pg, pg_irrep, parent_records, case)
        _transform, inverse = self.transform_irrep2complex_matrices(gid)
        kept_rows: np.ndarray | None = None
        selected_indices: list[int] = []
        selected_blocks: list[np.ndarray] = []
        for index, block in enumerate(candidates.rows):
            if np.max(np.abs(block)) <= tol:
                continue
            parent_dim = block.shape[0]
            half_dim = parent_dim // 2
            first = block.T @ inverse[:parent_dim, :half_dim]
            second = block.T @ inverse[:parent_dim, half_dim:parent_dim]
            trial_rows = np.vstack([first.reshape(1, -1), second.reshape(1, -1)])
            trial = trial_rows if kept_rows is None else np.vstack([kept_rows, trial_rows])
            reduced = self.zrowop2_like(trial, tol=tol)
            independent = bool(
                np.max(np.abs(reduced[-1, :].real) + np.abs(reduced[-1, :].imag)) > tol
            )
            if not independent:
                continue
            kept_rows = reduced
            selected_indices.append(index)
            selected_blocks.append(block)
            if len(selected_blocks) == project_count:
                break
        if len(selected_blocks) != project_count:
            raise RuntimeError(
                f"project_ type-2 selection found {len(selected_blocks)} blocks, expected {project_count}"
            )
        return ProjectSelection(
            parent_ops=candidates.parent_ops,
            site_ops=candidates.site_ops,
            pg_irrep=pg_irrep,
            selected_indices=tuple(selected_indices),
            real_blocks=np.array(selected_blocks, dtype=complex),
            reduced_rows=np.array(kept_rows, dtype=complex),
        )
