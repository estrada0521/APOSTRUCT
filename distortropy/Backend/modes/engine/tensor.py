"""Source-only rank-2 tensor projector for the Distortropy mode kernel."""

from __future__ import annotations

import math
from typing import Any

import numpy as np


RANK2_COMPONENTS = ("xx", "xy", "xz", "yy", "yz", "zz")
_RANK2_INDICES = ((0, 0), (0, 1), (0, 2), (1, 1), (1, 2), (2, 2))


def _cartesian_basis_for_sg(sg: int) -> np.ndarray:
    if 143 <= int(sg) <= 194:
        return np.array(
            [[1.0, -0.5, 0.0], [0.0, math.sqrt(3.0) / 2.0, 0.0], [0.0, 0.0, 1.0]],
            dtype=float,
        )
    return np.eye(3, dtype=float)


def _cartesian_point_matrix(decoder: Any, sg: int, point_op: int) -> np.ndarray:
    # tensor_basis_ changes each generated operation from PML to the default
    # International setting before indexing the Cartesian tensor vrep.  The
    # irrep matrix remains attached to the original PML operation record.
    lattice = int(decoder.iso.space["ispace_lattice"][int(sg) - 1])
    code = int(
        decoder.iso.space["ipoint_op_psettings"][(lattice - 1) * 72 + int(point_op) - 1]
    )
    pml_entries: list[float] = []
    for _ in range(9):
        pml_entries.append(float(code % 3 - 1))
        code //= 3
    pml = np.asarray(pml_entries, dtype=float).reshape(3, 3)
    pml_to_cinter = np.asarray(decoder.pml_to_cinter_matrix(int(sg)), dtype=float)
    # The setting matrices act on row-vector coordinates.  Transpose the
    # conjugated operation for the column-vector Cartesian tensor action.
    fractional = (
        np.linalg.inv(pml_to_cinter) @ pml @ pml_to_cinter
    ).T
    basis = _cartesian_basis_for_sg(int(sg))
    return basis @ fractional @ np.linalg.inv(basis)


def strain_normfactor(values: tuple[float, ...] | list[float]) -> float:
    """Return the reciprocal norm for one six-component strain row."""

    xx, xy, xz, yy, yz, zz = (float(value) for value in values)
    norm = math.sqrt(xx * xx + yy * yy + zz * zz + 0.5 * (xy * xy + xz * xz + yz * yz))
    return 0.0 if norm <= 1e-15 else 1.0 / norm


def strain_tensor_components(values: tuple[float, ...] | list[float]) -> dict[str, float]:
    """Convert ISO's Voigt shear coordinates to a symmetric 3x3 tensor."""

    xx, xy, xz, yy, yz, zz = (float(value) for value in values)
    return {
        "xx": xx,
        "xy": xy / 2.0,
        "xz": xz / 2.0,
        "yy": yy,
        "yz": yz / 2.0,
        "zz": zz,
    }


def _rank2_representation(rotation: np.ndarray) -> np.ndarray:
    columns: list[list[float]] = []
    for row, col in _RANK2_INDICES:
        tensor = np.zeros((3, 3), dtype=float)
        tensor[row, col] = 1.0
        tensor[col, row] = 1.0
        transformed = rotation @ tensor @ rotation.T
        columns.append([float(transformed[i, j]) for i, j in _RANK2_INDICES])
    return np.asarray(columns, dtype=float).T


def _canonical_copy(values: np.ndarray, *, tol: float = 1e-8) -> np.ndarray:
    out = np.array(values, dtype=float, copy=True)
    out[np.abs(out) <= tol] = 0.0
    nonzero = np.abs(out[np.abs(out) > tol])
    if not len(nonzero):
        return out
    out /= float(np.min(nonzero))
    out[np.abs(out) <= tol] = 0.0
    first = next(float(value) for value in out.ravel(order="C") if abs(float(value)) > tol)
    if first < 0:
        out *= -1.0
    return out


def _is_independent(candidate: np.ndarray, accepted: list[np.ndarray], *, tol: float = 1e-8) -> bool:
    vector = candidate.ravel(order="C")
    if not accepted:
        return np.linalg.norm(vector) > tol
    before = np.linalg.matrix_rank(np.vstack([item.ravel(order="C") for item in accepted]), tol=tol)
    after = np.linalg.matrix_rank(
        np.vstack([*(item.ravel(order="C") for item in accepted), vector]),
        tol=tol,
    )
    return int(after) > int(before)


def _hexagonal_source_carrier_gauge(anchor: np.ndarray, *, tol: float = 1e-8) -> np.ndarray:
    """Recover the Source vrep carrier after the hex2 Cartesian transform."""

    row_count = int(anchor.shape[0])
    pivot_columns: list[int] = []
    current = np.zeros((row_count, 0), dtype=float)
    current_rank = 0
    for column_index in range(int(anchor.shape[1])):
        candidate = np.column_stack((current, anchor[:, column_index]))
        candidate_rank = int(np.linalg.matrix_rank(candidate, tol=tol))
        if candidate_rank <= current_rank:
            continue
        pivot_columns.append(int(column_index))
        current = candidate
        current_rank = candidate_rank
        if current_rank == row_count:
            break
    if current_rank != row_count:
        return np.eye(row_count, dtype=float)
    pivot = anchor[:, pivot_columns]
    support = np.abs(pivot) > tol
    if np.all(np.sum(support, axis=0) == 1) and np.all(np.sum(support, axis=1) == 1):
        return np.eye(row_count, dtype=float)
    return np.linalg.inv(pivot)


def _rank2_tensor_copies(decoder: Any, sg: int, gid: int) -> tuple[np.ndarray, ...]:
    dimension = int(decoder.little_record_by_gid(int(gid)).full_dim)
    operations = decoder.generate_space_group_records(int(sg))
    representations: list[tuple[np.ndarray, np.ndarray]] = []
    for record in operations:
        rotation = _cartesian_point_matrix(decoder, int(sg), int(record[4]))
        tensor = _rank2_representation(rotation)
        phases = decoder.operation_record_phases(int(gid), record)
        irrep = np.real(decoder.little_phase_matrix_by_gid_record(int(gid), record, phases))
        representations.append((tensor, np.asarray(irrep, dtype=float)))

    accepted: list[np.ndarray] = []
    accepted_rank = 0
    order = float(len(representations))
    # tensor_basis_ scans the packed Cartesian tensor seeds first and the
    # Source-irrep columns second.  xrowop_ only tests whether the candidate
    # rows enlarge the accumulated tensor row space; it does not transform
    # the emitted carrier rows.
    for tensor_component in range(len(RANK2_COMPONENTS)):
        for source_component in range(dimension):
            copy = np.column_stack(
                [
                    sum(
                        irrep[target_component, source_component] * tensor[:, tensor_component]
                        for tensor, irrep in representations
                    )
                    * dimension
                    / order
                    for target_component in range(dimension)
                ]
            ).T
            copy = _canonical_copy(copy)
            if 143 <= int(sg) <= 194:
                copy[:, 1] *= 2.0
            combined = np.vstack((*accepted, copy)) if accepted else copy
            combined_rank = int(np.linalg.matrix_rank(combined, tol=1e-8))
            if combined_rank > accepted_rank:
                accepted.append(copy)
                accepted_rank = combined_rank
    irrep_type = int(decoder.little_record_by_gid(int(gid)).irrep_type)
    if accepted and irrep_type == 1 and 143 <= int(sg) <= 194:
        gauge = _hexagonal_source_carrier_gauge(accepted[0])
        accepted = [gauge @ copy for copy in accepted]
        for copy in accepted:
            copy[np.abs(copy) <= 1e-8] = 0.0
    return tuple(accepted)


def selected_rank2_tensor_rows(
    decoder: Any,
    sg: int,
    gid: int,
    direction_matrix: list[list[float]] | tuple[tuple[float, ...], ...],
) -> tuple[tuple[float, ...], ...]:
    """Apply an invariant OP direction to the faithful rank-2 tensor copies."""

    direction = np.asarray(direction_matrix, dtype=float)
    rows: list[tuple[float, ...]] = []
    for copy in _rank2_tensor_copies(decoder, int(sg), int(gid)):
        selected = direction.T @ copy
        for row in selected:
            values = np.asarray(row, dtype=float)
            values[np.abs(values) <= 1e-8] = 0.0
            if np.any(np.abs(values) > 1e-8):
                rows.append(tuple(float(value) for value in values))
    return tuple(rows)


def selected_rank2_tensor_row_groups(
    decoder: Any,
    sg: int,
    gid: int,
    direction_matrix: list[list[float]] | tuple[tuple[float, ...], ...],
    *,
    tolerance: float = 1e-8,
) -> tuple[tuple[tuple[float, ...], ...], ...]:
    """Return Source-ordered tensor-copy groups that add a new strain subspace."""

    direction = np.asarray(direction_matrix, dtype=float)
    accepted: list[np.ndarray] = []
    groups: list[tuple[tuple[float, ...], ...]] = []
    for copy in _rank2_tensor_copies(decoder, int(sg), int(gid)):
        group: list[tuple[float, ...]] = []
        for row in direction.T @ copy:
            values = np.asarray(row, dtype=float)
            values[np.abs(values) <= tolerance] = 0.0
            if not _is_independent(values, accepted, tol=tolerance):
                continue
            accepted.append(values)
            group.append(tuple(float(value) for value in values))
        if group:
            groups.append(tuple(group))
    return tuple(groups)


def totally_symmetric_rank2_rows(decoder: Any, sg: int, gid: int) -> tuple[tuple[float, ...], ...]:
    """Return selected GM totally-symmetric tensor copies in selection order."""

    dimension = int(decoder.little_record_by_gid(int(gid)).full_dim)
    return selected_rank2_tensor_rows(
        decoder,
        int(sg),
        int(gid),
        [[1.0 if row == col else 0.0 for col in range(dimension)] for row in range(dimension)],
    )
