"""``zrowop2_``, ``xrowop2_``, and subgroup constraint row logic."""

from __future__ import annotations

from typing import Iterable

import numpy as np


class RowOpsMixin:
    @staticmethod
    def zrowop2_like(rows: np.ndarray, tol: float = 1e-10) -> np.ndarray:
        """Return the complex reduced row form used by mode kernel ``zrowop2_``.

        This is ordinary complex Gauss-Jordan elimination: choose the first
        nonzero pivot in the current column, swap it into the current row,
        divide by the pivot, then remove that column from all other rows. This
        helper returns a copy rather than mutating the caller's array.
        """

        out = np.array(rows, dtype=complex, copy=True)
        if out.ndim != 2:
            raise ValueError("zrowop2_like expects a 2D matrix")
        pivot_row = 0
        nrow, ncol = out.shape
        for col in range(ncol):
            if pivot_row >= nrow:
                break
            pivot = None
            for row in range(pivot_row, nrow):
                if abs(out[row, col].real) + abs(out[row, col].imag) > tol:
                    pivot = row
                    break
            if pivot is None:
                continue
            if pivot != pivot_row:
                out[[pivot_row, pivot], :] = out[[pivot, pivot_row], :]
            out[pivot_row, :] = out[pivot_row, :] / out[pivot_row, col]
            for row in range(nrow):
                if row == pivot_row:
                    continue
                factor = out[row, col]
                if abs(factor.real) + abs(factor.imag) > tol:
                    out[row, :] = out[row, :] - factor * out[pivot_row, :]
            pivot_row += 1
        out[np.abs(out) < tol] = 0
        return out

    @staticmethod
    def xrowop2_like(rows: np.ndarray, tol: float = 1e-10) -> np.ndarray:
        """Return the real reduced row form used by mode kernel ``xrowop2_``."""

        out = np.array(rows, dtype=float, copy=True)
        if out.ndim != 2:
            raise ValueError("xrowop2_like expects a 2D matrix")
        pivot_row = 0
        nrow, ncol = out.shape
        for col in range(ncol):
            if pivot_row >= nrow:
                break
            pivot = None
            for row in range(pivot_row, nrow):
                if abs(out[row, col]) > tol:
                    pivot = row
                    break
            if pivot is None:
                continue
            if pivot != pivot_row:
                out[[pivot_row, pivot], :] = out[[pivot, pivot_row], :]
            out[pivot_row, :] = out[pivot_row, :] / out[pivot_row, col]
            for row in range(nrow):
                if row == pivot_row:
                    continue
                factor = out[row, col]
                if abs(factor) > tol:
                    out[row, :] = out[row, :] - factor * out[pivot_row, :]
            pivot_row += 1
        out[np.abs(out) < tol] = 0
        return out

    @staticmethod
    def subgroup_direction_from_reduced_matrix(
        reduced: np.ndarray,
        *,
        direction_stride: int = 48,
        tol: float = 1e-10,
    ) -> tuple[int, np.ndarray]:
        """Build `subgroup_to_orderparam_` direction vectors from an RREF matrix.

        `subgroup_to_orderparam_` scans the RREF rows produced by `xrowop2_`,
        records the first nonzero column of each nonzero row as a pivot, then
        emits one direction vector for each non-pivot column.  The output buffer
        uses a fixed 48-double stride per direction.
        """

        matrix = np.asarray(reduced, dtype=float)
        if matrix.ndim != 2:
            raise ValueError("subgroup_direction_from_reduced_matrix expects a 2D matrix")
        dim = matrix.shape[1]
        nonzero_rows = [row for row in matrix if np.max(np.abs(row)) > tol]
        rank = len(nonzero_rows)
        count = dim - rank
        directions = np.zeros((max(count, 0), direction_stride), dtype=float)
        pivot_cols: list[int] = []
        for row in nonzero_rows:
            nonzero = np.where(np.abs(row) > tol)[0]
            if len(nonzero) == 0:
                continue
            pivot_cols.append(int(nonzero[0]))

        free_index = 0
        # The routine walks matrix columns from right to left while writing
        # output coordinates from left to right.
        for write_col, matrix_col in enumerate(range(dim - 1, -1, -1)):
            if matrix_col in pivot_cols:
                continue
            if free_index >= count:
                break
            directions[free_index, write_col] = 1.0
            for row_index, pivot_col in enumerate(pivot_cols):
                pivot_value = nonzero_rows[row_index][pivot_col]
                if abs(pivot_value) <= tol:
                    continue
                directions[free_index, dim - pivot_col - 1] = (
                    -nonzero_rows[row_index][matrix_col] / pivot_value
                )
            free_index += 1
        directions[np.abs(directions) < tol] = 0
        return count, directions

    @staticmethod
    def subgroup_constraint_matrix_from_irrep_matrices(
        matrices: Iterable[np.ndarray],
        *,
        tol: float = 1e-12,
    ) -> np.ndarray:
        """Build the pre-RREF matrix used by ``subgroup_to_orderparam_``.

        After each internal ``get_irreps_`` call, the mode kernel appends
        ``dim`` real rows to the work matrix later passed to ``xrowop2_``.  The
        stored block is ``D - I`` with the columns written in reverse order;
        this matches the routine's later right-to-left free-column convention.
        """

        blocks: list[np.ndarray] = []
        dim: int | None = None
        for matrix in matrices:
            arr = np.asarray(matrix, dtype=float)
            if arr.ndim != 2 or arr.shape[0] != arr.shape[1]:
                raise ValueError(f"expected square irrep matrix, got shape {arr.shape}")
            if dim is None:
                dim = int(arr.shape[0])
            elif arr.shape[0] != dim:
                raise ValueError(f"mixed matrix dimensions: {dim} and {arr.shape[0]}")
            blocks.append((arr - np.eye(arr.shape[0]))[:, ::-1])
        if dim is None:
            return np.zeros((0, 0), dtype=float)
        out = np.vstack(blocks) if blocks else np.zeros((0, dim), dtype=float)
        out[np.abs(out) < tol] = 0.0
        return out
