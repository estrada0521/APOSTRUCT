"""Real-to-complex little-irrep transforms."""

from __future__ import annotations

import numpy as np


class TransformIrrep2ComplexMixin:
    def little_transform_block_count(self, gid: int) -> int:
        """Return the block count read by ``transform_irrep2complex_``.

        The block-count contract does not use the full-irrep
        ``little_irr_lif`` field. It indexes the Source small-k table by
        ``(space_group - 1) * 27 + kslot``; SG226/W maps to 3 blocks and
        SG219/L maps to 4.
        """

        sg = int(self.iso.little["little_irr_space_group"][gid - 1])
        kslot = int(self.iso.little["little_irr_k"][gid - 1])
        index = (sg - 1) * 27 + kslot
        return int(self.iso.little["little_k_star_count"][index - 1])

    def little_transform_uses_permutation(self, gid: int) -> bool:
        """Return the small-k flag used by ``transform_irrep2complex_``.

        Source stores this flag in the T/F-encoded
        ``little_k_star_minusk`` table adjacent to ``little_k_star_count``.
        It controls the initial real matrix for type-2 irreps.
        """

        sg = int(self.iso.little["little_irr_space_group"][gid - 1])
        kslot = int(self.iso.little["little_irr_k"][gid - 1])
        index = (sg - 1) * 27 + kslot - 1
        flags = self.iso.little.get("little_k_star_minusk", [])
        return bool(flags and flags[index])

    def little_real2_triples(self, gid: int) -> tuple[tuple[int, int, int], ...]:
        """Return the ``little_irr_real2`` triples used by type-2 conversion.

        The Source pointer is a Fortran 1-based triple index, so the flat
        Python start is ``(pointer - 1) * 3``.
        """

        little = self.little_record_by_gid(gid)
        if not little.real2_pointer:
            return ()
        half = little.full_dim // 2
        start = (int(little.real2_pointer) - 1) * 3
        values = self.iso.little["little_irr_real2"]
        triples = []
        for offset in range(half):
            base = start + offset * 3
            triples.append(tuple(int(x) for x in values[base:base + 3]))
        return tuple(triples)  # type: ignore[return-value]

    @staticmethod
    def _zmatmlt_fortran(
        left: np.ndarray,
        right: np.ndarray,
        *,
        active_dim: int,
        leading_dim: int,
    ) -> np.ndarray:
        """Multiply the active square block and clear the padded region."""

        out = np.array(right, dtype=complex, copy=True)
        out[:active_dim, :active_dim] = left[:active_dim, :active_dim] @ right[:active_dim, :active_dim]
        if active_dim < leading_dim:
            out[active_dim:, :] = 0
            out[:, active_dim:] = 0
        out[np.abs(out) < 1e-12] = 0
        return out

    def transform_irrep2complex_matrices(
        self,
        gid: int,
        *,
        leading_dim: int = 48,
    ) -> tuple[np.ndarray, np.ndarray]:
        """Return the two matrices emitted by ``transform_irrep2complex_``.

        The first matrix is constructed from the little-irrep dimensions and
        type/permutation metadata.  Type-2 rows with a ``real2`` pointer also
        apply the corresponding Source ``little_irr_real2`` triples.  The
        second matrix is its conjugate transpose, as consumed by projection.
        """

        little = self.little_record_by_gid(gid)
        active_dim = int(little.full_dim)
        dim = int(leading_dim)
        matrix = np.zeros((dim, dim), dtype=complex)
        matrix[:active_dim, :active_dim] = np.eye(active_dim, dtype=complex)
        if little.irrep_type == 1:
            return matrix, matrix.conj().T

        block_count = self.little_transform_block_count(gid)
        if block_count <= 0:
            raise ValueError(f"invalid transform block count for gid={gid}: {block_count}")
        block_dim = active_dim // block_count
        half = block_dim // 2
        root2 = np.sqrt(2.0)

        if (
            not little.real2_pointer
            and self.little_transform_uses_permutation(gid)
            and block_dim == 4
            and block_count % 2 == 1
        ):
            matrix = np.zeros((dim, dim), dtype=complex)
            pair_count = active_dim // 2
            for pair in range(pair_count):
                top = pair
                bottom = pair + pair_count
                col0 = (pair // 2) * 4 + (pair % 2)
                col1 = col0 + 2
                if pair % 2 == 0:
                    matrix[top, col0] = 1.0 / root2
                    matrix[bottom, col0] = 1.0 / root2
                    matrix[top, col1] = -1.0j / root2
                    matrix[bottom, col1] = 1.0j / root2
                else:
                    matrix[top, col0] = -1.0j / root2
                    matrix[bottom, col0] = 1.0j / root2
                    matrix[top, col1] = 1.0 / root2
                    matrix[bottom, col1] = 1.0 / root2
            inverse = matrix.conj().T
            return matrix, inverse

        if little.real2_pointer and self.little_transform_uses_permutation(gid):
            matrix = np.zeros((dim, dim), dtype=complex)
            perm = np.array(
                [
                    [1, 0, 0, 0],
                    [0, 0, 0, 1],
                    [0, 0, 1, 0],
                    [0, 1, 0, 0],
                ],
                dtype=complex,
            )
            if block_dim > 3:
                repeat = block_dim // 4
                for block in range(block_count):
                    base = block * block_dim
                    for outer in range(4):
                        row_start = base + (outer * block_dim) // 4
                        for inner in range(4):
                            col_start = base + (inner * block_dim) // 4
                            for offset in range(repeat):
                                matrix[row_start + offset, col_start + offset] = perm[outer, inner]

        first = np.zeros((dim, dim), dtype=complex)
        base = 0
        # Build the first complex transform when no real type-2 block is available.
        for _block in range(block_count):
            for outer in range(2):
                row_start = base + (outer * block_dim) // 2
                for inner in range(2):
                    col_start = base + (inner * block_dim) // 2
                    value = (1.0, -1.0j, 1.0, 1.0j)[outer * 2 + inner] / root2
                    for offset in range(half):
                        first[row_start + offset, col_start + offset] = value
            base += block_dim
        matrix = self._zmatmlt_fortran(first, matrix, active_dim=active_dim, leading_dim=dim)

        second = np.zeros((dim, dim), dtype=complex)
        wrapped = 1
        base = 0
        for _step in range(block_count * 2):
            row_start = base // 2
            col_start = ((wrapped - 1) * block_dim) // 2
            for offset in range(half):
                second[row_start + offset, col_start + offset] = 1.0
            wrapped += 2
            if block_count * 2 < wrapped:
                wrapped = 2
            base += block_dim
        matrix = self._zmatmlt_fortran(second, matrix, active_dim=active_dim, leading_dim=dim)

        if little.irrep_type == 2 and little.real2_pointer:
            third = np.zeros((dim, dim), dtype=complex)
            final_half = active_dim // 2
            for offset in range(final_half):
                third[offset, offset] = 1.0
            for column_offset, (row_code, real_code, imag_code) in enumerate(self.little_real2_triples(gid)):
                row = final_half + int(row_code) - 1
                col = final_half + column_offset
                third[row, col] = complex(float(self.iso.const[int(real_code)]), float(self.iso.const[int(imag_code)]))
            matrix = self._zmatmlt_fortran(third, matrix, active_dim=active_dim, leading_dim=dim)

        inverse = matrix.conj().T
        return matrix, inverse
