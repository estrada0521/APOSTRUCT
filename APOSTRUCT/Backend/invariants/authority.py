"""APOSTRUCT Source authority consumed by the invariant kernel."""

from __future__ import annotations

from fractions import Fraction
from functools import lru_cache
import math
from typing import Sequence

import numpy as np

from APOSTRUCT.Backend.isotropy.engine.dynamic_isotropy_file import (
    sort_dynamic_rows_for_file,
)
from APOSTRUCT.Backend.isotropy.engine.get_isotropy import (
    generate_dynamic_isotropy_rows,
)
from APOSTRUCT.Backend.invariants.domains import (
    direction_domain_from_subgroup,
)
from APOSTRUCT.Backend.isotropy.engine.source_data import SourceData
from APOSTRUCT.Backend.modes.engine.decoder import ModeDataDecoder
from APOSTRUCT.Backend.source.tables import SOURCE, source_tables


class InvariantSource:
    """Small Source facade required by DISPLAY INVARIANT calculations."""

    def __init__(self) -> None:
        tables = source_tables()
        self.iso = tables.iso
        self.decoder = ModeDataDecoder(SOURCE, tables=tables)
        self.source_data = SourceData(SOURCE, tables=tables.iso)
        self.space = tables.iso.space
        self.little = tables.iso.little
        self.irreps = tables.iso.irreps
        self.isotropy = tables.iso.isotropy
        self.const = tables.iso.const

    def little_record_by_gid(self, gid: int):
        return self.decoder.little_record_by_gid(int(gid))

    def generate_space_group_records(self, sg: int):
        return self.source_data.generate_space_group_records(int(sg))

    def get_irreps_matrix_for_record(
        self,
        gid: int,
        record: tuple[int, int, int, int, int],
        kparam: Sequence[int] | None = None,
    ):
        return self.source_data.get_irreps_matrix_for_record(
            int(gid),
            tuple(int(value) for value in record),
            tuple(int(value) for value in (kparam or ())),
        )

    def list_irreps(self, sg: int) -> list[tuple[int, str, str, str]]:
        parent = int(sg)
        lattice = int(self.space["ispace_lattice"][parent - 1])
        out: list[tuple[int, str, str, str]] = []
        for gid, row_sg in enumerate(self.little["little_irr_space_group"], start=1):
            if int(row_sg) != parent:
                continue
            kslot = int(self.little["little_irr_k"][gid - 1])
            k_index = (lattice - 1) * 27 + kslot - 1
            old_id = int(self.little["little_irr_old"][gid - 1])
            k_dim = int(self.little["little_k_dim"][k_index])
            out.append(
                (
                    gid,
                    str(self.little["little_irr_full_label"][gid - 1]).strip(),
                    str(self.little["little_k_label"][k_index]).strip(),
                    "parametric" if k_dim > 0 or old_id <= 0 else "fixed",
                )
            )
        return out

    def k_parameter_dimension_by_gid(self, gid: int) -> int:
        gid = int(gid)
        sg = int(self.little["little_irr_space_group"][gid - 1])
        kslot = int(self.little["little_irr_k"][gid - 1])
        lattice = int(self.space["ispace_lattice"][sg - 1])
        return int(self.little["little_k_dim"][(lattice - 1) * 27 + kslot - 1])

    @staticmethod
    def _kparam_record(values: Sequence[float | Fraction | int]) -> tuple[int, int, int, int]:
        raw = tuple(values)
        if len(raw) == 4 and int(raw[3]):
            return tuple(int(value) for value in raw)  # type: ignore[return-value]
        fractions = tuple(Fraction(value) for value in raw[:3])
        denominator = 1
        for value in fractions:
            denominator = math.lcm(denominator, value.denominator)
        padded = fractions + (Fraction(0),) * (3 - len(fractions))
        return (
            int(padded[0] * denominator),
            int(padded[1] * denominator),
            int(padded[2] * denominator),
            denominator,
        )

    def source_kparam_for_gid(
        self,
        gid: int,
        kparam: Sequence[float | Fraction | int] | None = None,
    ) -> tuple[int, int, int, int]:
        if kparam is not None:
            return self._kparam_record(kparam)
        vectors = self.decoder.little_k_vectors_by_gid(int(gid)).vectors
        vector = vectors[0] if vectors else (Fraction(0), Fraction(0), Fraction(0))
        return self._kparam_record(vector)

    def _direction_rows(self, gid: int) -> range:
        old_id = int(self.little_record_by_gid(int(gid)).old_id)
        if old_id <= 0:
            return range(0)
        pointers = self.isotropy["isotropy_irrep_pointer"]
        return range(int(pointers[old_id - 1]), int(pointers[old_id]))

    def isotropy_direction_row(self, gid: int, direction: str) -> int | None:
        wanted = str(direction).strip()
        for row_id in self._direction_rows(int(gid)):
            if str(self.isotropy["isotropy_orderparam_label"][row_id - 1]).strip() == wanted:
                return row_id
        return None

    def list_direction_labels(self, gid: int) -> list[str]:
        labels = [
            str(self.isotropy["isotropy_orderparam_label"][row_id - 1]).strip()
            for row_id in self._direction_rows(int(gid))
        ]
        if labels:
            return labels
        dimension = int(self.little_record_by_gid(int(gid)).full_dim)
        return [f"{dimension}D1" if dimension > 1 else "P1"]

    def _static_direction_matrix(self, row_id: int) -> tuple[tuple[float, ...], ...]:
        index = int(row_id) - 1
        dimension = int(self.isotropy["isotropy_orderparam_dim"][index])
        free = int(self.isotropy["isotropy_orderparam_freeparam"][index])
        pointer = int(self.isotropy["isotropy_orderparam_pointer"][index]) - 1
        codes = self.isotropy["isotropy_orderparam"][pointer:pointer + dimension * free]
        values = [float(self.const[int(code)]) for code in codes]
        return tuple(
            tuple(values[column * dimension + row] for column in range(free))
            for row in range(dimension)
        )

    @lru_cache(maxsize=4096)
    def _dynamic_rows(
        self,
        gid: int,
        kparam: tuple[int, int, int, int],
    ) -> tuple[object, ...]:
        return tuple(
            sort_dynamic_rows_for_file(
                generate_dynamic_isotropy_rows(
                    self.source_data,
                    gid=int(gid),
                    kparam=kparam,
                )
            )
        )

    def _dynamic_row_for_gid_kparam_direction(
        self,
        gid: int,
        kparam: Sequence[float | Fraction | int],
        direction: str,
    ):
        wanted = str(direction).strip()
        return next(
            (
                row
                for row in self._dynamic_rows(int(gid), self._kparam_record(kparam))
                if str(row.direction).strip() == wanted
            ),
            None,
        )

    def direction_matrix_by_gid_label(
        self,
        gid: int,
        direction: str,
        *,
        kparam: Sequence[float | Fraction | int] | None = None,
    ) -> tuple[int | None, tuple[tuple[float, ...], ...]]:
        little = self.little_record_by_gid(int(gid))
        if int(little.old_id) > 0:
            row_id = self.isotropy_direction_row(int(gid), direction)
            if row_id is None:
                raise KeyError(f"unknown Source OPD for gid={gid}: {direction}")
            return row_id, self._static_direction_matrix(row_id)
        if kparam is None:
            raise ValueError("parametric Source OPD requires k parameters")
        row = self._dynamic_row_for_gid_kparam_direction(int(gid), kparam, direction)
        if row is None or not row.matrix:
            raise KeyError(f"unknown dynamic Source OPD for gid={gid}: {direction}")
        free = len(row.matrix)
        dimension = len(row.matrix[0])
        return None, tuple(
            tuple(float(row.matrix[column][component]) for column in range(free))
            for component in range(dimension)
        )

    @staticmethod
    def _same_direction_family(
        left: Sequence[Sequence[float]],
        right: Sequence[Sequence[float]],
        *,
        tol: float = 1e-7,
    ) -> bool:
        left_matrix = np.asarray(left, dtype=float)
        right_matrix = np.asarray(right, dtype=float)
        if left_matrix.ndim != 2 or right_matrix.ndim != 2:
            return False
        left_matrix = left_matrix[np.linalg.norm(left_matrix, axis=1) > 1e-10]
        right_matrix = right_matrix[np.linalg.norm(right_matrix, axis=1) > 1e-10]
        if (
            left_matrix.shape[1] != right_matrix.shape[1]
            or len(left_matrix) != len(right_matrix)
        ):
            return False
        for source, target in ((left_matrix, right_matrix), (right_matrix, left_matrix)):
            for vector in source:
                coefficients, *_ = np.linalg.lstsq(target.T, vector, rcond=None)
                projected = target.T @ coefficients
                residual = np.linalg.norm(projected - vector) / max(np.linalg.norm(vector), 1e-12)
                if residual > tol:
                    return False
        return True

    def direction_label_for_matrix(
        self,
        gid: int,
        matrix: Sequence[Sequence[float]],
        *,
        kparam: Sequence[float | Fraction | int] | None = None,
    ) -> str | None:
        """Return the Source OPD family label, independent of its selected domain."""

        little = self.little_record_by_gid(int(gid))
        candidates: list[tuple[str, tuple[tuple[float, ...], ...]]] = []
        if int(little.old_id) > 0:
            candidates.extend(
                (
                    str(self.isotropy["isotropy_orderparam_label"][row_id - 1]).strip(),
                    self._static_direction_matrix(row_id),
                )
                for row_id in self._direction_rows(int(gid))
            )
        elif kparam is not None:
            for row in self._dynamic_rows(int(gid), self._kparam_record(kparam)):
                if not row.matrix:
                    continue
                free = len(row.matrix)
                dimension = len(row.matrix[0])
                candidates.append(
                    (
                        str(row.direction).strip(),
                        tuple(
                            tuple(float(row.matrix[column][component]) for column in range(free))
                            for component in range(dimension)
                        ),
                    )
                )
        return next(
            (
                label
                for label, candidate in candidates
                if self._same_direction_family(matrix, candidate)
            ),
            None,
        )

    def direction_label_domain_for_matrix(
        self,
        sg: int,
        gid: int,
        matrix: Sequence[Sequence[float]],
        *,
        kparam: Sequence[float | Fraction | int] | None = None,
    ) -> tuple[str, int] | None:
        """Return the Source OPD label and exact domain for a carrier subspace."""

        little = self.little_record_by_gid(int(gid))
        label = self.direction_label_for_matrix(int(gid), matrix, kparam=kparam)
        if label is None:
            return None
        candidates = []
        if int(little.old_id) > 0:
            for row_id in self._direction_rows(int(gid)):
                index = row_id - 1
                if str(self.isotropy["isotropy_orderparam_label"][index]).strip() != label:
                    continue
                candidates.append(
                    (
                        str(self.isotropy["isotropy_orderparam_label"][index]).strip(),
                        self._static_direction_matrix(row_id),
                        int(self.isotropy["isotropy_subgroup"][index]),
                        tuple(
                            int(value)
                            for value in self.isotropy["isotropy_basis"][index * 9:(index + 1) * 9]
                        ),
                        tuple(
                            int(value)
                            for value in self.isotropy["isotropy_origin"][index * 4:(index + 1) * 4]
                        ),
                    )
                )
            source_kparam: tuple[int, ...] = ()
            tolerance = 1e-7
        else:
            if kparam is None:
                raise ValueError("parametric Source OPD requires k parameters")
            source_kparam = self._kparam_record(kparam)
            for row in self._dynamic_rows(int(gid), source_kparam):
                if not row.matrix or str(row.direction).strip() != label:
                    continue
                free = len(row.matrix)
                dimension = len(row.matrix[0])
                candidates.append(
                    (
                        str(row.direction).strip(),
                        tuple(
                            tuple(float(row.matrix[column][component]) for column in range(free))
                            for component in range(dimension)
                        ),
                        int(row.subgroup_number),
                        tuple(int(value) for value in row.basis_values),
                        tuple(int(value) for value in row.origin_values),
                    )
                )
            tolerance = 1e-6

        matches = []
        for _label, candidate, child_sg, basis, origin in candidates:
            domain = direction_domain_from_subgroup(
                self,
                sg=int(sg),
                gid=int(gid),
                kparam=source_kparam,
                target=matrix,
                candidate=candidate,
                child_sg=child_sg,
                basis=basis,
                origin=origin,
                tolerance=tolerance,
            )
            if domain is not None:
                matches.append(domain)
        if not matches:
            return None
        return label, min(matches)

    def pml_to_cinter_matrix(self, sg: int):
        return self.decoder.pml_to_cinter_matrix(int(sg))


@lru_cache(maxsize=1)
def invariant_source() -> InvariantSource:
    return InvariantSource()
