"""Assembled Source authority consumed by the invariant kernel."""

from __future__ import annotations

from fractions import Fraction
from functools import lru_cache
import math
from typing import Sequence

from ISODISTORT.Assembled.Backend.isotropy.engine.dynamic_isotropy_file import (
    sort_dynamic_rows_for_file,
)
from ISODISTORT.Assembled.Backend.isotropy.engine.get_isotropy import (
    generate_dynamic_isotropy_rows,
)
from ISODISTORT.Assembled.Backend.isotropy.engine.source_data import SourceData
from ISODISTORT.Assembled.Backend.modes.engine.decoder import ModeDataDecoder
from ISODISTORT.Assembled.Backend.source.tables import SOURCE, source_tables


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

    def pml_to_cinter_matrix(self, sg: int):
        return self.decoder.pml_to_cinter_matrix(int(sg))


@lru_cache(maxsize=1)
def invariant_source() -> InvariantSource:
    return InvariantSource()
