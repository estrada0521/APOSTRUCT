from __future__ import annotations

from fractions import Fraction
from typing import Iterable

import numpy as np

from verification.backend.source.representation import _translation_mod_one
from verification.backend.modes.engine.records import (
    LittleKVectorRecord,
)


class IrrepMixin:
    def generated_space_point_op_positions(self, sg: int) -> dict[int, int]:
        """Return the point-op -> internal position map used by `get_irrep4_ssgn_`.

        `get_irrep4_ssgn_` calls `generate_space_group_` and builds a `mapop`
        array from the generated operation records.  The sparse
        `little_irr_full_matrices_*` tables are indexed by that internal
        position, not by the user-facing operation order used by
        `irrep_matrices_with_ops`.
        """

        pointer = int(self.iso.space["ispace_elements_pointer"][sg - 1])
        count = len(self.iso.generate_space_group_records(sg))
        raw = self.iso.space["ispace_elements"][(pointer - 1) * 5:(pointer - 1 + count) * 5]
        out: dict[int, int] = {}
        for index in range(count):
            record = raw[5 * index:5 * index + 5]
            if len(record) == 5:
                # The Fortran map is overwritten if a point op appears more
                # than once; mirror that behavior.
                out[int(record[4])] = index + 1
        return out

    def generated_space_operation_position(self, sg: int, record: tuple[int, int, int, int, int]) -> int:
        """Return the generated operation position matching a full 5-int record.

        Parent space groups can have several operations with the same point-op
        label but different centering/screw/glide translations.  The parent
        `get_irrep4_ssgn_` matrix lookup is keyed by the generated operation
        position, so point-op-only lookup is insufficient for centered groups.
        """

        tables = self.__dict__.setdefault("_generated_position_tables", {})
        table = tables.get(int(sg))
        if table is None:
            table = tables[int(sg)] = {}
            pointer = int(self.iso.space["ispace_elements_pointer"][sg - 1])
            for index in range(len(self.iso.generate_space_group_records(sg))):
                x, y, z, den, op = (
                    int(value)
                    for value in self.iso.space["ispace_elements"][
                        5 * (pointer - 1 + index):5 * (pointer - 1 + index + 1)
                    ]
                )
                table.setdefault((op, _translation_mod_one(x, y, z, den)), index + 1)
        key = (int(record[4]), _translation_mod_one(*(int(value) for value in record[:4])))
        position = table.get(key)
        if position is None:
            raise KeyError(f"operation record {record} not found in generated SG{sg} records")
        return position

    def space_operation_records(self, sg: int) -> tuple[tuple[int, int, int, int, int], ...]:
        """Return generated space-operation records as native 5-int tuples."""

        cache = self.__dict__.setdefault("_space_operation_records_cache", {})
        if sg not in cache:
            cache[sg] = self._space_operation_records(int(sg))
        return cache[sg]

    def _space_operation_records(self, sg: int) -> tuple[tuple[int, int, int, int, int], ...]:
        pointer = int(self.iso.space["ispace_elements_pointer"][sg - 1])
        count = len(self.iso.generate_space_group_records(sg))
        out: list[tuple[int, int, int, int, int]] = []
        for index in range(count):
            raw = self.iso.space["ispace_elements"][
                5 * (pointer - 1 + index):5 * (pointer - 1 + index + 1)
            ]
            if len(raw) == 5:
                out.append(tuple(int(x) for x in raw))  # type: ignore[arg-type]
        return tuple(out)

    def generated_space_tau(self, sg: int, point_op: int) -> tuple[Fraction, Fraction, Fraction]:
        """Return the generated operation translation for a point operation.

        This mirrors the helper used by the parent `get_irrep4_` phase
        verifier.  If a point operation appears more than once, the generated
        table's first matching operation is the one used for the sparse matrix
        representative; extra translation is carried as a Bloch phase.
        """

        for x, y, z, den, op in self.space_operation_records(sg):
            if int(op) == int(point_op):
                return (Fraction(x, den), Fraction(y, den), Fraction(z, den))
        raise KeyError(f"point op {point_op} not found in generated SG{sg} records")

    def site_pg_element_settings(self, site_pg: int) -> list[tuple[int, ...]]:
        """Return candidate site point-group operation settings.

        Each setting is a row from ``iwyckoff_pg_elements`` with trailing zeros
        removed.  In ``project_`` the selected row is converted into 5-int
        operation records ``(0,0,0,1, point_op)`` before calling
        ``get_irreps_`` for the site point-group irrep.
        """

        return list(self.iso.site_pg_element_settings(int(site_pg)))

    def little_k_vectors_by_gid(self, gid: int) -> LittleKVectorRecord:
        """Return the k-vector records copied by ``get_irrep4_``.

        For fixed special-k rows, the relevant static table is ``little_k_star_special``, not the later
        ``little_kvec`` section.  ``read_little`` loads it as 2100 records of
        four integers ``(h, k, l, denominator)``.  ``get_irrep4_`` addresses it
        with the 1-based ``little_k_star_ml_pointer`` and copies ``nmod``
        consecutive records into the first column of the runtime ``kvec.1``
        common block.

        Parametric rows require runtime ``Case.k_params`` and are handled by
        ``little_k_star_vectors_for_case`` using the 16-int ``little_k_star``
        records.
        """

        gid = int(gid)
        sg = int(self.iso.little["little_irr_space_group"][gid - 1])
        kslot = int(self.iso.little["little_irr_k"][gid - 1])
        slot = (sg - 1) * 27 + kslot - 1
        nmod = int(self.iso.little["little_k_star_count"][slot])
        pointer = int(self.iso.little["little_k_star_ml_pointer"][slot])
        if nmod < 0:
            raise ValueError(f"negative k-star count for gid={gid}: {nmod}")
        if nmod and pointer <= 0:
            raise ValueError(f"missing little_k_star_ml_pointer for gid={gid} SG{sg} kslot={kslot}")

        raw = self.iso.little["little_k_star_special"]
        records: list[tuple[int, int, int, int]] = []
        vectors: list[tuple[Fraction, Fraction, Fraction]] = []
        for offset in range(nmod):
            start = 4 * (pointer - 1 + offset)
            item = tuple(int(x) for x in raw[start:start + 4])
            if len(item) != 4:
                raise ValueError(f"k-vector pointer out of range for gid={gid}: pointer={pointer} nmod={nmod}")
            den = item[3]
            if den == 0:
                raise ValueError(f"zero k-vector denominator for gid={gid}: pointer={pointer} offset={offset}")
            records.append(item)  # type: ignore[arg-type]
            vectors.append((Fraction(item[0], den), Fraction(item[1], den), Fraction(item[2], den)))

        return LittleKVectorRecord(
            gid=gid,
            sg=sg,
            kslot=kslot,
            nmod=nmod,
            pointer=pointer,
            records=tuple(records),
            vectors=tuple(vectors),
        )
