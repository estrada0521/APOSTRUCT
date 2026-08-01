from __future__ import annotations

from fractions import Fraction
from typing import Iterable

from distortropy.Backend.modes.engine.records import WyckoffRow
from distortropy.Backend.source.iso_data import ISOData

class WyckoffGeometryMixin:
    def wyckoff_representative(
        self,
        row: WyckoffRow,
        site_params: Iterable[object] | None = None,
    ) -> tuple[Fraction, Fraction, Fraction]:
        """Return the input-dependent representative for a Wyckoff row.

        `iwyckoff_fract` stores four rational 3-vectors:
        `base + x*v_x + y*v_y + z*v_z`.  mode-kernel atom input uses those x/y/z
        slots directly, including zero placeholders for unused earlier slots.
        Constant rows therefore work with `site_params=None`, while parametric
        rows must pass the parsed atom-line parameters to recover the same site
        stabilizer that `project_` receives.
        """

        raw = tuple(int(x) for x in self.iso.wyckoff["iwyckoff_fract"][(row.row_id - 1) * 16:row.row_id * 16])
        den = int(raw[3])
        if den == 0:
            raise ValueError(f"zero Wyckoff representative denominator for row {row.row_id}")
        coords = [Fraction(raw[i], den) for i in range(3)]
        params = [Fraction(str(value)) for value in (site_params or ())]
        for slot in range(3):
            start = 4 + slot * 4
            slot_den = int(raw[start + 3])
            if slot_den == 0:
                continue
            value = params[slot] if slot < len(params) else Fraction(0, 1)
            if value == 0:
                continue
            for axis in range(3):
                coords[axis] += value * Fraction(raw[start + axis], slot_den)
        return tuple(coord % 1 for coord in coords)  # type: ignore[return-value]

    @staticmethod
    def _vmlt_fraction(
        matrix: Iterable[int],
        vector: tuple[Fraction, Fraction, Fraction],
    ) -> tuple[Fraction, Fraction, Fraction]:
        return ISOData.vmlt_fraction(matrix, vector)

    def wyc_pg_cosets_records(
        self,
        sg: int,
        site_records: Iterable[tuple[int, int, int, int, int]],
    ) -> tuple[tuple[int, int, int, int, int], ...]:
        records = tuple(tuple(int(value) for value in record) for record in site_records)
        return self.iso.wyc_pg_cosets_records(int(sg), records)
