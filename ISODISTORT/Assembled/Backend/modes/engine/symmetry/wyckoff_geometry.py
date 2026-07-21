from __future__ import annotations

from fractions import Fraction
import math
from typing import Iterable

from ISODISTORT.Assembled.Backend.modes.engine.records import (
    LittleIrrepRecord,
    WyckoffRow,
)
from ISODISTORT.Assembled.Backend.source.iso_data import ISOData


class WyckoffGeometryMixin:
    def project_parent_point_ops(self, sg: int, label: str, site_pg: int) -> tuple[int, ...]:
        """Select the parent-operation row that `project_` uses for a site PG.

        `little_ops` gives the little-group point operations.  `project_`
        then chooses the Wyckoff site-PG setting whose point operations are
        contained in that little group.  For SG205/R/a/PG17 this returns
        `(1, 9, 5, 25, 33, 29)`, matching the GDB `param_4` dump.
        """

        little_ops = set(self.little_parent_point_ops(sg, label))
        for setting in self.site_pg_element_settings(site_pg):
            if setting and all(op in little_ops for op in setting):
                return setting
        raise ValueError(f"no project parent-op setting for SG{sg} {label} PG{site_pg}")

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

    def wyckoff_fraction_vectors(self, row: WyckoffRow) -> tuple[tuple[Fraction, Fraction, Fraction], ...]:
        return self.iso.wyckoff_fraction_vectors(int(row.row_id))

    @staticmethod
    def _vmlt_fraction(
        matrix: Iterable[int],
        vector: tuple[Fraction, Fraction, Fraction],
    ) -> tuple[Fraction, Fraction, Fraction]:
        return ISOData.vmlt_fraction(matrix, vector)

    def vrot_fraction(
        self,
        sg: int,
        point_op: int,
        vector: tuple[Fraction, Fraction, Fraction],
    ) -> tuple[Fraction, Fraction, Fraction]:
        return self.iso.vrot_fraction(int(sg), int(point_op), vector)

    def wyc_pg_elements_records(self, sg: int, row: WyckoffRow) -> tuple[tuple[int, int, int, int, int], ...]:
        return self.iso.wyc_pg_elements_records(int(sg), int(row.row_id), int(row.site_pg))

    def wyc_pg_cosets_records(
        self,
        sg: int,
        site_records: Iterable[tuple[int, int, int, int, int]],
    ) -> tuple[tuple[int, int, int, int, int], ...]:
        records = tuple(tuple(int(value) for value in record) for record in site_records)
        return self.iso.wyc_pg_cosets_records(int(sg), records)

    def stabilizer_parent_point_ops(
        self,
        sg: int,
        row: WyckoffRow,
        site_params: Iterable[object] | None = None,
    ) -> tuple[int, ...]:
        """Return parent point operations that fix the Wyckoff representative."""

        x = self.wyckoff_representative(row, site_params)
        out: list[int] = []
        for op, tau in self.iso.space_ops(sg):
            R = self.iso.point_ops[op - 1]
            delta = tuple(
                (
                    sum(Fraction(int(R[i, j])) * x[j] for j in range(3))
                    + tau[i]
                    - x[i]
                )
                % 1
                for i in range(3)
            )
            if all(value == 0 for value in delta):
                out.append(int(op))
        return tuple(out)

    def stabilizer_parent_operation_records(
        self,
        sg: int,
        row: WyckoffRow,
        site_params: Iterable[object] | None = None,
    ) -> tuple[tuple[int, int, int, int, int], ...]:
        """Return full parent operation records that fix the representative.

        `project_` does not reduce the stabilizer translation modulo the unit
        cell before passing it to `get_irrep4_`.  If a generated operation maps
        the representative back only after an integer lattice shift, that
        integer shift remains in the 5-int record and later contributes the
        Bloch phase at non-Gamma k.
        """

        x = self.wyckoff_representative(row, site_params)
        out: list[tuple[int, int, int, int, int]] = []
        for record in self.space_operation_records(sg):
            op = int(record[4])
            generated_tau = tuple(Fraction(int(record[i]), int(record[3])) for i in range(3))
            R = self.iso.point_ops[op - 1]
            rotated = tuple(sum(Fraction(int(R[i, j])) * x[j] for j in range(3)) for i in range(3))
            delta = tuple((rotated[i] + generated_tau[i] - x[i]) % 1 for i in range(3))
            if all(value == 0 for value in delta):
                exact_tau = tuple(x[i] - rotated[i] for i in range(3))
                den = 1
                for value in exact_tau:
                    den = math.lcm(den, value.denominator)
                out.append(tuple(int(value * den) for value in exact_tau) + (den, op))  # type: ignore[arg-type]
        return tuple(out)

    def project_parent_point_ops_for_row(
        self,
        sg: int,
        label: str,
        row: WyckoffRow,
        site_params: Iterable[object] | None = None,
    ) -> tuple[int, ...]:
        """Select the parent operation setting for a specific Wyckoff row.

        The robust upstream rule is row-specific.  First compute the stabilizer
        of the Wyckoff representative in the full parent space group and match
        its point-operation set to one of the site-PG settings, preserving the
        site-PG setting order that the mode kernel uses for `param_4`.  This fixes cases
        where the little group is smaller than the site stabilizer.

        If the representative-only stabilizer does not have the declared
        site-PG order, fall back to the older little-group subset heuristic;
        those fallback failures remain visible in the broad verifier.
        """

        settings = self.site_pg_element_settings(row.site_pg)
        expected_order = len(settings[0]) if settings else 0
        stabilizer = tuple(record[4] for record in self.stabilizer_parent_operation_records(sg, row, site_params))
        if len(stabilizer) == expected_order:
            stabilizer_set = set(stabilizer)
            for setting in settings:
                if set(setting) == stabilizer_set:
                    return setting
        elif len(stabilizer) > expected_order:
            stabilizer_set = set(stabilizer)
            for setting in settings:
                if setting and set(setting).issubset(stabilizer_set):
                    return setting
        elif 0 < len(stabilizer) < expected_order:
            stabilizer_set = set(stabilizer)
            for setting in settings:
                if setting and stabilizer_set.issubset(set(setting)):
                    return setting
        return self.project_parent_point_ops(sg, label, row.site_pg)

    def project_parent_point_ops_for_little_row(
        self,
        sg: int,
        little: LittleIrrepRecord,
        row: WyckoffRow,
        site_params: Iterable[object] | None = None,
    ) -> tuple[int, ...]:
        """Row-specific parent operation setting using a native little record."""

        settings = self.site_pg_element_settings(row.site_pg)
        expected_order = len(settings[0]) if settings else 0
        stabilizer = tuple(record[4] for record in self.stabilizer_parent_operation_records(sg, row, site_params))
        if len(stabilizer) == expected_order:
            stabilizer_set = set(stabilizer)
            for setting in settings:
                if set(setting) == stabilizer_set:
                    return setting
        elif len(stabilizer) > expected_order:
            stabilizer_set = set(stabilizer)
            for setting in settings:
                if setting and set(setting).issubset(stabilizer_set):
                    return setting
        elif 0 < len(stabilizer) < expected_order:
            stabilizer_set = set(stabilizer)
            for setting in settings:
                if setting and stabilizer_set.issubset(set(setting)):
                    return setting

        little_ops = set(self.little_parent_point_ops_by_gid(sg, little.gid))
        for setting in settings:
            if setting and all(op in little_ops for op in setting):
                return setting
        if settings:
            # The caller's `param_4` is a Wyckoff-site affine setting, not
            # always a subset of the little-group point-op list.  In such
            # cases the mode kernel falls back to the default site-PG setting and
            # carries the required affine translations in the operation
            # records.
            return settings[0]
        raise ValueError(f"no project parent-op setting for SG{sg} {little.label} PG{row.site_pg}")

    def affine_record_for_wyckoff_point_op(
        self,
        row: WyckoffRow,
        point_op: int,
        site_params: Iterable[object] | None = None,
    ) -> tuple[int, int, int, int, int]:
        """Return the exact affine record `{R | x - R x}` for a site operation."""

        x = self.wyckoff_representative(row, site_params)
        R = self.iso.point_ops[int(point_op) - 1]
        tau = tuple(
            x[i] - sum(Fraction(int(R[i, j])) * x[j] for j in range(3))
            for i in range(3)
        )
        den = 1
        for value in tau:
            den = math.lcm(den, value.denominator)
        return tuple(int(value * den) for value in tau) + (den, int(point_op))  # type: ignore[return-value]
