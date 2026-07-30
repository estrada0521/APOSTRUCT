"""Typed access to the distributed ISOTROPY Source tables."""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from functools import lru_cache
import math
from numbers import Integral
from pathlib import Path
import re
import sys
from typing import Any

import numpy as np


ROOT = Path(__file__).resolve().parents[4]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
from ISODISTORT.Assembled.Backend.exactmath import (  # noqa: E402
    fraction_matrix_inverse3,
    fraction_matrix_multiply3,
    integer_determinant3,
)
from ISODISTORT.Assembled.Backend.source.iso_data import (  # noqa: E402
    ISOData,
    pml_to_cml_matrix_from_table,
    setting_to_cinter_affine_from_table,
)

# Distributed data_ text is the single source of truth (shared upstream Source/).
SOURCE = ROOT / "Source"
TYPE_DIVISOR = {1: 1, 2: 4, 3: 2}
DEFAULT_SPACE_GROUP_PREFERENCES = (
    "monoclinic axes a(b)c, monoclinic cell choice 1, "
    "orthorhombic axes abc, origin choice 2, hexagonal axes, "
    "SSG standard setting"
)
SPACE_SETTINGS_CINTER_BASE = 29408


def _fmt_frac(value: Fraction) -> str:
    return str(value.numerator) if value.denominator == 1 else f"{value.numerator}/{value.denominator}"


@dataclass(frozen=True)
class LittleIrrep:
    gid: int
    old_id: int
    label: str
    full_dim: int
    irrep_type: int
    lif: int
    real_pointer: int
    real2_pointer: int


@dataclass(frozen=True)
class WyckoffRow:
    offset0: int
    row_id: int
    label: str
    site_pg: int


@dataclass(frozen=True)
class IsotropySubduction:
    """One entry in data_isotropy:isotropy_subduce_* for an isotropy row."""

    irrep_old_id: int
    irrep_label: str
    sg: int
    frequency: int
    subgroup_row_id: int
    domain: int
    domain_old: int


class SourceTables:
    """Local reader for the Source records shared by every backend domain.

    The distributed `data_*.txt` files under the shared root `Source/` are the
    source of truth.  No decoded pickle cache is used here.
    """

    def __init__(self, source: Path = SOURCE):
        self.source = Path(source)
        self.iso = ISOData(self.source)
        self.space = self.iso.space
        self.irreps = self.iso.irreps
        self.little = self.iso.little
        self.isotropy = self.iso.isotropy
        self.wyckoff = self.iso.wyckoff
        self.images = self.iso.images
        self.const = self.iso.const

    def space_symbol(self, sg: int) -> str:
        return str(self.space["space_label_bc"][int(sg) - 1]).strip()

    def default_setting_space_symbol(self, sg: int) -> str:
        """Return the public default-setting symbol for a space group.

        ``space_label_bc`` is still useful for the decoded little-group/k-vector
        tables.  ISODISTORT's OPD subgroup table, however, displays subgroup
        symbols in the default inter-setting listed by ``ispace_inter_choice``.
        Keeping this separate prevents the parent/k-vector path from silently
        changing convention.
        """

        try:
            return str(self.default_inter_setting_record(sg)["label_short"]).strip()
        except (IndexError, KeyError, ValueError):
            return self.space_symbol(sg)

    def default_space_group_preferences(self) -> str:
        """Return the public ISODISTORT default setting-preference sentence.

        The uploaded parent can still be in a non-default setting.  ISODISTORT
        then adds a separate "Space-group preferences for parent" line; the
        default line itself is a UI preference, not proof that the parent uses
        that setting.
        """

        return DEFAULT_SPACE_GROUP_PREFERENCES

    def inter_setting_record(self, setting_id: int) -> dict[str, Any]:
        setting = int(setting_id)
        index = setting - 1
        if index < 0 or index >= len(self.space["ispace_inter_number"]):
            raise IndexError(f"inter setting {setting} out of range")
        return {
            "id": setting,
            "space_group": int(self.space["ispace_inter_number"][index]),
            "version": int(self.space["ispace_inter_version"][index]),
            "origin": int(self.space["ispace_inter_origin"][index]),
            "cell": int(self.space["ispace_inter_cell"][index]),
            "axis": str(self.space["space_inter_axis"][index]).strip(),
            "abc": str(self.space["space_inter_abc"][index]).strip(),
            "label_full": str(self.space["space_inter_label_full"][index]).strip(),
            "label_short": str(self.space["space_inter_label_short"][index]).strip(),
        }

    def default_inter_setting_record(self, sg: int) -> dict[str, Any]:
        return self.inter_setting_record(int(self.space["ispace_inter_choice"][int(sg) - 1]))

    def inter_setting_ids_for_space_group(self, sg: int) -> list[int]:
        sg = int(sg)
        return [
            index + 1
            for index, number in enumerate(self.space["ispace_inter_number"])
            if int(number) == sg
        ]

    def cml_to_cinter_matrix(
        self,
        sg: int,
        setting_id: int | None = None,
    ) -> tuple[tuple[Fraction, Fraction, Fraction], ...]:
        return setting_to_cinter_affine_from_table(self.space, sg, "cml", setting_id)[0]

    def cml_to_cinter_origin(
        self,
        sg: int,
        setting_id: int | None = None,
    ) -> tuple[Fraction, Fraction, Fraction]:
        return setting_to_cinter_affine_from_table(self.space, sg, "cml", setting_id)[1]

    def pml_to_cml_matrix(self, sg: int) -> tuple[tuple[Fraction, Fraction, Fraction], ...]:
        return pml_to_cml_matrix_from_table(self.space, sg)

    def pml_to_cinter_matrix(
        self,
        sg: int,
        setting_id: int | None = None,
    ) -> tuple[tuple[Fraction, Fraction, Fraction], ...]:
        return setting_to_cinter_affine_from_table(self.space, sg, "pml", setting_id)[0]

    @staticmethod
    def transform_basis_rows(
        basis: tuple[int, ...],
        matrix: tuple[tuple[Fraction, Fraction, Fraction], ...],
    ) -> tuple[tuple[Fraction, Fraction, Fraction], ...]:
        rows = tuple(tuple(Fraction(int(x)) for x in basis[3 * row:3 * row + 3]) for row in range(3))
        return fraction_matrix_multiply3(rows, matrix)

    @staticmethod
    def _fraction_matrix_to_int_rows(
        matrix: tuple[tuple[Fraction, Fraction, Fraction], ...],
    ) -> list[list[int]] | None:
        out: list[list[int]] = []
        for row in matrix:
            if any(value.denominator != 1 for value in row):
                return None
            out.append([int(value) for value in row])
        return out

    @staticmethod
    def _matrix4(raw: list[int] | tuple[int, ...]) -> list[list[int]]:
        if len(raw) != 16:
            raise ValueError(f"expected 16 matrix entries, got {len(raw)}")
        return [list(map(int, raw[row * 4:row * 4 + 4])) for row in range(4)]

    @staticmethod
    def _mat4_mul(left: list[list[int]], right: list[list[int]]) -> list[list[int]]:
        return [
            [sum(left[row][k] * right[k][col] for k in range(4)) for col in range(4)]
            for row in range(4)
        ]

    @classmethod
    def _trans4(cls, first: list[list[int]], middle: list[list[int]], last: list[list[int]]) -> list[list[int]]:
        """Apply ``trans4_``/``matmlt4_`` rules to row-major 4x4 integer matrices.

        This adapter implements ``matmlt4_(A, B)`` as ``B @ A``; therefore
        ``trans4_(A, B, C)`` returns ``C @ B @ A``.
        """

        return cls._mat4_mul(last, cls._mat4_mul(middle, first))

    @staticmethod
    def _mat4_inverse_fraction(matrix: list[list[Fraction]]) -> list[list[Fraction]]:
        size = 4
        aug: list[list[Fraction]] = [
            [matrix[row][col] for col in range(size)]
            + [Fraction(1 if row == col else 0) for col in range(size)]
            for row in range(size)
        ]
        for col in range(size):
            pivot = next((row for row in range(col, size) if aug[row][col] != 0), None)
            if pivot is None:
                raise ValueError(f"singular 4x4 setting matrix: {matrix}")
            if pivot != col:
                aug[col], aug[pivot] = aug[pivot], aug[col]
            scale = aug[col][col]
            aug[col] = [value / scale for value in aug[col]]
            for row in range(size):
                if row == col:
                    continue
                factor = aug[row][col]
                if factor:
                    aug[row] = [aug[row][idx] - factor * aug[col][idx] for idx in range(size * 2)]
        return [row[size:] for row in aug]

    @classmethod
    def _setting_matrix_inverse(cls, matrix: list[list[int]]) -> list[list[Fraction]]:
        den = int(matrix[3][3])
        if den == 0:
            raise ValueError(f"zero setting matrix denominator: {matrix}")
        normalized = [
            [Fraction(matrix[row][col], den) for col in range(4)]
            for row in range(4)
        ]
        normalized[0][3] = Fraction(0)
        normalized[1][3] = Fraction(0)
        normalized[2][3] = Fraction(0)
        normalized[3][3] = Fraction(1)
        return cls._mat4_inverse_fraction(normalized)

    @staticmethod
    def _fraction_values_to_int_den(values: list[Fraction]) -> tuple[list[int], int]:
        den = 1
        for value in values:
            den = math.lcm(den, value.denominator)
        ints = [int(value * den) for value in values]
        gcd = abs(den)
        for value in ints:
            gcd = math.gcd(gcd, abs(value))
        gcd = gcd or 1
        return [value // gcd for value in ints], den // gcd

    @classmethod
    def _reduce3_fraction(cls, matrix: list[list[Fraction]]) -> tuple[list[list[int]], int]:
        values = [matrix[row][col] for row in range(3) for col in range(3)]
        ints, den = cls._fraction_values_to_int_den(values)
        return [ints[row * 3:row * 3 + 3] for row in range(3)], den

    @classmethod
    def _reduce2_fraction(cls, record: list[Fraction]) -> tuple[int, int, int, int]:
        ints, den = cls._fraction_values_to_int_den(record[:3])
        return ints[0], ints[1], ints[2], den

    @staticmethod
    def _reduce3(matrix: list[list[int]], denominator: int) -> tuple[list[list[int]], int]:
        if denominator == 0:
            raise ValueError("zero denominator in subgroup basis")
        sign = -1 if denominator < 0 else 1
        values = [sign * int(value) for row in matrix[:3] for value in row[:3]]
        den = sign * int(denominator)
        gcd = abs(den)
        for value in values:
            gcd = math.gcd(gcd, abs(value))
        gcd = gcd or 1
        return (
            [values[row * 3:row * 3 + 3] for row in range(3)],
            den // gcd,
        ) if gcd == 1 else (
            [[value // gcd for value in values[row * 3:row * 3 + 3]] for row in range(3)],
            den // gcd,
        )

    @staticmethod
    def _reduce2(record: tuple[int, int, int, int]) -> tuple[int, int, int, int]:
        num = [int(record[0]), int(record[1]), int(record[2])]
        den = int(record[3])
        if den == 0:
            raise ValueError("zero denominator in subgroup origin")
        if den < 0:
            num = [-value for value in num]
            den = -den
        gcd = abs(den)
        for value in num:
            gcd = math.gcd(gcd, abs(value))
        gcd = gcd or 1
        return (num[0] // gcd, num[1] // gcd, num[2] // gcd, den // gcd)

    def _space_cinter_base_matrix(self, sg: int, half: int) -> list[list[int]]:
        """Return the SG-indexed conventional/cinter base matrix used by iso.

        The two SG-indexed halves occupy 32 integers per space group in
        ``data_space:ispace_settings`` starting at offset 29408. SG numbers
        index the table directly, so index 0 is a dummy identity block.
        """

        start = SPACE_SETTINGS_CINTER_BASE + int(sg) * 32 + int(half) * 16
        return self._matrix4(list(map(int, self.space["ispace_settings"][start:start + 16])))

    def _default_inter_matrix(self, sg: int, half: int) -> list[list[int]]:
        choice = int(self.space["ispace_inter_choice"][int(sg) - 1])
        return self._inter_matrix_by_id(choice, half)

    def _inter_matrix_by_id(self, setting_id: int, half: int) -> list[list[int]]:
        choice = int(setting_id)
        start = (choice - 1) * 32 + int(half) * 16
        return self._matrix4(list(map(int, self.space["ispace_settings_inter"][start:start + 16])))

    def subgroup_change_setting(
        self,
        parent_sg: int,
        subgroup: int,
        basis: tuple[int, ...],
        origin: tuple[int, int, int, int],
        *,
        setting: str = "cinter",
        parent_setting_id: int | None = None,
        subgroup_setting_id: int | None = None,
    ) -> dict[str, Any]:
        """Apply ``subgroup_change_setting_`` to OPD subgroup basis records.

        ``setting='cml'`` leaves the stored Miller-Love basis untouched, while
        non-``cml`` settings apply SG-level setting transforms. The public
        display path uses ``setting='cinter'``. The inverse transform supports
        parent/child-basis comparisons.
        """

        den = int(origin[3] or 1)
        matrix = [
            [int(basis[0]) * den, int(basis[1]) * den, int(basis[2]) * den, 0],
            [int(basis[3]) * den, int(basis[4]) * den, int(basis[5]) * den, 0],
            [int(basis[6]) * den, int(basis[7]) * den, int(basis[8]) * den, 0],
            [int(origin[0]), int(origin[1]), int(origin[2]), den],
        ]
        normalized_setting = setting.strip().lower()
        if normalized_setting != "cml":
            if normalized_setting != "cinter":
                raise NotImplementedError(f"subgroup_change_setting currently supports cml/cinter, not {setting!r}")
            matrix = self._trans4(
                self._space_cinter_base_matrix(int(parent_sg), 1),
                matrix,
                self._space_cinter_base_matrix(int(subgroup), 0),
            )
            matrix = self._trans4(
                self._inter_matrix_by_id(int(parent_setting_id), 0)
                if parent_setting_id is not None
                else self._default_inter_matrix(int(parent_sg), 0),
                matrix,
                self._inter_matrix_by_id(int(subgroup_setting_id), 1)
                if subgroup_setting_id is not None
                else self._default_inter_matrix(int(subgroup), 1),
            )
        rows, row_den = self._reduce3(matrix, int(matrix[3][3]))
        reduced_origin = self._reduce2((matrix[3][0], matrix[3][1], matrix[3][2], matrix[3][3]))
        inverse = self._setting_matrix_inverse(matrix)
        inverse_rows, inverse_den = self._reduce3_fraction(inverse)
        inverse_origin = self._reduce2_fraction(inverse[3])
        size = self.displayed_primitive_size(int(subgroup), rows)
        if row_den and row_den != 1 and size % abs(row_den) == 0:
            size //= abs(row_den)
        return {
            "basis": rows,
            "basis_denominator": row_den,
            "origin": reduced_origin,
            "inverse_basis": inverse_rows,
            "inverse_basis_denominator": inverse_den,
            "inverse_origin": inverse_origin,
            "size": size,
        }

    def subgroup_change_setting_cinter(
        self,
        parent_sg: int,
        subgroup: int,
        basis: tuple[int, ...],
        origin: tuple[int, int, int, int],
        *,
        parent_setting_id: int | None = None,
        subgroup_setting_id: int | None = None,
    ) -> dict[str, Any]:
        return self.subgroup_change_setting(
            parent_sg,
            subgroup,
            basis,
            origin,
            setting="cinter",
            parent_setting_id=parent_setting_id,
            subgroup_setting_id=subgroup_setting_id,
        )

    @staticmethod
    def basis_determinant(rows: list[list[int]]) -> int:
        flat = tuple(value for row in rows for value in row)
        return abs(integer_determinant3(flat))

    def centering_count_for_symbol(self, symbol: str) -> int:
        symbol = str(symbol).strip()
        if not symbol:
            return 1
        return {"P": 1, "A": 2, "B": 2, "C": 2, "I": 2, "F": 4, "R": 3}.get(symbol[0], 1)

    def displayed_primitive_size(self, subgroup: int, rows: list[list[int]]) -> int:
        det = self.basis_determinant(rows)
        centering = self.centering_count_for_symbol(self.default_setting_space_symbol(subgroup))
        return max(1, det // centering) if centering and det % centering == 0 else det

    @staticmethod
    def parent_preference_text(record: dict[str, Any]) -> str:
        parts: list[str] = []
        origin = int(record.get("origin") or 0)
        cell = int(record.get("cell") or 0)
        abc = str(record.get("abc") or "").strip()
        axis = str(record.get("axis") or "").strip()
        if axis and axis not in {" ", "0"}:
            parts.append(f"axis {axis}")
        if abc:
            parts.append(f"axes {abc}")
        if cell:
            parts.append(f"cell choice {cell}")
        if origin:
            parts.append(f"origin choice {origin}")
        return ", ".join(parts)

    def little_record_by_gid(self, gid: int) -> LittleIrrep:
        gid = int(gid)
        return LittleIrrep(
            gid=gid,
            old_id=int(self.little["little_irr_old"][gid - 1]),
            label=str(self.little["little_irr_full_label"][gid - 1]).strip(),
            full_dim=int(self.little["little_irr_full_dim"][gid - 1]),
            irrep_type=int(self.little["little_irr_type"][gid - 1]),
            lif=int(self.little["little_irr_lif"][gid - 1]),
            real_pointer=int(self.little["little_irr_real_pointer"][gid - 1]),
            real2_pointer=int(self.little["little_irr_real2_pointer"][gid - 1]),
        )

    def image_record(self, old_id: int) -> dict[str, Any]:
        image_id = int(self.irreps["irrep_image"][int(old_id) - 1])
        return {
            "id": image_id,
            "label": str(self.images["image_label"][image_id - 1]).strip(),
            "dimension": int(self.images["image_dimension"][image_id - 1]),
            "order": int(self.images["image_order"][image_id - 1]),
            "type": int(self.images["image_type"][image_id - 1])
            if image_id - 1 < len(self.images.get("image_type", []))
            else None,
            "landau": int(self.images["image_landau"][image_id - 1]),
            "direction_count": int(self.images["image_subgroup_count"][image_id - 1]),
            "subgroup_pointer": int(self.images["image_subgroup_pointer"][image_id - 1]),
        }

    def canonical_isotropy_row_id(self, rows: list[dict[str, Any]]) -> int | None:
        best_id: int | None = None
        best_key: tuple[int, int, int, int, int] | None = None
        for row in rows:
            subgroup = int(row["subgroup"]["number"])
            key = (
                int(row["free"]),
                int(row["det"]),
                -self.space_group_point_group_order(subgroup),
                0 if self.space_group_has_inversion(subgroup) else 1,
                int(row["row_id"]),
            )
            if best_key is None or key < best_key:
                best_key = key
                best_id = int(row["row_id"])
        return best_id

    def space_group_point_group_order(self, sg: int) -> int:
        point_group = int(self.iso.space["ispace_point_group"][int(sg) - 1])
        return int(self.iso.space["ipoint_group_order"][point_group - 1])

    def space_group_has_inversion(self, sg: int) -> bool:
        point_group = int(self.iso.space["ispace_point_group"][int(sg) - 1])
        count = int(self.iso.space["ipoint_group_gen_count"][point_group - 1])
        generators = self.iso.space["ipoint_group_gen"][(point_group - 1) * 5:(point_group - 1) * 5 + count]
        return any(int(generator) in {25, 61} for generator in generators)

    def isotropy_row_ids_for_old_irrep(self, old_id: int) -> range:
        pointers = self.isotropy["isotropy_irrep_pointer"]
        old_id = int(old_id)
        if old_id < 1 or old_id >= len(pointers):
            return range(0)
        return range(int(pointers[old_id - 1]), int(pointers[old_id]))

    def irrep_label_for_old_id(self, old_id: int) -> tuple[int, str]:
        old_id = int(old_id)
        return (
            int(self.iso.irreps["irrep_space_group"][old_id - 1]),
            str(self.iso.irreps["irrep_label"][old_id - 1]).strip(),
        )

    def isotropy_subductions_for_row(self, row_id: int) -> tuple[IsotropySubduction, ...]:
        """Return secondary/primary irreps allowed by one isotropy subgroup row.

        In ISODISTORT complete mode details, selecting an OPD fixes an
        isotropy subgroup.  The native ``isotropy_subduce_*`` CSR tables map
        that subgroup row back to the parent irreps/OPD rows that are allowed
        under the same subgroup, which is the Source-level hook for secondary
        mode blocks.
        """

        row_id = int(row_id)
        pointers = self.iso.isotropy["isotropy_subduce_pointer"]
        if row_id < 1 or row_id > len(pointers):
            return ()
        start = int(pointers[row_id - 1])
        end = (
            int(pointers[row_id])
            if row_id < len(pointers)
            else len(self.iso.isotropy["isotropy_subduce_irrep"]) + 1
        )
        out: list[IsotropySubduction] = []
        for index in range(start - 1, end - 1):
            old_id = int(self.iso.isotropy["isotropy_subduce_irrep"][index])
            sg, label = self.irrep_label_for_old_id(old_id)
            out.append(
                IsotropySubduction(
                    irrep_old_id=old_id,
                    irrep_label=label,
                    sg=sg,
                    frequency=int(self.iso.isotropy["isotropy_subduce_frequency"][index]),
                    subgroup_row_id=int(self.iso.isotropy["isotropy_subduce_subgroup"][index]),
                    domain=int(self.iso.isotropy["isotropy_subduce_domain"][index]),
                    domain_old=int(self.iso.isotropy["isotropy_subduce_domain_old"][index]),
                )
            )
        return tuple(out)

    def wyckoff_rows(self, sg: int) -> list[WyckoffRow]:
        sg = int(sg)
        start = int(self.wyckoff["iwyckoff_pointer"][sg - 1])
        count = int(self.wyckoff["iwyckoff_count"][sg - 1])
        return [
            WyckoffRow(
                offset0=offset0,
                row_id=start + offset0,
                label=str(self.wyckoff["wyckoff_label"][start + offset0 - 1]).strip(),
                site_pg=int(self.wyckoff["iwyckoff_pg"][start + offset0 - 1]),
            )
            for offset0 in range(count)
        ]

    def wyckoff_row_by_id(self, sg: int, row_id: int) -> WyckoffRow:
        for row in self.wyckoff_rows(sg):
            if int(row.row_id) == int(row_id):
                return row
        raise KeyError(f"unknown Wyckoff row {row_id} for SG{sg}")

    def wyckoff_row_by_label(self, sg: int, label: str) -> WyckoffRow:
        wanted = str(label).strip().casefold()
        for row in self.wyckoff_rows(sg):
            if row.label.casefold() == wanted:
                return row
        raise KeyError(f"unknown Wyckoff label {label!r} for SG{sg}")

    def wyckoff_fraction_vectors(self, row: WyckoffRow) -> tuple[tuple[Fraction, Fraction, Fraction], ...]:
        return self.iso.wyckoff_fraction_vectors(int(row.row_id))

    def site_pg_element_settings(self, site_pg: int) -> tuple[tuple[int, ...], ...]:
        return self.iso.site_pg_element_settings(int(site_pg))

    def generate_space_group_records(self, sg: int) -> tuple[tuple[int, int, int, int, int], ...]:
        return self.iso.generate_space_group_records(int(sg))

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
        site_records: tuple[tuple[int, int, int, int, int], ...],
    ) -> tuple[tuple[int, int, int, int, int], ...]:
        return self.iso.wyc_pg_cosets_records(int(sg), tuple(site_records))

    def _inter_wyckoff_fraction_vectors(
        self,
        sg: int,
        row: WyckoffRow,
        setting_id: int,
    ) -> tuple[tuple[Fraction, Fraction, Fraction], ...]:
        pointer = int(self.space["ispace_inter_wyckoff_pointer"][int(setting_id) - 1])
        raw_start = (pointer - 1 + int(row.offset0)) * 16
        raw = tuple(int(x) for x in self.space["ispace_inter_wyckoff_fract"][raw_start:raw_start + 16])
        vectors: list[tuple[Fraction, Fraction, Fraction]] = []
        for start in range(0, 16, 4):
            den = int(raw[start + 3])
            if den == 0:
                vectors.append((Fraction(0), Fraction(0), Fraction(0)))
            else:
                vectors.append(tuple(Fraction(raw[start + axis], den) for axis in range(3)))  # type: ignore[arg-type]
        matrix = self.pml_to_cinter_matrix(sg, setting_id)
        transformed = [
            tuple(sum(vector[col] * matrix[col][row_axis] for col in range(3)) for row_axis in range(3))
            for vector in vectors
        ]
        origin = self.cml_to_cinter_origin(sg, setting_id)
        transformed[0] = tuple(transformed[0][axis] + origin[axis] for axis in range(3))
        return tuple(transformed)  # type: ignore[return-value]

    def _wyckoff_formula_from_vectors(
        self,
        vectors: tuple[tuple[Fraction, Fraction, Fraction], ...],
    ) -> dict[str, Any]:
        raw_params = ("x", "y", "z")
        used_params = [
            raw_params[index - 1]
            for index in range(1, 4)
            if any(vectors[index])
        ]
        param_map = {name: name for name in used_params}
        components: list[str] = []
        for axis in range(3):
            variable_terms: list[str] = []
            base = vectors[0][axis]
            for index, raw_param in enumerate(raw_params, start=1):
                coef = vectors[index][axis]
                if not coef:
                    continue
                param = raw_param
                if coef == 1:
                    variable_terms.append(param)
                elif coef == -1:
                    variable_terms.append("-" + param)
                else:
                    variable_terms.append(f"{_fmt_frac(coef)}{param}")
            terms = list(variable_terms)
            if base:
                terms.append(_fmt_frac(base))
            if not terms:
                components.append("0")
                continue
            text = terms[0]
            for term in terms[1:]:
                text += term if term.startswith("-") else "+" + term
            components.append(text)
        return {"formula": "(" + ",".join(components) + ")", "free": used_params, "param_map": param_map}

    def wyckoff_formula(self, row: WyckoffRow) -> dict[str, Any]:
        return self._wyckoff_formula_from_vectors(self.wyckoff_fraction_vectors(row))

    def inter_wyckoff_formula(self, sg: int, row: WyckoffRow, setting_id: int) -> dict[str, Any]:
        return self._wyckoff_formula_from_vectors(self._inter_wyckoff_fraction_vectors(sg, row, setting_id))

    def inter_wyckoff_fraction_vectors(
        self,
        sg: int,
        row: WyckoffRow,
        setting_id: int | None = None,
    ) -> tuple[tuple[Fraction, Fraction, Fraction], ...]:
        return self._inter_wyckoff_fraction_vectors(
            int(sg),
            row,
            int(setting_id or self.space["ispace_inter_choice"][int(sg) - 1]),
        )

    def _lattice_centering_count(self, sg: int) -> int:
        symbol = self.space_symbol(sg).strip()
        if not symbol:
            return 1
        return self.centering_count_for_symbol(symbol)

    def wyckoff_multiplicity(self, sg: int, row: WyckoffRow) -> int:
        point_group = int(self.space["ispace_point_group"][int(sg) - 1])
        parent_order = int(self.space["ipoint_group_order"][point_group - 1])
        site_order = int(self.space["ipoint_group_order"][int(row.site_pg) - 1])
        return self._lattice_centering_count(sg) * parent_order // site_order

    @staticmethod
    def _mod_distance(value: float) -> float:
        return abs(value - round(value))

    @staticmethod
    @lru_cache(maxsize=None)
    def _smith_normal_decomp_int(
        rows: tuple[tuple[int, ...], ...],
    ) -> tuple[tuple[tuple[int, ...], ...], tuple[tuple[int, ...], ...], tuple[tuple[int, ...], ...]]:
        from sympy import Matrix, ZZ
        from sympy.matrices.normalforms import smith_normal_decomp

        diagonal, left, right = smith_normal_decomp(Matrix([list(row) for row in rows]), domain=ZZ)
        return (
            tuple(tuple(int(diagonal[row, col]) for col in range(diagonal.shape[1])) for row in range(diagonal.shape[0])),
            tuple(tuple(int(left[row, col]) for col in range(left.shape[1])) for row in range(left.shape[0])),
            tuple(tuple(int(right[row, col]) for col in range(right.shape[1])) for row in range(right.shape[0])),
        )

    @classmethod
    def _solve_integer_modular_system(
        cls,
        matrix: np.ndarray,
        rhs: np.ndarray,
        tol: float,
    ) -> np.ndarray | None:
        if matrix.size == 0:
            return np.array([], dtype=float)
        rounded = np.rint(matrix).astype(int)
        if not np.allclose(matrix, rounded, atol=tol):
            return None
        rows = tuple(tuple(int(rounded[row, col]) for col in range(rounded.shape[1])) for row in range(rounded.shape[0]))
        try:
            diagonal, left, right = cls._smith_normal_decomp_int(rows)
        except Exception:
            return None
        left_matrix = np.array(left, dtype=float)
        right_matrix = np.array(right, dtype=float)
        diag_matrix = np.array(diagonal, dtype=float)
        transformed = left_matrix @ rhs
        y = np.zeros(right_matrix.shape[0], dtype=float)
        for index in range(min(diag_matrix.shape)):
            diag = diag_matrix[index, index]
            value = float(transformed[index] % 1.0)
            if abs(value - 1.0) <= tol:
                value = 0.0
            if abs(diag) <= tol:
                if abs(value) > tol:
                    return None
                continue
            y[index] = value / diag
        for index in range(min(diag_matrix.shape), len(transformed)):
            value = float(transformed[index] % 1.0)
            if abs(value - 1.0) <= tol:
                value = 0.0
            if abs(value) > tol:
                return None
        solution = right_matrix @ y
        residual = matrix @ solution - rhs
        if all(cls._mod_distance(float(value)) <= tol for value in residual):
            return solution
        return None

    def _solve_wyckoff_params(
        self,
        row: WyckoffRow,
        fract: tuple[float, float, float],
        tol: float = 1e-5,
    ) -> dict[str, float] | None:
        return self._solve_wyckoff_params_from_vectors(self.wyckoff_fraction_vectors(row), fract, tol=tol)

    def _solve_wyckoff_params_from_vectors(
        self,
        vectors: tuple[tuple[Fraction, Fraction, Fraction], ...],
        fract: tuple[float, float, float],
        tol: float = 1e-5,
        *,
        modulo: bool = True,
    ) -> dict[str, float] | None:
        base = np.array([float(x) for x in vectors[0]], dtype=float)
        columns = [
            np.array([float(vectors[index][axis]) for axis in range(3)], dtype=float)
            for index in range(1, 4)
            if any(vectors[index])
        ]
        names = [
            name
            for index, name in enumerate(("x", "y", "z"), start=1)
            if any(vectors[index])
        ]
        target = np.array(fract, dtype=float)
        if not columns:
            delta = target - base
            return {} if all(self._mod_distance(float(x)) <= tol for x in delta) else None
        matrix = np.column_stack(columns)
        rhs = target - base
        solution = self._solve_integer_modular_system(matrix, rhs, tol)
        if solution is None:
            wrapped_rhs = np.array([float(value % 1.0) for value in rhs], dtype=float)
            wrapped_rhs[np.isclose(wrapped_rhs, 1.0, atol=tol)] = 0.0
            solution, *_ = np.linalg.lstsq(matrix, wrapped_rhs, rcond=None)
            residual = matrix @ solution - wrapped_rhs
            if float(np.linalg.norm(residual, ord=np.inf)) > tol:
                return None
        if solution is not None:
            return {
                name: float(value % 1.0) if modulo else float(value)
                for name, value in zip(names, solution)
            }
        return None

    @staticmethod
    def _source_param_value(value: float, tol: float) -> float:
        normalized = float(value)
        while abs(normalized - 1.0) < abs(normalized):
            normalized -= 1.0
        while abs(normalized + 1.0) < abs(normalized):
            normalized += 1.0
        return 0.0 if abs(normalized) <= tol else normalized

    def _solve_inter_wyckoff_params_source_branch(
        self,
        sg: int,
        row: WyckoffRow,
        fract: tuple[float, float, float],
        tol: float = 1e-5,
    ) -> dict[str, float] | None:
        source_setting_id = int(self.space["ispace_inter_choice"][int(sg) - 1])
        matrix = self.pml_to_cinter_matrix(sg, source_setting_id)
        inverse = fraction_matrix_inverse3(matrix)
        origin = self.cml_to_cinter_origin(sg, source_setting_id)
        cinter = tuple(Fraction(str(value)) for value in fract)
        target = tuple(
            sum((cinter[col] - origin[col]) * inverse[col][axis] for col in range(3))
            for axis in range(3)
        )
        base, *parameters = self.wyckoff_fraction_vectors(row)
        site_records = self.wyc_pg_elements_records(sg, row)
        cosets = self.wyc_pg_cosets_records(sg, site_records)
        branch_vectors: list[tuple[tuple[Fraction, Fraction, Fraction], ...]] = [(base, *parameters)]
        for coset in cosets[1:]:
            op = int(coset[4])
            den = int(coset[3])
            tau = tuple(Fraction(int(coset[axis]), den) for axis in range(3))
            rotated_base = self.vrot_fraction(sg, op, base)
            branch_vectors.append((
                tuple(rotated_base[axis] + tau[axis] for axis in range(3)),
                *(self.vrot_fraction(sg, op, parameter) for parameter in parameters),
            ))
        for vectors in branch_vectors:
            params = self._solve_wyckoff_params_from_vectors(
                vectors,
                tuple(float(x) for x in target),
                tol=tol,
                modulo=False,
            )
            if params is None:
                continue
            return {
                name: self._source_param_value(value, tol)
                for name, value in params.items()
            }
        return None

    def match_wyckoff_site(self, sg: int, site: dict[str, Any], tol: float = 1e-5) -> dict[str, Any] | None:
        if not site.get("fract"):
            return None
        try:
            fract = tuple(float(x) for x in site["fract"])
        except (TypeError, ValueError):
            return None
        wanted_mult: int | None = None
        if site.get("multiplicity"):
            match = re.search(r"\d+", str(site["multiplicity"]))
            if match:
                wanted_mult = int(match.group(0))
        candidates: list[dict[str, Any]] = []
        for row in self.wyckoff_rows(sg):
            mult = self.wyckoff_multiplicity(sg, row)
            if wanted_mult is not None and mult != wanted_mult:
                continue
            params = self._solve_wyckoff_params(row, fract, tol=tol)
            if params is None:
                continue
            formula = self.wyckoff_formula(row)
            mapped_params = {
                formula["param_map"].get(name, name): value
                for name, value in params.items()
            }
            display_params = {
                name: mapped_params[name]
                for name in formula["free"]
                if name in mapped_params
            }
            candidates.append({
                "row_id": row.row_id,
                "label": row.label,
                "multiplicity": mult,
                "formula": formula["formula"],
                "free": formula["free"],
                "params": display_params,
            })
        if not candidates:
            return None
        return candidates[0]

    def match_inter_wyckoff_site(
        self,
        sg: int,
        site: dict[str, Any],
        setting_id: int,
        tol: float = 1e-5,
    ) -> dict[str, Any] | None:
        if not site.get("fract"):
            return None
        try:
            fract = tuple(float(x) for x in site["fract"])
            source_fract = tuple(float(x) for x in site.get("_source_fract", site["fract"]))
        except (TypeError, ValueError):
            return None
        wanted_mult: int | None = None
        if site.get("multiplicity"):
            match = re.search(r"\d+", str(site["multiplicity"]))
            if match:
                wanted_mult = int(match.group(0))
        candidates: list[dict[str, Any]] = []
        for row in self.wyckoff_rows(sg):
            mult = self.wyckoff_multiplicity(sg, row)
            if wanted_mult is not None and mult != wanted_mult:
                continue
            vectors = self._inter_wyckoff_fraction_vectors(sg, row, setting_id)
            params = self._solve_wyckoff_params_from_vectors(vectors, fract, tol=tol)
            if params is None:
                continue
            display_source_params = self._solve_inter_wyckoff_params_source_branch(
                sg,
                row,
                source_fract,
                tol=tol,
            )
            formula = self._wyckoff_formula_from_vectors(vectors)
            if display_source_params is not None:
                params = display_source_params
            representative_fract = tuple(
                float(vectors[0][axis])
                + sum(
                    float(vectors[index][axis]) * float(params.get(name, 0.0))
                    for index, name in enumerate(("x", "y", "z"), start=1)
                )
                for axis in range(3)
            )
            mapped_params = {
                formula["param_map"].get(name, name): value
                for name, value in params.items()
            }
            display_params = {
                name: mapped_params[name]
                for name in formula["free"]
                if name in mapped_params
            }
            candidates.append({
                "row_id": row.row_id,
                "label": row.label,
                "multiplicity": mult,
                "formula": formula["formula"],
                "free": formula["free"],
                "params": display_params,
                "inter_setting_id": int(setting_id),
                "formula_representative_fract": representative_fract,
            })
        if not candidates:
            return None
        return candidates[0]

    def site_pg_irrep_old_id(self, site_pg: int, pg_irrep: int) -> int:
        start = (int(site_pg) - 1) * 12
        old_id = int(self.wyckoff["iwyckoff_pg_irrep"][start + int(pg_irrep) - 1])
        if old_id <= 0:
            raise KeyError(f"PG{site_pg} irrep {pg_irrep} is not defined")
        return old_id

    def site_pg_irrep_label(self, site_pg: int, pg_irrep: int) -> str:
        """Return the Source point-group label attached to a site irrep."""

        for field, value in (("site point group", site_pg), ("site irrep", pg_irrep)):
            if isinstance(value, bool) or not isinstance(value, Integral):
                raise TypeError(f"{field} must be an exact positive integer")
            if value <= 0:
                raise ValueError(f"{field} must be positive, got {value}")
        site_pg = int(site_pg)
        pg_irrep = int(pg_irrep)
        wyckoff = self.iso.wyckoff
        counts = wyckoff["iwyckoff_pg_irrep_count"]
        if site_pg > len(counts):
            raise KeyError(f"site point group out of range: {site_pg}")
        raw_count = counts[site_pg - 1]
        if isinstance(raw_count, bool) or not isinstance(raw_count, Integral):
            raise TypeError(
                f"PG{site_pg} Source irrep count must be an exact integer"
            )
        count = int(raw_count)
        if count <= 0 or count > 12:
            raise ValueError(f"PG{site_pg} Source irrep count is invalid: {count}")
        if pg_irrep > count:
            raise KeyError(f"PG{site_pg} irrep {pg_irrep} is not defined")
        raw_old_id = wyckoff["iwyckoff_pg_irrep"][
            (site_pg - 1) * 12 + pg_irrep - 1
        ]
        if isinstance(raw_old_id, bool) or not isinstance(raw_old_id, Integral):
            raise TypeError(
                f"PG{site_pg} irrep {pg_irrep} Source old id must be an exact integer"
            )
        old_id = int(raw_old_id)
        labels = self.iso.irreps["irrep_label_pg"]
        if old_id <= 0 or old_id > len(labels):
            raise KeyError(
                f"PG{site_pg} irrep {pg_irrep} Source old id is invalid: {old_id}"
            )
        raw_label = labels[old_id - 1]
        if not isinstance(raw_label, str) or not raw_label.strip():
            raise KeyError(
                f"PG{site_pg} irrep {pg_irrep} has no Source point-group label"
            )
        return raw_label.strip()

    def site_vector_reps(self, site_pg: int) -> tuple[int, ...]:
        values = self.wyckoff["iwyckoff_pg_vector_reps"]
        start = (int(site_pg) - 1) * 6
        return tuple(int(x) for x in values[start:start + 6])

    def wyckoff_subduction_pairs_for_little_gid(
        self,
        sg: int,
        row_id: int,
        gid: int,
    ) -> tuple[tuple[int, int], ...]:
        """Return site-PG irrep/frequency pairs for one Wyckoff row and little irrep.

        This is the Assembled backend equivalent of the upstream SMODES
        pre-project count path.  Fixed-k irreps use ``iwyckoff_*`` through
        ``little_irr_old``; parametric little rows use ``little_subduce_*``.
        """

        little = self.little_record_by_gid(gid)
        row = self.wyckoff_row_by_id(sg, row_id)
        if little.old_id > 0:
            base = int(self.wyckoff["iwyckoff_irrep_pointer"][little.old_id - 1])
            entry = base + row.offset0
            pointer = int(self.wyckoff["iwyckoff_subduce_pointer"][entry - 1])
            count = int(self.wyckoff["iwyckoff_subduce_count"][entry - 1])
            pg_values = self.wyckoff["iwyckoff_subduce_pg_irrep"]
            freq_values = self.wyckoff["iwyckoff_subduce_frequency"]
        else:
            base = int(self.little["little_subduce_irr_pointer"][gid - 1])
            if base <= 0:
                return ()
            entry = base + row.offset0
            pointer = int(self.little["little_subduce_pointer"][entry - 1])
            count = int(self.little["little_subduce_count"][entry - 1])
            pg_values = self.little["little_subduce_pg_irrep"]
            freq_values = self.little["little_subduce_frequency"]
        if pointer <= 0 or count <= 0:
            return ()
        return tuple(
            (int(pg_values[pointer - 1 + index]), int(freq_values[pointer - 1 + index]))
            for index in range(count)
        )

    def _vector_mode_total_for_row_id(self, sg: int, row_id: int, gid: int, vector_slice: slice) -> int:
        """Return vector-mode count for one site/irrep pair.

        ``vector_slice`` selects polar displacement reps (0:3) or axial
        magnetic-moment reps (3:6) from ``iwyckoff_pg_vector_reps``.
        """

        little = self.little_record_by_gid(gid)
        row = self.wyckoff_row_by_id(sg, row_id)
        divisor = TYPE_DIVISOR.get(little.irrep_type)
        if divisor is None:
            return 0
        project_counts: dict[int, int] = {}
        site_dims: dict[int, int] = {}
        for pg_irrep, frequency in self.wyckoff_subduction_pairs_for_little_gid(sg, row_id, gid):
            site_old_id = self.site_pg_irrep_old_id(row.site_pg, pg_irrep)
            site_image = self.image_record(site_old_id)
            site_dims[int(pg_irrep)] = int(site_image["dimension"])
            factor = 1 if int(site_image["type"] or 1) == 1 else 2
            numerator = int(frequency) * factor
            if numerator % divisor:
                return 0
            project_counts[int(pg_irrep)] = numerator // divisor
        project_total = Fraction(0, 1)
        for pg_irrep in self.site_vector_reps(row.site_pg)[vector_slice]:
            pg_irrep = int(pg_irrep)
            if pg_irrep in project_counts and site_dims.get(pg_irrep):
                project_total += Fraction(project_counts[pg_irrep], site_dims[pg_irrep])
        total = project_total * little.full_dim
        return int(total) if total.denominator == 1 else 0

    def displacive_mode_total_for_row_id(self, sg: int, row_id: int, gid: int) -> int:
        """Return total displacement-mode count for one site/irrep pair."""

        return self._vector_mode_total_for_row_id(sg, row_id, gid, slice(0, 3))

    def magnetic_mode_total_for_row_id(self, sg: int, row_id: int, gid: int) -> int:
        """Return total axial magnetic-moment mode count for one site/irrep pair."""

        return self._vector_mode_total_for_row_id(sg, row_id, gid, slice(3, 6))


@lru_cache(maxsize=1)
def source_tables() -> SourceTables:
    return SourceTables()
