"""Typed access to the distributed ISOTROPY Source tables."""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from functools import lru_cache
import math
from pathlib import Path
import re
from typing import Any, Sequence

from verification.backend import native
from verification.backend.exactmath import (
    integer_determinant3,
)
from verification.backend.cif_numbers import float_cif_number
from verification.backend.source.iso_data import (
    ISOData,
    pml_to_cml_matrix_from_table,
    setting_to_cinter_affine_from_table,
)


ROOT = Path(__file__).resolve().parents[3]
# Distributed data_ text is the single source of truth (shared upstream Source/).
SOURCE = ROOT / "source"
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


@dataclass(frozen=True)
class WyckoffRow:
    offset0: int
    row_id: int
    label: str
    site_pg: int


class SourceTables:
    """Local reader for the Source records shared by every backend domain.

    The distributed `data_*.txt` files under the shared root `Source/` are the
    source of truth. The packaged runtime verifies those files against its
    manifest and loads individual decoded sections only when requested.
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
        tables. The public OPD table, however, displays subgroup
        symbols in the default inter-setting listed by ``ispace_inter_choice``.
        Keeping this separate prevents the parent/k-vector path from silently
        changing convention.
        """

        return str(self.default_inter_setting_record(sg)["label_short"]).strip()

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
    def _normalized_affine4(
        matrix: Sequence[Sequence[int | Fraction]],
    ) -> list[list[Fraction]]:
        denominator = Fraction(matrix[3][3])
        if denominator == 0:
            raise ValueError("zero affine denominator")
        normalized = [
            [Fraction(matrix[row][col]) / denominator for col in range(4)]
            for row in range(4)
        ]
        for row in range(3):
            normalized[row][3] = Fraction(0)
        normalized[3][3] = Fraction(1)
        return normalized

    def subgroup_embedding_from_cinter(
        self,
        parent_sg: int,
        subgroup: int,
        basis: Sequence[Sequence[Fraction]],
        origin: Sequence[Fraction],
        *,
        parent_setting_id: int | None = None,
        subgroup_setting_id: int | None = None,
    ) -> tuple[tuple[int, ...], tuple[int, int, int, int]]:
        """Convert a public International-setting embedding to Source PML."""

        if len(basis) != 3 or any(len(row) != 3 for row in basis):
            raise ValueError("subgroup basis must contain three 3-vectors")
        if len(origin) != 3:
            raise ValueError("subgroup origin must contain three values")
        public = [
            [Fraction(basis[row][col]) for col in range(3)] + [Fraction(0)]
            for row in range(3)
        ]
        public.append([Fraction(value) for value in origin] + [Fraction(1)])

        parent_inter = (
            self._inter_matrix_by_id(int(parent_setting_id), 0)
            if parent_setting_id is not None
            else self._default_inter_matrix(int(parent_sg), 0)
        )
        child_inter = (
            self._inter_matrix_by_id(int(subgroup_setting_id), 1)
            if subgroup_setting_id is not None
            else self._default_inter_matrix(int(subgroup), 1)
        )
        left = self._mat4_mul(
            self._normalized_affine4(child_inter),
            self._normalized_affine4(self._space_cinter_base_matrix(int(subgroup), 0)),
        )
        right = self._mat4_mul(
            self._normalized_affine4(self._space_cinter_base_matrix(int(parent_sg), 1)),
            self._normalized_affine4(parent_inter),
        )
        source = self._mat4_mul(
            self._mat4_mul(self._mat4_inverse_fraction(left), public),
            self._mat4_inverse_fraction(right),
        )
        source_basis = tuple(source[row][col] for row in range(3) for col in range(3))
        if any(value.denominator != 1 for value in source_basis):
            raise ValueError("subgroup embedding does not map to an integral Source basis")
        return tuple(int(value) for value in source_basis), self._reduce2_fraction(list(source[3]))

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

    def little_record_by_gid(self, gid: int) -> LittleIrrep:
        gid = int(gid)
        cache = self.__dict__.setdefault("_little_record_by_gid_cache", {})
        record = cache.get(gid)
        if record is None:
            record = cache[gid] = self._little_record_by_gid(gid)
        return record

    def _little_record_by_gid(self, gid: int) -> LittleIrrep:
        return LittleIrrep(
            gid=gid,
            old_id=int(self.little["little_irr_old"][gid - 1]),
            label=native.irrep_full_label(gid),
            full_dim=int(self.little["little_irr_full_dim"][gid - 1]),
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

    def point_group_order(self, point_group: int) -> int:
        """The number of operations generated by the point group's generators."""

        cache = self.__dict__.setdefault("_point_group_order_cache", {})
        if point_group not in cache:
            space = self.iso.space
            count = int(space["ipoint_group_gen_count"][point_group - 1])
            generators = [int(g) for g in space["ipoint_group_gen"][(point_group - 1) * 5:(point_group - 1) * 5 + count]]
            matrix = lambda op: tuple(int(x) for x in space["ipoint_op"][9 * (op - 1):9 * op])
            identity = (1, 0, 0, 0, 1, 0, 0, 0, 1)
            group = {identity} | {matrix(g) for g in generators}
            while True:
                products = {tuple(sum(a[3 * i + k] * b[3 * k + j] for k in range(3)) for i in range(3) for j in range(3))
                            for a in group for b in group}
                if products <= group:
                    break
                group |= products
            cache[point_group] = len(group)
        return cache[point_group]

    def space_group_point_group_order(self, sg: int) -> int:
        return self.point_group_order(int(self.iso.space["ispace_point_group"][int(sg) - 1]))

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
            native.old_irrep_space_group(old_id),
            native.irrep_old_label(old_id),
        )

    def wyckoff_rows(self, sg: int) -> list[WyckoffRow]:
        sg = int(sg)
        cache = self.__dict__.setdefault("_wyckoff_rows_cache", {})
        rows = cache.get(sg)
        if rows is None:
            rows = cache[sg] = tuple(self._wyckoff_rows(sg))
        return list(rows)

    def _wyckoff_rows(self, sg: int) -> list[WyckoffRow]:
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

    def _inter_wyckoff_pml_vectors(
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
        return tuple(vectors)

    def _inter_wyckoff_fraction_vectors(
        self,
        sg: int,
        row: WyckoffRow,
        setting_id: int,
    ) -> tuple[tuple[Fraction, Fraction, Fraction], ...]:
        cache = self.__dict__.setdefault("_inter_wyckoff_fraction_vectors_cache", {})
        key = (int(sg), row, int(setting_id))
        cached = cache.get(key)
        if cached is None:
            cached = cache[key] = self._inter_wyckoff_fraction_vectors_uncached(sg, row, setting_id)
        return cached

    def _inter_wyckoff_fraction_vectors_uncached(
        self,
        sg: int,
        row: WyckoffRow,
        setting_id: int,
    ) -> tuple[tuple[Fraction, Fraction, Fraction], ...]:
        vectors = self._inter_wyckoff_pml_vectors(sg, row, setting_id)
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

    def _lattice_centering_count(self, sg: int) -> int:
        symbol = self.space_symbol(sg).strip()
        if not symbol:
            return 1
        return self.centering_count_for_symbol(symbol)

    def wyckoff_multiplicity(self, sg: int, row: WyckoffRow) -> int:
        cache = self.__dict__.setdefault("_wyckoff_multiplicity_cache", {})
        key = (int(sg), row)
        if key not in cache:
            cache[key] = self._wyckoff_multiplicity(int(sg), row)
        return cache[key]

    def _wyckoff_multiplicity(self, sg: int, row: WyckoffRow) -> int:
        point_group = int(self.space["ispace_point_group"][int(sg) - 1])
        parent_order = self.point_group_order(point_group)
        site_order = self.point_group_order(int(row.site_pg))
        return self._lattice_centering_count(sg) * parent_order // site_order

    def match_inter_wyckoff_site(
        self,
        sg: int,
        site: dict[str, Any],
        setting_id: int,
        tol: float = 1e-5,
    ) -> dict[str, Any] | None:
        if not site.get("fract"):
            return None
        fract_values = tuple(float_cif_number(x) for x in site["fract"])
        source_values = tuple(
            float_cif_number(x)
            for x in site.get("_source_fract", site["fract"])
        )
        if any(value is None for value in (*fract_values, *source_values)):
            return None
        wanted_mult: int | None = None
        if site.get("multiplicity"):
            match = re.search(r"\d+", str(site["multiplicity"]))
            if match:
                wanted_mult = int(match.group(0))
        found = native.match_inter_wyckoff(
            int(sg),
            int(setting_id),
            tuple(float(value) for value in fract_values),
            tuple(float(value) for value in source_values),
            wanted_mult,
            tol,
        )
        if found is None:
            return None
        offset, values, mask, representative_fract = found
        row = self.wyckoff_rows(sg)[offset]
        formula = self._inter_wyckoff_formula(int(sg), row, int(setting_id))
        params = {name: value for name, value, used in zip(("x", "y", "z"), values, mask) if used}
        mapped_params = {formula["param_map"].get(name, name): value for name, value in params.items()}
        return {
            "row_id": row.row_id,
            "label": row.label,
            "multiplicity": self.wyckoff_multiplicity(sg, row),
            "formula": formula["formula"],
            "free": list(formula["free"]),
            "params": {name: mapped_params[name] for name in formula["free"] if name in mapped_params},
            "inter_setting_id": int(setting_id),
            "formula_representative_fract": representative_fract,
        }

    def _inter_wyckoff_formula(self, sg: int, row: WyckoffRow, setting_id: int) -> dict[str, Any]:
        cache = self.__dict__.setdefault("_inter_wyckoff_formula_cache", {})
        key = (sg, row, setting_id)
        formula = cache.get(key)
        if formula is None:
            formula = cache[key] = self._wyckoff_formula_from_vectors(
                self._inter_wyckoff_fraction_vectors(sg, row, setting_id)
            )
        return formula

    def site_pg_irrep_old_id(self, site_pg: int, pg_irrep: int) -> int:
        start = (int(site_pg) - 1) * 12
        old_id = int(self.wyckoff["iwyckoff_pg_irrep"][start + int(pg_irrep) - 1])
        if old_id <= 0:
            raise KeyError(f"PG{site_pg} irrep {pg_irrep} is not defined")
        return old_id

    def site_pg_irrep_label(self, site_pg: int, pg_irrep: int) -> str:
        """Return the Source point-group label attached to a site irrep."""

        old_id = self.site_pg_irrep_old_id(site_pg, pg_irrep)
        return str(self.iso.irreps["irrep_label_pg"][old_id - 1]).strip()

    def site_vector_reps(self, site_pg: int) -> tuple[int, ...]:
        values = self.wyckoff["iwyckoff_pg_vector_reps"]
        start = (int(site_pg) - 1) * 6
        return tuple(int(x) for x in values[start:start + 6])


@lru_cache(maxsize=1)
def source_tables() -> SourceTables:
    return SourceTables()
