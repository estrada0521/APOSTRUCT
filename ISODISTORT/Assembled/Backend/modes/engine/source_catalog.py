from __future__ import annotations

from pathlib import Path
from typing import Any

import numpy as np

from ISODISTORT.Assembled.Backend.exactmath import integer_determinant3
from ISODISTORT.Assembled.Backend.modes.engine.records import (
    ImageRecord,
    IsotropySubgroupRow,
    LittleIrrepRecord,
    WyckoffRow,
    WyckoffSubduction,
)
from ISODISTORT.Assembled.Backend.source.iso_data import ISOData
from ISODISTORT.Assembled.Backend.source.tables import SOURCE, SourceTables


class ModeSourceCatalogMixin:
    def __init__(self, data_dir: str | Path = SOURCE, *, iso: Any | None = None):
        self.data_dir = Path(data_dir)
        if iso is None:
            self.iso = ISOData(self.data_dir)
        else:
            table_dir = Path(getattr(iso, "data_dir", self.data_dir))
            if table_dir.resolve() != self.data_dir.resolve():
                raise ValueError(
                    f"mode table directory mismatch: {table_dir} != {self.data_dir}"
                )
            self.iso = iso
        self.images = self.iso.images

    def old_irrep_id(self, sg: int, label: str) -> int:
        return self.iso.old_irrep_id(sg, label)

    def little_record(self, sg: int, label: str) -> LittleIrrepRecord:
        old_id = self.old_irrep_id(sg, label)
        gids = [
            index + 1
            for index, value in enumerate(self.iso.little["little_irr_old"])
            if value == old_id and self.iso.little["little_irr_space_group"][index] == sg
        ]
        if len(gids) != 1:
            raise KeyError(f"expected one little irrep for SG{sg} {label}, got {gids}")
        gid = gids[0]
        real2_pointer = int(self.iso.little["little_irr_real2_pointer"][gid - 1])
        return LittleIrrepRecord(
            gid=gid,
            old_id=old_id,
            label=self.iso.little["little_irr_full_label"][gid - 1].strip(),
            full_dim=int(self.iso.little["little_irr_full_dim"][gid - 1]),
            irrep_type=int(self.iso.little["little_irr_type"][gid - 1]),
            lif=int(self.iso.little["little_irr_lif"][gid - 1]),
            real_pointer=int(self.iso.little["little_irr_real_pointer"][gid - 1]),
            real2_pointer=real2_pointer,
            real2_slice=self._pointed_slice(
                self.iso.little["little_irr_real2_pointer"],
                self.iso.little["little_irr_real2"],
                gid,
            ),
        )

    def little_record_by_gid(self, gid: int) -> LittleIrrepRecord:
        """Return little-irrep metadata by the native data_little row id.

        Some `little_irr_full_label` values are combined labels that are not
        directly present as `data_irreps:irrep_label`.  The upstream path
        already has the little row / old irrep id, so broad comparisons should
        not depend on label reverse lookup.
        """

        gid = int(gid)
        real2_pointer = int(self.iso.little["little_irr_real2_pointer"][gid - 1])
        return LittleIrrepRecord(
            gid=gid,
            old_id=int(self.iso.little["little_irr_old"][gid - 1]),
            label=self.iso.little["little_irr_full_label"][gid - 1].strip(),
            full_dim=int(self.iso.little["little_irr_full_dim"][gid - 1]),
            irrep_type=int(self.iso.little["little_irr_type"][gid - 1]),
            lif=int(self.iso.little["little_irr_lif"][gid - 1]),
            real_pointer=int(self.iso.little["little_irr_real_pointer"][gid - 1]),
            real2_pointer=real2_pointer,
            real2_slice=self._pointed_slice(
                self.iso.little["little_irr_real2_pointer"],
                self.iso.little["little_irr_real2"],
                gid,
            ),
        )

    def isotropy_rows_for_old_irrep(self, old_id: int) -> tuple[IsotropySubgroupRow, ...]:
        """Return `data_isotropy` candidates attached to an old irrep id.

        MAIN__ scans these rows when it chooses the transformation printed as
        "Vectors defining superlattice".  The exact row-selection rule is still
        being ported, so this helper exposes all candidates.
        """

        old_id = int(old_id)
        pointers = self.iso.isotropy["isotropy_irrep_pointer"]
        if old_id < 1 or old_id >= len(pointers):
            raise KeyError(f"old irrep id out of isotropy range: {old_id}")
        start = int(pointers[old_id - 1])
        end = int(pointers[old_id])
        rows: list[IsotropySubgroupRow] = []
        for row_id in range(start, end):
            idx = row_id - 1
            rows.append(
                IsotropySubgroupRow(
                    row_id=row_id,
                    subgroup=int(self.iso.isotropy["isotropy_subgroup"][idx]),
                    arms=int(self.iso.isotropy["isotropy_arms"][idx]),
                    direction=int(self.iso.isotropy["isotropy_direction"][idx]),
                    basis=tuple(
                        int(x) for x in self.iso.isotropy["isotropy_basis"][idx * 9:(idx + 1) * 9]
                    ),
                    origin=tuple(
                        int(x) for x in self.iso.isotropy["isotropy_origin"][idx * 4:(idx + 1) * 4]
                    ),  # type: ignore[arg-type]
                )
            )
        return tuple(rows)

    def isotropy_rows_for_little_gid(self, gid: int) -> tuple[IsotropySubgroupRow, ...]:
        little = self.little_record_by_gid(gid)
        if little.old_id <= 0:
            return ()
        return self.isotropy_rows_for_old_irrep(little.old_id)

    irrep_label_for_old_id = SourceTables.irrep_label_for_old_id
    isotropy_subductions_for_row = SourceTables.isotropy_subductions_for_row

    def isotropy_orderparam_matrix(self, row_id: int, full_dim: int) -> tuple[tuple[float, ...], ...]:
        """Return the real order-parameter matrix for one `data_isotropy` row.

        MAIN__ reads `isotropy_orderparam_pointer` and converts the stored
        constant codes with `constant_()` before building the atom/mode source
        matrix consumed by the `project_` -> `project_vector_` bridge.  The
        matrix has `isotropy_orderparam_dim[row]` rows and `full_dim` columns.
        """

        row_id = int(row_id)
        full_dim = int(full_dim)
        idx = row_id - 1
        if idx < 0 or idx >= len(self.iso.isotropy["isotropy_orderparam_dim"]):
            raise KeyError(f"isotropy row out of range: {row_id}")
        row_count = int(self.iso.isotropy["isotropy_orderparam_dim"][idx])
        pointer = int(self.iso.isotropy["isotropy_orderparam_pointer"][idx]) - 1
        values = self.iso.isotropy["isotropy_orderparam"]
        out: list[tuple[float, ...]] = []
        for row_index in range(row_count):
            row: list[float] = []
            for col_index in range(full_dim):
                code_index = pointer + row_index * full_dim + col_index
                row.append(float(self.iso.const[int(values[code_index])]))
            out.append(tuple(row))
        return tuple(out)

    def isotropy_orderparam_dim(self, row_id: int) -> int:
        row_id = int(row_id)
        return int(self.iso.isotropy["isotropy_orderparam_dim"][row_id - 1])

    def isotropy_orderparam_freeparam(self, row_id: int) -> int:
        row_id = int(row_id)
        return int(self.iso.isotropy["isotropy_orderparam_freeparam"][row_id - 1])

    def isotropy_basis_abs_det(self, row: IsotropySubgroupRow) -> int:
        return abs(integer_determinant3(row.basis))

    space_group_has_inversion = SourceTables.space_group_has_inversion
    space_group_point_group_order = SourceTables.space_group_point_group_order

    def canonical_isotropy_row(self, rows: tuple[IsotropySubgroupRow, ...]) -> IsotropySubgroupRow | None:
        """Return the fixed-k isotropy row selected by MAIN__ for mode sources.

        The fixed-k branch minimizes the number of free order-parameter
        parameters, then the determinant of the subgroup basis.  Ties prefer
        larger subgroup point-group order and finally the generator-list flag
        used by MAIN__ for point operations 25/61.
        """

        best: IsotropySubgroupRow | None = None
        best_key: tuple[int, int, int, int] | None = None
        for row in rows:
            key = (
                self.isotropy_orderparam_freeparam(row.row_id),
                self.isotropy_basis_abs_det(row),
                -self.space_group_point_group_order(row.subgroup),
                0 if self.space_group_has_inversion(row.subgroup) else 1,
                row.row_id,
            )
            if best_key is None or key < best_key:
                best = row
                best_key = key
        return best

    def image_record(self, old_id: int) -> ImageRecord:
        image_id = int(self.iso.irreps["irrep_image"][old_id - 1])

        def optional(name: str):
            values = self.images[name]
            return values[image_id - 1] if image_id - 1 < len(values) else None

        return ImageRecord(
            image_id=image_id,
            label=self.images["image_label"][image_id - 1].strip(),
            dimension=int(self.images["image_dimension"][image_id - 1]),
            order=int(self.images["image_order"][image_id - 1]),
            landau=int(self.images["image_landau"][image_id - 1]),
            subgroup_count=int(self.images["image_subgroup_count"][image_id - 1]),
            subgroup_pointer=int(self.images["image_subgroup_pointer"][image_id - 1]),
            label2=(optional("image_label2") or "").strip() or None,
            image_type=None if optional("image_type") is None else int(optional("image_type")),
            generator_count=None
            if optional("image_generator_count") is None
            else int(optional("image_generator_count")),
            diagonalize_pointer=None
            if optional("image_diagonalize_pointer") is None
            else int(optional("image_diagonalize_pointer")),
        )

    def wyckoff_rows(self, sg: int) -> list[WyckoffRow]:
        start = int(self.iso.wyckoff["iwyckoff_pointer"][sg - 1])
        count = int(self.iso.wyckoff["iwyckoff_count"][sg - 1])
        rows = []
        for offset0 in range(count):
            row_id = start + offset0
            rows.append(
                WyckoffRow(
                    offset0=offset0,
                    row_id=row_id,
                    label=self.iso.wyckoff["wyckoff_label"][row_id - 1].strip(),
                    site_pg=int(self.iso.wyckoff["iwyckoff_pg"][row_id - 1]),
                )
            )
        return rows

    def wyckoff_subductions(self, sg: int, label: str) -> list[WyckoffSubduction]:
        old_id = self.old_irrep_id(sg, label)
        return self.wyckoff_subductions_for_old(sg, old_id)

    def wyckoff_subductions_for_old(self, sg: int, old_id: int) -> list[WyckoffSubduction]:
        base = int(self.iso.wyckoff["iwyckoff_irrep_pointer"][old_id - 1])
        out = []
        for row in self.wyckoff_rows(sg):
            entry = base + row.offset0
            pointer = int(self.iso.wyckoff["iwyckoff_subduce_pointer"][entry - 1])
            count = int(self.iso.wyckoff["iwyckoff_subduce_count"][entry - 1])
            pairs = tuple(
                (
                    int(self.iso.wyckoff["iwyckoff_subduce_pg_irrep"][pointer - 1 + j]),
                    int(self.iso.wyckoff["iwyckoff_subduce_frequency"][pointer - 1 + j]),
                )
                for j in range(count)
            )
            out.append(WyckoffSubduction(row, entry, pointer, count, pairs))
        return out

    def wyckoff_subductions_for_little_gid(self, sg: int, gid: int) -> list[WyckoffSubduction]:
        """Return Wyckoff subductions for a native little-irrep row.

        Fixed special-point irreps point back to ``data_irreps`` via
        ``little_irr_old`` and use the historical ``iwyckoff_*`` tables.
        Parametric k irreps have ``little_irr_old == 0``; the mode kernel then uses the
        parallel ``little_subduce_*`` tables keyed directly by the little row.
        """

        little = self.little_record_by_gid(gid)
        if little.old_id > 0:
            return self.wyckoff_subductions_for_old(sg, little.old_id)

        base = int(self.iso.little["little_subduce_irr_pointer"][gid - 1])
        if base <= 0:
            raise KeyError(f"missing little_subduce pointer for gid={gid}")
        out = []
        for row in self.wyckoff_rows(sg):
            entry = base + row.offset0
            pointer = int(self.iso.little["little_subduce_pointer"][entry - 1])
            count = int(self.iso.little["little_subduce_count"][entry - 1])
            pairs = tuple(
                (
                    int(self.iso.little["little_subduce_pg_irrep"][pointer - 1 + j]),
                    int(self.iso.little["little_subduce_frequency"][pointer - 1 + j]),
                )
                for j in range(count)
            )
            out.append(WyckoffSubduction(row, entry, pointer, count, pairs))
        return out

    def site_vector_reps(self, site_pg: int) -> tuple[int, ...]:
        values = self.iso.wyckoff["iwyckoff_pg_vector_reps"]
        return tuple(int(x) for x in values[(site_pg - 1) * 6:site_pg * 6])

    def site_vector_basis_codes(self, site_pg: int) -> tuple[int, ...]:
        values = self.iso.wyckoff["iwyckoff_pg_vector_basis"]
        return tuple(int(x) for x in values[(site_pg - 1) * 18:site_pg * 18])

    def site_pg_irrep_old_id(self, site_pg: int, pg_irrep: int) -> int:
        start = (site_pg - 1) * 12
        old_id = int(self.iso.wyckoff["iwyckoff_pg_irrep"][start + pg_irrep - 1])
        if old_id == 0:
            raise KeyError(f"PG{site_pg} irrep {pg_irrep} is not defined")
        return old_id

    def site_pg_irrep_label(self, site_pg: int, pg_irrep: int) -> str:
        start = (site_pg - 1) * 12
        return self.iso.wyckoff["wyckoff_pg_irrep_label"][start + pg_irrep - 1].strip()

    def site_irrep_matrices(
        self,
        site_pg: int,
        pg_irrep: int,
    ) -> dict[int, np.ndarray]:
        """Return site-irrep matrices keyed by canonical site point-op index."""

        old_id = self.site_pg_irrep_old_id(site_pg, pg_irrep)
        sg = int(self.iso.irreps["irrep_space_group"][old_id - 1])
        label = self.iso.irreps["irrep_label"][old_id - 1].strip()
        elements = self.iso.irrep_matrices_with_ops(sg, label)
        return {int(element.op_index): element.D for element in elements[3:] if element.op_index is not None}

    def little_gid_for_old_id(self, old_id: int) -> int:
        """Return the unique little-row id for one exact Source old-irrep id."""

        index = getattr(self, "_little_gids_by_old_id", None)
        if index is None:
            grouped: dict[int, list[int]] = {}
            for gid, value in enumerate(
                self.iso.little["little_irr_old"], start=1
            ):
                grouped.setdefault(int(value), []).append(gid)
            index = {key: tuple(values) for key, values in grouped.items()}
            self._little_gids_by_old_id = index
        matches = index.get(int(old_id), ())
        if len(matches) != 1:
            raise KeyError(
                f"expected one little irrep row for old id {old_id}, "
                f"got {list(matches)}"
            )
        return int(matches[0])

    @staticmethod
    def _pointed_slice(pointers: list[int], values: list[int], owner_index1: int) -> tuple[int, ...]:
        pointer = int(pointers[owner_index1 - 1])
        if pointer == 0:
            return tuple()
        later = sorted({int(p) for p in pointers if int(p) > pointer})
        end = later[0] if later else len(values) + 1
        return tuple(int(x) for x in values[pointer - 1:end - 1])
