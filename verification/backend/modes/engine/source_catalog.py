from __future__ import annotations

from pathlib import Path

from verification.backend.modes.engine.records import (
    LittleIrrepRecord,
    WyckoffRow,
)
from verification.backend.source.tables import SOURCE, SourceTables


class ModeSourceCatalogMixin:
    def __init__(
        self,
        data_dir: str | Path = SOURCE,
        *,
        tables: SourceTables | None = None,
    ):
        self.data_dir = Path(data_dir)
        if tables is None:
            tables = SourceTables(self.data_dir)
        else:
            self.data_dir = Path(tables.source)
        self.tables = tables
        self.iso = tables.iso
        self.images = tables.images

    def old_irrep_id(self, sg: int, label: str) -> int:
        return self.iso.old_irrep_id(sg, label)

    def little_record_by_gid(self, gid: int) -> LittleIrrepRecord:
        """Return little-irrep metadata by its `data_little` row id.

        Some full little-irrep labels are combined labels that are not
        directly present as old (data_irreps) labels. Callers that already
        have the little-row id therefore need not depend on label reverse
        lookup.
        """

        return self.tables.little_record_by_gid(gid)

    irrep_label_for_old_id = SourceTables.irrep_label_for_old_id
    generate_space_group_records = SourceTables.generate_space_group_records
    wyckoff_fraction_vectors = SourceTables.wyckoff_fraction_vectors
    vrot_fraction = SourceTables.vrot_fraction
    wyc_pg_elements_records = SourceTables.wyc_pg_elements_records

    space_group_has_inversion = SourceTables.space_group_has_inversion
    point_group_order = SourceTables.point_group_order
    space_group_point_group_order = SourceTables.space_group_point_group_order

    def wyckoff_rows(self, sg: int) -> list[WyckoffRow]:
        return self.tables.wyckoff_rows(sg)

    def site_vector_reps(self, site_pg: int) -> tuple[int, ...]:
        return self.tables.site_vector_reps(site_pg)

    def site_pg_irrep_old_id(self, site_pg: int, pg_irrep: int) -> int:
        return self.tables.site_pg_irrep_old_id(site_pg, pg_irrep)

    site_pg_irrep_label = SourceTables.site_pg_irrep_label

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
