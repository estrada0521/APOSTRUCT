"""Set up parent and site operation records for mode projection."""

from __future__ import annotations

from typing import Iterable

from ISODISTORT.Assembled.Backend.modes.engine.records import LittleIrrepRecord, WyckoffRow


class ProjectSetupMixin:
    def project_parent_operation_records_for_little_row(
        self,
        sg: int,
        little: LittleIrrepRecord,
        row: WyckoffRow,
        site_params: Iterable[object] | None = None,
    ) -> tuple[tuple[int, int, int, int, int], ...]:
        """Return full parent operation records in `project_` setting order."""

        _ = little, site_params
        return self.wyc_pg_elements_records(sg, row)

    def map_parent_ops_to_site_project_records(
        self,
        site_pg: int,
        parent_point_ops: Iterable[int],
    ) -> tuple[tuple[int, int, int, int, int], ...]:
        """Map parent-space operations to site-irrep operations.

        This mirrors the `project_` setting loop.  It tries each candidate
        setting from `iwyckoff_pg_elements` until all setting point operations
        are found in the parent operation list.  The returned operation records
        keep the parent operation order but replace each matched parent point
        op by the operation in the canonical first site setting at the same
        position.

        SG205/R/Wyckoff-a example:

        parent `[1,9,5,25,33,29]` matches setting 2 and maps to
        `[49,51,53,61,63,65]`.
        """

        parent = [int(op) for op in parent_point_ops]
        settings = self.site_pg_element_settings(site_pg)
        if not settings:
            raise ValueError(f"site point group {site_pg} has no element settings")
        canonical = settings[0]
        for setting in settings:
            if len(setting) != len(canonical):
                continue
            mapped: list[tuple[int, int, int, int, int] | None] = [None] * len(parent)
            ok = True
            for setting_op, canonical_op in zip(setting, canonical):
                try:
                    parent_index = parent.index(setting_op)
                except ValueError:
                    ok = False
                    break
                mapped[parent_index] = (0, 0, 0, 1, canonical_op)
            if ok:
                return tuple(record for record in mapped if record is not None)
        raise ValueError(f"no site setting for PG{site_pg} matches parent ops {parent}")

    def vector_basis_id_for_site_records(
        self,
        site_pg: int,
        site_records: Iterable[tuple[int, int, int, int, int]],
    ) -> int:
        """Return the vector-basis row id used by project-vector expansion.

        The row id is aligned with the site point-group setting order in
        `iwyckoff_pg_elements`.  The projector uses it to index
        `iwyckoff_pg_vector_basis` in 18-int blocks.
        """

        ops = tuple(int(record[4]) for record in site_records)
        pointer = int(self.iso.wyckoff["iwyckoff_pg_elements_pointer"][int(site_pg) - 1])
        for index, setting in enumerate(self.site_pg_element_settings(site_pg)):
            if tuple(setting) == ops:
                return pointer + index
        for index, setting in enumerate(self.site_pg_element_settings(site_pg)):
            if set(setting) == set(ops):
                return pointer + index
        raise ValueError(f"no vector basis setting for PG{site_pg} ops={ops}")
