"""Set up parent and site operation records for mode projection."""

from __future__ import annotations

from typing import Iterable

from verification.backend.modes.engine.records import LittleIrrepRecord, WyckoffRow


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
