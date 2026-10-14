from __future__ import annotations

from typing import Iterable


class WyckoffGeometryMixin:
    def wyc_pg_cosets_records(
        self,
        sg: int,
        site_records: Iterable[tuple[int, int, int, int, int]],
    ) -> tuple[tuple[int, int, int, int, int], ...]:
        records = tuple(tuple(int(value) for value in record) for record in site_records)
        return self.iso.wyc_pg_cosets_records(int(sg), records)
