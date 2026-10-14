"""Mode-count helpers for the pre-``project_`` mode-kernel data path."""

from __future__ import annotations

from verification.backend.modes.engine.input import k_label_from_irrep
from verification.backend.modes.engine.records import LittleIrrepRecord


def little_records_for_k(decoder, sg: int, k_label: str) -> list[LittleIrrepRecord]:
    memo = decoder.__dict__.setdefault("_little_records_for_k", {})
    key = (int(sg), str(k_label))
    if key not in memo:
        memo[key] = tuple(_little_records_for_k(decoder, *key))
    return list(memo[key])


def _little_records_for_k(decoder, sg: int, k_label: str) -> list[LittleIrrepRecord]:
    by_sg = decoder.__dict__.get("_little_gids_by_sg")
    if by_sg is None:
        by_sg = {}
        for gid, row_sg in enumerate(decoder.iso.little["little_irr_space_group"], start=1):
            by_sg.setdefault(int(row_sg), []).append(gid)
        decoder.__dict__["_little_gids_by_sg"] = by_sg
    out: list[LittleIrrepRecord] = []
    seen: set[int] = set()
    for gid in by_sg.get(sg, ()):
        old_id = int(decoder.iso.little["little_irr_old"][gid - 1])
        seen_key = old_id if old_id > 0 else -gid
        if seen_key in seen:
            continue
        record = decoder.little_record_by_gid(gid)
        if k_label_from_irrep(record.label) == k_label:
            seen.add(seen_key)
            out.append(record)
    return out

