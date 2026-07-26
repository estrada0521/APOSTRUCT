"""Trace helpers for site-irrep projection inputs."""

from __future__ import annotations

from ISODISTORT.Assembled.Backend.modes.engine.input import Case, k_label_from_irrep
from ISODISTORT.Assembled.Backend.modes.engine.records import WyckoffRow
from ISODISTORT.Assembled.Backend.modes.engine.project.mode_counts import project_counts_for_little


def project_entry_trace(
    decoder,
    case: Case,
    row: WyckoffRow,
    *,
    vector_setting: int = 1,
) -> list[dict[str, int]]:
    start = (int(vector_setting) - 1) * 3
    vector_slots = set(int(slot) for slot in decoder.site_vector_reps(row.site_pg)[start:start + 3])
    out: list[dict[str, int]] = []
    seen_rows: set[int] = set()
    for gid, sg in enumerate(decoder.iso.little["little_irr_space_group"], start=1):
        if int(sg) != case.sg:
            continue
        old_id = int(decoder.iso.little["little_irr_old"][gid - 1])
        seen_key = old_id if old_id > 0 else -gid
        if seen_key in seen_rows:
            continue
        little = decoder.little_record_by_gid(gid)
        if k_label_from_irrep(little.label) != case.k_label:
            continue
        seen_rows.add(seen_key)
        subductions = [
            item for item in decoder.wyckoff_subductions_for_little_gid(case.sg, little.gid)
            if item.wyckoff.row_id == row.row_id
        ]
        if len(subductions) != 1:
            raise KeyError(f"expected one subduction for {case} gid={gid}")
        for pg_irrep, _frequency in subductions[0].pairs:
            if int(pg_irrep) in vector_slots:
                out.append({
                    "gid": gid,
                    "pg_irrep": int(pg_irrep),
                    "param5_second": row.offset0 + 1,
                })
    return out


def project_local408_trace(
    decoder,
    case: Case,
    row: WyckoffRow,
    project_entries: list[dict[str, int]],
) -> list[dict[str, object]]:
    out: list[dict[str, object]] = []
    for entry in project_entries:
        little = decoder.little_record_by_gid(int(entry["gid"]))
        sub_counts = project_counts_for_little(decoder, case.sg, row.row_id, little)
        pg_irrep = int(entry["pg_irrep"])
        site_old_id = decoder.site_pg_irrep_old_id(row.site_pg, pg_irrep)
        local_records = decoder.site_pg_project_records(row.site_pg, 0)
        local_ops = tuple(int(record[4]) for record in local_records)
        out.append({
            "site_old_id": int(site_old_id),
            "project_count": int(sub_counts[pg_irrep]),
            "site_op_count": len(local_ops),
            "local_ops": list(local_ops),
        })
    return out


def site_get_irreps_trace(project_local408: list[dict[str, object]]) -> dict[str, object]:
    groups = [list(project["local_ops"]) for project in project_local408]
    records = [
        {"op_record": [0, 0, 0, 1, int(op)]}
        for group in groups
        for op in group
    ]
    return {"records": records, "groups": groups}
