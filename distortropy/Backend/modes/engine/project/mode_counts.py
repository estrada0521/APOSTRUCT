"""Mode-count helpers for the pre-``project_`` mode-kernel data path."""

from __future__ import annotations

from fractions import Fraction

from distortropy.Backend.modes.engine.input import Case, k_label_from_irrep
from distortropy.Backend.modes.engine.records import LittleIrrepRecord, WyckoffRow


def little_records_for_k(decoder, sg: int, k_label: str) -> list[LittleIrrepRecord]:
    out: list[LittleIrrepRecord] = []
    seen: set[int] = set()
    for gid, row_sg in enumerate(decoder.iso.little["little_irr_space_group"], start=1):
        if int(row_sg) != sg:
            continue
        old_id = int(decoder.iso.little["little_irr_old"][gid - 1])
        seen_key = old_id if old_id > 0 else -gid
        if seen_key in seen:
            continue
        record = decoder.little_record_by_gid(gid)
        if k_label_from_irrep(record.label) == k_label:
            seen.add(seen_key)
            out.append(record)
    return out


def mode_total_for_row(decoder, sg: int, row: WyckoffRow, little: LittleIrrepRecord) -> int:
    """Current decoded upstream mode-count layer for one Wyckoff/little-irrep row."""

    subductions = [
        item for item in decoder.wyckoff_subductions_for_little_gid(sg, little.gid)
        if item.wyckoff.row_id == row.row_id
    ]
    if len(subductions) != 1:
        raise KeyError(f"expected one subduction for SG{sg} row {row.row_id} old {little.old_id}")
    divisor = decoder.TYPE_DIVISOR[little.irrep_type]
    project_counts: dict[int, int] = {}
    site_dims: dict[int, int] = {}
    for pg_irrep, frequency in subductions[0].pairs:
        site_old = decoder.site_pg_irrep_old_id(row.site_pg, int(pg_irrep))
        site_image = decoder.image_record(site_old)
        site_dims[int(pg_irrep)] = int(site_image.dimension)
        factor = 1 if site_image.image_type == 1 else 2
        numerator = int(frequency) * factor
        if numerator % divisor:
            raise ValueError(
                f"nonintegral count for SG{sg} {row.label} {little.label}: "
                f"frequency={frequency} factor={factor} divisor={divisor}"
            )
        project_counts[int(pg_irrep)] = numerator // divisor

    vector_slots = decoder.site_vector_reps(row.site_pg)[:3]
    project_total = Fraction(0, 1)
    for pg_irrep in vector_slots:
        pg_irrep = int(pg_irrep)
        if pg_irrep in project_counts:
            project_total += Fraction(project_counts[pg_irrep], site_dims[pg_irrep])
    total = project_total * little.full_dim
    if total.denominator != 1:
        raise ValueError(f"nonintegral vector-projected total for SG{sg} {row.label} {little.label}: {total}")
    return int(total)


def mode_totals(decoder, case: Case) -> dict[str, int]:
    """Return mode totals for one atom-displacement input."""
    rows = [row for row in decoder.wyckoff_rows(case.sg) if row.label == case.wyckoff]
    if len(rows) != 1:
        raise KeyError(f"expected one Wyckoff row for SG{case.sg} {case.wyckoff}, got {len(rows)}")
    row = rows[0]
    expected: dict[str, int] = {}
    for little in little_records_for_k(decoder, case.sg, case.k_label):
        total = mode_total_for_row(decoder, case.sg, row, little)
        if total > 0:
            expected[little.label] = total
    return expected


def project_counts_for_little(
    decoder,
    sg: int,
    row_id: int,
    little: LittleIrrepRecord,
) -> dict[int, int]:
    """Return project counts keyed by site PG irrep for a native little row."""

    rows = [
        item for item in decoder.wyckoff_subductions_for_little_gid(sg, little.gid)
        if item.wyckoff.row_id == row_id
    ]
    if len(rows) != 1:
        raise KeyError(f"expected one subduction for SG{sg} row_id={row_id} gid={little.gid}")
    divisor = decoder.TYPE_DIVISOR[little.irrep_type]
    out: dict[int, int] = {}
    for pg_irrep, frequency in rows[0].pairs:
        site_old_id = decoder.site_pg_irrep_old_id(rows[0].wyckoff.site_pg, int(pg_irrep))
        site_image = decoder.image_record(site_old_id)
        factor = 1 if site_image.image_type == 1 else 2
        numerator = int(frequency) * factor
        if numerator % divisor:
            raise ValueError(
                f"nonintegral project count: SG{sg} gid={little.gid} "
                f"pg={pg_irrep} frequency={frequency} factor={factor} divisor={divisor}"
            )
        out[int(pg_irrep)] = numerator // divisor
    return out
