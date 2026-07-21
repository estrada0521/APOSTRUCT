"""Subgroup-first state for ISODISTORT displacive Method 2.

This module is intentionally upstream of any mode-vector synthesis.  The
selected OPD fixes an isotropy subgroup; that subgroup is the state from which
complete mode details must be generated.
"""

from __future__ import annotations

import re
from typing import Any

from ISODISTORT.Assembled.Backend.source.tables import SourceTables


def _isotropy(selected_opd: dict[str, Any] | None) -> dict[str, Any] | None:
    if not isinstance(selected_opd, dict):
        return None
    value = selected_opd.get("isotropy")
    return value if isinstance(value, dict) else None


def _k_label_from_irrep_label(label: str) -> str:
    match = re.match(r"^[A-Za-z]+", str(label).strip())
    return match.group(0) if match else str(label).strip()


def _as_int(value: Any) -> int | None:
    try:
        return int(value)
    except (TypeError, ValueError):
        return None


def _little_gids_for_old_id(data: SourceTables, parent_sg: int, old_id: int) -> tuple[int, ...]:
    return tuple(
        gid
        for gid, value in enumerate(data.little["little_irr_old"], start=1)
        if int(value) == int(old_id)
        and int(data.little["little_irr_space_group"][gid - 1]) == int(parent_sg)
    )


def _displacive_mode_total(
    data: SourceTables,
    parent_sg: int,
    old_id: int,
    displacive_row_ids: tuple[int, ...],
) -> int | None:
    if not displacive_row_ids:
        return None
    total = 0
    for gid in _little_gids_for_old_id(data, parent_sg, old_id):
        for row_id in displacive_row_ids:
            total += data.displacive_mode_total_for_row_id(parent_sg, row_id, gid)
    return total


def build_subgroup_state(
    *,
    data: SourceTables,
    parent_sg: int,
    selected_irrep: dict[str, Any],
    selected_opd: dict[str, Any] | None,
    displacive_row_ids: list[int] | tuple[int, ...] | None = None,
) -> dict[str, Any]:
    """Return the selected OPD's isotropy-subgroup state.

    The returned object is the intended handoff point for the future
    complete-mode-detail kernel:

    parent/IR/OPD -> isotropy subgroup -> primary/secondary sources.
    """

    iso = _isotropy(selected_opd)
    primary_old_id = _as_int(selected_irrep.get("old_id")) or 0
    primary_label = str(selected_irrep.get("symbol") or selected_irrep.get("label") or "")
    if not iso:
        return {
            "status": "missing",
            "reason": "selected OPD has no concrete isotropy row",
            "primary": {
                "old_id": primary_old_id,
                "label": primary_label,
                "k_label": _k_label_from_irrep_label(primary_label),
            },
            "sources": [],
        }

    row_id = _as_int(iso.get("row_id"))
    if row_id is None:
        return {
            "status": "missing",
            "reason": "selected OPD isotropy row has no row_id",
            "primary": {
                "old_id": primary_old_id,
                "label": primary_label,
                "k_label": _k_label_from_irrep_label(primary_label),
            },
            "isotropy": iso,
            "sources": [],
        }

    row_ids = tuple(sorted({int(row_id) for row_id in displacive_row_ids or []}))
    subgroup = iso.get("subgroup") if isinstance(iso.get("subgroup"), dict) else {}
    if iso.get("dynamic") or primary_old_id <= 0:
        return {
            "status": "ok",
            "row_id": row_id,
            "subgroup": subgroup,
            "basis": iso.get("basis"),
            "origin": iso.get("origin"),
            "arms": iso.get("arms"),
            "det": iso.get("det"),
            "direction": selected_opd.get("direction") if isinstance(selected_opd, dict) else None,
            "opd": {
                "label": iso.get("opd_label"),
                "display": iso.get("display_opd"),
                "stored": iso.get("stored_opd"),
                "free": iso.get("free"),
                "op_rows": iso.get("op_rows"),
            },
            "primary": {
                "old_id": primary_old_id,
                "label": primary_label,
                "k_label": _k_label_from_irrep_label(primary_label),
            },
            "sources": [],
            "displacive_row_ids": list(row_ids),
            "notes": [
                "Dynamic parametric-k OPD row; static data_isotropy row identity is not available.",
                "Static data_isotropy secondary-source rows are not available for this dynamic row_id.",
            ],
        }

    sources: list[dict[str, Any]] = []
    seen: set[tuple[int, int]] = set()
    for order, entry in enumerate(data.isotropy_subductions_for_row(row_id), start=1):
        if int(entry.sg) != int(parent_sg):
            continue
        mode_total = _displacive_mode_total(data, parent_sg, int(entry.irrep_old_id), row_ids)
        if mode_total is not None and mode_total <= 0:
            continue
        key = (int(entry.irrep_old_id), int(entry.subgroup_row_id))
        if key in seen:
            continue
        seen.add(key)
        is_primary = (
            int(entry.irrep_old_id) == int(primary_old_id)
            and int(entry.subgroup_row_id) == int(row_id)
        )
        sources.append(
            {
                "source_order": order,
                "role": "primary" if is_primary else "secondary",
                "old_id": int(entry.irrep_old_id),
                "label": str(entry.irrep_label),
                "k_label": _k_label_from_irrep_label(str(entry.irrep_label)),
                "frequency": int(entry.frequency),
                "isotropy_row_id": int(entry.subgroup_row_id),
                "domain": int(entry.domain),
                "domain_old": int(entry.domain_old),
                "displacive_mode_total": mode_total,
            }
        )

    if primary_old_id > 0 and not any(item["role"] == "primary" for item in sources):
        sources.insert(
            0,
            {
                "source_order": 0,
                "role": "primary",
                "old_id": primary_old_id,
                "label": primary_label,
                "k_label": _k_label_from_irrep_label(primary_label),
                "frequency": 1,
                "isotropy_row_id": row_id,
                "domain": None,
                "domain_old": None,
                "displacive_mode_total": _displacive_mode_total(data, parent_sg, primary_old_id, row_ids),
                "fallback": True,
            },
        )

    return {
        "status": "ok",
        "row_id": row_id,
        "subgroup": subgroup,
        "basis": iso.get("basis"),
        "origin": iso.get("origin"),
        "arms": iso.get("arms"),
        "det": iso.get("det"),
        "direction": selected_opd.get("direction") if isinstance(selected_opd, dict) else None,
        "opd": {
            "label": iso.get("opd_label"),
            "display": iso.get("display_opd"),
            "stored": iso.get("stored_opd"),
            "free": iso.get("free"),
            "op_rows": iso.get("op_rows"),
        },
        "primary": {
            "old_id": primary_old_id,
            "label": primary_label,
            "k_label": _k_label_from_irrep_label(primary_label),
        },
        "sources": sources,
        "displacive_row_ids": list(row_ids),
        "notes": [
            "This is the subgroup-first state. Mode vectors must be generated in this selected subgroup basis/origin.",
            "Sources come from data_isotropy:isotropy_subduce_* and are filtered by the selected displacive Wyckoff rows when those rows are known.",
        ],
    }
