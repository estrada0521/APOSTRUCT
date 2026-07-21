"""Enumerate order-parameter directions and select an isotropy row."""

from __future__ import annotations

import re
from typing import Any

from ISODISTORT.Assembled.Backend.isotropy import catalog as isotropy_catalog
from ISODISTORT.Assembled.Backend.parent import parent_cell_tuple
from ISODISTORT.Assembled.Backend.reciprocal.k_points import normalize_k_params
from ISODISTORT.Assembled.Backend.reciprocal.irreps import select_irrep
from ISODISTORT.Assembled.Backend.isotropy.coupled import coupled_opd_rows
from ISODISTORT.Assembled.Backend.isotropy.state import build_subgroup_state
from ISODISTORT.Assembled.Backend.source.tables import source_tables


_COUPLED_OPD_PART_RE = re.compile(
    r"(?P<label>(?:P|C|S|(?:[4-9]|[1-9][0-9])D)[1-9][0-9]*)"
    r"\((?P<domain>[1-9][0-9]*)\)"
)


def _coupled_selected_labels(
    orderparam: int | str | None,
    *,
    slot_count: int,
    selected_opd_only: bool,
) -> tuple[str, ...] | None:
    """Parse exact per-slot Source labels from a complete OPD heading.

    B: the heading is a conditioned input.  It only bounds Source OPD
    discovery; no Web row values are accepted as computed output.
    """

    if not selected_opd_only or not isinstance(orderparam, str) or slot_count < 2:
        return None
    text = orderparam.strip()
    parts = tuple(_COUPLED_OPD_PART_RE.finditer(text))
    if len(parts) != int(slot_count):
        return None
    if "".join(part.group(0) for part in parts) != text:
        return None
    return tuple(part.group("label") for part in parts)


def _rows_for_selected_label(
    rows: list[dict[str, Any]],
    label: str | None,
) -> list[dict[str, Any]]:
    if label is None:
        return rows
    return [
        row
        for row in rows
        if label
        in {
            str((row.get("isotropy") or {}).get("opd_label") or "").strip(),
            str((row.get("direction") or {}).get("label") or "").strip(),
        }
    ]


def _select_orderparam(
    opd_rows: list[dict[str, Any]],
    *,
    orderparam: int | str | None = None,
    orderparam_index: int = 1,
) -> dict[str, Any] | None:
    if isinstance(orderparam, str) and orderparam.strip():
        label = orderparam.strip()
        selected = next(
            (
                row
                for row in opd_rows
                if label
                in {
                    str((row.get("isotropy") or {}).get("opd_label") or "").strip(),
                    str((row.get("direction") or {}).get("label") or "").strip(),
                }
            ),
            None,
        )
        if selected is not None:
            return selected

    orderparam_row_id = None
    if isinstance(orderparam, int) and not isinstance(orderparam, bool):
        orderparam_row_id = int(orderparam)
    elif isinstance(orderparam, str) and re.fullmatch(r"\s*\d+\s*", orderparam):
        orderparam_row_id = int(orderparam)
    if orderparam_row_id is not None:
        selected = next(
            (
                row
                for row in opd_rows
                if int(((row.get("isotropy") or {}).get("row_id") or 0)) == orderparam_row_id
            ),
            None,
        )
        if selected is not None:
            return selected
    if 1 <= int(orderparam_index) <= len(opd_rows):
        return opd_rows[int(orderparam_index) - 1]
    return None


def select_order_parameter_direction(
    reciprocal_state: dict[str, Any],
    *,
    k_params: dict[str, str] | None = None,
    display_k_params: dict[str, str] | None = None,
    irrep: str | None = None,
    irrep_index: int = 1,
    orderparam: int | str | None = None,
    orderparam_index: int = 1,
    selections: list[dict[str, Any]] | None = None,
    opd_row_limit: int | None = None,
    selected_opd_only: bool = False,
) -> dict[str, Any]:
    sg = int(reciprocal_state["input"]["parent"]["number"])
    if selections and len(selections) > 1:
        return _select_coupled_order_parameter(
            reciprocal_state,
            selections=selections,
            orderparam=orderparam,
            orderparam_index=orderparam_index,
            opd_row_limit=opd_row_limit,
            selected_opd_only=selected_opd_only,
        )
    k_params = normalize_k_params(k_params)
    display_k_params = normalize_k_params(display_k_params) or k_params
    local_data = source_tables()
    selection = (selections or [{}])[0]
    selected_irrep = select_irrep(
        reciprocal_state["selected"]["irreps"],
        selection.get("irrep", irrep),
        int(selection.get("irrep_index") or irrep_index),
    )
    old_id = int(selected_irrep["old_id"])
    if selected_irrep.get("magnetic"):
        opd_rows = isotropy_catalog.magnetic_opd_rows(
            old_id,
            parent_sg=sg,
            gid=int(selected_irrep["gid"]),
            k_params=k_params,
            display_k_params=display_k_params,
            selected_orderparam=orderparam if selected_opd_only else None,
            parent_setting_id=int(reciprocal_state["_internal"]["parent_inter_setting_id"]),
            parent_cell=parent_cell_tuple(reciprocal_state["input"]),
        )
    else:
        opd_rows = isotropy_catalog.opd_rows(
            old_id,
            parent_sg=sg,
            gid=int(selected_irrep["gid"]),
            k_params=k_params,
            display_k_params=display_k_params,
            selected_orderparam=orderparam if selected_opd_only else None,
            parent_setting_id=int(reciprocal_state["_internal"]["parent_inter_setting_id"]),
            parent_cell=parent_cell_tuple(reciprocal_state["input"]),
        )
    selected_opd = _select_orderparam(
        opd_rows,
        orderparam=orderparam,
        orderparam_index=orderparam_index,
    )
    if selected_irrep.get("magnetic"):
        subgroup_state = {
            "status": "unsupported",
            "reason": "magnetic OPD subgroup_state is not wired yet",
        }
    else:
        subgroup_state = build_subgroup_state(
            data=local_data,
            parent_sg=sg,
            selected_irrep=selected_irrep,
            selected_opd=selected_opd,
            displacive_row_ids=reciprocal_state["_internal"]["displacive_row_ids"],
        )
    return {
        **reciprocal_state,
        "schema": "isodistort.assembled.opd_direction.v1",
        "selected": {
            **reciprocal_state["selected"],
            "irrep": selected_irrep,
            "selection_slots": [
                {
                    "slot": 1,
                    "kpoint": reciprocal_state["selected"].get("kpoint"),
                    "irrep": selected_irrep,
                    "k_params": k_params,
                    "display_k_params": display_k_params,
                }
            ],
            "image": None if old_id <= 0 else local_data.image_record(old_id),
            "opd_rows": opd_rows,
            "orderparam": selected_opd,
            "subgroup_state": subgroup_state,
        },
    }


def _select_coupled_order_parameter(
    reciprocal_state: dict[str, Any],
    *,
    selections: list[dict[str, Any]],
    orderparam: int | str | None,
    orderparam_index: int,
    opd_row_limit: int | None,
    selected_opd_only: bool,
) -> dict[str, Any]:
    sg = int(reciprocal_state["input"]["parent"]["number"])
    data = source_tables()
    irrep_slots = list((reciprocal_state.get("selected") or {}).get("irrep_slots") or [])
    if len(irrep_slots) != len(selections):
        raise ValueError("K/IR slot count changed before coupled OPD generation")
    selected_labels = _coupled_selected_labels(
        orderparam,
        slot_count=len(selections),
        selected_opd_only=selected_opd_only,
    )
    slots: list[dict[str, Any]] = []
    for slot_index, (slot_state, selection) in enumerate(zip(irrep_slots, selections), start=1):
        selected_irrep = select_irrep(
            slot_state["irreps"],
            selection.get("irrep"),
            int(selection.get("irrep_index") or 1),
        )
        params = normalize_k_params(selection.get("k_params"))
        display_params = normalize_k_params(selection.get("display_k_params")) or params
        old_id = int(selected_irrep["old_id"])
        selected_label = (
            selected_labels[slot_index - 1]
            if selected_labels is not None
            else None
        )
        if selected_irrep.get("magnetic"):
            rows = isotropy_catalog.magnetic_opd_rows(
                old_id,
                parent_sg=sg,
                gid=int(selected_irrep["gid"]),
                k_params=params,
                display_k_params=display_params,
                selected_orderparam=selected_label,
                parent_setting_id=int(reciprocal_state["_internal"]["parent_inter_setting_id"]),
                parent_cell=parent_cell_tuple(reciprocal_state["input"]),
            )
        else:
            rows = isotropy_catalog.opd_rows(
                old_id,
                parent_sg=sg,
                gid=int(selected_irrep["gid"]),
                k_params=params,
                display_k_params=display_params,
                selected_orderparam=selected_label,
                parent_setting_id=int(reciprocal_state["_internal"]["parent_inter_setting_id"]),
                parent_cell=parent_cell_tuple(reciprocal_state["input"]),
            )
        rows = _rows_for_selected_label(rows, selected_label)
        source_kparam = isotropy_catalog._dynamic_kparam(  # noqa: SLF001
            data,
            gid=int(selected_irrep["gid"]),
            k_params=params,
        ) or (0, 0, 0, 1)
        slots.append({
            "slot": slot_index,
            "kpoint": slot_state["kpoint"],
            "irrep": selected_irrep,
            "k_params": params,
            "display_k_params": display_params,
            "source_kparam": source_kparam,
            "opd_rows": rows,
        })
    rows = coupled_opd_rows(
        isotropy_catalog._subgroup_core_data(),  # noqa: SLF001
        display_data=data,
        parent_sg=sg,
        slots=slots,
        parent_setting_id=int(reciprocal_state["_internal"]["parent_inter_setting_id"]),
        parent_cell=parent_cell_tuple(reciprocal_state["input"]),
        row_limit=opd_row_limit,
    )
    selected_opd = _select_orderparam(rows, orderparam=orderparam, orderparam_index=orderparam_index)
    selected_iso = (selected_opd or {}).get("isotropy") or {}
    subgroup_state = {
        "status": "ok" if selected_opd else "missing",
        "coupled": True,
        "subgroup": selected_iso.get("subgroup"),
        "basis": selected_iso.get("basis"),
        "origin": selected_iso.get("origin"),
        "primaries": [
            {
                "slot": slot["slot"],
                "old_id": int(slot["irrep"].get("old_id") or 0),
                "label": slot["irrep"].get("symbol"),
                "k_label": slot["kpoint"].get("label"),
            }
            for slot in slots
        ],
        "sources": [],
    }
    return {
        **reciprocal_state,
        "schema": "isodistort.assembled.opd_direction.multi.v1",
        "selected": {
            **reciprocal_state["selected"],
            "irrep": slots[0]["irrep"],
            "irreps_selected": [slot["irrep"] for slot in slots],
            "selection_slots": slots,
            "opd_rows": rows,
            "orderparam": selected_opd,
            "subgroup_state": subgroup_state,
        },
    }
