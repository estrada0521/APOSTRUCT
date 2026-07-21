"""Enumerate irreps for the selected reciprocal-space k vector."""

from __future__ import annotations

import copy
from typing import Any

from ISODISTORT.Assembled.Backend.reciprocal.catalog import kpoints
from ISODISTORT.Assembled.Backend.source.magnetic import data as magnetic_data
from ISODISTORT.Assembled.Backend.source.tables import SourceTables, source_tables


def select_irrep(irreps: list[dict[str, Any]], label: str | None, index: int) -> dict[str, Any]:
    if label:
        wanted = label.casefold()
        for irrep in irreps:
            if str(irrep["symbol"]).casefold() == wanted:
                return irrep
        raise KeyError(f"irrep {label!r} not found")
    if index < 1 or index > len(irreps):
        raise IndexError(f"irrep index {index} out of range 1..{len(irreps)}")
    return irreps[index - 1]


def _displacive_visible(
    *,
    local_data: SourceTables,
    sg: int,
    displacive_row_ids: list[int],
    candidate: dict[str, Any],
) -> int:
    return sum(
        local_data.displacive_mode_total_for_row_id(sg, row_id, int(candidate["gid"]))
        for row_id in displacive_row_ids
    )


def _strain_visible(local_data: SourceTables, kpoint: dict[str, Any], candidate: dict[str, Any]) -> bool:
    """Return whether a little irrep can carry macroscopic strain.

    This is intentionally limited to IR visibility.  Strain OPD row generation
    still uses the ordinary isotropy machinery downstream.
    """

    if str(kpoint.get("label")) != "GM":
        return False
    label = str(candidate.get("symbol") or "")
    if label in {"GM1", "GM1+"}:
        return True
    old_id = int(candidate.get("old_id") or 0)
    if old_id <= 0:
        return False
    for row_id in local_data.isotropy_row_ids_for_old_irrep(old_id):
        if int(local_data.isotropy["isotropy_ferroelastic"][row_id - 1]) == 1:
            return True
    return False


def _apply_order_parameter_visibility_filter(
    *,
    local_data: SourceTables,
    sg: int,
    k_catalog: dict[str, Any],
    displacive_row_ids: list[int],
    include_strain: bool,
) -> None:
    for kpoint in k_catalog["kpoints"]:
        kpoint["_ordinary_irreps_unfiltered"] = list(kpoint["irreps"])
        filtered_irreps = []
        for candidate in kpoint["irreps"]:
            total = _displacive_visible(
                local_data=local_data,
                sg=sg,
                displacive_row_ids=displacive_row_ids,
                candidate=candidate,
            )
            strain = include_strain and _strain_visible(local_data, kpoint, candidate)
            if total > 0 or strain:
                item = dict(candidate)
                if total > 0:
                    item["displacive_mode_total"] = total
                if strain:
                    item["strain_visible"] = True
                filtered_irreps.append(item)
        kpoint["irreps"] = filtered_irreps
        kpoint["n_irreps"] = len(filtered_irreps)


def _magnetic_variant(candidate: dict[str, Any]) -> dict[str, Any]:
    item = dict(candidate)
    item["symbol"] = "m" + str(candidate["symbol"])
    kov_label = str(candidate.get("kov_label") or "")
    item["kov_label"] = "m" + kov_label if kov_label and not kov_label.startswith("m") else kov_label
    item["magnetic"] = True
    item["ordinary_symbol"] = candidate.get("symbol")
    item["order_parameter_type"] = "magnetic"
    return item


def _has_magnetic_variant(candidate: dict[str, Any]) -> bool:
    old_id = int(candidate.get("old_id") or 0)
    if old_id <= 0:
        # Parametric/composite rows such as LD1LE1 have no single old_id, but
        # Web still exposes the matching m-prefixed option in all-mode.
        return True
    return magnetic_data().irrep_has_magnetic_isotropy(old_id)


def _magnetic_visible(
    *,
    local_data: SourceTables,
    sg: int,
    displacive_row_ids: list[int],
    candidate: dict[str, Any],
) -> int:
    return sum(
        local_data.magnetic_mode_total_for_row_id(sg, row_id, int(candidate["gid"]))
        for row_id in displacive_row_ids
    )


def _append_magnetic_irreps(
    *,
    local_data: SourceTables,
    sg: int,
    k_catalog: dict[str, Any],
    displacive_row_ids: list[int],
) -> None:
    for kpoint in k_catalog["kpoints"]:
        ordinary = list(kpoint["irreps"])
        magnetic_source = list(kpoint.get("_ordinary_irreps_unfiltered") or ordinary)
        magnetic = []
        for candidate in magnetic_source:
            if not _has_magnetic_variant(candidate):
                continue
            total = _magnetic_visible(
                local_data=local_data,
                sg=sg,
                displacive_row_ids=displacive_row_ids,
                candidate=candidate,
            )
            if total <= 0:
                continue
            item = _magnetic_variant(candidate)
            item["magnetic_mode_total"] = total
            magnetic.append(item)
        if not magnetic:
            continue
        kpoint["irreps"] = ordinary + magnetic
        kpoint["n_irreps"] = len(kpoint["irreps"])


def select_irreps(reciprocal_state: dict[str, Any]) -> dict[str, Any]:
    sg = int(reciprocal_state["input"]["parent"]["number"])
    local_data = source_tables()
    k_catalog = copy.deepcopy((reciprocal_state.get("_internal") or {}).get("k_catalog"))
    if not k_catalog:
        k_catalog = kpoints(sg)
    mode_set = str((reciprocal_state.get("_internal") or {}).get("mode_set") or reciprocal_state.get("mode_set") or "displacive")
    internal = reciprocal_state.get("_internal") or {}
    selection = internal.get("distortion_selection") or {}
    displacive_row_ids = list(internal.get("displacive_row_ids") or [])
    magnetic_row_ids = list(internal.get("magnetic_row_ids") or [])
    include_strain = bool(selection.get("strain", mode_set in {"strain_displacive", "all"}))
    _apply_order_parameter_visibility_filter(
        local_data=local_data,
        sg=sg,
        k_catalog=k_catalog,
        displacive_row_ids=displacive_row_ids,
        include_strain=include_strain,
    )
    if magnetic_row_ids:
        _append_magnetic_irreps(
            local_data=local_data,
            sg=sg,
            k_catalog=k_catalog,
            displacive_row_ids=magnetic_row_ids,
        )
    selected_kpoints = list((reciprocal_state.get("selected") or {}).get("kpoints") or [])
    if not selected_kpoints:
        selected_kpoints = [{"slot": 1, **((reciprocal_state.get("selected") or {}).get("kpoint") or {})}]
    irrep_slots: list[dict[str, Any]] = []
    for slot, selected_kpoint in enumerate(selected_kpoints, start=1):
        selected_k = next(
            (k for k in k_catalog["kpoints"] if str(k["label"]) == str(selected_kpoint.get("label"))),
            k_catalog["kpoints"][0],
        )
        irrep_slots.append({
            "slot": int(selected_kpoint.get("slot") or slot),
            "kpoint": selected_kpoint,
            "irreps": selected_k["irreps"],
        })
    selected_k = next(
        (k for k in k_catalog["kpoints"] if str(k["label"]) == str(selected_kpoints[0].get("label"))),
        k_catalog["kpoints"][0],
    )
    return {
        **reciprocal_state,
        "schema": "isodistort.assembled.irreps.v1",
        "mode_set": mode_set,
        "irrep_catalog": [
            {
                "kslot": k["kslot"],
                "kid": k["kid"],
                "k_label": k["label"],
                "kvector": k["kvector"],
                "isodistort_kvector": k["isodistort_kvector"],
                "dimension": k["dimension"],
                "star_size": k["star_size"],
                "little_order": k["little_order"],
                "irreps": k["irreps"],
            }
            for k in k_catalog["kpoints"]
        ],
        "selected": {
            **(reciprocal_state.get("selected") or {}),
            "irreps": selected_k["irreps"],
            "irrep_slots": irrep_slots,
        },
        "_internal": {
            **(reciprocal_state.get("_internal") or {}),
            "k_catalog": k_catalog,
            "mode_set": mode_set,
        },
    }
