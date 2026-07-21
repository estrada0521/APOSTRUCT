"""Enumerate and select reciprocal-space k vectors."""

from __future__ import annotations

from typing import Any

from ISODISTORT.Assembled.Backend.reciprocal.catalog import kpoints


def normalize_k_params(value: Any) -> dict[str, str]:
    if not value:
        return {}
    if isinstance(value, dict):
        return {str(key): str(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        out: dict[str, str] = {}
        for item in value:
            text = str(item)
            if "=" not in text:
                continue
            key, param_value = text.split("=", 1)
            out[key.strip()] = param_value.strip()
        return out
    return {}


def select_k(kpoints: list[dict[str, Any]], label: str | None, index: int) -> dict[str, Any]:
    if label:
        wanted = label.casefold()
        for kpoint in kpoints:
            if str(kpoint["label"]).casefold() == wanted:
                return kpoint
        raise KeyError(f"k label {label!r} not found")
    if index < 1 or index > len(kpoints):
        raise IndexError(f"k index {index} out of range 1..{len(kpoints)}")
    return kpoints[index - 1]


def select_k_points(
    parent_state: dict[str, Any],
    *,
    k_label: str | None = "GM",
    k_index: int = 1,
    k_params: dict[str, str] | None = None,
    display_k_params: dict[str, str] | None = None,
    selections: list[dict[str, Any]] | None = None,
) -> dict[str, Any]:
    sg = int(parent_state["input"]["parent"]["number"])
    k_params = normalize_k_params(k_params)
    display_k_params = normalize_k_params(display_k_params) or k_params
    k_catalog = kpoints(sg)
    selection_specs = list(selections or [{
        "k_label": k_label,
        "k_index": k_index,
        "k_params": k_params,
        "display_k_params": display_k_params,
    }])
    selected_kpoints: list[dict[str, Any]] = []
    for slot, spec in enumerate(selection_specs, start=1):
        selected = select_k(
            k_catalog["kpoints"],
            spec.get("k_label", k_label),
            int(spec.get("k_index") or k_index),
        )
        slot_params = normalize_k_params(spec.get("k_params"))
        slot_display_params = normalize_k_params(spec.get("display_k_params")) or slot_params
        selected_kpoints.append({
            "slot": slot,
            "kslot": selected["kslot"],
            "kid": selected["kid"],
            "label": selected["label"],
            "display": selected["display"],
            "kvector": selected["kvector"],
            "isodistort_kvector": selected["isodistort_kvector"],
            "parameters": slot_params,
            "display_parameters": slot_display_params,
            "dimension": selected["dimension"],
            "star": selected["star"],
            "little_order": selected["little_order"],
        })
    selected_k = selected_kpoints[0]
    return {
        **parent_state,
        "schema": "isodistort.assembled.k_points.v1",
        "space_group": k_catalog["space_group"],
        "kpoints": [
            {
                "kslot": k["kslot"],
                "kid": k["kid"],
                "label": k["label"],
                "display": k["display"],
                "kvector": k["kvector"],
                "isodistort_kvector": k["isodistort_kvector"],
                "dimension": k["dimension"],
                "star_size": k["star_size"],
                "little_order": k["little_order"],
                "n_irreps": k["n_irreps"],
            }
            for k in k_catalog["kpoints"]
        ],
        "selected": {
            "kpoint": {key: value for key, value in selected_k.items() if key != "slot"},
            "kpoints": selected_kpoints,
        },
        "_internal": {
            **(parent_state.get("_internal") or {}),
            "k_catalog": k_catalog,
        },
    }
