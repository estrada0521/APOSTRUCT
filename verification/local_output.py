"""The comparator's view of an apo output: the mode details with every vector definition expanded into per-atom rows
(the compact form keeps one atom table and one vector list per definition)."""

from __future__ import annotations

from typing import Any, Mapping

VECTOR_DEFINITION_FIELDS = ("displacive_definitions", "magnetic_definitions")


def definition_rows(mode_details: Mapping[str, Any], definition: Mapping[str, Any]) -> list[dict[str, Any]]:
    """The per-atom rows of one definition, from either form."""

    if "rows" in definition:
        return list(definition["rows"] or [])
    atoms = mode_details.get("atoms") or []
    block = (mode_details.get("atom_blocks") or [])[int(definition["atom_block"])]
    return [
        {
            "atom": atoms[index]["atom"],
            "atom_id": atoms[index]["atom_id"],
            "child_site": atoms[index]["child_site"],
            "wyckoff_site": atoms[index]["wyckoff_site"],
            "xyz": list(atoms[index]["xyz"]),
            "dxyz": list(vector),
        }
        for index, vector in zip(block, definition["vectors"], strict=True)
    ]


def with_definition_rows(mode_details: Mapping[str, Any]) -> dict[str, Any]:
    """A shallow copy whose vector definitions carry their per-atom rows."""

    out = dict(mode_details)
    for field in VECTOR_DEFINITION_FIELDS:
        out[field] = [
            {**definition, "rows": definition_rows(mode_details, definition)} if isinstance(definition, Mapping) else definition
            for definition in mode_details.get(field) or []
        ]
    return out


def comparison_preview(preview: Any) -> Any:
    """The preview the comparator reads: the selected mode details with their per-atom rows."""

    if not isinstance(preview, dict):
        return preview
    selected = preview.get("selected")
    if not isinstance(selected, dict) or not isinstance(selected.get("mode_details"), dict):
        return preview
    return {**preview, "selected": {**selected, "mode_details": with_definition_rows(selected["mode_details"])}}


def request_selections(input_payload: Mapping[str, Any]) -> list[dict[str, Any]]:
    """The k selections of a case in the comparator's form."""

    return [
        {"k_label": slot["label"], "k_params": dict(slot.get("params") or {}), "irrep": slot["ir"]}
        for slot in input_payload["k"]
    ]


__all__ = ["VECTOR_DEFINITION_FIELDS", "comparison_preview", "definition_rows", "request_selections", "with_definition_rows"]
