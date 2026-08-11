"""Compact, stable JSON projections for command-line consumers."""

from __future__ import annotations

from collections.abc import Sequence
from typing import Any

from APOSTRUCT.case_input import parameter_names_for_kpoint


def _selected_fields(row: dict[str, Any], fields: Sequence[str]) -> dict[str, Any]:
    return {field: row[field] for field in fields if field in row}


def _kpoint_result(row: dict[str, Any]) -> dict[str, Any]:
    result = _selected_fields(
        row,
        (
            "slot",
            "label",
            "display",
            "kid",
            "kslot",
            "dimension",
            "parameters",
            "display_parameters",
            "miller_love_parameters",
            "sg_specific",
            "star_size",
            "little_order",
            "n_irreps",
            "star",
        ),
    )
    if row.get("display_kvector") is not None:
        result["kvector"] = row["display_kvector"]
    if row.get("kvector") is not None:
        result["miller_love_kvector"] = row["kvector"]
    result["parameter_names"] = list(parameter_names_for_kpoint(row))
    return result


def _irrep_result(row: dict[str, Any]) -> dict[str, Any]:
    result = _selected_fields(
        row,
        (
            "symbol",
            "gid",
            "k_label",
            "kov_label",
            "type",
            "type_label",
            "little_dim",
            "full_dim",
            "displacive_mode_total",
            "magnetic_mode_total",
        ),
    )
    result["strain_visible"] = bool(row.get("strain_visible"))
    return result


def _subgroup_result(raw: Any) -> dict[str, Any]:
    subgroup = dict(raw or {})
    magnetic = subgroup.get("ordinary_number") is not None
    ordinary_number = subgroup.get("ordinary_number")
    if ordinary_number is None:
        ordinary_number = subgroup.get("number")
    result = {
        "display_label": str(
            subgroup.get("display_label") or ordinary_number or ""
        ),
        "ordinary_number": ordinary_number,
        "symbol": subgroup.get("symbol"),
        "magnetic": magnetic,
    }
    if magnetic:
        result["magnetic_id"] = subgroup.get("number")
    return result


def compact_opd_result(row: dict[str, Any]) -> dict[str, Any]:
    direction = dict(row.get("direction") or {})
    isotropy = dict(row.get("isotropy") or {})
    result = {
        "label": direction.get("label"),
        "opd": direction.get("opd"),
        "subgroup": _subgroup_result(isotropy.get("subgroup")),
        "basis": isotropy.get("basis_text"),
        "origin": isotropy.get("origin"),
        "index": isotropy.get("i"),
        "cell_index": isotropy.get("s"),
        "parameter_count": isotropy.get("free"),
        "k_active": isotropy.get("k_active_vectors") or isotropy.get("k_active"),
        "ferroic_classified": bool(isotropy.get("ferroic_classified")),
        "ferroic_properties": isotropy.get("ferroic_properties") or [],
    }
    if direction.get("components") is not None:
        result["components"] = direction["components"]
    if isotropy.get("embedding_selected") is True:
        result["embedding_selected"] = True
        result["direction_row"] = isotropy.get("embedding_direction_row")
    return result


def _public_value(value: Any) -> Any:
    if isinstance(value, dict):
        return {
            key: _public_value(item)
            for key, item in value.items()
            if not str(key).startswith("_")
        }
    if isinstance(value, list):
        return [_public_value(item) for item in value]
    if isinstance(value, tuple):
        return [_public_value(item) for item in value]
    return value


def _attach_definition_ids(mode_details: dict[str, Any]) -> None:
    for kind, field in (
        ("displacive", "displacive_definitions"),
        ("magnetic", "magnetic_definitions"),
    ):
        for index, definition in enumerate(mode_details.get(field) or [], start=1):
            if isinstance(definition, dict):
                definition["definition_id"] = f"{kind}-{index}"


def _invariant_factor_result(row: dict[str, Any]) -> dict[str, Any]:
    return _public_value(_selected_fields(
        row,
        (
            "slot",
            "role",
            "gid",
            "label",
            "magnetic",
            "k_parameters",
            "source_k_parameters",
            "opd",
            "domain",
            "domain_display",
            "parameter_count",
            "parameter_offset",
            "domain_presentation",
            "direction_matrix",
        ),
    ))


def compact_result(value: dict[str, Any], *, stage: str) -> dict[str, Any]:
    """Project one pipeline state to the public data needed at ``stage``."""

    input_state = dict(value.get("input") or {})
    selected = dict(value.get("selected") or {})
    space_group = dict(value.get("space_group") or input_state.get("parent") or {})
    input_symbol = str((input_state.get("parent") or {}).get("symbol") or "").strip()
    if input_symbol:
        space_group["symbol"] = input_symbol
    common = {
        "schema": f"APOSTRUCT.cli.{stage}",
        "space_group": space_group,
    }
    if stage == "info":
        common.update({
            "formula": input_state.get("formula"),
            "cell": input_state.get("lattice"),
            "sites": [
                _selected_fields(
                    dict(site),
                    (
                        "type",
                        "label",
                        "wyckoff_multiplicity",
                        "wyckoff",
                        "wyckoff_formula",
                        "wyckoff_parameters",
                        "wyckoff_params",
                        "coordinates_realized",
                        "fract",
                        "occupancy",
                        "wyckoff_mapping_error",
                    ),
                )
                for site in input_state.get("atom_sites") or []
            ],
        })
        return common
    if stage == "kpoints":
        common["kpoints"] = [
            _kpoint_result(dict(row)) for row in value.get("kpoints") or []
        ]
        return common
    if stage == "irreps":
        common["distortion_selection"] = value.get("distortion_selection") or {}
        common["selections"] = [
            {
                "slot": slot.get("slot"),
                "kpoint": _kpoint_result(dict(slot.get("kpoint") or {})),
                "irreps": [
                    _irrep_result(dict(row)) for row in slot.get("irreps") or []
                ],
            }
            for slot in selected.get("irrep_slots") or []
        ]
        return common
    if stage == "opds":
        common["distortion_selection"] = value.get("distortion_selection") or {}
        common["selections"] = [
            {
                "slot": slot.get("slot"),
                "kpoint": _kpoint_result(dict(slot.get("kpoint") or {})),
                "irrep": _irrep_result(dict(slot.get("irrep") or {})),
            }
            for slot in selected.get("selection_slots") or []
        ]
        common["opds"] = [
            compact_opd_result(dict(row))
            for row in selected.get("opd_rows") or []
        ]
        return common
    if stage == "modes":
        common["distortion_selection"] = value.get("distortion_selection") or {}
        common["selections"] = [
            {
                "slot": slot.get("slot"),
                "kpoint": _kpoint_result(dict(slot.get("kpoint") or {})),
                "irrep": _irrep_result(dict(slot.get("irrep") or {})),
            }
            for slot in selected.get("selection_slots") or []
        ]
        common["opd"] = compact_opd_result(
            dict(selected.get("orderparam") or {})
        )
        invariant_factors = selected.get("invariant_factors")
        if isinstance(invariant_factors, list):
            common["invariant_factors"] = [
                _invariant_factor_result(dict(factor))
                for factor in invariant_factors
            ]
        mode_details = _public_value(selected.get("mode_details") or {})
        _attach_definition_ids(mode_details)
        subgroup_details = mode_details.get("subgroup_details")
        if isinstance(subgroup_details, dict):
            subgroup_details.pop("subgroup", None)
            subgroup_details.pop("magnetic_subgroup", None)
            subgroup_details["selected_subgroup"] = common["opd"]["subgroup"]
        common["mode_details"] = mode_details
        return common
    return value


__all__ = ["compact_opd_result", "compact_result"]
