"""Local CIF-to-mode-details implementation surface.

The assembled pipeline follows the ISODISTORT-native path:

    CIF -> k vector -> irrep -> OPD -> mode details

The backend supports displacive, strain+displacive, and magnetic-inclusive
domain composition.
"""

from __future__ import annotations

from collections.abc import Mapping
from fractions import Fraction
import math
from typing import Any

from ISODISTORT.Assembled.Backend.pipeline import build_state_from_text
from ISODISTORT.Assembled.Backend.modes.mode_detail_text import render_mode_detail_text


TOOL = {
    "id": "isodistort",
    "category": "Representation Analysis",
    "title": "ISODISTORT Assembled",
    "frontend_dir": "Frontend",
}

_STATE_FIELDS = {
    "cif",
    "mode_set",
    "distortion_selection",
    "k_label",
    "k_index",
    "k_params",
    "irrep",
    "irrep_index",
    "orderparam",
    "orderparam_index",
    "selections",
    "include_opd",
    "include_mode_details",
}
_SELECTION_FIELDS = {
    "k_label",
    "k_index",
    "k_params",
    "display_k_params",
    "irrep",
    "irrep_index",
}
_DISTORTION_FIELDS = {"strain", "displacive_sites", "magnetic_sites"}


def _reject_unknown_fields(
    value: Mapping[Any, Any], *, allowed: set[str], field: str
) -> None:
    if any(not isinstance(key, str) for key in value):
        raise ValueError(f"{field} field names must be strings")
    unknown = set(value) - allowed
    if unknown:
        raise ValueError(f"unknown {field} fields: {sorted(unknown)!r}")


def _exact_text(value: Any, *, field: str) -> str:
    if not isinstance(value, str) or not value or value != value.strip():
        raise ValueError(f"{field} must be a non-empty string without outer whitespace")
    return value


def _optional_text(value: Any, *, field: str) -> str | None:
    if value is None:
        return None
    return _exact_text(value, field=field)


def _positive_index(value: Any, *, field: str) -> int:
    if type(value) is not int or value < 1:
        raise ValueError(f"{field} must be a positive integer")
    return value


def _boolean(value: Any, *, field: str) -> bool:
    if type(value) is not bool:
        raise ValueError(f"{field} must be boolean")
    return value


def _rational_params(value: Any, *, field: str) -> dict[str, str]:
    if value is None:
        return {}
    if not isinstance(value, Mapping):
        raise ValueError(f"{field} must be an object")
    out: dict[str, str] = {}
    for raw_key, raw_value in value.items():
        key = _exact_text(raw_key, field=f"{field} key")
        if isinstance(raw_value, str):
            text = raw_value.strip()
            if not text:
                continue
        elif type(raw_value) is int:
            text = str(raw_value)
        elif type(raw_value) is float and math.isfinite(raw_value):
            text = str(raw_value)
        else:
            raise ValueError(f"{field}.{key} must be a finite rational value")
        try:
            Fraction(text)
        except (ValueError, ZeroDivisionError) as exc:
            raise ValueError(f"{field}.{key} must be a finite rational value") from exc
        out[key] = text
    return out


def _site_indexes(value: Any, *, field: str) -> str | list[int]:
    if isinstance(value, str):
        if value not in {"all", "none"}:
            raise ValueError(f"{field} must be 'all', 'none', or a list of indexes")
        return value
    if not isinstance(value, list):
        raise ValueError(f"{field} must be 'all', 'none', or a list of indexes")
    indexes = [_positive_index(item, field=f"{field} item") for item in value]
    if len(indexes) != len(set(indexes)):
        raise ValueError(f"{field} must not contain duplicate indexes")
    return indexes


def _distortion_selection(value: Any) -> dict[str, Any] | None:
    if value is None:
        return None
    if not isinstance(value, Mapping):
        raise ValueError("distortion_selection must be an object")
    _reject_unknown_fields(
        value, allowed=_DISTORTION_FIELDS, field="distortion_selection"
    )
    out: dict[str, Any] = {}
    if "strain" in value:
        out["strain"] = _boolean(value["strain"], field="distortion_selection.strain")
    for name in ("displacive_sites", "magnetic_sites"):
        if name in value:
            out[name] = _site_indexes(value[name], field=f"distortion_selection.{name}")
    return out


def _selection(value: Any, *, ordinal: int) -> dict[str, Any]:
    field = f"selections[{ordinal}]"
    if not isinstance(value, Mapping):
        raise ValueError(f"{field} must be an object")
    _reject_unknown_fields(value, allowed=_SELECTION_FIELDS, field=field)
    if "k_label" not in value and "k_index" not in value:
        raise ValueError(f"{field} requires k_label or k_index")
    out: dict[str, Any] = {}
    if "k_label" in value:
        out["k_label"] = _exact_text(value["k_label"], field=f"{field}.k_label")
    if "k_index" in value:
        out["k_index"] = _positive_index(value["k_index"], field=f"{field}.k_index")
    for name in ("k_params", "display_k_params"):
        if name in value:
            out[name] = _rational_params(value[name], field=f"{field}.{name}")
    if "irrep" in value:
        out["irrep"] = _optional_text(value["irrep"], field=f"{field}.irrep")
    if "irrep_index" in value:
        out["irrep_index"] = _positive_index(
            value["irrep_index"], field=f"{field}.irrep_index"
        )
    return out


def _state_request(payload: Any) -> dict[str, Any]:
    if not isinstance(payload, Mapping):
        raise ValueError("request body must be a JSON object")
    _reject_unknown_fields(payload, allowed=_STATE_FIELDS, field="state")
    cif_text = payload.get("cif")
    if not isinstance(cif_text, str) or not cif_text.strip():
        raise ValueError("cif is required and must be text")
    mode_set = payload.get("mode_set", "strain_displacive")
    mode_set = _exact_text(mode_set, field="mode_set")
    if mode_set not in {"displacive", "strain_displacive", "magnetic", "all"}:
        raise ValueError(f"unknown mode_set: {mode_set}")
    selections = payload.get("selections")
    if selections is not None:
        if not isinstance(selections, list) or not 1 <= len(selections) <= 9:
            raise ValueError("selections must contain one to nine ordered objects")
        selections = [
            _selection(value, ordinal=ordinal)
            for ordinal, value in enumerate(selections, 1)
        ]
    orderparam = payload.get("orderparam")
    if orderparam is not None:
        if type(orderparam) is int:
            orderparam = _positive_index(orderparam, field="orderparam")
        else:
            orderparam = _exact_text(orderparam, field="orderparam")
    include_opd = _boolean(payload.get("include_opd", True), field="include_opd")
    include_mode_details = _boolean(
        payload.get("include_mode_details", include_opd),
        field="include_mode_details",
    )
    if include_mode_details and not include_opd:
        raise ValueError("include_mode_details requires include_opd")
    return {
        "cif_text": cif_text,
        "mode_set": mode_set,
        "distortion_selection": _distortion_selection(
            payload.get("distortion_selection")
        ),
        "k_label": _exact_text(payload.get("k_label", "GM"), field="k_label"),
        "k_index": _positive_index(payload.get("k_index", 1), field="k_index"),
        "k_params": _rational_params(payload.get("k_params"), field="k_params"),
        "irrep": _optional_text(payload.get("irrep"), field="irrep"),
        "irrep_index": _positive_index(
            payload.get("irrep_index", 1), field="irrep_index"
        ),
        "orderparam": orderparam,
        "orderparam_index": _positive_index(
            payload.get("orderparam_index", 1), field="orderparam_index"
        ),
        "selections": selections,
        "include_opd": include_opd,
        "include_mode_details": include_mode_details,
    }


def handle_api(endpoint: str, query: dict[str, list[str]]):
    raise ValueError(f"unknown endpoint: {endpoint}")


def handle_api_post(endpoint: str, query: dict[str, list[str]], payload: dict[str, Any]):
    if endpoint == "state":
        request = _state_request(payload)
        state = build_state_from_text(
            request.pop("cif_text"),
            **request,
        )
        if request["include_mode_details"]:
            state["selected"] = {
                **(state.get("selected") or {}),
                "mode_detail_text": render_mode_detail_text(state),
            }
        return state
    raise ValueError(f"unknown endpoint: {endpoint}")
