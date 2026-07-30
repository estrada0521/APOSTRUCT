"""Local CIF-to-mode-details implementation surface.

The assembled pipeline follows the ISODISTORT-native path:

    CIF -> k vector -> irrep -> OPD -> mode details

The backend supports displacive, strain+displacive, and magnetic-inclusive
domain composition.
"""

from __future__ import annotations

from collections.abc import Mapping
from fractions import Fraction
from functools import lru_cache
import math
from typing import Any

from ISODISTORT.Assembled.Backend.pipeline import build_state_from_text
from ISODISTORT.Assembled.Backend.modes.mode_detail_text import render_mode_detail_text


TOOL = {
    "id": "isodistort",
    "category": "Representation Analysis",
    "title": "Distortropy",
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
_INVARIANT_FIELDS = {
    "space_group",
    "factors",
    "minimum_degree",
    "maximum_degree",
}
_INVARIANT_FACTOR_FIELDS = {
    "slot",
    "gid",
    "label",
    "magnetic",
    "k_parameters",
    "opd",
    "domain",
    "parameter_offset",
}


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


def _rational_sequence(value: Any, *, field: str) -> list[str]:
    if not isinstance(value, list):
        raise ValueError(f"{field} must be an array")
    out: list[str] = []
    for index, raw_value in enumerate(value):
        item_field = f"{field}[{index}]"
        if isinstance(raw_value, str):
            text = raw_value.strip()
            if not text:
                raise ValueError(f"{item_field} must be a finite rational value")
        elif type(raw_value) is int:
            text = str(raw_value)
        elif type(raw_value) is float and math.isfinite(raw_value):
            text = str(raw_value)
        else:
            raise ValueError(f"{item_field} must be a finite rational value")
        try:
            Fraction(text)
        except (ValueError, ZeroDivisionError) as exc:
            raise ValueError(f"{item_field} must be a finite rational value") from exc
        out.append(text)
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


def _invariant_service():
    from ISODISTORT.Assembled.Backend.invariants.service import invariant_service

    return invariant_service()


@lru_cache(maxsize=4096)
def _invariant_domain_option(
    sg: int,
    gid: int,
    magnetic: bool,
    k_parameters: tuple[str, ...],
    opd: str,
    domain: int,
) -> tuple[str, str, str, str, str, int, list[str]]:
    options = _invariant_service().domains(
        {
            "space_group": int(sg),
            "gid": int(gid),
            "magnetic": bool(magnetic),
            "k_parameters": list(k_parameters),
            "direction": str(opd),
            "parameter_offset": 0,
        }
    )["domains"]
    item = next(item for item in options if item["number"] == domain)
    return (
        str(item["display_opd"]),
        str(item["subgroup"]),
        str(item["basis"]),
        str(item["origin"]),
        str(item["display"]),
        int(item["parameter_count"]),
        [str(value) for value in item.get("ferroic_properties") or []],
    )


def _complete_invariant_factors(state: dict[str, Any]) -> None:
    """Finish offsets and classify only unnamed induced subspaces."""
    selected = state.get("selected") or {}
    mode_details = selected.get("mode_details") or {}
    factors = mode_details.get("secondary_opd_factors") or []
    sg = int((state.get("space_group") or state.get("input", {}).get("parent"))["number"])
    parameter_offset = 0
    for factor in selected.get("opd_factors") or []:
        parameter_count = factor.get("parameter_count")
        if type(parameter_count) is not int or parameter_count < 1:
            raise ValueError("selected primary OPD lost its Source parameter count")
        factor["parameter_offset"] = parameter_offset
        parameter_offset += parameter_count
    completed_factors: list[dict[str, Any]] = []
    for raw_factor in factors:
        factor = dict(raw_factor)
        direction_matrix = factor.pop("direction_matrix", None)
        if direction_matrix is not None:
            resolved = _invariant_service().match_subspace(
                {
                    "space_group": sg,
                    "gid": int(factor["gid"]),
                    "magnetic": bool(factor.get("magnetic")),
                    "k_parameters": list(factor.get("k_parameters") or []),
                    "direction_matrix": direction_matrix,
                    "include_domain_details": True,
                }
            )
            if resolved is None:
                continue
            factor.update(resolved)
        if not factor.get("domain_display"):
            option = _invariant_domain_option(
                sg,
                int(factor["gid"]),
                bool(factor.get("magnetic")),
                tuple(str(value) for value in factor.get("k_parameters") or []),
                str(factor["opd"]),
                int(factor["domain"]),
            )
            factor["domain_display"] = option[4]
            factor["parameter_count"] = option[5]
            factor["domain_presentation"] = {
                "display_opd": option[0],
                "subgroup": option[1],
                "basis": option[2],
                "origin": option[3],
            }
            factor["ferroic_properties"] = option[6]
        parameter_count = factor.get("parameter_count")
        if type(parameter_count) is not int or parameter_count < 1:
            raise ValueError("secondary OPD lost its Source parameter count")
        factor["parameter_offset"] = parameter_offset
        parameter_offset += parameter_count
        completed_factors.append(factor)
    mode_details["secondary_opd_factors"] = completed_factors


def _invariant_request(payload: Any) -> dict[str, Any]:
    if not isinstance(payload, Mapping):
        raise ValueError("request body must be a JSON object")
    _reject_unknown_fields(payload, allowed=_INVARIANT_FIELDS, field="invariants")
    factors = payload.get("factors")
    if not isinstance(factors, list) or not factors:
        raise ValueError("invariants.factors must be a non-empty array")
    space_group = _positive_index(
        payload.get("space_group"), field="invariants.space_group"
    )
    minimum_degree = _positive_index(
        payload.get("minimum_degree"), field="invariants.minimum_degree"
    )
    maximum_degree = _positive_index(
        payload.get("maximum_degree"), field="invariants.maximum_degree"
    )
    service = _invariant_service()
    resolved: list[dict[str, Any]] = []
    for ordinal, factor in enumerate(factors, 1):
        field = f"invariants.factors[{ordinal}]"
        if not isinstance(factor, Mapping):
            raise ValueError(f"{field} must be an object")
        _reject_unknown_fields(
            factor,
            allowed=_INVARIANT_FACTOR_FIELDS,
            field=field,
        )
        gid = _positive_index(factor.get("gid"), field=f"{field}.gid")
        parameter_count = service.projection.k_parameter_dimension_by_gid(gid)
        parameters = _rational_sequence(
            factor.get("k_parameters"), field=f"{field}.k_parameters"
        )
        if len(parameters) != parameter_count:
            raise ValueError(
                f"{field}.k_parameters must contain exactly {parameter_count} values"
            )
        resolved_factor = {
            "gid": gid,
            "k_parameters": parameters,
            "magnetic": _boolean(
                factor.get("magnetic", False), field=f"{field}.magnetic"
            ),
            "direction": _exact_text(factor.get("opd"), field=f"{field}.opd"),
            "domain": _positive_index(
                factor.get("domain"), field=f"{field}.domain"
            ),
        }
        if "parameter_offset" in factor:
            parameter_offset = factor["parameter_offset"]
            if type(parameter_offset) is not int or parameter_offset < 0:
                raise ValueError(
                    f"{field}.parameter_offset must be a nonnegative integer"
                )
            resolved_factor["parameter_offset"] = parameter_offset
        resolved.append(resolved_factor)
    return {
        "space_group": space_group,
        "mode": "opd",
        "factors": resolved,
        "minimum_degree": minimum_degree,
        "maximum_degree": maximum_degree,
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
            _complete_invariant_factors(state)
        return state
    if endpoint == "invariants":
        if query:
            raise ValueError("invariants does not accept query parameters")
        return _invariant_service().compute(_invariant_request(payload))
    raise ValueError(f"unknown endpoint: {endpoint}")
