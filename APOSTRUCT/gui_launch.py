"""Send one existing APOSTRUCT case to a running graphical interface."""

from __future__ import annotations

from copy import deepcopy
from fractions import Fraction
import json
import threading
from typing import Any, Callable, Mapping, Sequence
from urllib.error import HTTPError, URLError
from urllib.request import Request, urlopen

from APOSTRUCT.case_input import APOSTRUCTCase, pipeline_request
from APOSTRUCT.cli_mode_combination import parse_definition_weights


GUI_LAUNCH_SCHEMA = "APOSTRUCT.gui.launch.v1"
GUI_LAUNCH_STAGES = ("kpoints", "irreps", "opds", "modes")
_STATE_REQUEST_FIELDS = (
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
)


def _case_stage(case: APOSTRUCTCase) -> str:
    if case.opd is not None:
        return "modes"
    if case.k and all(slot.irrep for slot in case.k):
        return "opds"
    if case.k:
        return "irreps"
    return "kpoints"


def _generic_parent(request: Mapping[str, Any], case: APOSTRUCTCase) -> dict[str, Any]:
    parent_state = request.get("parent_state")
    if not isinstance(parent_state, Mapping):
        raise ValueError("generic case returned no parent state")
    parent_input = parent_state.get("input")
    if not isinstance(parent_input, Mapping):
        raise ValueError("generic case returned no parent input")
    return {
        "space_group": int(case.space_group or 0),
        "sites": [
            {
                "wyckoff": str(site.get("wyckoff") or ""),
                "parameters": dict(site.get("wyckoff_params") or {}),
            }
            for site in parent_input.get("atom_sites") or []
        ],
    }


def gui_launch_from_case(
    case: APOSTRUCTCase,
    *,
    raw_amplitudes: Sequence[str] = (),
) -> dict[str, Any]:
    """Build the GUI request for one already-normalized CLI case."""

    stage = _case_stage(case)
    request = pipeline_request(
        case,
        require_irreps=stage in {"opds", "modes"},
        allow_empty_distortion=stage == "kpoints",
    )
    state_request = {
        key: request[key]
        for key in _STATE_REQUEST_FIELDS
        if key in request and request[key] is not None
    }
    if case.structure is None:
        state_request["generic_parent"] = _generic_parent(request, case)
        source_name = f"SG {case.space_group}"
    else:
        try:
            state_request["cif"] = case.structure.read_text(encoding="utf-8")
        except (OSError, UnicodeError) as exc:
            raise ValueError(f"could not read structure CIF: {case.structure}") from exc
        source_name = case.structure.name
    state_request["include_opd"] = stage in {"opds", "modes"}
    state_request["include_mode_details"] = stage == "modes"

    if raw_amplitudes and stage != "modes":
        raise ValueError("--amplitude requires a case with opd")
    amplitudes = [
        {"definition_id": definition_id, "value": str(value)}
        for definition_id, value in (
            parse_definition_weights(raw_amplitudes) if raw_amplitudes else ()
        )
    ]
    return {
        "schema": GUI_LAUNCH_SCHEMA,
        "stage": stage,
        "source_name": source_name,
        "state_request": state_request,
        "amplitudes": amplitudes,
    }


def validate_gui_launch(value: Any) -> dict[str, Any]:
    """Validate the small transport envelope accepted by the local server."""

    if not isinstance(value, Mapping) or set(value) != {
        "schema",
        "stage",
        "source_name",
        "state_request",
        "amplitudes",
    }:
        raise ValueError("GUI launch has an invalid shape")
    if value.get("schema") != GUI_LAUNCH_SCHEMA:
        raise ValueError("GUI launch has an unknown schema")
    stage = value.get("stage")
    if stage not in GUI_LAUNCH_STAGES:
        raise ValueError("GUI launch has an unknown stage")
    source_name = value.get("source_name")
    if not isinstance(source_name, str) or not source_name.strip():
        raise ValueError("GUI launch source_name must be non-empty text")
    state_request = value.get("state_request")
    if not isinstance(state_request, Mapping):
        raise ValueError("GUI launch state_request must be an object")
    amplitudes = value.get("amplitudes")
    if not isinstance(amplitudes, list):
        raise ValueError("GUI launch amplitudes must be an array")
    normalized_amplitudes: list[dict[str, str]] = []
    seen: set[str] = set()
    for row in amplitudes:
        if not isinstance(row, Mapping) or set(row) != {"definition_id", "value"}:
            raise ValueError("GUI launch amplitude has an invalid shape")
        definition_id = row.get("definition_id")
        if not isinstance(definition_id, str) or not definition_id:
            raise ValueError("GUI launch amplitude definition_id must be non-empty text")
        if definition_id in seen:
            raise ValueError(f"duplicate amplitude for definition {definition_id}")
        seen.add(definition_id)
        try:
            amplitude = Fraction(str(row.get("value")))
        except (ValueError, ZeroDivisionError) as exc:
            raise ValueError(
                f"amplitude for definition {definition_id} must be an exact number"
            ) from exc
        normalized_amplitudes.append({
            "definition_id": definition_id,
            "value": str(amplitude),
        })
    if normalized_amplitudes and stage != "modes":
        raise ValueError("GUI launch amplitudes require the modes stage")
    return {
        "schema": GUI_LAUNCH_SCHEMA,
        "stage": stage,
        "source_name": source_name.strip(),
        "state_request": dict(state_request),
        "amplitudes": normalized_amplitudes,
    }


class GuiLaunchMailbox:
    """Keep only the latest launch for frontends attached to this process."""

    def __init__(self) -> None:
        self._lock = threading.Lock()
        self._revision = 0
        self._launch: dict[str, Any] | None = None

    def publish(self, launch: Mapping[str, Any]) -> int:
        with self._lock:
            self._revision += 1
            self._launch = deepcopy(dict(launch))
            return self._revision

    def read(self, after: int) -> dict[str, Any]:
        with self._lock:
            launch = (
                deepcopy(self._launch)
                if self._launch is not None and self._revision > after
                else None
            )
            return {"revision": self._revision, "launch": launch}


GUI_LAUNCH_MAILBOX = GuiLaunchMailbox()


def publish_gui_launch(
    launch: Mapping[str, Any],
    *,
    compute_state: Callable[[dict[str, Any], bool], dict[str, Any]],
) -> dict[str, Any]:
    normalized = validate_gui_launch(launch)
    state = compute_state(
        normalized["state_request"],
        normalized["stage"] == "kpoints",
    )
    if normalized["amplitudes"]:
        selected = state.get("selected") if isinstance(state, Mapping) else None
        mode_details = selected.get("mode_details") if isinstance(selected, Mapping) else None
        if not isinstance(mode_details, Mapping) or mode_details.get("status") != "ok":
            raise ValueError("GUI launch returned no mode definitions for amplitudes")
        available = {
            f"{kind}-{index}"
            for kind, field in (
                ("displacive", "displacive_definitions"),
                ("magnetic", "magnetic_definitions"),
                ("strain", "strain_definitions"),
            )
            for index, _definition in enumerate(mode_details.get(field) or [], start=1)
        }
        unknown = [
            row["definition_id"]
            for row in normalized["amplitudes"]
            if row["definition_id"] not in available
        ]
        if unknown:
            raise ValueError(
                f"mode definitions are unavailable in this case: {unknown}; "
                f"available IDs: {sorted(available)}"
            )
    prepared = {
        **normalized,
        "state": state,
    }
    revision = GUI_LAUNCH_MAILBOX.publish(prepared)
    return {"revision": revision, "stage": normalized["stage"]}


def send_gui_launch(launch: Mapping[str, Any], *, server: str) -> dict[str, Any]:
    base = str(server).strip().rstrip("/")
    if not base:
        raise ValueError("--server must be a non-empty URL")
    endpoint = f"{base}/api/APOSTRUCT/gui_launch"
    request = Request(
        endpoint,
        data=json.dumps(dict(launch), ensure_ascii=False).encode("utf-8"),
        headers={"Content-Type": "application/json"},
        method="POST",
    )
    try:
        with urlopen(request) as response:
            result = json.loads(response.read().decode("utf-8"))
    except HTTPError as exc:
        try:
            error = json.loads(exc.read().decode("utf-8")).get("error")
        except (UnicodeDecodeError, json.JSONDecodeError, AttributeError):
            error = None
        raise ValueError(str(error or f"GUI server returned HTTP {exc.code}")) from exc
    except (OSError, URLError) as exc:
        raise ValueError(f"could not reach APOSTRUCT GUI server at {base}") from exc
    if not isinstance(result, dict) or type(result.get("revision")) is not int:
        raise ValueError("GUI server returned an invalid launch response")
    return {**result, "server": base}


__all__ = [
    "GUI_LAUNCH_MAILBOX",
    "GUI_LAUNCH_SCHEMA",
    "GuiLaunchMailbox",
    "gui_launch_from_case",
    "publish_gui_launch",
    "send_gui_launch",
    "validate_gui_launch",
]
