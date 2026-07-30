"""Parse CLI, JSON, and line-oriented ISODISTORT case inputs."""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path
from fractions import Fraction
import json
import re
import shlex
from typing import Any

from ISODISTORT.Assembled.Backend.parent import read_cif_summary
from ISODISTORT.Assembled.Backend.reciprocal.catalog import (
    ISODISTORT_K_PARAMS,
    kpoints,
)


@dataclass(frozen=True)
class CaseKSlot:
    label: str
    values: tuple[str, ...] = ()
    params: tuple[tuple[str, str], ...] = ()
    irrep: str | None = None


@dataclass(frozen=True)
class AssembledCase:
    structure: Path
    k: tuple[CaseKSlot, ...] = ()
    sites: tuple[tuple[str, bool, bool], ...] = ()
    strain: bool = False
    opd: str | None = None


def _text(value: Any, *, field: str) -> str:
    if not isinstance(value, str) or not value.strip():
        raise ValueError(f"{field} must be a non-empty string")
    return value.strip()


def _rational(value: Any, *, field: str) -> str:
    if isinstance(value, bool):
        raise ValueError(f"{field} must be rational")
    try:
        return str(Fraction(str(value).strip()))
    except (ValueError, ZeroDivisionError) as exc:
        raise ValueError(f"{field} must be rational") from exc


def _structure_path(value: Any, *, base: Path) -> Path:
    raw = Path(_text(value, field="structure")).expanduser()
    path = raw if raw.is_absolute() else base / raw
    path = path.resolve()
    if not path.is_file():
        raise ValueError(f"structure CIF does not exist: {path}")
    return path


def _json_slot(raw: Any, ordinal: int) -> CaseKSlot:
    if not isinstance(raw, Mapping):
        raise ValueError(f"k[{ordinal}] must be an object")
    allowed = {"label", "params", "ir"}
    if not {"label"} <= set(raw) <= allowed:
        raise ValueError(f"k[{ordinal}] fields must be drawn from {sorted(allowed)!r}")
    raw_params = raw.get("params")
    if raw_params is None:
        raw_params = {}
    elif not isinstance(raw_params, Mapping):
        raise ValueError(f"k[{ordinal}].params must be an object")
    params = tuple(
        (str(name), _rational(value, field=f"k[{ordinal}].params.{name}"))
        for name, value in sorted(raw_params.items(), key=lambda item: str(item[0]))
    )
    irrep = raw.get("ir")
    return CaseKSlot(
        label=_text(raw.get("label"), field=f"k[{ordinal}].label"),
        params=params,
        irrep=None if irrep is None else _text(irrep, field=f"k[{ordinal}].ir"),
    )


def case_from_json(payload: Any, *, base: Path) -> AssembledCase:
    if not isinstance(payload, Mapping):
        raise ValueError("case JSON must be an object")
    allowed = {"structure", "k", "sites", "strain", "opd"}
    if not {"structure", "k"} <= set(payload) <= allowed:
        raise ValueError(f"case fields must be drawn from {sorted(allowed)!r}")
    raw_k = payload.get("k")
    if (
        not isinstance(raw_k, Sequence)
        or isinstance(raw_k, (str, bytes))
        or not 1 <= len(raw_k) <= 4
    ):
        raise ValueError("k must contain one to four ordered slots")

    raw_sites = payload.get("sites")
    if raw_sites is None:
        raw_sites = {}
    elif not isinstance(raw_sites, Mapping):
        raise ValueError("sites must be an object")
    sites: list[tuple[str, bool, bool]] = []
    for name, flags in raw_sites.items():
        element = _text(name, field="sites element")
        if not isinstance(flags, Mapping) or set(flags) != {"displacive", "magnetic"}:
            raise ValueError(
                f"sites[{element!r}] must contain exactly displacive and magnetic"
            )
        displacive = flags.get("displacive")
        magnetic = flags.get("magnetic")
        if not isinstance(displacive, bool) or not isinstance(magnetic, bool):
            raise ValueError(f"sites[{element!r}] flags must be boolean")
        sites.append((element, displacive, magnetic))

    strain = payload.get("strain", False)
    if not isinstance(strain, bool):
        raise ValueError("strain must be boolean")
    opd = payload.get("opd")
    return AssembledCase(
        structure=_structure_path(payload.get("structure"), base=base),
        k=tuple(_json_slot(raw, ordinal) for ordinal, raw in enumerate(raw_k, 1)),
        sites=tuple(sites),
        strain=strain,
        opd=None if opd is None else _text(opd, field="opd"),
    )


def _parse_in_lines(text: str) -> dict[str, Any]:
    structure: str | None = None
    slots: list[dict[str, Any]] = []
    displacive: list[str] = []
    magnetic: list[str] = []
    strain = False
    opd: str | None = None
    for lineno, line in enumerate(text.splitlines(), 1):
        try:
            tokens = shlex.split(line, comments=True, posix=True)
        except ValueError as exc:
            raise ValueError(f"line {lineno}: {exc}") from exc
        if not tokens:
            continue
        keyword = tokens[0].upper()
        values = tokens[1:]
        if keyword == "CIF":
            if len(values) != 1 or structure is not None:
                raise ValueError(f"line {lineno}: CIF requires one unique path")
            structure = values[0]
        elif keyword == "K":
            if not values or len(slots) >= 4:
                raise ValueError(f"line {lineno}: K requires a label and at most four slots")
            slots.append({"label": values[0], "values": values[1:], "irrep": None})
        elif keyword == "IR":
            if len(values) != 1 or not slots or slots[-1]["irrep"] is not None:
                raise ValueError(f"line {lineno}: IR must follow one K exactly once")
            slots[-1]["irrep"] = values[0]
        elif keyword == "STRAIN":
            if values:
                raise ValueError(f"line {lineno}: STRAIN takes no values")
            strain = True
        elif keyword == "DISPLACIVE":
            if not values:
                raise ValueError(f"line {lineno}: DISPLACIVE requires site names")
            displacive.extend(values)
        elif keyword == "MAGNETIC":
            if not values:
                raise ValueError(f"line {lineno}: MAGNETIC requires site names")
            magnetic.extend(values)
        elif keyword == "OPD":
            if not values or opd is not None:
                raise ValueError(f"line {lineno}: OPD requires one unique label")
            opd = " ".join(values)
        else:
            raise ValueError(f"line {lineno}: unknown directive {tokens[0]!r}")
    if structure is None:
        raise ValueError("case input requires CIF")
    if not slots:
        raise ValueError("case input requires at least one K")
    sites = {
        element: {
            "displacive": element in displacive,
            "magnetic": element in magnetic,
        }
        for element in dict.fromkeys([*displacive, *magnetic])
    }
    return {
        "structure": structure,
        "k": slots,
        "sites": sites,
        "strain": strain,
        "opd": opd,
    }


def case_from_in(text: str, *, base: Path) -> AssembledCase:
    raw = _parse_in_lines(text)
    return AssembledCase(
        structure=_structure_path(raw["structure"], base=base),
        k=tuple(
            CaseKSlot(
                label=_text(slot["label"], field=f"k[{ordinal}].label"),
                values=tuple(
                    _rational(value, field=f"k[{ordinal}].values")
                    for value in slot["values"]
                ),
                irrep=(
                    None
                    if slot["irrep"] is None
                    else _text(slot["irrep"], field=f"k[{ordinal}].ir")
                ),
            )
            for ordinal, slot in enumerate(raw["k"], 1)
        ),
        sites=tuple(
            (name, bool(flags["displacive"]), bool(flags["magnetic"]))
            for name, flags in raw["sites"].items()
        ),
        strain=bool(raw["strain"]),
        opd=raw["opd"],
    )


def load_case(path: Path) -> AssembledCase:
    resolved = Path(path).expanduser().resolve()
    if not resolved.is_file():
        raise ValueError(f"case file does not exist: {resolved}")
    text = resolved.read_text(encoding="utf-8")
    if resolved.suffix.casefold() == ".json":
        try:
            payload = json.loads(text)
        except json.JSONDecodeError as exc:
            raise ValueError(f"invalid case JSON: {exc}") from exc
        return case_from_json(payload, base=resolved.parent)
    return case_from_in(text, base=resolved.parent)


def case_from_cli(
    structure: Path,
    *,
    k: Sequence[Sequence[str]],
    irreps: Sequence[str],
    displacive: Sequence[str],
    magnetic: Sequence[str],
    strain: bool,
    opd: str | None,
) -> AssembledCase:
    if not 0 <= len(k) <= 4:
        raise ValueError("--k must contain at most four ordered slots")
    if irreps and len(irreps) != len(k):
        raise ValueError("each --k must have one corresponding --irrep")
    if not isinstance(strain, bool):
        raise ValueError("strain must be boolean")
    slots = tuple(
        CaseKSlot(
            label=_text(values[0] if values else None, field=f"k[{ordinal}].label"),
            values=tuple(
                _rational(value, field=f"k[{ordinal}].values")
                for value in values[1:]
            ),
            irrep=irreps[ordinal - 1] if irreps else None,
        )
        for ordinal, values in enumerate(k, 1)
    )
    sites = tuple(
        (element, element in displacive, element in magnetic)
        for element in dict.fromkeys([*displacive, *magnetic])
    )
    return AssembledCase(
        structure=Path(structure).expanduser().resolve(),
        k=slots,
        sites=sites,
        strain=strain,
        opd=opd,
    )


def _parameter_names(kpoint: Mapping[str, Any]) -> tuple[str, ...]:
    expression = str(kpoint.get("isodistort_kvector") or "")
    present = set(re.findall(r"(?<![A-Za-z])[abg](?![A-Za-z])", expression))
    return tuple(name for name in ISODISTORT_K_PARAMS if name in present)


def _slot_params(
    slot: CaseKSlot,
    *,
    catalog: Sequence[Mapping[str, Any]],
    ordinal: int,
) -> dict[str, str]:
    matches = [
        item for item in catalog
        if str(item.get("label") or "").casefold() == slot.label.casefold()
    ]
    if len(matches) != 1:
        raise ValueError(f"k[{ordinal}] label {slot.label!r} is not uniquely available")
    names = _parameter_names(matches[0])
    if slot.params:
        params = dict(slot.params)
        if set(params) != set(names):
            raise ValueError(
                f"k[{ordinal}] {slot.label} parameters must be exactly {list(names)!r}"
            )
        return params
    if len(slot.values) != len(names):
        raise ValueError(
            f"k[{ordinal}] {slot.label} requires {len(names)} parameters {list(names)!r}"
        )
    return dict(zip(names, slot.values, strict=True))


def pipeline_request(
    case: AssembledCase,
    *,
    require_irreps: bool,
    default_displacive: bool = False,
    allow_empty_distortion: bool = False,
) -> dict[str, Any]:
    if not case.structure.is_file():
        raise ValueError(f"structure CIF does not exist: {case.structure}")
    return pipeline_request_from_cif_info(
        case,
        read_cif_summary(case.structure),
        require_irreps=require_irreps,
        default_displacive=default_displacive,
        allow_empty_distortion=allow_empty_distortion,
    )


def pipeline_request_from_cif_info(
    case: AssembledCase,
    cif_info: Mapping[str, Any],
    *,
    require_irreps: bool,
    default_displacive: bool = False,
    allow_empty_distortion: bool = False,
) -> dict[str, Any]:
    cif_info = dict(cif_info)
    site_types: list[str] = []
    for index, atom in enumerate(cif_info.get("atom_sites") or [], 1):
        name = str(atom.get("type") or atom.get("label") or index)
        if name not in site_types:
            site_types.append(name)
    requested_sites = {name: (displacive, magnetic) for name, displacive, magnetic in case.sites}
    unknown = sorted(set(requested_sites) - set(site_types))
    if unknown:
        raise ValueError(f"unknown CIF site types: {unknown}")
    displacive = [
        index for index, name in enumerate(site_types, 1)
        if requested_sites.get(name, (False, False))[0]
    ]
    magnetic = [
        index for index, name in enumerate(site_types, 1)
        if requested_sites.get(name, (False, False))[1]
    ]
    if not case.sites and default_displacive:
        displacive = list(range(1, len(site_types) + 1))
    if (
        not allow_empty_distortion
        and not case.strain
        and not displacive
        and not magnetic
    ):
        raise ValueError("case must include strain, displacive, or magnetic distortion")

    catalog = kpoints(int(cif_info["parent"]["number"]))["kpoints"]
    selections: list[dict[str, Any]] = []
    for ordinal, slot in enumerate(case.k, 1):
        if require_irreps and not slot.irrep:
            raise ValueError(f"k[{ordinal}] requires an irrep")
        params = _slot_params(slot, catalog=catalog, ordinal=ordinal)
        selections.append({
            "k_label": slot.label,
            "k_index": 1,
            "k_params": params,
            "display_k_params": params,
            "irrep": slot.irrep,
            "irrep_index": 1,
        })
    mode_set = "all" if magnetic else "strain_displacive" if case.strain else "displacive"
    first = selections[0] if selections else {
        "k_label": "GM",
        "k_index": 1,
        "k_params": {},
        "display_k_params": {},
        "irrep": None,
        "irrep_index": 1,
    }
    return {
        "cif_info": cif_info,
        "mode_set": mode_set,
        "distortion_selection": {
            "strain": case.strain,
            "displacive_sites": displacive,
            "magnetic_sites": magnetic,
        },
        "k_label": first["k_label"],
        "k_index": first["k_index"],
        "k_params": first["k_params"],
        "display_k_params": first["display_k_params"],
        "irrep": first["irrep"],
        "irrep_index": first["irrep_index"],
        "orderparam": case.opd,
        "selections": selections or None,
    }


__all__ = [
    "AssembledCase",
    "CaseKSlot",
    "case_from_cli",
    "case_from_in",
    "case_from_json",
    "load_case",
    "pipeline_request",
    "pipeline_request_from_cif_info",
]
