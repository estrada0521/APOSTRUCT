"""Normalize CLI and JSON inputs for the APOSTRUCT pipeline."""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path
from fractions import Fraction
import json
import re
from typing import Any

from APOSTRUCT.Backend.parent import (
    build_generic_parent_state,
    read_cif_summary,
)
from APOSTRUCT.Backend.reciprocal.catalog import (
    DISPLAY_K_PARAMS,
    kpoints,
)


@dataclass(frozen=True)
class CaseKSlot:
    label: str
    values: tuple[str, ...] = ()
    params: tuple[tuple[str, str], ...] = ()
    irrep: str | None = None


@dataclass(frozen=True)
class APOSTRUCTCase:
    structure: Path | None
    k: tuple[CaseKSlot, ...] = ()
    sites: tuple[tuple[str, bool, bool], ...] = ()
    strain: bool = False
    opd: str | None = None
    space_group: int | None = None
    wyckoff: tuple[Any, ...] = ()


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


def _cif_asset_path(value: Any, *, base: Path) -> Path:
    cif_id = _text(value, field="cif").casefold()
    if re.fullmatch(r"[0-9a-f]{12}", cif_id) is None:
        raise ValueError("cif must be a 12-character content ID")
    roots: list[Path] = []
    for start in (base.resolve(), Path.cwd().resolve(), Path(__file__).resolve()):
        for root in (start, *start.parents):
            if root not in roots:
                roots.append(root)
    for root in roots:
        candidate = root / "Assets" / "cif" / f"{cif_id}.cif"
        if candidate.is_file():
            return candidate.resolve()
    raise ValueError(f"CIF asset does not exist: Assets/cif/{cif_id}.cif")


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


def case_from_json(payload: Any, *, base: Path) -> APOSTRUCTCase:
    if not isinstance(payload, Mapping):
        raise ValueError("case JSON must be an object")
    allowed = {"cif", "structure", "sg", "wyckoff", "k", "sites", "strain", "opd"}
    if not set(payload) <= allowed:
        raise ValueError(f"case fields must be drawn from {sorted(allowed)!r}")
    structure_fields = set(payload) & {"cif", "structure", "sg"}
    if len(structure_fields) != 1:
        raise ValueError("case must contain exactly one of cif, structure, or sg")
    if "wyckoff" in payload and "sg" not in payload:
        raise ValueError("wyckoff requires sg")

    space_group: int | None = None
    wyckoff: tuple[str, ...] = ()
    if "sg" in payload:
        raw_sg = payload.get("sg")
        if isinstance(raw_sg, bool) or not isinstance(raw_sg, int):
            raise ValueError("sg must be an integer from 1 through 230")
        space_group = int(raw_sg)
        if not 1 <= space_group <= 230:
            raise ValueError("sg must be an integer from 1 through 230")
        raw_wyckoff = payload.get("wyckoff", ())
        if not isinstance(raw_wyckoff, Sequence) or isinstance(
            raw_wyckoff, (str, bytes)
        ):
            raise ValueError("wyckoff must be an array of positions")
        normalized_wyckoff: list[Any] = []
        for index, value in enumerate(raw_wyckoff, 1):
            field = f"wyckoff[{index}]"
            if isinstance(value, Mapping):
                if set(value) - {"wyckoff", "parameters"}:
                    raise ValueError(
                        f"{field} fields must be drawn from ['parameters', 'wyckoff']"
                    )
                raw_parameters = value.get("parameters", {})
                if not isinstance(raw_parameters, Mapping):
                    raise ValueError(f"{field}.parameters must be an object")
                normalized_wyckoff.append({
                    "wyckoff": _text(value.get("wyckoff"), field=f"{field}.wyckoff"),
                    "parameters": {
                        _text(name, field=f"{field}.parameters key"): _rational(
                            parameter, field=f"{field}.parameters.{name}"
                        )
                        for name, parameter in raw_parameters.items()
                    },
                })
            else:
                normalized_wyckoff.append(_text(value, field=field))
        wyckoff = tuple(normalized_wyckoff)
    raw_k = payload.get("k", ())
    if (
        not isinstance(raw_k, Sequence)
        or isinstance(raw_k, (str, bytes))
        or len(raw_k) > 4
    ):
        raise ValueError("k must contain zero to four ordered slots")

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
        if displacive or magnetic:
            sites.append((element, displacive, magnetic))
    if space_group is not None and sites and not wyckoff:
        raise ValueError("wyckoff is required when a space-group case selects sites")

    strain = payload.get("strain", False)
    if not isinstance(strain, bool):
        raise ValueError("strain must be boolean")
    opd = payload.get("opd")
    return APOSTRUCTCase(
        structure=(
            _cif_asset_path(payload.get("cif"), base=base)
            if "cif" in payload
            else (
                _structure_path(payload.get("structure"), base=base)
                if "structure" in payload
                else None
            )
        ),
        k=tuple(_json_slot(raw, ordinal) for ordinal, raw in enumerate(raw_k, 1)),
        sites=tuple(sites),
        strain=strain,
        opd=None if opd is None else _text(opd, field="opd"),
        space_group=space_group,
        wyckoff=wyckoff,
    )


def load_case(path: Path) -> APOSTRUCTCase:
    resolved = Path(path).expanduser().resolve()
    if not resolved.is_file():
        raise ValueError(f"case file does not exist: {resolved}")
    text = resolved.read_text(encoding="utf-8")
    if resolved.suffix.casefold() != ".json":
        raise ValueError("case input must be JSON")
    try:
        payload = json.loads(text)
    except json.JSONDecodeError as exc:
        raise ValueError(f"invalid case JSON: {exc}") from exc
    return case_from_json(payload, base=resolved.parent)


def case_from_cli(
    structure: Path,
    *,
    k: Sequence[Sequence[str]],
    irreps: Sequence[str],
    displacive: Sequence[str],
    magnetic: Sequence[str],
    strain: bool,
    opd: str | None,
) -> APOSTRUCTCase:
    if not 0 <= len(k) <= 4:
        raise ValueError("--k must contain at most four ordered slots")
    if irreps and len(irreps) != len(k):
        raise ValueError("each --k must have one corresponding --irrep")
    if not isinstance(strain, bool):
        raise ValueError("strain must be boolean")
    slots: list[CaseKSlot] = []
    for ordinal, values in enumerate(k, 1):
        label = _text(values[0] if values else None, field=f"k[{ordinal}].label")
        raw_values = tuple(values[1:])
        named = tuple("=" in str(value) for value in raw_values)
        if any(named) and not all(named):
            raise ValueError(
                f"k[{ordinal}] parameters must be entirely positional or NAME=VALUE"
            )
        params: list[tuple[str, str]] = []
        positional: list[str] = []
        if raw_values and all(named):
            seen: set[str] = set()
            for value in raw_values:
                raw_name, raw_parameter = str(value).split("=", 1)
                name = _text(raw_name, field=f"k[{ordinal}] parameter name")
                if name in seen:
                    raise ValueError(f"k[{ordinal}] parameter {name!r} is duplicated")
                seen.add(name)
                params.append(
                    (name, _rational(raw_parameter, field=f"k[{ordinal}].params.{name}"))
                )
        else:
            positional.extend(
                _rational(value, field=f"k[{ordinal}].values")
                for value in raw_values
            )
        slots.append(
            CaseKSlot(
                label=label,
                values=tuple(positional),
                params=tuple(params),
                irrep=irreps[ordinal - 1] if irreps else None,
            )
        )
    sites = tuple(
        (element, element in displacive, element in magnetic)
        for element in dict.fromkeys([*displacive, *magnetic])
    )
    return APOSTRUCTCase(
        structure=Path(structure).expanduser().resolve(),
        k=tuple(slots),
        sites=sites,
        strain=strain,
        opd=opd,
    )


def parameter_names_for_kpoint(kpoint: Mapping[str, Any]) -> tuple[str, ...]:
    expression = str(kpoint.get("display_kvector") or "")
    present = set(re.findall(r"(?<![A-Za-z])[abg](?![A-Za-z])", expression))
    return tuple(name for name in DISPLAY_K_PARAMS if name in present)


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
    names = parameter_names_for_kpoint(matches[0])
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
    case: APOSTRUCTCase,
    *,
    require_irreps: bool,
    default_displacive: bool = False,
    allow_empty_distortion: bool = False,
) -> dict[str, Any]:
    if case.structure is None:
        if case.space_group is None:
            raise ValueError("case has no parent structure authority")
        parent_state = build_generic_parent_state(case.space_group, case.wyckoff)
        return {
            **pipeline_request_from_cif_info(
                case,
                parent_state["input"],
                require_irreps=require_irreps,
                default_displacive=default_displacive,
                allow_empty_distortion=allow_empty_distortion,
            ),
            "parent_state": parent_state,
        }
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
    case: APOSTRUCTCase,
    cif_info: Mapping[str, Any],
    *,
    require_irreps: bool,
    default_displacive: bool = False,
    allow_empty_distortion: bool = False,
) -> dict[str, Any]:
    cif_info = dict(cif_info)
    atom_sites = list(cif_info.get("atom_sites") or [])
    atom_types = [
        str(atom.get("type") or atom.get("label") or index)
        for index, atom in enumerate(atom_sites, 1)
    ]
    atom_labels = [
        str(atom.get("label") or atom.get("type") or index)
        for index, atom in enumerate(atom_sites, 1)
    ]
    site_types = list(dict.fromkeys(atom_types))
    site_selectors = list(dict.fromkeys([*site_types, *atom_labels]))
    requested_sites = {
        name: (displacive, magnetic)
        for name, displacive, magnetic in case.sites
        if displacive or magnetic
    }
    unknown = sorted(set(requested_sites) - set(site_selectors))
    if unknown:
        raise ValueError(
            f"unknown CIF site selectors: {unknown}; available selectors: {site_selectors}"
        )
    displacive = [
        index
        for index, (site_type, site_label) in enumerate(
            zip(atom_types, atom_labels, strict=True), 1
        )
        if requested_sites.get(site_type, (False, False))[0]
        or requested_sites.get(site_label, (False, False))[0]
    ]
    magnetic = [
        index
        for index, (site_type, site_label) in enumerate(
            zip(atom_types, atom_labels, strict=True), 1
        )
        if requested_sites.get(site_type, (False, False))[1]
        or requested_sites.get(site_label, (False, False))[1]
    ]
    if not requested_sites and default_displacive:
        displacive = list(range(1, len(atom_sites) + 1))
    if (
        not allow_empty_distortion
        and not case.strain
        and not displacive
        and not magnetic
    ):
        raise ValueError("case must include strain, displacive, or magnetic distortion")

    pure_strain_irreps = (
        not require_irreps
        and not allow_empty_distortion
        and case.strain
        and not displacive
        and not magnetic
    )
    if pure_strain_irreps and any(
        slot.label.casefold() != "gm" for slot in case.k
    ):
        raise ValueError("strain-only irreps require k label GM or no k")
    if not case.k and not allow_empty_distortion and not pure_strain_irreps:
        raise ValueError("case requires one to four ordered k slots")

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
    "APOSTRUCTCase",
    "CaseKSlot",
    "case_from_cli",
    "case_from_json",
    "load_case",
    "parameter_names_for_kpoint",
    "pipeline_request",
    "pipeline_request_from_cif_info",
]
