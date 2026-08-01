"""Symbolic space-group/Wyckoff parent structures.

This input route establishes the crystallographic identities needed through
OPD and invariant analysis without inventing coordinates for free Wyckoff
parameters.  Coordinates are realized only when every required parameter is
supplied for mode-detail calculations.
"""

from __future__ import annotations

import copy
from collections.abc import Mapping, Sequence
from fractions import Fraction
import re
from typing import Any

from distortropy.Backend.source.tables import source_tables


def _mode_realization_cell(sg: int) -> dict[str, str]:
    if sg <= 2:
        values = ("4", "5", "6", "70", "80", "75")
    elif sg <= 15:
        values = ("4", "5", "6", "90", "105", "90")
    elif sg <= 74:
        values = ("4", "5", "6", "90", "90", "90")
    elif sg <= 142:
        values = ("4", "4", "6", "90", "90", "90")
    elif sg <= 194:
        values = ("4", "4", "6", "90", "90", "120")
    else:
        values = ("4", "4", "4", "90", "90", "90")
    return dict(zip(("a", "b", "c", "alpha", "beta", "gamma"), values, strict=True))


def _space_group(sg: Any) -> int:
    if type(sg) is not int or not 1 <= sg <= 230:
        raise ValueError("space_group must be an integer from 1 through 230")
    return sg


def _position(value: Any) -> tuple[str, int | None]:
    match = re.fullmatch(r"(?:(\d+))?([A-Za-z]+)", str(value).strip())
    if match is None:
        raise ValueError(f"invalid Wyckoff position {value!r}")
    multiplicity, label = match.groups()
    return label.casefold(), None if multiplicity is None else int(multiplicity)


def _parameters(value: Any, *, field: str) -> dict[str, Fraction]:
    if value is None:
        return {}
    if not isinstance(value, Mapping):
        raise ValueError(f"{field} must be an object")
    out: dict[str, Fraction] = {}
    for raw_name, raw_value in value.items():
        name = str(raw_name).strip().casefold()
        if not name:
            raise ValueError(f"{field} contains an empty parameter name")
        try:
            out[name] = Fraction(str(raw_value).strip())
        except (ValueError, ZeroDivisionError) as exc:
            raise ValueError(f"{field}.{name} must be rational") from exc
    return out


def _site_spec(value: Any, ordinal: int) -> tuple[str, int | None, dict[str, Fraction]]:
    if isinstance(value, str):
        position, separator, assignments = value.partition(":")
        params = None
        if separator:
            params = {}
            for assignment in assignments.split(","):
                name, equals, raw_value = assignment.partition("=")
                if not equals:
                    raise ValueError(
                        f"sites[{ordinal}] parameters must use name=value"
                    )
                params[name.strip()] = raw_value.strip()
    elif isinstance(value, Mapping):
        unknown = set(value) - {"wyckoff", "parameters"}
        if unknown:
            raise ValueError(f"sites[{ordinal}] has unknown fields: {sorted(unknown)!r}")
        position = value.get("wyckoff")
        params = value.get("parameters")
    else:
        raise ValueError(f"sites[{ordinal}] must be a Wyckoff string or object")
    label, multiplicity = _position(position)
    return label, multiplicity, _parameters(params, field=f"sites[{ordinal}].parameters")


def generic_parent_catalog(sg: int | None = None) -> dict[str, Any]:
    data = source_tables()
    if sg is None:
        return {
            "space_groups": [
                {
                    "number": number,
                    "symbol": str(data.default_inter_setting_record(number)["label_short"]),
                }
                for number in range(1, 231)
            ]
        }
    sg = _space_group(sg)
    setting = data.default_inter_setting_record(sg)
    setting_id = int(setting["id"])
    return {
        "space_group": {"number": sg, "symbol": str(setting["label_short"])},
        "wyckoff": [
            {
                "label": row.label,
                "multiplicity": data.wyckoff_multiplicity(sg, row),
                "row": int(row.row_id),
                "site_pg": data.site_point_group_hm_label(row.site_pg),
                "formula": formula["formula"],
                "parameters": list(formula["free"]),
            }
            for row in data.wyckoff_rows(sg)
            for formula in (data.inter_wyckoff_formula(sg, row, setting_id),)
        ],
    }


def build_generic_parent_state(
    space_group: int,
    sites: Sequence[Any] = (),
) -> dict[str, Any]:
    sg = _space_group(space_group)
    data = source_tables()
    setting = data.default_inter_setting_record(sg)
    setting_id = int(setting["id"])
    atom_sites: list[dict[str, Any]] = []
    for ordinal, raw_site in enumerate(sites, 1):
        label, supplied_multiplicity, parameters = _site_spec(raw_site, ordinal)
        try:
            row = data.wyckoff_row_by_label(sg, label)
        except KeyError as exc:
            raise ValueError(str(exc)) from exc
        multiplicity = data.wyckoff_multiplicity(sg, row)
        if supplied_multiplicity is not None and supplied_multiplicity != multiplicity:
            raise ValueError(
                f"Wyckoff position {label!r} has multiplicity {multiplicity} in SG{sg}"
            )
        formula = data.inter_wyckoff_formula(sg, row, setting_id)
        free = tuple(str(name) for name in formula["free"])
        unknown = sorted(set(parameters) - set(free))
        if unknown:
            raise ValueError(
                f"SG{sg} Wyckoff {label} has no parameters {unknown}; expected {list(free)}"
            )
        site: dict[str, Any] = {
            "type": label,
            "label": f"{label}{ordinal}",
            "multiplicity": str(multiplicity),
            "occupancy": "1",
            "wyckoff": label,
            "wyckoff_row_id": int(row.row_id),
            "site_pg": data.site_point_group_hm_label(row.site_pg),
            "wyckoff_formula": formula["formula"],
            "wyckoff_parameters": list(free),
            "wyckoff_params": {name: str(value) for name, value in parameters.items()},
            "wyckoff_multiplicity": multiplicity,
            "coordinates_realized": set(parameters) == set(free),
        }
        if site["coordinates_realized"]:
            vectors = data.inter_wyckoff_fraction_vectors(sg, row, setting_id)
            point = tuple(
                (
                    vectors[0][axis]
                    + sum(
                        parameters[name] * vectors[index + 1][axis]
                        for index, name in enumerate(("x", "y", "z"))
                        if name in parameters
                    )
                )
                % 1
                for axis in range(3)
            )
            site["fract"] = [str(float(value)) for value in point]
        atom_sites.append(site)

    lattice = int(data.space["ispace_lattice"][sg - 1])
    parent_inter_setting = data.inter_setting_record(setting_id)
    cif_info = {
        "path": None,
        "block": f"space_group_{sg}",
        "formula": None,
        "parent": {"number": sg, "symbol": str(setting["label_short"])},
        "explicit_parent_preference": None,
        "source": "space_group_wyckoff",
        "lattice": {
            key: None for key in ("a", "b", "c", "alpha", "beta", "gamma")
        },
        "atom_sites": atom_sites,
        "symmetry_operations": [],
    }
    return {
        "schema": "distortropy.parent.v1",
        "input": cif_info,
        "space_group": {
            "number": sg,
            "symbol": str(setting["label_short"]),
            "lattice": lattice,
            "lattice_type": str(data.space["lattice_label"][lattice - 1]).strip(),
        },
        "setting": {
            "default_preferences": data.default_space_group_preferences(),
            "default_inter_setting": setting,
            "parent_preferences": data.parent_preference_text(parent_inter_setting) or None,
            "parent_inter_setting": parent_inter_setting,
        },
        "_internal": {
            "parent_inter_setting_id": setting_id,
            "displacive_row_ids": sorted(
                {int(site["wyckoff_row_id"]) for site in atom_sites}
            ),
            "generic_parent": True,
        },
    }


def realize_generic_cif_info(cif_info: dict[str, Any]) -> dict[str, Any]:
    if cif_info.get("source") != "space_group_wyckoff":
        return cif_info
    missing = [
        f"{site.get('label')}:{name}"
        for site in cif_info.get("atom_sites") or []
        for name in site.get("wyckoff_parameters") or []
        if name not in (site.get("wyckoff_params") or {})
    ]
    if missing:
        raise ValueError(
            "mode details require Wyckoff parameters: " + ", ".join(missing)
        )
    state = build_generic_parent_state(
        int(cif_info["parent"]["number"]),
        [
            {
                "wyckoff": site["wyckoff"],
                "parameters": site.get("wyckoff_params") or {},
            }
            for site in cif_info.get("atom_sites") or []
        ],
    )
    realized = copy.deepcopy(state["input"])
    realized["lattice"] = _mode_realization_cell(int(cif_info["parent"]["number"]))
    realized["lattice_source"] = "generic_mode_geometry"
    for site in realized.get("atom_sites") or []:
        site["wyckoff_params"] = {
            name: str(float(Fraction(value)))
            for name, value in (site.get("wyckoff_params") or {}).items()
        }
    return realized


__all__ = [
    "build_generic_parent_state",
    "generic_parent_catalog",
    "realize_generic_cif_info",
]
