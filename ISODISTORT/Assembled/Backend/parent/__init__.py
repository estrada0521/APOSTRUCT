"""Parse the parent structure and establish its displayed setting."""

from __future__ import annotations

import copy
from fractions import Fraction
import math
import re
from pathlib import Path
from typing import Any

import gemmi

from ISODISTORT.Assembled.Backend.source.tables import SourceTables, source_tables


def float_cif_number(value: Any) -> float | None:
    text = str(value or "").strip()
    if not text or text in {".", "?"}:
        return None
    text = re.sub(r"\([^)]*\)", "", text)
    try:
        return float(text)
    except ValueError:
        return None


def parent_cell_tuple(cif_info: dict[str, Any]) -> tuple[float, float, float, float, float, float] | None:
    lattice = cif_info.get("lattice") or {}
    values = tuple(float_cif_number(lattice.get(key)) for key in ("a", "b", "c", "alpha", "beta", "gamma"))
    if any(value is None for value in values):
        return None
    return values  # type: ignore[return-value]


def _strip_cif(value: str | None) -> str | None:
    if value is None:
        return None
    value = str(value).strip()
    if len(value) >= 2 and value[0] == value[-1] and value[0] in {"'", '"'}:
        value = value[1:-1]
    return value


def _cif_value(block: gemmi.cif.Block, tags: list[str]) -> str | None:
    for tag in tags:
        value = _strip_cif(block.find_value(tag))
        if value:
            return value
    return None


def _parse_parent_space_group(block: gemmi.cif.Block) -> dict[str, Any]:
    number_text = _cif_value(
        block,
        [
            "_symmetry_Int_Tables_number",
            "_space_group_IT_number",
            "_space_group.it_number",
            "_space_group_IT_number",
        ],
    )
    if number_text is None:
        raise ValueError("CIF does not contain an explicit IT space-group number")
    match = re.search(r"\d+", number_text)
    if not match:
        raise ValueError(f"cannot parse IT space-group number from {number_text!r}")
    number = int(match.group(0))
    symbol = _cif_value(
        block,
        [
            "_symmetry_space_group_name_H-M",
            "_space_group_name_H-M_alt",
            "_space_group.name_H-M_alt",
        ],
    )
    return {"number": number, "symbol": symbol or ""}


def _parse_explicit_parent_preference(block: gemmi.cif.Block) -> str | None:
    for value in (
        _cif_value(block, ["_space_group.reference_setting"]),
        _cif_value(block, ["_space_group_reference_setting"]),
        _cif_value(block, ["_symmetry_space_group_setting"]),
    ):
        if value:
            return value
    for tag in (
        "_symmetry_space_group_name_H-M",
        "_space_group_name_H-M_alt",
        "_space_group.name-H-M_alt",
        "_space_group.name_H-M_alt",
    ):
        value = _cif_value(block, [tag])
        if not value:
            continue
        match = re.search(r"\(([^)]*(?:origin|choice|axis|setting)[^)]*)\)", value, flags=re.I)
        if match:
            return match.group(1).strip()
    return None


def _chemical_formula(block: gemmi.cif.Block) -> str | None:
    formula = _cif_value(
        block,
        [
            "_chemical_formula_structural",
            "_chemical_formula_sum",
            "_chemical_formula_moiety",
            "_chemical_formula.structural",
            "_chemical_formula.sum",
            "_chemical_formula.moiety",
        ],
    )
    return None if formula in {None, ".", "?"} else formula


def _element_symbol(value: Any) -> str | None:
    match = re.match(r"\s*([A-Z][a-z]?)", str(value or ""))
    return None if match is None else match.group(1)


def _formula_number(value: Any, *, default: Fraction | None = None) -> Fraction | None:
    text = str(value or "").strip()
    if not text or text in {".", "?"}:
        return default
    text = re.sub(r"\([^)]*\)$", "", text)
    try:
        return Fraction(text)
    except (ValueError, ZeroDivisionError):
        return None


def _formula_from_sites(sites: list[dict[str, Any]]) -> str | None:
    totals: dict[str, Fraction] = {}
    order: list[str] = []
    for site in sites:
        element = _element_symbol(site.get("type") or site.get("label"))
        multiplicity = _formula_number(
            site.get("wyckoff_multiplicity") or site.get("multiplicity")
        )
        occupancy = _formula_number(site.get("occupancy"), default=Fraction(1))
        if element is None or multiplicity is None or occupancy is None:
            return None
        if element not in totals:
            totals[element] = Fraction(0)
            order.append(element)
        totals[element] += multiplicity * occupancy
    if not totals or any(value <= 0 for value in totals.values()):
        return None
    denominator = math.lcm(*(value.denominator for value in totals.values()))
    counts = [int(totals[element] * denominator) for element in order]
    divisor = math.gcd(*counts)
    counts = [count // divisor for count in counts]
    return " ".join(
        element if count == 1 else f"{element}{count}"
        for element, count in zip(order, counts, strict=True)
    )


def _atom_sites(block: gemmi.cif.Block) -> list[dict[str, Any]]:
    label_col = block.find_loop("_atom_site_label")
    if not label_col:
        return []
    loop = label_col.get_loop()
    tags = [tag.lower() for tag in loop.tags]

    def col(name: str) -> int | None:
        lname = name.lower()
        return tags.index(lname) if lname in tags else None

    idx_type = col("_atom_site_type_symbol")
    idx_label = col("_atom_site_label")
    idx_mult = col("_atom_site_symmetry_multiplicity")
    idx_x = col("_atom_site_fract_x")
    idx_y = col("_atom_site_fract_y")
    idx_z = col("_atom_site_fract_z")
    idx_occ = col("_atom_site_occupancy")
    idx_wyckoff = None
    for tag in (
        "_atom_site_Wyckoff_label",
        "_atom_site_wyckoff_label",
        "_atom_site_Wyckoff_symbol",
        "_atom_site_wyckoff_symbol",
    ):
        idx_wyckoff = col(tag)
        if idx_wyckoff is not None:
            break
    out: list[dict[str, Any]] = []
    for row in range(loop.length()):
        item: dict[str, Any] = {}
        if idx_type is not None:
            item["type"] = _strip_cif(loop[row, idx_type])
        if idx_label is not None:
            item["label"] = _strip_cif(loop[row, idx_label])
        if idx_mult is not None:
            item["multiplicity"] = _strip_cif(loop[row, idx_mult])
        if idx_wyckoff is not None:
            item["wyckoff"] = _strip_cif(loop[row, idx_wyckoff])
        if idx_x is not None and idx_y is not None and idx_z is not None:
            item["fract"] = [
                _strip_cif(loop[row, idx_x]),
                _strip_cif(loop[row, idx_y]),
                _strip_cif(loop[row, idx_z]),
            ]
        if idx_occ is not None:
            item["occupancy"] = _strip_cif(loop[row, idx_occ])
        out.append(item)
    return out


def _symmetry_operations(block: gemmi.cif.Block) -> list[str]:
    for tag in (
        "_symmetry_equiv_pos_as_xyz",
        "_space_group_symop_operation_xyz",
        "_space_group_symop.operation_xyz",
    ):
        col = block.find_loop(tag)
        if not col:
            continue
        loop = col.get_loop()
        tags = [item.lower() for item in loop.tags]
        idx = tags.index(tag.lower())
        return [
            str(_strip_cif(loop[row, idx]) or "").strip()
            for row in range(loop.length())
            if str(_strip_cif(loop[row, idx]) or "").strip()
        ]
    return []


def _lattice_parameters(block: gemmi.cif.Block) -> dict[str, str | None]:
    return {
        "a": _cif_value(block, ["_cell_length_a", "_cell.length_a"]),
        "b": _cif_value(block, ["_cell_length_b", "_cell.length_b"]),
        "c": _cif_value(block, ["_cell_length_c", "_cell.length_c"]),
        "alpha": _cif_value(block, ["_cell_angle_alpha", "_cell.angle_alpha"]),
        "beta": _cif_value(block, ["_cell_angle_beta", "_cell.angle_beta"]),
        "gamma": _cif_value(block, ["_cell_angle_gamma", "_cell.angle_gamma"]),
    }


def read_cif_summary_from_block(block: gemmi.cif.Block, path: str | None = None) -> dict[str, Any]:
    parent = _parse_parent_space_group(block)
    return {
        "path": path,
        "block": block.name,
        "formula": _chemical_formula(block),
        "parent": parent,
        "explicit_parent_preference": _parse_explicit_parent_preference(block),
        "lattice": _lattice_parameters(block),
        "atom_sites": _atom_sites(block),
        "symmetry_operations": _symmetry_operations(block),
    }


def read_cif_summary(cif: Path) -> dict[str, Any]:
    block = gemmi.cif.read_file(str(cif)).sole_block()
    return read_cif_summary_from_block(block, str(cif))


def read_cif_summary_from_text(cif_text: str) -> dict[str, Any]:
    block = gemmi.cif.read_string(cif_text).sole_block()
    return read_cif_summary_from_block(block, None)


def _equivalent_site_images(site: dict[str, Any], operations: list[str]) -> list[dict[str, Any]]:
    if not site.get("fract"):
        return []
    try:
        fract = tuple(float(x) % 1.0 for x in site["fract"])
    except (TypeError, ValueError):
        return []
    out = [{
        "fract": fract,
        "operation": "x,y,z",
        "operation_index": None,
    }]
    seen = {tuple(round(x, 10) for x in fract)}
    for operation_index, text in enumerate(operations):
        try:
            transformed = gemmi.Op(text).apply_to_xyz(fract)
        except Exception:
            continue
        wrapped = tuple(float(x) % 1.0 for x in transformed)
        key = tuple(round(x, 10) for x in wrapped)
        if key not in seen:
            seen.add(key)
            out.append({
                "fract": wrapped,
                "operation": str(text),
                "operation_index": operation_index,
            })
    return out


def _with_representative_image(
    match: dict[str, Any],
    site: dict[str, Any],
    image: dict[str, Any],
) -> dict[str, Any]:
    enriched = dict(match)
    enriched["representative_source_fract"] = list(site.get("fract") or [])
    enriched["representative_matched_fract"] = list(image["fract"])
    enriched["representative_operation"] = image["operation"]
    enriched["representative_operation_index"] = image["operation_index"]
    return enriched


def _canonical_representative_mapping(
    site: dict[str, Any],
    match: dict[str, Any],
    operations: list[str],
    tol: float = 2e-5,
) -> dict[str, Any] | None:
    canonical = match.get("formula_representative_fract")
    if not site.get("fract") or canonical is None:
        return None
    try:
        source = tuple(float(value) for value in site["fract"])
        target = tuple(float(value) for value in canonical)
    except (TypeError, ValueError):
        return None
    candidates = [(None, "x,y,z"), *enumerate(operations)]
    for operation_index, text in candidates:
        try:
            applied = tuple(float(value) for value in gemmi.Op(text).apply_to_xyz(source))
        except Exception:
            continue
        translation = tuple(applied[axis] - target[axis] for axis in range(3))
        integer_translation = tuple(int(round(value)) for value in translation)
        if all(abs(value - integer) <= tol for value, integer in zip(translation, integer_translation)):
            return {
                "representative_operation": str(text),
                "representative_operation_index": operation_index,
                "representative_source_fract": list(site["fract"]),
                "representative_matched_fract": list(target),
                "representative_lattice_translation": list(integer_translation),
            }
    return None


def _match_wyckoff_site_with_symmetry(
    local_data: SourceTables,
    sg: int,
    site: dict[str, Any],
    operations: list[str],
) -> dict[str, Any] | None:
    direct = local_data.match_wyckoff_site(sg, site)
    if direct:
        return _with_representative_image(direct, site, _equivalent_site_images(site, [])[0])
    for image in _equivalent_site_images(site, operations)[1:]:
        fract = image["fract"]
        candidate = dict(site)
        candidate["fract"] = [str(value) for value in fract]
        match = local_data.match_wyckoff_site(sg, candidate)
        if match:
            return _with_representative_image(match, site, image)
    lattice = int(local_data.space["ispace_lattice"][int(sg) - 1])
    if lattice == 2 and site.get("fract"):
        try:
            x, y, z = (float(value) % 1.0 for value in site["fract"])
        except (TypeError, ValueError):
            return None
        candidate = dict(site)
        candidate["fract"] = [str(x), str(z), str(y)]
        match = local_data.match_wyckoff_site(sg, candidate)
        if match:
            return _with_representative_image(match, site, {
                "fract": (x, z, y),
                "operation": "x,z,y",
                "operation_index": None,
            })
    return None


def _match_inter_wyckoff_site_with_symmetry(
    local_data: SourceTables,
    sg: int,
    site: dict[str, Any],
    operations: list[str],
    setting_id: int,
) -> dict[str, Any] | None:
    direct = local_data.match_inter_wyckoff_site(sg, site, setting_id)
    if direct:
        return _with_representative_image(direct, site, _equivalent_site_images(site, [])[0])
    for image in _equivalent_site_images(site, operations)[1:]:
        fract = image["fract"]
        candidate = dict(site)
        candidate["_source_fract"] = site.get("fract")
        candidate["fract"] = [str(value) for value in fract]
        match = local_data.match_inter_wyckoff_site(sg, candidate, setting_id)
        if match:
            return _with_representative_image(match, site, image)
    return None


def _match_parent_inter_setting(
    local_data: SourceTables,
    sg: int,
    atom_sites: list[dict[str, Any]],
    operations: list[str],
    setting_id: int,
) -> tuple[
    tuple[int, int, int, int, int],
    tuple[dict[str, Any] | None, ...],
]:
    matched = 0
    direct = 0
    constrained = 0
    matches: list[dict[str, Any] | None] = []
    for site in atom_sites:
        if not site.get("fract"):
            matches.append(None)
            continue
        match = local_data.match_inter_wyckoff_site(sg, site, setting_id)
        if match:
            direct += 1
        else:
            match = _match_inter_wyckoff_site_with_symmetry(
                local_data, sg, site, operations, setting_id
            )
        matches.append(match)
        if match:
            matched += 1
            constrained += max(0, 3 - len(match.get("free") or []))
    default = local_data.default_inter_setting_record(sg)
    default_id = int(default["id"])
    default_bonus = 1 if int(setting_id) == default_id else 0
    setting = local_data.inter_setting_record(setting_id)
    same_orientation = 1 if all(
        str(setting.get(key, "")) == str(default.get(key, ""))
        for key in ("cell", "axis", "abc")
    ) else 0
    return (
        matched,
        constrained,
        same_orientation,
        direct,
        default_bonus,
    ), tuple(matches)


def _select_parent_inter_setting(
    local_data: SourceTables,
    sg: int,
    atom_sites: list[dict[str, Any]],
    operations: list[str],
) -> tuple[int, tuple[dict[str, Any] | None, ...]]:
    candidates = local_data.inter_setting_ids_for_space_group(sg)
    if not candidates:
        candidates = (int(local_data.default_inter_setting_record(sg)["id"]),)
    scored = []
    for setting_id in candidates:
        score, matches = _match_parent_inter_setting(
            local_data, sg, atom_sites, operations, setting_id
        )
        scored.append((score, int(setting_id), matches))
    scored.sort(key=lambda item: (item[0], -int(item[1])), reverse=True)
    return int(scored[0][1]), scored[0][2]


def build_parent_state_from_cif_info(cif_info: dict[str, Any]) -> dict[str, Any]:
    sg = int(cif_info["parent"]["number"])
    local_data = source_tables()
    cif_info = copy.deepcopy(cif_info)
    atom_sites = [site for site in cif_info.get("atom_sites", []) if isinstance(site, dict)]
    parent_inter_setting_id, parent_site_matches = _select_parent_inter_setting(
        local_data,
        sg,
        atom_sites,
        cif_info.get("symmetry_operations", []),
    )
    for site, selected_match in zip(
        cif_info.get("atom_sites", []), parent_site_matches, strict=True
    ):
        # CIFs variously report a Wyckoff letter ("b") or a composed label
        # ("4b").  Neither spelling is the site identity used downstream:
        # establish that identity from the selected-setting coordinates.
        site.pop("wyckoff", None)
        site.pop("wyckoff_row_id", None)
        site.pop("wyckoff_formula", None)
        site.pop("wyckoff_params", None)
        site.pop("wyckoff_multiplicity", None)
        match = selected_match
        if not match:
            match = _match_wyckoff_site_with_symmetry(
                local_data,
                sg,
                site,
                cif_info.get("symmetry_operations", []),
            )
        if match:
            canonical_mapping = _canonical_representative_mapping(
                site,
                match,
                cif_info.get("symmetry_operations", []),
            )
            if canonical_mapping:
                match.update(canonical_mapping)
            site["wyckoff"] = match["label"]
            site["wyckoff_row_id"] = match["row_id"]
            site["wyckoff_formula"] = match["formula"]
            site["wyckoff_params"] = match["params"]
            site["wyckoff_multiplicity"] = match["multiplicity"]
            for key in (
                "representative_source_fract",
                "representative_matched_fract",
                "representative_operation",
                "representative_operation_index",
                "representative_lattice_translation",
            ):
                if key in match:
                    site[key] = match[key]

    if not str(cif_info.get("formula") or "").strip():
        cif_info["formula"] = _formula_from_sites(cif_info.get("atom_sites", []))

    default_setting = local_data.default_inter_setting_record(sg)
    parent_inter_setting = local_data.inter_setting_record(parent_inter_setting_id)
    parent_preference = str(cif_info.get("explicit_parent_preference") or "").strip()
    if not parent_preference:
        parent_preference = local_data.parent_preference_text(parent_inter_setting)
    setting = {
        "default_preferences": local_data.default_space_group_preferences(),
        "default_inter_setting": default_setting,
        "parent_preferences": parent_preference or None,
        "parent_inter_setting": parent_inter_setting,
    }
    matched_sites = [site for site in cif_info.get("atom_sites", []) if site.get("wyckoff_row_id")]
    displacive_row_ids = sorted({
        int(site["wyckoff_row_id"])
        for site in matched_sites
    }) if len(matched_sites) == len(cif_info.get("atom_sites", [])) else []
    lattice = int(local_data.space["ispace_lattice"][sg - 1])
    return {
        "schema": "isodistort.assembled.parent.v1",
        "input": cif_info,
        "space_group": {
            "number": sg,
            "symbol": cif_info.get("parent", {}).get("symbol", ""),
            "lattice": lattice,
            "lattice_type": str(local_data.space["lattice_label"][lattice - 1]).strip(),
        },
        "setting": setting,
        "_internal": {
            "parent_inter_setting_id": int(parent_inter_setting_id),
            "displacive_row_ids": displacive_row_ids,
        },
    }


def build_parent_state(cif: Path) -> dict[str, Any]:
    return build_parent_state_from_cif_info(read_cif_summary(cif))


def build_parent_state_from_text(cif_text: str) -> dict[str, Any]:
    return build_parent_state_from_cif_info(read_cif_summary_from_text(cif_text))
