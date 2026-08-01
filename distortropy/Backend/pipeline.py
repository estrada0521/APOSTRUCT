"""Compose parent, reciprocal-space, isotropy, and mode-domain services.

Frontend code calls this boundary and renders the returned state. It does not
recompute domain state.
"""

from __future__ import annotations

from typing import Any

from distortropy.Backend.modes import (  # noqa: E402
    build_mode_details,
)
from distortropy.Backend.reciprocal.k_points import normalize_k_params  # noqa: E402
from distortropy.Backend.reciprocal.irreps import select_irreps  # noqa: E402
from distortropy.Backend.reciprocal.k_points import select_k_points  # noqa: E402
from distortropy.Backend.isotropy.selection import (  # noqa: E402
    select_embedding_direction,
    select_order_parameter_direction,
)
from distortropy.Backend.parent import (  # noqa: E402
    build_generic_parent_state,
    build_parent_state_from_cif_info,
    read_cif_summary_from_text,
    realize_generic_cif_info,
)


def _site_selection(
    requested: list[int] | str | None,
    *,
    site_count: int,
    default_all: bool,
) -> list[int]:
    if requested is None:
        return list(range(1, site_count + 1)) if default_all else []
    if requested == "all":
        return list(range(1, site_count + 1))
    if requested == "none":
        return []
    selected = sorted({int(value) for value in requested})
    invalid = [value for value in selected if value < 1 or value > site_count]
    if invalid:
        raise ValueError(f"site indexes out of range 1..{site_count}: {invalid}")
    return selected


def _apply_distortion_selection(
    parent_state: dict[str, Any],
    *,
    mode_set: str,
    distortion_selection: dict[str, Any] | None,
) -> dict[str, Any]:
    atoms = list((parent_state.get("input") or {}).get("atom_sites") or [])
    site_types = [
        str(atom.get("type") or atom.get("label") or index)
        for index, atom in enumerate(atoms, start=1)
    ]
    site_labels = [
        str(atom.get("label") or atom.get("type") or index)
        for index, atom in enumerate(atoms, start=1)
    ]
    selectable_sites = list(range(1, len(atoms) + 1))
    requested = distortion_selection or {}
    include_strain = bool(requested.get("strain", mode_set in {"strain_displacive", "all"}))
    displacive_sites = (
        list(selectable_sites)
        if requested.get("displacive_sites") is None or requested.get("displacive_sites") == "all"
        else _site_selection(
            requested.get("displacive_sites"),
            site_count=len(atoms),
            default_all=True,
        )
    )
    magnetic_sites = (
        list(selectable_sites) if requested.get("magnetic_sites") == "all"
        else list(selectable_sites) if requested.get("magnetic_sites") is None and mode_set in {"magnetic", "all"}
        else [] if requested.get("magnetic_sites") is None
        else _site_selection(
            requested.get("magnetic_sites"),
            site_count=len(atoms),
            default_all=False,
        )
    )
    for name, selected in (("displacive", displacive_sites), ("magnetic", magnetic_sites)):
        invalid = [index for index in selected if index not in selectable_sites]
        if invalid:
            raise ValueError(f"{name} site indexes are not selectable atom sites: {invalid}")
    if not include_strain and not displacive_sites and not magnetic_sites:
        raise ValueError("You must include some kind of distortion")

    def row_ids(site_indexes: list[int]) -> list[int]:
        return sorted({
            int(atoms[index - 1]["wyckoff_row_id"])
            for index in site_indexes
            if atoms[index - 1].get("wyckoff_row_id")
        })

    selection = {
        "strain": include_strain,
        "displacive_sites": displacive_sites,
        "magnetic_sites": magnetic_sites,
        "selectable_sites": selectable_sites,
        "site_types": site_types,
        "site_labels": site_labels,
    }
    return {
        **parent_state,
        "mode_set": mode_set,
        "distortion_selection": selection,
        "_internal": {
            **(parent_state.get("_internal") or {}),
            "mode_set": mode_set,
            "distortion_selection": selection,
            "displacive_row_ids": row_ids(displacive_sites),
            "magnetic_row_ids": row_ids(magnetic_sites),
            "displacive_site_indexes": displacive_sites,
            "magnetic_site_indexes": magnetic_sites,
        },
    }


def attach_mode_details(
    state: dict[str, Any],
    *,
    k_params: dict[str, str],
    include_mode_structure: bool = True,
) -> dict[str, Any]:
    selected = state.get("selected") or {}
    mode_cif_info = state["input"]
    try:
        mode_cif_info = realize_generic_cif_info(mode_cif_info)
        mode_details = build_mode_details(
            cif_info=mode_cif_info,
            selected_k=selected.get("kpoint"),
            selected_irrep=selected.get("irrep"),
            selected_opd=selected.get("orderparam"),
            selected_slots=selected.get("selection_slots"),
            selected_factors=selected.get("opd_factors"),
            k_params=k_params,
            include_strain=bool((state.get("distortion_selection") or {}).get("strain")),
            include_magnetic=bool((state.get("distortion_selection") or {}).get("magnetic_sites")),
            displacive_row_ids=(state.get("_internal") or {}).get("displacive_row_ids"),
            magnetic_row_ids=(state.get("_internal") or {}).get("magnetic_row_ids"),
            displacive_site_indexes=(state.get("_internal") or {}).get("displacive_site_indexes"),
            magnetic_site_indexes=(state.get("_internal") or {}).get("magnetic_site_indexes"),
            parent_inter_setting_id=(state.get("_internal") or {}).get("parent_inter_setting_id"),
            include_structure=include_mode_structure,
        )
    except Exception as exc:
        mode_details = {
            "status": "error",
            "reason": f"{type(exc).__name__}: {exc}",
        }
    out = {**state, "input": mode_cif_info}
    out["selected"] = {**selected, "mode_details": mode_details}
    return out


def build_state_from_cif_info(
    cif_info: dict[str, Any],
    *,
    mode_set: str = "strain_displacive",
    distortion_selection: dict[str, Any] | None = None,
    k_label: str | None = "GM",
    k_index: int = 1,
    k_params: dict[str, str] | None = None,
    display_k_params: dict[str, str] | None = None,
    irrep: str | None = None,
    irrep_index: int = 1,
    orderparam: int | str | None = None,
    orderparam_index: int = 1,
    selections: list[dict[str, Any]] | None = None,
    include_opd: bool = True,
    include_mode_details: bool = False,
    include_mode_structure: bool = True,
    opd_row_limit: int | None = None,
    selected_opd_only: bool = False,
    direction_selection: tuple[dict[str, Any], int] | None = None,
) -> dict[str, Any]:
    k_params = normalize_k_params(k_params)
    display_k_params = normalize_k_params(display_k_params) or k_params
    return build_state_from_parent_state(
        build_parent_state_from_cif_info(cif_info),
        mode_set=mode_set,
        distortion_selection=distortion_selection,
        k_label=k_label,
        k_index=k_index,
        k_params=k_params,
        display_k_params=display_k_params,
        irrep=irrep,
        irrep_index=irrep_index,
        orderparam=orderparam,
        orderparam_index=orderparam_index,
        selections=selections,
        include_opd=include_opd,
        include_mode_details=include_mode_details,
        include_mode_structure=include_mode_structure,
        opd_row_limit=opd_row_limit,
        selected_opd_only=selected_opd_only,
        direction_selection=direction_selection,
    )


def build_state_from_generic(
    space_group: int,
    sites: list[Any] | tuple[Any, ...] = (),
    **options: Any,
) -> dict[str, Any]:
    """Run the shared pipeline from symbolic space-group/Wyckoff input."""

    return build_state_from_parent_state(
        build_generic_parent_state(space_group, sites),
        **options,
    )


def build_state_from_parent_state(
    parent_state: dict[str, Any],
    *,
    mode_set: str = "strain_displacive",
    distortion_selection: dict[str, Any] | None = None,
    k_label: str | None = "GM",
    k_index: int = 1,
    k_params: dict[str, str] | None = None,
    display_k_params: dict[str, str] | None = None,
    irrep: str | None = None,
    irrep_index: int = 1,
    orderparam: int | str | None = None,
    orderparam_index: int = 1,
    selections: list[dict[str, Any]] | None = None,
    include_opd: bool = True,
    include_mode_details: bool = False,
    include_mode_structure: bool = True,
    opd_row_limit: int | None = None,
    selected_opd_only: bool = False,
    direction_selection: tuple[dict[str, Any], int] | None = None,
) -> dict[str, Any]:
    k_params = normalize_k_params(k_params)
    display_k_params = normalize_k_params(display_k_params) or k_params
    parent_state = _apply_distortion_selection(
        parent_state,
        mode_set=mode_set,
        distortion_selection=distortion_selection,
    )
    k_state = select_k_points(
        parent_state,
        k_label=k_label,
        k_index=k_index,
        k_params=k_params,
        display_k_params=display_k_params,
        selections=selections,
    )
    reciprocal_state = select_irreps(k_state)
    if not include_opd:
        return reciprocal_state
    isotropy_state = select_order_parameter_direction(
        reciprocal_state,
        k_params=k_params,
        display_k_params=display_k_params,
        irrep=irrep,
        irrep_index=irrep_index,
        orderparam=orderparam,
        orderparam_index=orderparam_index,
        selections=selections,
        opd_row_limit=opd_row_limit,
        selected_opd_only=selected_opd_only,
    )
    if direction_selection is not None:
        isotropy_state = select_embedding_direction(
            isotropy_state,
            direction_selection[0],
            direction_selection[1],
        )
    if include_mode_details:
        return attach_mode_details(
            isotropy_state,
            k_params=k_params,
            include_mode_structure=include_mode_structure,
        )
    return isotropy_state


def build_state_from_text(
    cif_text: str,
    *,
    mode_set: str = "strain_displacive",
    distortion_selection: dict[str, Any] | None = None,
    k_label: str | None = "GM",
    k_index: int = 1,
    k_params: dict[str, str] | None = None,
    display_k_params: dict[str, str] | None = None,
    irrep: str | None = None,
    irrep_index: int = 1,
    orderparam: int | str | None = None,
    orderparam_index: int = 1,
    selections: list[dict[str, Any]] | None = None,
    include_opd: bool = True,
    include_mode_details: bool = False,
    include_mode_structure: bool = True,
    opd_row_limit: int | None = None,
    selected_opd_only: bool = False,
    direction_selection: tuple[dict[str, Any], int] | None = None,
) -> dict[str, Any]:
    return build_state_from_cif_info(
        read_cif_summary_from_text(cif_text),
        mode_set=mode_set,
        distortion_selection=distortion_selection,
        k_label=k_label,
        k_index=k_index,
        k_params=k_params,
        display_k_params=display_k_params,
        irrep=irrep,
        irrep_index=irrep_index,
        orderparam=orderparam,
        orderparam_index=orderparam_index,
        selections=selections,
        include_opd=include_opd,
        include_mode_details=include_mode_details,
        include_mode_structure=include_mode_structure,
        opd_row_limit=opd_row_limit,
        selected_opd_only=selected_opd_only,
        direction_selection=direction_selection,
    )
