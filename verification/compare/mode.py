"""Assemble the formal complete-mode comparison verdict."""

from __future__ import annotations

from pathlib import Path
from typing import Any

from verification.compare.structure import (
    _atom_payload,
    _basis_equivalent_atom_payload,
    _basis_equivalent_lattice_payload,
    _lattice_payload,
    _parent_atom_payload,
    _parent_lattice_payload,
    _spacegroup_payload,
    _subgroup_payload,
)
from verification.compare.basis_mode import basis_equivalent_mode_frame
from verification.compare.definitions import Definitions
from verification.compare.structure import (
    _exact_integer_value,
    _validated_basis_equivalent_atom_result,
    _validated_basis_equivalent_lattice_result,
)
from verification.web.text import CompleteModeDetails


def _input_standardization_payload(preview: dict[str, Any]) -> dict[str, Any]:
    sites = ((preview.get("input") or {}).get("atom_sites") or [])
    rows = []
    for site in sites:
        if not isinstance(site, dict):
            continue
        rows.append(
            {
                "label": site.get("label"),
                "type": site.get("type"),
                "fract": site.get("fract"),
                "multiplicity": site.get("multiplicity"),
                "wyckoff": site.get("wyckoff"),
                "wyckoff_row_id": site.get("wyckoff_row_id"),
            }
        )
    mapped = sum(1 for row in rows if row.get("wyckoff_row_id"))
    return {
        "status": "ok" if mapped == len(rows) else "partial",
        "mapped": mapped,
        "total": len(rows),
        "rows": rows,
    }


def _has_complete_mode_payload(web: CompleteModeDetails, definitions: Definitions) -> bool:
    return any(
        (
            web.parent is not None,
            web.subgroup is not None,
            web.parent_lattice is not None,
            web.undistorted_lattice is not None,
            bool(web.parent_atoms),
            bool(web.undistorted_atoms),
            any(definitions.web_count(kind) for kind in ("dsp", "mag", "occ", "strain")),
        )
    )


_SECTION_NAMES = {"dsp": "displacive", "mag": "magnetic", "occ": "occupancy", "strain": "strain"}


def _mode_frame(
    *,
    selected_state: object,
    subgroup: dict[str, Any],
    lattice: dict[str, Any],
    atoms: dict[str, Any],
):
    """The request child frame of a closed equivalent-frame proof, or None."""

    proof_result = {
        "selected_state": selected_state,
        "subgroup": subgroup,
        "undistorted_lattice": lattice,
    }
    if not _validated_basis_equivalent_lattice_result(proof_result) or not isinstance(selected_state, dict):
        return None
    request = selected_state.get("request")
    local_selected = selected_state.get("local")
    if (
        not isinstance(request, (tuple, list))
        or not isinstance(local_selected, (tuple, list))
        or len(request) != 4
        or len(local_selected) != 4
        or request[0] != local_selected[0]
        or (
            request[3] != local_selected[3]
            and selected_state.get("status") != "equivalent_embedding"
        )
    ):
        return None
    request_child = _exact_integer_value(request[1])
    local_child = _exact_integer_value(local_selected[1])
    if request_child is None or local_child is None:
        return None
    if not _validated_basis_equivalent_atom_result(
        atoms,
        basis_change=lattice.get("basis_change"),
        child_space_group=request_child,
    ):
        return None
    return basis_equivalent_mode_frame(
        basis_change=lattice.get("basis_change"),
        request_origin=request[3],
        local_origin=local_selected[3],
        request_child_space_group=request_child,
        local_child_space_group=local_child,
        request_basis=request[2],
    )


def compare_mode_details(
    *,
    case: str,
    cif: Path,
    web: CompleteModeDetails,
    preview: dict[str, Any],
    definitions: Definitions,
    tol: float,
    meta: dict[str, Any] | None = None,
) -> dict[str, Any]:
    input_standardization = _input_standardization_payload(preview)
    if not _has_complete_mode_payload(web, definitions):
        return {
            "case": case,
            "status": "context_mismatch",
            "reason": "saved Web FINAL contains no complete-mode payload",
            "cif": str(cif),
            "input_standardization": input_standardization,
            **(meta or {}),
        }
    local = ((preview.get("selected") or {}).get("mode_details") or {})
    if local.get("status") in {"unsupported", "missing", "error"}:
        return {
            "case": case,
            "status": "unsupported",
            "reason": local.get("reason"),
            "mode_detail_status": local.get("status"),
            "cif": str(cif),
            "input_standardization": input_standardization,
            **(meta or {}),
        }
    parent = _spacegroup_payload(web, preview)
    parent_lattice = _parent_lattice_payload(web, preview, tol)
    parent_atoms = _parent_atom_payload(web, preview, tol)
    opd_comparison = (meta or {}).get("opd_comparison")
    operation_set_equivalent = (
        isinstance(opd_comparison, dict)
        and opd_comparison.get("origin_equivalence") == "operation_set"
    )
    subgroup = _subgroup_payload(
        web,
        preview,
        operation_set_equivalent=operation_set_equivalent,
    )
    selected_state = (meta or {}).get("selected_state")
    lattice = _lattice_payload(web, local, tol)
    transported_lattice = _basis_equivalent_lattice_payload(
        web,
        local,
        parent_lattice,
        selected_state,
        subgroup,
        tol,
    )
    if transported_lattice is not None:
        lattice = transported_lattice
    atoms = _atom_payload(web, local, tol)
    accepted_atom_statuses = {"ok", "site_gauge"}
    basis_transport_context = (
        isinstance(selected_state, dict)
        and selected_state.get("status")
        in {"equivalent_basis", "equivalent_embedding"}
        and subgroup.get("status")
        in {"basis_equivalent", "embedding_equivalent"}
        and lattice.get("status") == "basis_equivalent"
    )
    transported_atoms = None
    if atoms.get("status") not in accepted_atom_statuses or basis_transport_context:
        transported_atoms = _basis_equivalent_atom_payload(
            web,
            local,
            selected_state,
            subgroup,
            lattice,
            tol,
        )
        if (
            atoms.get("status") not in accepted_atom_statuses
            and transported_atoms is not None
        ):
            atoms = transported_atoms
    mode_transport_atoms = transported_atoms or atoms
    frame = None
    for kind in ("dsp", "mag"):
        if definitions.accepted(kind):
            continue
        if frame is None:
            frame = _mode_frame(
                selected_state=selected_state,
                subgroup=subgroup,
                lattice=lattice,
                atoms=mode_transport_atoms,
            ) or ()
        if frame:
            definitions.transport(kind, *frame)
    if isinstance(selected_state, dict) and selected_state.get("status") == "diff":
        status = "context_mismatch"
    elif input_standardization["status"] != "ok":
        status = "input_standardization"
    elif parent["status"] != "ok":
        status = "parent"
    elif parent_lattice["status"] != "ok":
        status = "parent_lattice"
    elif parent_atoms["status"] != "ok":
        status = "parent_atoms"
    elif subgroup["status"] not in {
        "ok",
        "basis_equivalent",
        "embedding_equivalent",
    }:
        status = "subgroup"
    elif lattice["status"] not in {"ok", "basis_equivalent"}:
        status = "undistorted_lattice"
    elif atoms["status"] not in accepted_atom_statuses:
        status = "undistorted_atoms"
    else:
        failed_kind = next(
            (kind for kind in ("dsp", "mag", "occ", "strain") if definitions.payload(kind)["strict_status"] != "ok"),
            None,
        )
        status = (
            "ok"
            if failed_kind is None
            else f"definitions_{_SECTION_NAMES[failed_kind]}_{definitions.payload(failed_kind)['strict_status']}"
        )
    result = {
        "case": case,
        "status": status,
        "cif": str(cif),
        "input_standardization": input_standardization,
        "parent": parent,
        "parent_lattice": parent_lattice,
        "parent_atoms": parent_atoms,
        "subgroup": subgroup,
        "undistorted_lattice": lattice,
        "undistorted_atoms": atoms,
        **(meta or {}),
    }
    result["definitions"] = definitions
    return result


__all__ = ["compare_mode_details"]
