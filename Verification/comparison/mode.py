"""Assemble the formal complete-mode comparison verdict."""

from __future__ import annotations

from pathlib import Path
from typing import Any

from Verification.comparison.structure import (
    _atom_payload,
    _basis_equivalent_atom_payload,
    _basis_equivalent_lattice_payload,
    _lattice_payload,
    _parent_atom_payload,
    _parent_lattice_payload,
    _spacegroup_payload,
    _subgroup_payload,
)
from Verification.comparison.vectors import (
    _basis_equivalent_mode_payload,
    _strain_modes_payload,
    _vector_modes_probe,
    _vector_modes_payload,
)
from Verification.parsers.complete_mode_text import CompleteModeDetails


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


def _has_complete_mode_payload(web: CompleteModeDetails) -> bool:
    return any(
        (
            web.parent is not None,
            web.subgroup is not None,
            web.parent_lattice is not None,
            web.undistorted_lattice is not None,
            bool(web.parent_atoms),
            bool(web.undistorted_atoms),
            bool(web.displacive_definitions),
            bool(web.magnetic_definitions),
            bool(web.strain_definitions),
        )
    )


def _local_modes_in_exact_web_order(
    web_modes: tuple[Any, ...],
    local_modes: list[dict[str, Any]],
) -> list[dict[str, Any]]:
    """Use complete unique labels as the row-comparison identity."""

    if len(web_modes) != len(local_modes):
        return local_modes
    web_labels = tuple(getattr(mode, "label", None) for mode in web_modes)
    local_labels = tuple(mode.get("label") for mode in local_modes)
    if (
        any(not isinstance(label, str) or not label for label in web_labels)
        or any(not isinstance(label, str) or not label for label in local_labels)
        or len(set(web_labels)) != len(web_labels)
        or len(set(local_labels)) != len(local_labels)
        or set(web_labels) != set(local_labels)
    ):
        return local_modes
    by_label = dict(zip(local_labels, local_modes, strict=True))
    return [by_label[label] for label in web_labels]


def compare_mode_details(
    *,
    case: str,
    cif: Path,
    web: CompleteModeDetails,
    preview: dict[str, Any],
    tol: float,
    meta: dict[str, Any] | None = None,
    structure_only: bool = False,
) -> dict[str, Any]:
    input_standardization = _input_standardization_payload(preview)
    if not _has_complete_mode_payload(web):
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
    if atoms.get("status") not in accepted_atom_statuses or (
        not structure_only and basis_transport_context
    ):
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
    accepted_mode_statuses = {"ok", "sign", "signed_permutation", "family_basis"}
    local_displacive = _local_modes_in_exact_web_order(
        web.displacive_definitions,
        list(local.get("displacive_definitions") or []),
    )
    local_magnetic = _local_modes_in_exact_web_order(
        web.magnetic_definitions,
        list(local.get("magnetic_definitions") or []),
    )
    use_transport_probe = (
        not structure_only
        and basis_transport_context
    )
    modes_probe = (
        _vector_modes_probe(
            web.displacive_definitions,
            local_displacive,
            tol,
            vector_kind="dsp",
        )
        if use_transport_probe
        else None
    )
    magnetic_probe = (
        _vector_modes_probe(
            web.magnetic_definitions,
            local_magnetic,
            tol,
            vector_kind="mag",
        )
        if use_transport_probe
        else None
    )
    modes = None if structure_only else (
        modes_probe
        if modes_probe is not None
        and modes_probe.get("status") not in accepted_mode_statuses
        else _vector_modes_payload(
            web.displacive_definitions,
            local_displacive,
            tol,
            vector_kind="dsp",
        )
    )
    magnetic_modes = None if structure_only else (
        magnetic_probe
        if magnetic_probe is not None
        and magnetic_probe.get("status") not in accepted_mode_statuses
        else _vector_modes_payload(
            web.magnetic_definitions,
            local_magnetic,
            tol,
            vector_kind="mag",
        )
    )
    if modes is not None and modes.get("status") not in accepted_mode_statuses:
        transported_modes = _basis_equivalent_mode_payload(
            web.displacive_definitions,
            local_displacive,
            vector_kind="dsp",
            selected_state=selected_state,
            subgroup=subgroup,
            lattice=lattice,
            atoms=mode_transport_atoms,
            tol=tol,
        )
        if transported_modes is not None:
            modes = transported_modes
        elif modes_probe is not None:
            modes = _vector_modes_payload(
                web.displacive_definitions,
                local_displacive,
                tol,
                vector_kind="dsp",
            )
    if (
        magnetic_modes is not None
        and magnetic_modes.get("status") not in accepted_mode_statuses
    ):
        transported_magnetic_modes = _basis_equivalent_mode_payload(
            web.magnetic_definitions,
            local_magnetic,
            vector_kind="mag",
            selected_state=selected_state,
            subgroup=subgroup,
            lattice=lattice,
            atoms=mode_transport_atoms,
            tol=tol,
        )
        if transported_magnetic_modes is not None:
            magnetic_modes = transported_magnetic_modes
        elif magnetic_probe is not None:
            magnetic_modes = _vector_modes_payload(
                web.magnetic_definitions,
                local_magnetic,
                tol,
                vector_kind="mag",
            )
    strain_modes = None if structure_only else _strain_modes_payload(
        web.strain_definitions,
        list(local.get("strain_definitions") or []),
        tol,
    )
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
    elif structure_only:
        status = "ok"
    else:
        assert modes is not None and magnetic_modes is not None and strain_modes is not None
        kind_statuses = (
            ("displacive", modes),
            ("magnetic", magnetic_modes),
            ("strain", strain_modes),
        )
        failed_kind = next((kind for kind, payload in kind_statuses if payload["strict_status"] != "ok"), None)
        status = "ok" if failed_kind is None else f"definitions_{failed_kind}_{dict(kind_statuses)[failed_kind]['strict_status']}"
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
    if modes is not None:
        result["definitions"] = modes
        result["magnetic_definitions"] = magnetic_modes
        result["strain_definitions"] = strain_modes
    return result


__all__ = ["compare_mode_details"]
