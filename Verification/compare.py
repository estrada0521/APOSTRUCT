"""The actual Local-vs-Web comparison, independent of any content-addressed store.

Extracted out of Validation/core/validate.py:judge_local_output() so there is
exactly one implementation of "how do we turn a Local preview and a Web text
export into an opd/comparison verdict pair." judge_local_output() unwraps
this monorepo's own provenance-checked ValidationCase/PairedOracleCase/
LocalArtifact objects and calls this; Verification/verify_comparison.py builds
the same plain arguments directly from an `apo` CLI output and a saved Web
text export and calls this too. Neither caller reimplements the comparison
itself — this is the one place it lives. It depends only on the standalone
Verification package and the explicitly public product kernels used by its
physical-equivalence checks.
"""

from __future__ import annotations

from typing import Any, Mapping, Sequence

from Verification.comparison.mode import compare_mode_details
from Verification.comparison.opd import compare_opd, local_opd
from Verification.comparison.selected_state import compare_selected_state
from Verification.parsers.complete_mode_text import (
    parse_complete_mode_details_from_text,
)


def compare_local_to_web(
    *,
    opd_exists: bool,
    preview: Mapping[str, Any] | None,
    mode_text: str,
    opd_label: str,
    k_selections: Sequence[Mapping[str, Any]],
    case_label: str,
    cif_label: str,
    comparison_source: str,
    tol: float,
) -> tuple[dict[str, Any], dict[str, Any] | None]:
    """Return (opd, comparison) — exactly what ComputedCase.opd/.comparison hold."""

    if opd_exists is False:
        if preview is not None:
            raise ValueError("absent local OPD unexpectedly has a preview")
        return (
            {
                "exist": False,
                "level": "absent",
                "source": "compute_pairing",
                "label": opd_label,
                "physical": None,
                "strict": None,
            },
            None,
        )
    if opd_exists is not True or not isinstance(preview, Mapping):
        raise ValueError("present local OPD has no structured preview")
    selected = (preview.get("selected") or {}).get("orderparam")
    if not isinstance(selected, Mapping):
        raise ValueError("local output has no selected OPD payload")
    web = parse_complete_mode_details_from_text(mode_text)
    if web.subgroup is None:
        raise ValueError("Web mode-details output has no subgroup details")
    subgroup = web.subgroup
    web_opd = {
        "label": opd_label,
        "child_sg": subgroup.display_label or subgroup.number,
        "symbol": subgroup.symbol,
        "basis": subgroup.basis,
        "origin": subgroup.origin,
        "s": subgroup.s,
        "i": subgroup.i,
        "k_active": subgroup.k_active,
    }
    magnetic = "." in str(web_opd["child_sg"] or "")
    isotropy = selected.get("isotropy")
    if not isinstance(isotropy, Mapping):
        isotropy = {}
    local_subgroup = isotropy.get("subgroup")
    if not isinstance(local_subgroup, Mapping):
        local_subgroup = {}
    selected_payload = preview.get("selected")
    if not isinstance(selected_payload, Mapping):
        selected_payload = {}
    mode_details = selected_payload.get("mode_details")
    if not isinstance(mode_details, Mapping):
        mode_details = {}
    subgroup_details = mode_details.get("subgroup_details")
    if not isinstance(subgroup_details, Mapping):
        subgroup_details = {}
    opd = compare_opd(
        web_opd,
        local_opd(selected, magnetic=magnetic),
        source=comparison_source,
        direction_matched_by_label=True,
        parent_sg=(web.parent.number if web.parent is not None else None),
        parent_setting_id=subgroup_details.get("parent_inter_setting_id"),
        subgroup_setting_id=None,
        child_ordinary_sg=(
            local_subgroup.get("ordinary_number") or local_subgroup.get("number")
        ),
        child_magnetic_group=(local_subgroup.get("number") if magnetic else None),
        local_operation_records=isotropy.get("source_operation_records"),
    )
    selected_state = compare_selected_state(
        preview,
        selections=k_selections,
        web_subgroup={
            "number": subgroup.number,
            "basis": subgroup.basis,
            "origin": subgroup.origin,
        },
        operation_set_equivalent=opd.get("origin_equivalence") == "operation_set",
    )
    comparison = compare_mode_details(
        case=case_label,
        cif=cif_label,
        web=web,
        preview=preview,
        tol=tol,
        meta={"opd": opd_label, "opd_comparison": opd, "selected_state": selected_state},
    )
    return opd, comparison


__all__ = ["compare_local_to_web"]
