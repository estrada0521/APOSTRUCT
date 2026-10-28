"""Project raw comparison payloads onto the formal validation vocabulary."""

from __future__ import annotations

from typing import Any, Mapping


COMPARATOR_VERSION = "validation.v10"
DEFAULT_TOLERANCE = 5e-5


def _accepted(payload: object) -> bool:
    return isinstance(payload, Mapping) and payload.get("status") in {
        "ok",
        "basis_equivalent",
    }


def _accepted_undistorted_atoms(payload: object) -> bool:
    return isinstance(payload, Mapping) and payload.get("status") in {
        "ok",
        "basis_equivalent",
        "site_gauge",
    }


def _failure_detail(payload: object) -> dict[str, Any]:
    if not isinstance(payload, Mapping):
        return {"status": "missing"}
    result: dict[str, Any] = {"status": str(payload.get("status") or "missing")}
    matches = payload.get("matches")
    if isinstance(matches, Mapping):
        failed = sorted(str(field) for field, value in matches.items() if value is not True)
        if failed:
            result["fields"] = failed
    for field in ("web_count", "local_count", "missing", "extra"):
        value = payload.get(field)
        if value not in (None, [], {}):
            result[field] = value
    return result


def _mode_result(
    comparison: Mapping[str, Any],
) -> tuple[dict[str, Any], str | None, bool, bool]:
    parent_pass = (
        _accepted(comparison.get("parent"))
        and _accepted(comparison.get("parent_lattice"))
        and (comparison.get("input_standardization") or {}).get("status") == "ok"
    )
    atom_pass = _accepted(comparison.get("parent_atoms"))
    undistorted_pass = _accepted(
        comparison.get("undistorted_lattice")
    ) and _accepted_undistorted_atoms(comparison.get("undistorted_atoms"))
    parent_structure: dict[str, Any] = {"pass": parent_pass}
    if not parent_pass:
        failures: dict[str, Any] = {}
        if not _accepted(comparison.get("parent")):
            failures["sg"] = _failure_detail(comparison.get("parent"))
        if not _accepted(comparison.get("parent_lattice")):
            failures["lattice"] = _failure_detail(comparison.get("parent_lattice"))
        standardization = comparison.get("input_standardization")
        if not isinstance(standardization, Mapping) or standardization.get("status") != "ok":
            failures["standardization"] = _failure_detail(standardization)
        parent_structure["failures"] = failures

    atom_site: dict[str, Any] = {"pass": atom_pass}
    if not atom_pass:
        atom_site["failure"] = _failure_detail(comparison.get("parent_atoms"))

    atom_site_gauge = (
        isinstance(comparison.get("undistorted_atoms"), Mapping)
        and comparison["undistorted_atoms"].get("status") == "site_gauge"
    )
    undistorted_structure: dict[str, Any] = {"pass": undistorted_pass}
    if atom_site_gauge:
        undistorted_structure["strict_pass"] = False
        undistorted_structure["gauge"] = ["wyckoff_label"]
    if not undistorted_pass:
        failures = {}
        if not _accepted(comparison.get("undistorted_lattice")):
            failures["lattice"] = _failure_detail(
                comparison.get("undistorted_lattice")
            )
        if not _accepted_undistorted_atoms(comparison.get("undistorted_atoms")):
            failures["atoms"] = _failure_detail(comparison.get("undistorted_atoms"))
        undistorted_structure["failures"] = failures

    mode: dict[str, Any] = {
        "parent_structure": parent_structure,
        "atom_site": atom_site,
        "undistorted_structure": undistorted_structure,
        "mode_definition": None,
    }
    if not parent_pass:
        return mode, "mode.parent_structure", False, False
    if not atom_pass:
        return mode, "mode.atom_site", False, False
    if not undistorted_pass:
        return mode, "mode.undistorted_structure", False, False
    definition, divergence, physical, strict = comparison["definitions"].mode_definition()
    mode["mode_definition"] = definition
    if atom_site_gauge and physical:
        return (
            mode,
            "mode.undistorted_structure.site_gauge",
            physical,
            False,
        )
    if atom_site_gauge:
        strict = False
    return mode, divergence, physical, strict


def _opd_divergence(opd: Mapping[str, Any]) -> str | None:
    if opd.get("exist") is not True:
        return "opd.absent"
    for field in (
        "child_sg",
        "k_active",
        "dir",
        "sublattice",
        "origin_coset",
        "s",
        "i",
    ):
        if (opd.get("physical") or {}).get(field) is not True:
            return f"opd.physical.{field}"
    for field in ("basis", "origin", "k_active"):
        if (opd.get("strict") or {}).get(field) is not True:
            return f"opd.strict.{field}"
    return None


def _opd_result(opd: Mapping[str, Any]) -> dict[str, Any]:
    exists = opd.get("exist") is True
    result: dict[str, Any] = {
        "exist": exists,
        "level": str(opd.get("level") or ("present" if exists else "absent")),
        "source": opd.get("source"),
        "physical": None,
        "strict": None,
    }
    if exists:
        result["physical"] = dict(opd.get("physical") or {})
        result["strict"] = dict(opd.get("strict") or {})
    return result


VERDICTS = ("identical", "equivalent", "unproven")
_FORMER_VERDICTS = {"strict_pass": "identical", "physical_pass": "equivalent", "physical_fail": "unproven"}


def canonical_verdict(verdict: object) -> str:
    """The verdict name, reading records written before the rename; a missing verdict is an error."""

    if verdict is None:
        return "error"
    return _FORMER_VERDICTS.get(str(verdict), str(verdict))


def result_for_comparison(
    opd: Mapping[str, Any],
    comparison: Mapping[str, Any] | None,
    *,
    tol: float,
) -> dict[str, Any]:
    """Project one explicit comparison onto the scientific verdict vocabulary."""

    base = {
        "comparator_version": COMPARATOR_VERSION,
        "status": "ok",
        "verdict": None,
        "first_divergence": {"opd": _opd_divergence(opd), "mode": None},
        "error": None,
        "opd_result": _opd_result(opd),
        "mode": None,
    }
    if opd.get("exist") is not True:
        return {**base, "status": "local_opd_absent"}
    if not isinstance(comparison, Mapping):
        raise ValueError("present OPD has no mode comparison payload")
    if comparison.get("status") in {"unsupported", "context_mismatch", "error", "timeout"}:
        return {
            **base,
            "status": str(comparison.get("status")),
            "error": {
                "stage": "local_mode_details",
                "code": str(comparison.get("status")),
                "detail": str(comparison.get("reason") or ""),
            },
        }
    if opd.get("level") not in {"physical", "strict"}:
        return {**base, "verdict": "unproven"}
    mode, mode_divergence, mode_physical, mode_strict = _mode_result(comparison)
    opd_physical = opd.get("level") in {"physical", "strict"}
    opd_strict = opd.get("level") == "strict"
    if opd_physical and mode_physical:
        verdict = "identical" if opd_strict and mode_strict else "equivalent"
    else:
        verdict = "unproven"
    return {
        **base,
        "verdict": verdict,
        "first_divergence": {
            "opd": _opd_divergence(opd),
            "mode": mode_divergence,
        },
        "mode": mode,
    }


__all__ = [
    "COMPARATOR_VERSION",
    "DEFAULT_TOLERANCE",
    "VERDICTS",
    "canonical_verdict",
    "result_for_comparison",
]
