"""Project raw comparison payloads onto the formal validation vocabulary."""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from fractions import Fraction
from functools import lru_cache
import math
import re
from typing import Any, Iterable, Mapping


COMPARATOR_VERSION = "validation.v9"
DEFAULT_TOLERANCE = 5e-5
KINDS = ("dsp", "mag", "strain")


@lru_cache(maxsize=4096)
def physical_k_label(value: str) -> str:
    """Fold a numeric three-component K prefix modulo reciprocal vectors."""

    match = re.fullmatch(
        r"(?P<prefix>[^\[\]]*)"
        r"\[(?P<k>[^\[\]]+)]"
        r"(?P<label>[^\[\]]+)",
        value,
    )
    if match is None:
        return value
    raw_components = match.group("k").split(",")
    if len(raw_components) != 3:
        return value
    try:
        components = tuple(Fraction(component.strip()) % 1 for component in raw_components)
    except (ValueError, ZeroDivisionError):
        return value

    def text(component: Fraction) -> str:
        if component.denominator == 1:
            return str(component.numerator)
        return f"{component.numerator}/{component.denominator}"

    return (
        f"{match.group('prefix')}"
        f"[{','.join(text(component) for component in components)}]"
        f"{match.group('label')}"
    )


@dataclass(frozen=True)
class ModeIdentity:
    kind: str
    irrep: str
    physical_irrep: str
    site: str | None
    site_irrep: str | None
    ordinal: int
    irrep_direction: str
    site_direction: str

    @property
    def site_key(self) -> tuple[str, str, str | None]:
        return self.kind, self.irrep, self.site

    @property
    def site_irrep_key(self) -> tuple[str, str, str | None, str | None]:
        return self.kind, self.irrep, self.site, self.site_irrep

    @property
    def family_key(self) -> tuple[Any, ...]:
        return (*self.site_irrep_key, self.ordinal)

def _split_final_direction(value: str, *, field: str) -> tuple[str, str]:
    match = re.fullmatch(r"(?P<head>.*?)(?P<direction>\([^()]*\))", value)
    if match is None or not match.group("head"):
        raise ValueError(f"cannot parse {field} from mode label {value!r}")
    return match.group("head"), match.group("direction")


def _base_and_ordinal(value: str) -> tuple[str, int]:
    match = re.fullmatch(r"(?P<base>.*?)(?:_(?P<ordinal>\d+))?", value)
    if match is None or not match.group("base"):
        raise ValueError(f"invalid site-irrep family {value!r}")
    return match.group("base"), int(match.group("ordinal") or 1)


def parse_mode_identity(label: object, *, expected_kind: str) -> ModeIdentity:
    """Parse the visible family label without using hidden Web identifiers."""

    if not isinstance(label, str) or not label:
        raise ValueError("mode label must be non-empty text")
    first_open = label.find("[")
    first_close = label.find("]", first_open + 1)
    if first_open < 1 or first_close < 0:
        raise ValueError(f"mode label has no K vector: {label!r}")
    k_vector = re.sub(r"\s+", "", label[first_open + 1 : first_close])
    tail = label[first_close + 1 :]
    second_open = tail.find("[")
    if expected_kind == "strain":
        match = re.fullmatch(
            r"(?P<irrep>.*?)(?P<irrep_direction>\([^()]*\))"
            r"strain(?:_(?P<ordinal>\d+))?(?P<site_direction>\([^()]*\))",
            tail,
        )
        if match is None:
            raise ValueError(f"cannot parse strain family label {label!r}")
        irrep_identity = f"[{k_vector}]{match.group('irrep')}"
        return ModeIdentity(
            kind="strain",
            irrep=irrep_identity,
            physical_irrep=physical_k_label(irrep_identity),
            site=None,
            site_irrep=None,
            ordinal=int(match.group("ordinal") or 1),
            irrep_direction=match.group("irrep_direction"),
            site_direction=match.group("site_direction"),
        )
    if second_open < 0:
        raise ValueError(f"mode label has no site identity: {label!r}")
    second_close = tail.find("]", second_open + 1)
    if second_close < 0:
        raise ValueError(f"mode label has an unterminated site identity: {label!r}")
    irrep, irrep_direction = _split_final_direction(
        tail[:second_open], field="parent irrep"
    )
    site_parts = tail[second_open + 1 : second_close].split(":")
    if len(site_parts) != 3:
        raise ValueError(f"invalid mode site identity in {label!r}")
    atom, wyckoff, kind = site_parts
    if kind != expected_kind or kind not in {"dsp", "mag"}:
        raise ValueError(
            f"mode label kind {kind!r} does not match section {expected_kind!r}"
        )
    site_tail = tail[second_close + 1 :]
    placeholder = re.fullmatch(r"mode(?:_(?P<ordinal>\d+))?", site_tail)
    if placeholder is not None:
        site_irrep = None
        ordinal = int(placeholder.group("ordinal") or 1)
        site_direction = ""
    else:
        site_family, site_direction = _split_final_direction(
            site_tail, field="site irrep"
        )
        site_irrep, ordinal = _base_and_ordinal(site_family)
    irrep_identity = f"[{k_vector}]{irrep}"
    return ModeIdentity(
        kind=kind,
        irrep=irrep_identity,
        physical_irrep=physical_k_label(irrep_identity),
        site=f"{atom}:{wyckoff}",
        site_irrep=site_irrep,
        ordinal=ordinal,
        irrep_direction=irrep_direction,
        site_direction=site_direction,
    )


def _payloads(comparison: Mapping[str, Any]) -> dict[str, Mapping[str, Any]]:
    return {
        "dsp": comparison.get("definitions") or {},
        "mag": comparison.get("magnetic_definitions") or {},
        "strain": comparison.get("strain_definitions") or {},
    }


def _mode_labels(payload: Mapping[str, Any], *, side: str, kind: str) -> list[str]:
    labels: list[str] = []
    for mode in payload.get("modes") or []:
        if not isinstance(mode, Mapping):
            continue
        if kind == "strain":
            value = (mode.get(side) or {}).get("label") if isinstance(mode.get(side), Mapping) else None
        else:
            value = mode.get(f"{side}_label")
        if value is not None:
            labels.append(str(value))
    tail_key = f"unpaired_{side}_tail_labels"
    labels.extend(str(value) for value in payload.get(tail_key) or [])
    return labels


def _identities(
    payloads: Mapping[str, Mapping[str, Any]], *, side: str
) -> tuple[ModeIdentity, ...]:
    out: list[ModeIdentity] = []
    for kind in KINDS:
        out.extend(
            parse_mode_identity(label, expected_kind=kind)
            for label in _mode_labels(payloads[kind], side=side, kind=kind)
        )
    return tuple(out)


def _diff_layer(web: Iterable[str], local: Iterable[str]) -> dict[str, Any]:
    web_set = set(web)
    local_set = set(local)
    missing = sorted(web_set - local_set)
    extra = sorted(local_set - web_set)
    result: dict[str, Any] = {"pass": not missing and not extra}
    if missing:
        result["missing"] = missing
    if extra:
        result["extra"] = extra
    return result


def _identity_token(
    identity: ModeIdentity, level: str, *, physical: bool = True
) -> str:
    parts = [identity.kind, identity.physical_irrep if physical else identity.irrep]
    if level in {"site", "site_irrep", "family"}:
        parts.append(identity.site or "-")
    if level in {"site_irrep", "family"}:
        parts.append(identity.site_irrep or "-")
    if level == "family":
        parts.append(str(identity.ordinal))
    return ":".join(parts)


def _count_layer(
    web: Iterable[ModeIdentity], local: Iterable[ModeIdentity]
) -> dict[str, Any]:
    def counts(items: Iterable[ModeIdentity]) -> Counter[str]:
        return Counter(_identity_token(item, "site_irrep") for item in items)

    web_counts = counts(web)
    local_counts = counts(local)
    mismatches = {
        key: {"web": web_counts.get(key, 0), "local": local_counts.get(key, 0)}
        for key in sorted(set(web_counts) | set(local_counts))
        if web_counts.get(key, 0) != local_counts.get(key, 0)
    }
    result: dict[str, Any] = {"pass": not mismatches}
    if mismatches:
        result["mismatches"] = mismatches
    return result


def _normfactor_layer(
    payloads: Mapping[str, Mapping[str, Any]], tol: float
) -> dict[str, Any]:
    # Mode comparison already owns the full-label or sequence pairing authority.
    mismatches: list[str] = []
    for kind in KINDS:
        for index, mode in enumerate(payloads[kind].get("modes") or [], start=1):
            if not isinstance(mode, Mapping):
                continue
            if kind == "strain":
                web = mode.get("web")
                local = mode.get("local")
                left = web.get("normfactor") if isinstance(web, Mapping) else None
                right = (
                    local.get("normfactor") if isinstance(local, Mapping) else None
                )
            else:
                left = mode.get("web_normfactor")
                right = mode.get("local_normfactor")
            if left is None or right is None:
                equal = left is None and right is None
            else:
                equal = math.isclose(
                    float(left), float(right), rel_tol=0.0, abs_tol=tol
                )
            if not equal:
                mismatches.append(f"{kind}:{index}")
    result: dict[str, Any] = {"pass": not mismatches}
    if mismatches:
        result["mismatches"] = mismatches
    return result


def _complete_mode_presentation_proof(payload: Mapping[str, Any]) -> bool:
    """Require a complete physical assignment before ignoring zip-order failures."""
    if (
        payload.get("count_match") is not True
        or payload.get("unpaired_web_tail_labels") not in (None, [], ())
        or payload.get("unpaired_local_tail_labels") not in (None, [], ())
    ):
        return False
    modes = payload.get("modes")
    if not isinstance(modes, (list, tuple)) or not modes:
        return False
    web_count = payload.get("web_count")
    local_count = payload.get("local_count")
    if (
        isinstance(web_count, bool)
        or isinstance(local_count, bool)
        or not isinstance(web_count, int)
        or not isinstance(local_count, int)
        or web_count != len(modes)
        or local_count != len(modes)
        or any(not isinstance(mode, Mapping) for mode in modes)
    ):
        return False
    status = payload.get("status")
    if status == "signed_permutation":
        return payload.get("strict_status") == "signed_permutation"
    if status != "family_basis" or payload.get("strict_status") != "family_basis":
        return False
    proof = payload.get("physical_family_basis")
    if not isinstance(proof, Mapping):
        return False

    counts = []
    for field in (
        "block_count",
        "signed_permutation_block_count",
        "family_basis_block_count",
        "parent_direction_gauge_block_count",
    ):
        value = proof.get(field)
        if isinstance(value, bool) or not isinstance(value, int) or value < 0:
            return False
        counts.append(value)
    block_count, signed_count, family_count, direction_count = counts
    if (
        block_count <= 0
        or family_count <= 0
        or signed_count + family_count != block_count
        or direction_count > block_count
    ):
        return False
    numeric_proof: dict[str, float] = {}
    for field in (
        "max_condition",
        "max_residual",
        "visible_component_rounding_bound",
        "max_rounding_ratio",
    ):
        value = proof.get(field)
        if isinstance(value, bool) or not isinstance(value, (int, float)):
            return False
        numeric = float(value)
        if not math.isfinite(numeric) or numeric < 0.0:
            return False
        if field == "visible_component_rounding_bound" and numeric <= 0.0:
            return False
        numeric_proof[field] = numeric
    residual = numeric_proof["max_residual"]
    bound = numeric_proof["visible_component_rounding_bound"]
    ratio = numeric_proof["max_rounding_ratio"]
    return (
        residual <= bound
        and ratio <= 1.0
        and math.isclose(ratio, residual / bound, rel_tol=1.0e-12, abs_tol=1.0e-15)
    )


def _row_atom_layer(payloads: Mapping[str, Mapping[str, Any]]) -> dict[str, Any]:
    exact_mismatches = 0
    gauge_mismatches = 0
    real_mismatches = 0
    positions = 0
    presentation_assignments = 0
    for kind in ("dsp", "mag"):
        payload = payloads[kind]
        has_presentation_proof = _complete_mode_presentation_proof(payload)
        for mode in payload.get("modes") or []:
            if not isinstance(mode, Mapping):
                continue
            web_identity = parse_mode_identity(
                mode.get("web_label"), expected_kind=kind
            )
            local_identity = parse_mode_identity(
                mode.get("local_label"), expected_kind=kind
            )
            web_mode_parent = (
                web_identity.site.split(":", 1)[0] if web_identity.site else None
            )
            local_mode_parent = (
                local_identity.site.split(":", 1)[0] if local_identity.site else None
            )
            if (
                mode.get("count_match") is not True
                or mode.get("status") == "position"
            ):
                if has_presentation_proof:
                    presentation_assignments += 1
                else:
                    positions += 1
                continue
            for row in mode.get("rows") or []:
                if not isinstance(row, Mapping):
                    continue
                web = row.get("web") if isinstance(row.get("web"), Mapping) else {}
                local = row.get("local") if isinstance(row.get("local"), Mapping) else {}
                web_atom = web.get("atom")
                local_atom = local.get("atom")
                if web_atom == local_atom and web_atom is not None:
                    continue
                exact_mismatches += 1
                web_parent = (
                    re.sub(r"_\d+$", "", web_atom)
                    if isinstance(web_atom, str)
                    else web_mode_parent
                )
                local_parent = (
                    re.sub(r"_\d+$", "", local_atom)
                    if isinstance(local_atom, str)
                    else local_mode_parent
                )
                if web_parent is not None and web_parent == local_parent:
                    gauge_mismatches += 1
                else:
                    real_mismatches += 1
    strict_pass = (
        positions == 0
        and exact_mismatches == 0
        and presentation_assignments == 0
    )
    result: dict[str, Any] = {"pass": positions == 0 and real_mismatches == 0}
    if not strict_pass:
        result["strict_pass"] = False
    gauges = []
    if gauge_mismatches:
        gauges.append("orbit_ordinal")
        result["label_gauge_mismatches"] = gauge_mismatches
    if presentation_assignments:
        gauges.append("mode_family_assignment")
        result["presentation_assignment_modes"] = presentation_assignments
    if gauges:
        result["gauge"] = gauges
    if real_mismatches:
        result["label_real_mismatches"] = real_mismatches
    if positions:
        result["position_failures"] = positions
    return result


def _vector_layer(payloads: Mapping[str, Mapping[str, Any]]) -> dict[str, Any]:
    gauge: set[str] = set()
    failures: dict[str, str] = {}
    for kind in KINDS:
        payload = payloads[kind]
        status = str(payload.get("status") or "missing")
        if status == "ok":
            pass
        elif status == "sign":
            gauge.add("sign")
        elif status == "signed_permutation":
            gauge.add("signed_permutation")
        elif status == "family_basis":
            gauge.add("family_basis")
        elif int(payload.get("web_count") or 0) == 0 and int(payload.get("local_count") or 0) == 0:
            pass
        else:
            failures[kind] = status
        if payload.get("validated") is True and payload.get("basis_change") is not None:
            gauge.add("basis_frame")
    result: dict[str, Any] = {"pass": not failures}
    if gauge:
        result["gauge"] = sorted(gauge)
    if failures:
        result["failures"] = failures
    return result


def _mode_definition(
    comparison: Mapping[str, Any], *, tol: float
) -> tuple[dict[str, Any], str | None, bool, bool]:
    payloads = _payloads(comparison)
    web_ids = _identities(payloads, side="web")
    local_ids = _identities(payloads, side="local")
    layers: dict[str, Any] = {
        "kind_set": None,
        "irrep_set": None,
        "site_set": None,
        "site_irrep": None,
        "mode_count": None,
        "normfactor": None,
        "row_atom": None,
        "vector": None,
    }
    web_kinds = {
        kind for kind, payload in payloads.items() if int(payload.get("web_count") or 0) > 0
    }
    local_kinds = {
        kind for kind, payload in payloads.items() if int(payload.get("local_count") or 0) > 0
    }
    kinds = _diff_layer(web_kinds, local_kinds)
    layers["kind_set"] = kinds
    if not kinds["pass"]:
        return layers, "mode.mode_definition.kind_set", False, False
    irreps = _diff_layer(
        (_identity_token(item, "irrep") for item in web_ids),
        (_identity_token(item, "irrep") for item in local_ids),
    )
    presented_irreps = _diff_layer(
        (_identity_token(item, "irrep", physical=False) for item in web_ids),
        (_identity_token(item, "irrep", physical=False) for item in local_ids),
    )
    if not presented_irreps["pass"]:
        irreps["strict_pass"] = False
        irreps["gauge"] = ["reciprocal_coordinate"]
        if presented_irreps.get("missing"):
            irreps["presentation_missing"] = presented_irreps["missing"]
        if presented_irreps.get("extra"):
            irreps["presentation_extra"] = presented_irreps["extra"]
    layers["irrep_set"] = irreps
    if not irreps["pass"]:
        return layers, "mode.mode_definition.irrep_set", False, False
    sites = _diff_layer(
        (_identity_token(item, "site") for item in web_ids),
        (_identity_token(item, "site") for item in local_ids),
    )
    layers["site_set"] = sites
    if not sites["pass"]:
        return layers, "mode.mode_definition.site_set", False, False
    site_irreps = _diff_layer(
        (_identity_token(item, "site_irrep") for item in web_ids),
        (_identity_token(item, "site_irrep") for item in local_ids),
    )
    layers["site_irrep"] = site_irreps
    if not site_irreps["pass"]:
        return layers, "mode.mode_definition.site_irrep", False, False
    section_count_mismatches = {
        kind: {
            "web": int(payload.get("web_count") or 0),
            "local": int(payload.get("local_count") or 0),
        }
        for kind, payload in payloads.items()
        if payload.get("count_match") is not True
    }
    counts = _count_layer(web_ids, local_ids)
    if section_count_mismatches:
        counts["pass"] = False
        counts["section_mismatches"] = section_count_mismatches
    layers["mode_count"] = counts
    if not counts["pass"]:
        return layers, "mode.mode_definition.mode_count", False, False

    norms = _normfactor_layer(payloads, tol)
    layers["normfactor"] = norms
    rows = _row_atom_layer(payloads)
    layers["row_atom"] = rows
    if not rows["pass"]:
        return layers, "mode.mode_definition.row_atom", False, False
    vectors = _vector_layer(payloads)
    layers["vector"] = vectors
    physical = bool(vectors["pass"])
    strict = (
        physical
        and irreps.get("strict_pass") is not False
        and not vectors.get("gauge")
        and rows.get("strict_pass") is not False
        and bool(norms["pass"])
    )
    if not vectors["pass"]:
        divergence = "mode.mode_definition.vector"
    elif irreps.get("strict_pass") is False:
        divergence = "mode.mode_definition.irrep_set.gauge"
    elif not norms["pass"]:
        divergence = "mode.mode_definition.normfactor"
    elif rows.get("strict_pass") is False:
        divergence = "mode.mode_definition.row_atom.gauge"
    elif vectors.get("gauge"):
        divergence = "mode.mode_definition.vector.gauge"
    else:
        divergence = None
    return layers, divergence, physical, strict


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
    comparison: Mapping[str, Any], *, tol: float
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
    definition, divergence, physical, strict = _mode_definition(comparison, tol=tol)
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
        return {**base, "verdict": "physical_fail"}
    if opd.get("level") not in {"physical", "strict"}:
        return {**base, "verdict": "physical_fail"}
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
    mode, mode_divergence, mode_physical, mode_strict = _mode_result(
        comparison, tol=tol
    )
    opd_physical = opd.get("level") in {"physical", "strict"}
    opd_strict = opd.get("level") == "strict"
    if opd_physical and mode_physical:
        verdict = "strict_pass" if opd_strict and mode_strict else "physical_pass"
    else:
        verdict = "physical_fail"
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
    "ModeIdentity",
    "parse_mode_identity",
    "physical_k_label",
    "result_for_comparison",
]
