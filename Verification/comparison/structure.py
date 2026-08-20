"""Compare parent, subgroup, lattice, and atom structure payloads."""

from __future__ import annotations

from collections import Counter
from fractions import Fraction
import math
import re
from typing import Any, Callable

import gemmi

from APOSTRUCT.Backend.modes.mode_detail_text import parent_display_xyz
from APOSTRUCT.Backend.modes.structure.magnetic_wyckoff import magnetic_orbit_points
from Verification.comparison.basis_atom import (
    basis_equivalent_atom_frame_transport,
    basis_equivalent_ordinary_atom_frame_transport,
    exact_unimodular_matrix3,
)
from Verification.comparison.basis_lattice import (
    basis_equivalent_lattice_transport,
    unimodular_basis_change as _unimodular_basis_change,
)
from Verification.comparison.magnetic_frame import (
    presentation_to_standard_child_cinter,
    selected_magnetic_frame,
)
from Verification.comparison.opd import _k_active_relation
from Verification.comparison.selected_state import _rational_tokens
from Verification.parsers.complete_mode_text import CompleteModeDetails, StructureAtom


def _close(a: float, b: float, tol: float) -> bool:
    return math.isfinite(a) and math.isfinite(b) and abs(a - b) <= tol


def _vec_close(a: tuple[float, ...] | list[float], b: tuple[float, ...] | list[float], tol: float) -> bool:
    return len(a) == len(b) and all(_close(float(x), float(y), tol) for x, y in zip(a, b, strict=True))


def _frac_close(a: float, b: float, tol: float) -> bool:
    diff = abs((float(a) - float(b)) % 1.0)
    return min(diff, abs(diff - 1.0)) <= tol


def _frac_vec_close(a: tuple[float, ...] | list[float], b: tuple[float, ...] | list[float], tol: float) -> bool:
    return len(a) == len(b) and all(_frac_close(float(x), float(y), tol) for x, y in zip(a, b, strict=True))


def _periodic_orbit_candidate_indices(
    web_items: list[tuple[str, str, tuple[float, float, float]]],
    local_items: list[tuple[str, str, Any, Any]],
    local_orbits: list[list[tuple[float, float, float]]],
    tol: float,
    *,
    match_site: bool,
    transform: Callable[
        [tuple[float, float, float]], tuple[float, float, float]
    ]
    | None = None,
) -> list[list[int]] | None:
    """Index orbit points while retaining the legacy exact predicate and row order."""

    if len(local_items) != len(local_orbits) or not math.isfinite(tol) or tol <= 0.0:
        return None
    bucket_width = 2.0 * tol
    inverse_width = 1.0 / bucket_width
    if not math.isfinite(inverse_width):
        return None
    bucket_count = max(1, int(math.floor(inverse_width)))
    buckets: dict[tuple[str, str | None, int, int, int], list[int]] = {}
    available_groups: set[tuple[str, str | None]] = set()

    def bucket_key(
        kind: str,
        site: str,
        point: tuple[float, float, float],
    ) -> tuple[str, str | None, int, int, int] | None:
        if len(point) != 3:
            return None
        try:
            folded = tuple(float(value) % 1.0 for value in point)
        except (OverflowError, TypeError, ValueError):
            return None
        if not all(math.isfinite(value) for value in folded):
            return None
        cell = tuple(
            min(bucket_count - 1, int(math.floor(value / bucket_width)))
            for value in folded
        )
        return (kind, site if match_site else None, *cell)

    for local_index, ((kind, site, _xyz, _explicit), orbit) in enumerate(
        zip(local_items, local_orbits, strict=True)
    ):
        group = (kind, site if match_site else None)
        available_groups.add(group)
        local_keys: set[tuple[str, str | None, int, int, int]] = set()
        for point in orbit:
            key = bucket_key(kind, site, point)
            if key is None:
                return None
            local_keys.add(key)
        for key in sorted(local_keys):
            buckets.setdefault(key, []).append(local_index)

    result: list[list[int]] = []
    for kind, site, xyz in web_items:
        group = (kind, site if match_site else None)
        if group not in available_groups:
            result.append([])
            continue
        target = transform(xyz) if transform is not None else xyz
        center_key = bucket_key(kind, site, target)
        if center_key is None:
            return None
        center = center_key[-3:]
        candidates: set[int] = set()
        for first in (-1, 0, 1):
            for second in (-1, 0, 1):
                for third in (-1, 0, 1):
                    candidates.update(
                        buckets.get(
                            (
                                group[0],
                                group[1],
                                (center[0] + first) % bucket_count,
                                (center[1] + second) % bucket_count,
                                (center[2] + third) % bucket_count,
                            ),
                            (),
                        )
                    )
        result.append(
            [
                local_index
                for local_index in sorted(candidates)
                if any(
                    len(point) == 3
                    and all(math.isfinite(value) for value in (*target, *point))
                    and all(
                        abs(((left - right + 0.5) % 1.0) - 0.5) <= tol
                        for left, right in zip(target, point, strict=True)
                    )
                    for point in local_orbits[local_index]
                )
            ]
        )
    return result


def _spacegroup_payload(web: CompleteModeDetails, preview: dict[str, Any]) -> dict[str, Any]:
    local = preview.get("space_group") or ((preview.get("input") or {}).get("parent") or {})
    web_parent = web.parent
    status = (
        "ok"
        if web_parent is not None
        and int(web_parent.number or 0) == int(local.get("number") or 0)
        else "fail"
    )
    return {
        "status": status,
        "web": None if web_parent is None else {"number": web_parent.number, "symbol": web_parent.symbol},
        "local": {"number": local.get("number"), "symbol": local.get("symbol")},
    }


def _lattice_payload(web: CompleteModeDetails, local: dict[str, Any], tol: float) -> dict[str, Any]:
    local_lattice = local.get("lattice") or {}
    keys = ("a", "b", "c", "alpha", "beta", "gamma")
    matches = {}
    if web.undistorted_lattice is None:
        return {"status": "missing_web", "matches": matches}
    for key in keys:
        matches[key] = _close(float(getattr(web.undistorted_lattice, key)), float(local_lattice.get(key, float("nan"))), tol)
    return {
        "status": "ok" if all(matches.values()) else "fail",
        "matches": matches,
        "web": {key: getattr(web.undistorted_lattice, key) for key in keys},
        "local": {key: local_lattice.get(key) for key in keys},
    }


def _basis_equivalent_lattice_payload(
    web: CompleteModeDetails,
    local: dict[str, Any],
    parent_lattice: dict[str, Any],
    selected_state: object,
    subgroup: dict[str, Any],
    tol: float,
) -> dict[str, Any] | None:
    """Validate a child lattice printed in an equivalent selected basis."""

    if (
        web.undistorted_lattice is None
        or not isinstance(selected_state, dict)
        or selected_state.get("status")
        not in {"equivalent_basis", "equivalent_embedding"}
        or subgroup.get("status")
        not in {"basis_equivalent", "embedding_equivalent"}
        or parent_lattice.get("status") != "ok"
    ):
        return None
    selected_request = selected_state.get("request")
    selected_local = selected_state.get("local")
    subgroup_matches = subgroup.get("matches")
    parent_matches = parent_lattice.get("matches")
    if (
        not isinstance(selected_request, (tuple, list))
        or not isinstance(selected_local, (tuple, list))
        or len(selected_request) != 4
        or len(selected_local) != 4
        or selected_request[0] != selected_local[0]
        or selected_request[1] != selected_local[1]
        or (
            selected_request[3] != selected_local[3]
            and selected_state.get("status") != "equivalent_embedding"
        )
        or not isinstance(subgroup_matches, dict)
        or not isinstance(parent_matches, dict)
        or any(parent_matches.get(key) is not True for key in ("a", "b", "c", "alpha", "beta", "gamma"))
        or (
            subgroup_matches.get("basis") is not False
            and selected_state.get("status") != "equivalent_embedding"
        )
        or any(
            subgroup_matches.get(key) is not True
            for key in ("number", "display_label", "s", "i", "k_active")
        )
        or (
            subgroup_matches.get("origin") is not True
            and selected_state.get("status") != "equivalent_embedding"
        )
    ):
        return None
    subgroup_web = subgroup.get("web")
    subgroup_local = subgroup.get("local")
    if not isinstance(subgroup_web, dict) or not isinstance(subgroup_local, dict):
        return None
    try:
        web_number = int(subgroup_web.get("number") or 0)
        local_number = int(subgroup_local.get("number") or 0)
        selected_request_number = int(selected_request[1] or 0)
        selected_local_number = int(selected_local[1] or 0)
    except (OverflowError, TypeError, ValueError):
        return None
    if web_number != selected_request_number or local_number != selected_local_number:
        return None
    selected_change = selected_state.get("basis_change")
    subgroup_change = subgroup.get("basis_change")
    try:
        selected_change_key = tuple(tuple(int(value) for value in row) for row in selected_change)
        subgroup_change_key = tuple(tuple(int(value) for value in row) for row in subgroup_change)
    except (TypeError, ValueError):
        return None
    if selected_change_key != subgroup_change_key:
        return None

    keys = ("a", "b", "c", "alpha", "beta", "gamma")
    parent_local = parent_lattice.get("local")
    local_lattice = local.get("lattice")
    if not isinstance(parent_local, dict) or not isinstance(local_lattice, dict):
        return None
    try:
        parent_cell = tuple(float(parent_local[key]) for key in keys)
    except (KeyError, OverflowError, TypeError, ValueError):
        return None
    transport = basis_equivalent_lattice_transport(
        parent_cell,
        selected_request[2],
        selected_local[2],
    )
    if transport is None or transport.basis_change != selected_change_key:
        return None
    try:
        local_matches = {
            key: _close(transport.local_cell[index], float(local_lattice[key]), tol)
            for index, key in enumerate(keys)
        }
        web_matches = {
            key: _close(
                transport.request_cell[index],
                float(getattr(web.undistorted_lattice, key)),
                tol,
            )
            for index, key in enumerate(keys)
        }
    except (KeyError, OverflowError, TypeError, ValueError):
        return None
    if not all(local_matches.values()) or not all(web_matches.values()):
        return None
    return {
        "status": "basis_equivalent",
        "validated": True,
        "provenance": "M:exact_gl3z_metric_transport",
        "matches": web_matches,
        "local_matches": local_matches,
        "basis_change": transport.basis_change,
        "metric_residual": transport.max_metric_residual,
        "condition_number": transport.condition_number,
        "web": {key: getattr(web.undistorted_lattice, key) for key in keys},
        "local": {key: local_lattice.get(key) for key in keys},
        "transported": {
            key: transport.request_cell[index]
            for index, key in enumerate(keys)
        },
    }


def _parent_lattice_payload(web: CompleteModeDetails, preview: dict[str, Any], tol: float) -> dict[str, Any]:
    local_lattice = (preview.get("input") or {}).get("lattice") or {}
    keys = ("a", "b", "c", "alpha", "beta", "gamma")
    matches = {}
    if web.parent_lattice is None:
        return {"status": "missing_web", "matches": matches}
    for key in keys:
        matches[key] = _close(float(getattr(web.parent_lattice, key)), float(local_lattice.get(key, float("nan"))), tol)
    return {
        "status": "ok" if all(matches.values()) else "fail",
        "matches": matches,
        "web": {key: getattr(web.parent_lattice, key) for key in keys},
        "local": {key: local_lattice.get(key) for key in keys},
    }


def _atom_kind(label: str | None) -> str:
    text = str(label or "")
    return re.sub(r"_\d+$", "", text)


def _atom_key_from_web(row: StructureAtom, tol: float) -> tuple[str, str, tuple[int, int, int]]:
    scale = max(1, int(round(1 / tol)))
    xyz = tuple(int(round((float(value) % 1.0) * scale)) % scale for value in row.xyz)
    return (_atom_kind(row.label), str(row.site), xyz)


def _atom_key_from_local(row: dict[str, Any], tol: float) -> tuple[str, str, tuple[int, int, int]] | None:
    if row.get("xyz"):
        xyz = row.get("xyz")
    elif row.get("fract"):
        xyz = parent_display_xyz(row)
    else:
        return None
    scale = max(1, int(round(1 / tol)))
    key_xyz = tuple(int(round((float(value) % 1.0) * scale)) % scale for value in xyz)
    site = str(row.get("site") or "")
    if not site and row.get("wyckoff"):
        mult = str(row.get("wyckoff_multiplicity") or row.get("multiplicity") or "")
        wyckoff = str(row.get("wyckoff") or "")
        site = wyckoff if wyckoff[:1].isdigit() else f"{mult}{wyckoff}"
    return (_atom_kind(row.get("label")), site, key_xyz)


def _quantized_xyz(xyz: tuple[float, float, float] | list[float], tol: float) -> tuple[int, int, int]:
    scale = max(1, int(round(1 / tol)))
    return tuple(int(round((float(value) % 1.0) * scale)) % scale for value in xyz)  # type: ignore[return-value]


def _orbit_xyz_keys(
    sg: int | None,
    xyz_key: tuple[int, int, int],
    tol: float,
) -> set[tuple[int, int, int]]:
    if sg is None:
        return {xyz_key}
    scale = max(1, int(round(1 / tol)))
    xyz = [value / scale for value in xyz_key]
    try:
        ops = gemmi.find_spacegroup_by_number(int(sg)).operations()
    except Exception:
        return {xyz_key}
    out = {xyz_key}
    for op in ops:
        try:
            point = op.apply_to_xyz(xyz)
        except Exception:
            continue
        out.add(_quantized_xyz(point, tol))
    return out


def _magnetic_orbit_xyz_keys(
    magnetic_group: int | None,
    xyz_key: tuple[int, int, int],
    tol: float,
    transform: Callable[[tuple[float, float, float]], tuple[float, float, float]] | None,
) -> set[tuple[int, int, int]]:
    if magnetic_group is None or transform is None:
        return {xyz_key}
    scale = max(1, int(round(1 / tol)))
    xyz = transform(tuple(value / scale for value in xyz_key))
    try:
        return {
            _quantized_xyz(tuple(float(value) for value in point), tol)
            for point in magnetic_orbit_points(int(magnetic_group), xyz)
        }
    except (IndexError, KeyError, TypeError, ValueError, ZeroDivisionError):
        return {xyz_key}


def _orbit_payload(
    web_counter,
    local_counter,
    *,
    orbit_sg: int | None,
    magnetic_group: int | None,
    magnetic_transform: Callable[[tuple[float, float, float]], tuple[float, float, float]] | None,
    tol: float,
) -> dict[str, Any]:
    local_items: list[tuple[str, str, tuple[int, int, int], set[tuple[int, int, int]]]] = []
    for kind, site, xyz in local_counter.elements():
        orbit = (
            _magnetic_orbit_xyz_keys(magnetic_group, xyz, tol, magnetic_transform)
            if magnetic_group is not None and magnetic_transform is not None
            else _orbit_xyz_keys(orbit_sg, xyz, tol)
        )
        local_items.append((kind, site, xyz, orbit))
    used = [False] * len(local_items)
    missing: list[list[Any]] = []
    matched = 0
    for kind, site, xyz in web_counter.elements():
        target_xyz = xyz
        if magnetic_group is not None and magnetic_transform is not None:
            scale = max(1, int(round(1 / tol)))
            target_xyz = _quantized_xyz(
                magnetic_transform(tuple(value / scale for value in xyz)), tol
            )
        hit = None
        for index, (lkind, lsite, _lxyz, orbit) in enumerate(local_items):
            if used[index] or lkind != kind or lsite != site:
                continue
            if target_xyz in orbit:
                hit = index
                break
        if hit is None:
            missing.append([kind, site, list(xyz)])
        else:
            used[hit] = True
            matched += 1
    extra = [
        [kind, site, list(xyz)]
        for used_flag, (kind, site, xyz, _orbit) in zip(used, local_items, strict=True)
        if not used_flag
    ]
    return {
        "status": "ok" if not missing and not extra else "fail",
        "orbit_sg": orbit_sg,
        "magnetic_group": magnetic_group,
        "matched": matched,
        "web_count": sum(web_counter.values()),
        "local_count": sum(local_counter.values()),
        "missing": missing[:20],
        "extra": extra[:20],
    }


def _structure_atom_payload(
    web_atoms: tuple[StructureAtom, ...],
    local_atoms: list[dict[str, Any]],
    tol: float,
    *,
    orbit_sg: int | None = None,
    magnetic_group: int | None = None,
    magnetic_transform: Callable[[tuple[float, float, float]], tuple[float, float, float]] | None = None,
) -> dict[str, Any]:
    def finite_xyz(value: object) -> tuple[float, float, float]:
        if (
            not isinstance(value, (list, tuple))
            or len(value) != 3
            or any(isinstance(item, bool) for item in value)
        ):
            raise ValueError("atom coordinates must contain exactly three numbers")
        xyz = tuple(float(item) % 1.0 for item in value)
        if not all(math.isfinite(item) for item in xyz):
            raise ValueError("atom coordinates must be finite")
        return xyz  # type: ignore[return-value]

    web_items = []
    invalid_web_rows: list[int] = []
    for index, row in enumerate(web_atoms):
        try:
            xyz = finite_xyz(row.xyz)
        except (OverflowError, TypeError, ValueError):
            invalid_web_rows.append(index)
            continue
        web_items.append((_atom_kind(row.label), str(row.site), xyz))

    local_items = []
    invalid_local_rows: list[int] = []
    for index, row in enumerate(local_atoms):
        if not isinstance(row, dict):
            invalid_local_rows.append(index)
            continue
        try:
            if row.get("xyz"):
                xyz = finite_xyz(row["xyz"])
            elif row.get("fract"):
                xyz = finite_xyz(parent_display_xyz(row))
            else:
                raise ValueError("atom row has no coordinates")
        except (IndexError, KeyError, OverflowError, TypeError, ValueError, ZeroDivisionError):
            invalid_local_rows.append(index)
            continue
        site = str(row.get("site") or "")
        if not site and row.get("wyckoff"):
            mult = str(row.get("wyckoff_multiplicity") or row.get("multiplicity") or "")
            wyckoff = str(row.get("wyckoff") or "")
            site = wyckoff if wyckoff[:1].isdigit() else f"{mult}{wyckoff}"
        explicit_orbit = row.get("_presentation_orbit_points")
        if isinstance(explicit_orbit, (list, tuple)) and explicit_orbit:
            try:
                explicit_points = tuple(finite_xyz(point) for point in explicit_orbit)
            except (OverflowError, TypeError, ValueError):
                invalid_local_rows.append(index)
                continue
        elif explicit_orbit not in (None, [], ()):
            invalid_local_rows.append(index)
            continue
        else:
            explicit_points = ()
        local_items.append(
            (
                _atom_kind(row.get("label")),
                site,
                xyz,
                explicit_points,
            )
        )

    if invalid_web_rows or invalid_local_rows:
        return {
            "status": "fail",
            "reason": "malformed_atom_rows",
            "invalid_web_rows": invalid_web_rows[:20],
            "invalid_local_rows": invalid_local_rows[:20],
            "web_count": len(web_atoms),
            "local_count": len(local_atoms),
            "count_match": len(web_atoms) == len(local_atoms),
        }

    def close(left, right) -> bool:
        return (
            len(left) == 3
            and len(right) == 3
            and all(math.isfinite(value) for value in (*left, *right))
            and all(
                abs(((a - b + 0.5) % 1.0) - 0.5) <= tol
                for a, b in zip(left, right, strict=True)
            )
        )

    def orbit_points(*, symmetry: bool) -> list[list[tuple[float, float, float]]]:
        operations = []
        if symmetry and magnetic_group is None and orbit_sg is not None:
            try:
                operations = list(gemmi.find_spacegroup_by_number(int(orbit_sg)).operations())
            except Exception:
                operations = []
        local_orbits: list[list[tuple[float, float, float]]] = []
        for _kind, _site, xyz, explicit_points in local_items:
            if symmetry and explicit_points:
                points = list(explicit_points)
            elif symmetry and magnetic_group is not None and magnetic_transform is not None:
                try:
                    standard_xyz = magnetic_transform(xyz)
                    points = [
                        tuple(float(value) % 1.0 for value in point)
                        for point in magnetic_orbit_points(int(magnetic_group), standard_xyz)
                    ]
                except (IndexError, KeyError, TypeError, ValueError, ZeroDivisionError):
                    points = [xyz]
            else:
                points = [xyz]
                points.extend(
                    tuple(float(value) % 1.0 for value in operation.apply_to_xyz(xyz))
                    for operation in operations
                )
            local_orbits.append(points)
        return local_orbits

    def maximum_match(
        local_orbits: list[list[tuple[float, float, float]]],
        *,
        symmetry: bool,
        match_site: bool,
    ) -> tuple[set[int], set[int]]:
        standard_transform = (
            magnetic_transform
            if symmetry and magnetic_group is not None and magnetic_transform is not None
            else None
        )
        explicit_frames = [bool(explicit) for _kind, _site, _xyz, explicit in local_items]
        # Explicit orbits share Web presentation coordinates; generated
        # magnetic orbits use the standard child frame.
        mixed_frames = (
            standard_transform is not None
            and any(explicit_frames)
            and not all(explicit_frames)
        )
        indexed_transform = (
            None
            if standard_transform is not None and any(explicit_frames)
            else standard_transform
        )
        edges = None if mixed_frames else _periodic_orbit_candidate_indices(
            web_items,
            local_items,
            local_orbits,
            tol,
            match_site=match_site,
            transform=indexed_transform,
        )
        if edges is None:
            edges = [
                [
                    local_index
                    for local_index, (
                        (local_kind, local_site, _xyz, explicit),
                        orbit,
                    ) in enumerate(zip(local_items, local_orbits, strict=True))
                    if local_kind == web_kind
                    and (not match_site or local_site == web_site)
                    and any(
                        close(
                            standard_transform(web_xyz)
                            if standard_transform is not None and not explicit
                            else web_xyz,
                            point,
                        )
                        for point in orbit
                    )
                ]
                for web_kind, web_site, web_xyz in web_items
            ]
        local_to_web: dict[int, int] = {}

        def augment(web_index: int, visited: set[int]) -> bool:
            for local_index in edges[web_index]:
                if local_index in visited:
                    continue
                visited.add(local_index)
                prior = local_to_web.get(local_index)
                if prior is None or augment(prior, visited):
                    local_to_web[local_index] = web_index
                    return True
            return False

        for web_index in range(len(web_items)):
            augment(web_index, set())
        return set(local_to_web), set(local_to_web.values())

    web_counter = Counter(
        (kind, site, _quantized_xyz(xyz, tol))
        for kind, site, xyz in web_items
    )
    local_counter = Counter(
        (kind, site, _quantized_xyz(xyz, tol))
        for kind, site, xyz, _explicit in local_items
    )
    site_web = Counter((kind, site) for kind, site, _xyz in web_counter.elements())
    site_local = Counter((kind, site) for kind, site, _xyz in local_counter.elements())
    direct_orbits = orbit_points(symmetry=False)
    symmetry_orbits = orbit_points(symmetry=True)
    direct_local, direct_web = maximum_match(
        direct_orbits,
        symmetry=False,
        match_site=True,
    )
    orbit_local, orbit_web = maximum_match(
        symmetry_orbits,
        symmetry=True,
        match_site=True,
    )
    kind_orbit_local, kind_orbit_web = maximum_match(
        symmetry_orbits,
        symmetry=True,
        match_site=False,
    )
    direct_match = len(direct_web) == len(web_items) and len(direct_local) == len(local_items)
    orbit_match = len(orbit_web) == len(web_items) and len(orbit_local) == len(local_items)
    kind_orbit_match = (
        len(kind_orbit_web) == len(web_items)
        and len(kind_orbit_local) == len(local_items)
    )
    payload = {
        "status": (
            "ok"
            if direct_match or orbit_match
            else ("site_gauge" if kind_orbit_match else "fail")
        ),
        "direct_match": direct_match,
        "orbit_match": orbit_match,
        "kind_orbit_match": kind_orbit_match,
        "orbit_kind": "magnetic" if magnetic_group is not None else "ordinary",
        "count_match": sum(web_counter.values()) == sum(local_counter.values()),
        "site_kind_match": site_web == site_local,
        "web_count": sum(web_counter.values()),
        "local_count": sum(local_counter.values()),
        "web_site_counts": {f"{kind}:{site}": count for (kind, site), count in sorted(site_web.items())},
        "local_site_counts": {f"{kind}:{site}": count for (kind, site), count in sorted(site_local.items())},
        "missing": [list(key) for key in sorted((web_counter - local_counter).elements())[:20]],
        "extra": [list(key) for key in sorted((local_counter - web_counter).elements())[:20]],
    }
    if not orbit_match and site_web == site_local:
        payload["orbit"] = _orbit_payload(
            web_counter,
            local_counter,
            orbit_sg=orbit_sg,
            magnetic_group=magnetic_group,
            magnetic_transform=magnetic_transform,
            tol=tol,
        )
    return payload


def _parent_atom_payload(web: CompleteModeDetails, preview: dict[str, Any], tol: float) -> dict[str, Any]:
    parent_sg = web.parent.number if web.parent is not None else None
    return _structure_atom_payload(
        web.parent_atoms,
        (preview.get("input") or {}).get("atom_sites") or [],
        tol,
        orbit_sg=parent_sg,
    )


def _compact_text(value: object) -> str:
    return re.sub(r"\s+", "", str(value or ""))


def _subgroup_payload(
    web: CompleteModeDetails,
    preview: dict[str, Any],
    *,
    operation_set_equivalent: bool = False,
) -> dict[str, Any]:
    selected = preview.get("selected") or {}
    opd = selected.get("orderparam") or {}
    iso = opd.get("isotropy") if isinstance(opd, dict) and isinstance(opd.get("isotropy"), dict) else opd
    if not isinstance(iso, dict) or web.subgroup is None:
        return {"status": "missing", "web": None if web.subgroup is None else web.subgroup, "local": iso}
    local_sg = iso.get("subgroup") or {}
    magnetic_display = bool(web.subgroup.display_label and "." in web.subgroup.display_label)
    local = {
        "number": local_sg.get("ordinary_number") if magnetic_display else local_sg.get("number"),
        "display_label": local_sg.get("display_label"),
        "symbol": local_sg.get("symbol"),
        "basis": iso.get("basis_text"),
        "origin": iso.get("origin"),
        "s": iso.get("s", iso.get("arms")),
        "i": iso.get("i"),
        "k_active": iso.get("k_active"),
    }
    physical_k_active, _ = _k_active_relation(
        web.subgroup.k_active, local.get("k_active")
    )
    matches = {
        "number": int(web.subgroup.number or 0) == int(local.get("number") or 0),
        "display_label": (
            _compact_text(web.subgroup.display_label) == _compact_text(local.get("display_label"))
            if magnetic_display
            else True
        ),
        "basis": _compact_text(web.subgroup.basis) == _compact_text(local.get("basis")),
        "origin": _compact_text(web.subgroup.origin) == _compact_text(local.get("origin")),
        "s": int(web.subgroup.s or 0) == int(local.get("s") or 0),
        "i": int(web.subgroup.i or 0) == int(local.get("i") or 0),
        "k_active": physical_k_active,
    }
    basis_change = _unimodular_basis_change(
        _rational_tokens(web.subgroup.basis),
        _rational_tokens(local.get("basis")),
    )
    basis_equivalent = (
        not matches["basis"]
        and basis_change is not None
        and all(matches[key] for key in matches if key != "basis")
    )
    embedding_equivalent = (
        operation_set_equivalent
        and basis_change is not None
        and all(
            matches[key]
            for key in matches
            if key not in {"basis", "origin"}
        )
    )
    return {
        "status": (
            "ok"
            if all(matches.values())
            else "basis_equivalent"
            if basis_equivalent
            else "embedding_equivalent"
            if embedding_equivalent
            else "fail"
        ),
        "matches": matches,
        **(
            {"basis_change": basis_change}
            if basis_equivalent or embedding_equivalent
            else {}
        ),
        "web": {
            "number": web.subgroup.number,
            "display_label": web.subgroup.display_label,
            "symbol": web.subgroup.symbol,
            "basis": web.subgroup.basis,
            "origin": web.subgroup.origin,
            "s": web.subgroup.s,
            "i": web.subgroup.i,
            "k_active": web.subgroup.k_active,
        },
        "local": local,
    }


def _atom_payload(web: CompleteModeDetails, local: dict[str, Any], tol: float) -> dict[str, Any]:
    local_atoms = local.get("undistorted_atoms") or local.get("distorted_atoms") or []
    child_sg = web.subgroup.number if web.subgroup is not None else None
    subgroup_details = local.get("subgroup_details") or {}
    magnetic_group = None
    magnetic_transform = None
    if isinstance(subgroup_details, dict):
        frame = selected_magnetic_frame(subgroup_details)
        parent_sg = web.parent.number if web.parent is not None else None
        if frame is None or parent_sg is None or int(frame["child_sg"]) != int(child_sg or 0):
            frame = None
        else:
            magnetic_group = int(frame["magnetic_group"])

            def magnetic_transform(xyz):
                return presentation_to_standard_child_cinter(
                    xyz,
                    parent_sg=int(parent_sg),
                    child_sg=int(frame["child_sg"]),
                    presentation_basis=frame["presentation_basis"],
                    presentation_origin=frame["presentation_origin"],
                    source_basis=frame["source_basis"],
                    source_origin=frame["source_origin"],
                    parent_setting_id=frame.get("parent_setting_id"),
                )
    return _structure_atom_payload(
        web.undistorted_atoms,
        local_atoms,
        tol,
        orbit_sg=child_sg,
        magnetic_group=magnetic_group,
        magnetic_transform=magnetic_transform,
    )


def _complete_basis_equivalent_atom_orbits(payload: object) -> bool:
    if not isinstance(payload, dict) or payload.get("count_match") is not True:
        return False
    return (
        payload.get("status") == "ok"
        and payload.get("site_kind_match") is True
        and payload.get("orbit_match") is True
    ) or (
        payload.get("status") == "site_gauge"
        and payload.get("kind_orbit_match") is True
    )


def _basis_equivalent_atom_payload(
    web: CompleteModeDetails,
    local: dict[str, Any],
    selected_state: object,
    subgroup: dict[str, Any],
    lattice: dict[str, Any],
    tol: float,
) -> dict[str, Any] | None:
    """Adopt request-frame atoms only after the complete lattice proof closes."""

    proof_result = {
        "selected_state": selected_state,
        "subgroup": subgroup,
        "undistorted_lattice": lattice,
    }
    if not _validated_basis_equivalent_lattice_result(proof_result):
        return None
    if not isinstance(selected_state, dict):
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
    if request_child <= 0 or request_child != local_child:
        return None
    subgroup_web = subgroup.get("web")
    subgroup_local = subgroup.get("local")
    if not isinstance(subgroup_web, dict) or not isinstance(subgroup_local, dict):
        return None
    try:
        if int(subgroup_web.get("number")) != request_child or int(subgroup_local.get("number")) != local_child:
            return None
    except (OverflowError, TypeError, ValueError):
        return None

    subgroup_details = local.get("subgroup_details")
    if not isinstance(subgroup_details, dict):
        return None
    local_atoms = local.get("undistorted_atoms") or local.get("distorted_atoms")
    frame = selected_magnetic_frame(subgroup_details)
    magnetic_group: int | None = None
    if isinstance(frame, dict):
        try:
            frame_child = int(frame["child_sg"])
            magnetic_group = int(frame["magnetic_group"])
        except (KeyError, OverflowError, TypeError, ValueError):
            return None
        if frame_child != request_child:
            return None
        transport = basis_equivalent_atom_frame_transport(
            local_atoms,
            basis_change=lattice.get("basis_change"),
            request_origin=request[3],
            local_origin=local_selected[3],
            child_space_group=request_child,
            magnetic_group=magnetic_group,
            request_basis=request[2],
        )
        provenance = "M:exact_gl3z_magnetic_atom_transport"
    else:
        web_subgroup = web.subgroup
        web_display = getattr(web_subgroup, "display_label", None)
        marker = subgroup_details.get("magnetic_subgroup")
        if (
            web_subgroup is None
            or _exact_integer_value(getattr(web_subgroup, "number", None)) != request_child
            or not isinstance(web_display, str)
            or _exact_integer_value(web_display) != request_child
        ):
            return None
        if marker is not None:
            if not isinstance(marker, dict) or marker.get("ordinary_number") is not None:
                return None
            marker_display = marker.get("display_label")
            if (
                _exact_integer_value(marker.get("number")) != request_child
                or (
                    marker_display is not None
                    and (
                        not isinstance(marker_display, str)
                        or _exact_integer_value(marker_display) != request_child
                    )
                )
            ):
                return None
        transport = basis_equivalent_ordinary_atom_frame_transport(
            local_atoms,
            basis_change=lattice.get("basis_change"),
            request_origin=request[3],
            local_origin=local_selected[3],
            child_space_group=request_child,
            tol=tol,
            request_basis=request[2],
        )
        provenance = "M:exact_gl3z_ordinary_atom_transport"
    if transport is None:
        return None
    transformed_rows = [
        {
            "label": atom.label,
            "site": atom.site,
            "xyz": atom.xyz,
            "_presentation_orbit_points": atom.orbit_points,
        }
        for atom in transport.atoms
    ]
    payload = _structure_atom_payload(
        web.undistorted_atoms,
        transformed_rows,
        tol,
        orbit_sg=request_child,
        magnetic_group=magnetic_group,
        magnetic_transform=(lambda xyz: xyz) if magnetic_group is not None else None,
    )
    if (
        not _complete_basis_equivalent_atom_orbits(payload)
        or payload.get("web_count") != len(transport.atoms)
        or payload.get("local_count") != len(transport.atoms)
    ):
        return None
    return {
        **payload,
        "validated": True,
        "provenance": provenance,
        "basis_change": transport.basis_change,
        "child_space_group": transport.child_space_group,
        **({"magnetic_group": transport.magnetic_group} if magnetic_group is not None else {}),
    }


def _validated_basis_equivalent_atom_result(
    atoms: object,
    *,
    basis_change: object,
    child_space_group: int,
) -> bool:
    if (
        not isinstance(atoms, dict)
        or not _complete_basis_equivalent_atom_orbits(atoms)
        or atoms.get("validated") is not True
        or atoms.get("provenance")
        not in {
            "M:exact_gl3z_magnetic_atom_transport",
            "M:exact_gl3z_ordinary_atom_transport",
        }
    ):
        return False
    expected = exact_unimodular_matrix3(basis_change)
    actual = exact_unimodular_matrix3(atoms.get("basis_change"))
    atom_child = _exact_integer_value(atoms.get("child_space_group"))
    if atom_child is None:
        return False
    common = (
        expected is not None
        and actual == expected
        and atom_child == child_space_group
    )
    if atoms.get("provenance") == "M:exact_gl3z_ordinary_atom_transport":
        return common and atoms.get("magnetic_group") is None
    magnetic_group = _exact_integer_value(atoms.get("magnetic_group"))
    return common and magnetic_group is not None and magnetic_group > 0


def _exact_integer_value(value: object) -> int | None:
    if isinstance(value, bool):
        return None
    try:
        fraction = value if isinstance(value, Fraction) else Fraction(str(value))
    except (OverflowError, TypeError, ValueError, ZeroDivisionError):
        return None
    return int(fraction) if fraction.denominator == 1 else None


def _validated_basis_equivalent_lattice_result(result: dict[str, Any]) -> bool:
    selected_state = result.get("selected_state")
    subgroup = result.get("subgroup")
    lattice = result.get("undistorted_lattice")
    if (
        not isinstance(selected_state, dict)
        or selected_state.get("status")
        not in {"equivalent_basis", "equivalent_embedding"}
        or not isinstance(subgroup, dict)
        or subgroup.get("status")
        not in {"basis_equivalent", "embedding_equivalent"}
        or not isinstance(lattice, dict)
        or lattice.get("status") != "basis_equivalent"
        or lattice.get("validated") is not True
        or lattice.get("provenance") != "M:exact_gl3z_metric_transport"
    ):
        return False

    changes = (
        exact_unimodular_matrix3(selected_state.get("basis_change")),
        exact_unimodular_matrix3(subgroup.get("basis_change")),
        exact_unimodular_matrix3(lattice.get("basis_change")),
    )
    if changes[0] is None or changes[0] != changes[1] or changes[0] != changes[2]:
        return False
    keys = ("a", "b", "c", "alpha", "beta", "gamma")
    matches = lattice.get("matches")
    local_matches = lattice.get("local_matches")
    web_values = lattice.get("web")
    local_values = lattice.get("local")
    transported_values = lattice.get("transported")
    if (
        not isinstance(matches, dict)
        or not isinstance(local_matches, dict)
        or any(matches.get(key) is not True for key in keys)
        or any(local_matches.get(key) is not True for key in keys)
        or not isinstance(web_values, dict)
        or not isinstance(local_values, dict)
        or not isinstance(transported_values, dict)
    ):
        return False
    try:
        scalar_values = (
            float(lattice["metric_residual"]),
            float(lattice["condition_number"]),
            *(float(values[key]) for values in (web_values, local_values, transported_values) for key in keys),
        )
    except (KeyError, OverflowError, TypeError, ValueError):
        return False
    residual, condition = scalar_values[:2]
    return (
        all(math.isfinite(value) for value in scalar_values)
        and residual >= 0.0
        and 1.0 <= condition <= 1e10
    )
