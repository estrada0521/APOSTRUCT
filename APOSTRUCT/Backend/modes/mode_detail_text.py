"""Render the APOSTRUCT mode-detail text surface."""

from __future__ import annotations

import math
import re
from typing import Any

from APOSTRUCT.Backend.fraction_expression import (
    evaluate_fraction_expression,
    split_coordinate_expression3,
)
from APOSTRUCT.Backend.modes.common import _site_label


def _fmt(value: object, digits: int = 5) -> str:
    return f"{float(value):.{digits}f}"


def _lattice_lines(lattice: dict[str, Any], *, digits: int = 5) -> list[str]:
    return [
        f"a={_fmt(lattice.get('a', 0), digits)}, b={_fmt(lattice.get('b', 0), digits)}, c={_fmt(lattice.get('c', 0), digits)}",
        f"alpha={_fmt(lattice.get('alpha', 0), digits)}, beta={_fmt(lattice.get('beta', 0), digits)}, gamma={_fmt(lattice.get('gamma', 0), digits)}",
    ]


def _fmt_basis(value: object) -> str:
    if not isinstance(value, list):
        return str(value or "")
    rows = []
    for row in value:
        if isinstance(row, list):
            rows.append("(" + ",".join(str(item) for item in row) + ")")
    return "{" + ",".join(rows) + "}" if rows else str(value)


def _fmt_origin(value: object) -> str:
    if isinstance(value, list) and len(value) == 4:
        x, y, z, den = value
        try:
            den = int(den)
            vals = [
                str(int(item) // den) if int(item) % den == 0 else f"{int(item)}/{den}"
                for item in (x, y, z)
            ]
            return "(" + ",".join(vals) + ")"
        except (TypeError, ValueError, ZeroDivisionError):
            pass
    if isinstance(value, list):
        return "(" + ",".join(str(item) for item in value) + ")"
    return str(value or "")


def _safe_eval_fraction_expr(expr: str, params: dict[str, Any]) -> float:
    expr = re.sub(
        r"(?<![A-Za-z0-9_])([+-]?(?:\d+(?:/\d+)?|\d*\.\d+))([xyz])\b", r"\1*\2", expr
    )
    return float(
        evaluate_fraction_expression(
            expr,
            params,
            "unknown Wyckoff parameter {name}",
            f"unsupported Wyckoff expression {expr!r}",
        )
    )


def parent_display_xyz(site: dict[str, Any]) -> list[float]:
    """Return the Source-formula representative for a mapped CIF site.

    The CIF may contain any symmetry-equivalent representative, while
    the local display uses the representative implied by the matched Source
    Wyckoff formula. Use the stored formula and fitted parameters when
    available, then fall back to the raw CIF fraction.
    """

    parts = split_coordinate_expression3(site.get("wyckoff_formula"))
    params = site.get("wyckoff_params") or {}
    if parts and isinstance(params, dict):
        try:
            return [_safe_eval_fraction_expr(part, params) % 1.0 for part in parts]
        except Exception:
            pass
    xyz = site.get("fract") or (0, 0, 0)
    return [float(value) % 1.0 for value in xyz]


def _render_parent_atoms(
    lines: list[str], atom_sites: list[dict[str, Any]], *, digits: int = 5
) -> None:
    lines.append("atom site    x         y         z")
    for site in atom_sites:
        xyz = parent_display_xyz(site)
        lines.append(
            f"{str(site.get('label') or ''):<6} {_site_label(site):<5} "
            f"{float(xyz[0]):{digits + 4}.{digits}f} "
            f"{float(xyz[1]):{digits + 4}.{digits}f} "
            f"{float(xyz[2]):{digits + 4}.{digits}f}"
        )


def _mode_rows_grouped_by_child_orbit(
    rows: list[dict[str, Any]],
    atoms: list[dict[str, Any]],
) -> list[dict[str, Any]]:
    """Validate and render the child-site blocks already carried by mode rows."""

    row_atom_ids = tuple(row.get("atom_id") for row in rows)
    if (
        not rows
        or not all(isinstance(atom_id, str) and atom_id for atom_id in row_atom_ids)
        or len(set(row_atom_ids)) != len(row_atom_ids)
    ):
        raise ValueError("mode definition has no complete atom identity")
    row_atom_id_set = set(row_atom_ids)
    expected_sites: list[tuple[str, tuple[str, ...]]] = []
    for atom in atoms:
        child_site = atom.get("child_site")
        atom_ids = atom.get("atom_ids")
        multiplicity = atom.get("multiplicity")
        if (
            not isinstance(child_site, str)
            or not child_site
            or atom.get("label") != child_site
            or not isinstance(atom_ids, list)
            or not atom_ids
            or not all(isinstance(atom_id, str) and atom_id for atom_id in atom_ids)
            or isinstance(multiplicity, bool)
            or not isinstance(multiplicity, int)
            or multiplicity != len(atom_ids)
        ):
            raise ValueError("undistorted child site has no complete atom identity")
        site_atom_ids = tuple(atom_ids)
        matched = row_atom_id_set.intersection(site_atom_ids)
        if matched:
            if matched != set(site_atom_ids):
                raise ValueError("mode definition splits an undistorted child site")
            expected_sites.append((child_site, site_atom_ids))
    expected = [
        (child_site, atom_id, member_order)
        for child_site, atom_ids in expected_sites
        for member_order, atom_id in enumerate(atom_ids)
    ]
    if len(rows) != len(expected):
        raise ValueError(
            "mode definition does not cover the undistorted child atoms: "
            f"{len(rows)} != {len(expected)}"
        )
    for row, (child_site, atom_id, member_order) in zip(rows, expected, strict=True):
        if row.get("child_site") != child_site or row.get("atom_id") != atom_id:
            raise ValueError("mode definition child-site identity/order changed")
        expected_label = child_site if member_order == 0 else None
        if (row.get("atom") or None) != expected_label:
            raise ValueError("mode definition child-site label run is malformed")
    return [dict(row) for row in rows]


def _web_subgroup_line(row: object) -> str | None:
    """Format a selected OPD row as the web ``Subgroup details`` single line.

    Matches ``06_completemodesdetails`` exactly: ``<num> <symbol>,
    basis={...}, origin=(...), s=N, i=M, k-active= (...)``.  Per upstream, the
    leading ``P1 (a,...)`` direction/OPD prefix is intentionally NOT included.
    The selected OPD row (``backend_state.selected.orderparam``) is verified (236/236
    web captures) to be byte-identical to what 06 shows here.
    """
    if not isinstance(row, dict):
        return None
    # The OPD row nests the row fields under ``isotropy`` (sibling of ``direction``).
    iso = row.get("isotropy") if isinstance(row.get("isotropy"), dict) else row
    sg = iso.get("subgroup")
    if not isinstance(sg, dict):
        return None
    number = sg.get("display_label") or sg.get("number", "")
    symbol = sg.get("symbol", "")
    basis = iso.get("basis_text") or _fmt_basis(iso.get("basis"))
    origin = iso.get("origin")
    if not isinstance(origin, str):
        origin = _fmt_origin(origin)
    s = iso.get("s", iso.get("arms", ""))
    i = iso.get("i", "")
    k_active = iso.get("k_active")
    if not k_active:
        vectors = iso.get("k_active_vectors")
        k_active = ",".join(vectors) if isinstance(vectors, list) else ""
    return f"{number} {symbol}, basis={{{basis}}}, origin={origin}, s={s}, i={i}, k-active= {k_active}"


def _supercell_lattice(parent: object, basis: object) -> dict[str, float] | object:
    """Undistorted supercell cell = OPD basis applied to the parent lattice.

    ``supercell_vectors = basis @ parent_vectors``; verified to match the web
    ``Undistorted/Distorted superstructure`` (and the ``undistorted supercell``
    line on the distort page).  Derived from parent lattice params + OPD basis,
    so it auto-corrects when the upstream basis presentation is fixed.
    """
    if not basis or not isinstance(parent, dict):
        return parent
    try:
        a, b, c = float(parent["a"]), float(parent["b"]), float(parent["c"])
        al = math.radians(float(parent["alpha"]))
        be = math.radians(float(parent["beta"]))
        ga = math.radians(float(parent["gamma"]))
    except (KeyError, TypeError, ValueError):
        return parent
    cx = c * math.cos(be)
    cy = c * (math.cos(al) - math.cos(be) * math.cos(ga)) / math.sin(ga)
    cz = math.sqrt(max(c * c - cx * cx - cy * cy, 0.0))
    P = [(a, 0.0, 0.0), (b * math.cos(ga), b * math.sin(ga), 0.0), (cx, cy, cz)]
    rows = [
        [sum(float(row[k]) * P[k][axis] for k in range(3)) for axis in range(3)]
        for row in basis
    ]

    def norm(vector: list[float]) -> float:
        return math.sqrt(sum(value * value for value in vector))

    la, lb, lc = norm(rows[0]), norm(rows[1]), norm(rows[2])

    def angle(u: list[float], w: list[float], lu: float, lw: float) -> float:
        if not lu or not lw:
            return 0.0
        d = sum(u[k] * w[k] for k in range(3)) / (lu * lw)
        return math.degrees(math.acos(max(-1.0, min(1.0, d))))

    return {
        "a": la,
        "b": lb,
        "c": lc,
        "alpha": angle(rows[1], rows[2], lb, lc),
        "beta": angle(rows[0], rows[2], la, lc),
        "gamma": angle(rows[0], rows[1], la, lb),
    }


def render_mode_detail_text(
    backend_state: dict[str, Any], *, precision: int = 5
) -> str:
    """Render only the local comparison-relevant mode-detail sections.

    This deliberately does not try to reproduce the full public-web HTML/stdout.
    In particular, the upstream interactive amplitude table is omitted because
    it is interactive amplitude state rather than the mode-basis kernel.
    """

    input_info = backend_state.get("input") or {}
    parent = backend_state.get("space_group") or input_info.get("parent") or {}
    selected = backend_state.get("selected") or {}
    mode_details = selected.get("mode_details") or {}
    status = mode_details.get("status")
    if status not in (None, "ok"):
        reason = mode_details.get("reason")
        raise ValueError(
            str(reason) if reason else f"mode details returned status {status!r}"
        )
    parent_lattice = input_info.get("lattice") or {}
    _opd = selected.get("orderparam")
    _opd_iso = (
        _opd.get("isotropy")
        if isinstance(_opd, dict) and isinstance(_opd.get("isotropy"), dict)
        else (_opd or {})
    )
    supercell_lattice = _supercell_lattice(parent_lattice, _opd_iso.get("basis"))
    parent_sites = [
        site for site in input_info.get("atom_sites", []) if isinstance(site, dict)
    ]
    raw_atoms = mode_details.get("undistorted_atoms")
    if not isinstance(raw_atoms, list) or not all(
        isinstance(atom, dict) for atom in raw_atoms
    ):
        raise ValueError("mode details have no canonical undistorted child sites")
    atoms = [dict(atom) for atom in raw_atoms]
    definitions = mode_details.get("displacive_definitions") or []
    magnetic_definitions = mode_details.get("magnetic_definitions") or []
    strain_definitions = mode_details.get("strain_definitions") or []
    lines: list[str] = []
    lines.append("Parent structure")
    lines.append(f"{parent.get('number', '')} {parent.get('symbol', '')}".rstrip())
    if not 5 <= precision <= 17:
        raise ValueError("mode-detail text precision must be between 5 and 17")
    vector_precision = 4 if precision == 5 else precision
    lines.extend(_lattice_lines(parent_lattice, digits=precision))
    lines.append("")
    _render_parent_atoms(lines, parent_sites, digits=precision)
    lines.append("")
    lines.append("Subgroup details")
    web_line = _web_subgroup_line(selected.get("orderparam"))
    if web_line:
        lines.append(web_line)
    else:
        lines.append("(not available yet)")
    lines.append("")
    lines.append("Undistorted superstructure")
    lines.extend(_lattice_lines(supercell_lattice, digits=precision))
    lines.append("atom site    x         y         z")
    for atom in atoms:
        xyz = atom.get("xyz") or (0, 0, 0)
        site = atom.get("site") or ""
        lines.append(
            f"{str(atom.get('child_site') or ''):<6} {str(site):<5} "
            f"{float(xyz[0]):{precision + 4}.{precision}f} "
            f"{float(xyz[1]):{precision + 4}.{precision}f} "
            f"{float(xyz[2]):{precision + 4}.{precision}f}"
        )
    lines.append("")
    lines.append("Displacive mode definitions")
    lines.append("")
    lines.append("atom    x        y        z        dx      dy      dz")
    for mode in definitions:
        norm = mode.get("normfactor")
        norm_text = "?" if norm is None else _fmt(norm, precision)
        lines.append(f"{mode.get('label', '')} normfactor = {norm_text}")
        for row in _mode_rows_grouped_by_child_orbit(
            list(mode.get("rows") or []),
            atoms,
        ):
            xyz = row.get("xyz") or (0, 0, 0)
            dxyz = row.get("dxyz") or (0, 0, 0)
            lines.append(
                f"{str(row.get('atom') or ''):<6} "
                f"{float(xyz[0]):{precision + 3}.{precision}f} "
                f"{float(xyz[1]):{precision + 3}.{precision}f} "
                f"{float(xyz[2]):{precision + 3}.{precision}f} "
                f"{float(dxyz[0]):{vector_precision + 3}.{vector_precision}f} "
                f"{float(dxyz[1]):{vector_precision + 3}.{vector_precision}f} "
                f"{float(dxyz[2]):{vector_precision + 3}.{vector_precision}f}"
            )
        lines.append("")
    if magnetic_definitions:
        lines.append("Magnetic mode definitions")
        lines.append("")
        lines.append("atom    x        y        z        dmx     dmy     dmz")
        for mode in magnetic_definitions:
            norm = mode.get("normfactor")
            norm_text = "?" if norm is None else _fmt(norm, precision)
            lines.append(f"{mode.get('label', '')} normfactor = {norm_text}")
            for row in _mode_rows_grouped_by_child_orbit(
                list(mode.get("rows") or []),
                atoms,
            ):
                xyz = row.get("xyz") or (0, 0, 0)
                vector = row.get("dxyz") or (0, 0, 0)
                lines.append(
                    f"{str(row.get('atom') or ''):<6} "
                    f"{float(xyz[0]):{precision + 3}.{precision}f} "
                    f"{float(xyz[1]):{precision + 3}.{precision}f} "
                    f"{float(xyz[2]):{precision + 3}.{precision}f} "
                    f"{float(vector[0]):{vector_precision + 3}.{vector_precision}f} "
                    f"{float(vector[1]):{vector_precision + 3}.{vector_precision}f} "
                    f"{float(vector[2]):{vector_precision + 3}.{vector_precision}f}"
                )
            lines.append("")
    if strain_definitions:
        lines.append("Parent-cell strain mode definitions")
        lines.append("")
        lines.append("    e1      e2      e3      e4      e5      e6")
        for mode in strain_definitions:
            norm = mode.get("normfactor")
            norm_text = "?" if norm is None else _fmt(norm, precision)
            components = list(mode.get("components") or (0, 0, 0, 0, 0, 0))
            lines.append(f"{mode.get('label', '')} normfactor = {norm_text}")
            lines.append(
                " "
                + " ".join(
                    f"{float(value):{vector_precision + 3}.{vector_precision}f}"
                    for value in components[:6]
                )
            )
    return "\n".join(lines).rstrip() + "\n"
