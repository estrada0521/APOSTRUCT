"""Render the local ISODISTORT mode-detail comparison surface."""

from __future__ import annotations

import ast
from fractions import Fraction
import math
import re
from typing import Any

import gemmi

from ISODISTORT.Assembled.Backend.modes.structure_runtime import (
    _mode_rows_grouped_by_presentation_orbits,
)


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
            vals = [str(int(item) // den) if int(item) % den == 0 else f"{int(item)}/{den}" for item in (x, y, z)]
            return "(" + ",".join(vals) + ")"
        except (TypeError, ValueError, ZeroDivisionError):
            pass
    if isinstance(value, list):
        return "(" + ",".join(str(item) for item in value) + ")"
    return str(value or "")


def _site_label(site: dict[str, Any]) -> str:
    wyckoff = str(site.get("wyckoff") or "")
    multiplicity = str(site.get("wyckoff_multiplicity") or site.get("multiplicity") or "")
    if wyckoff and wyckoff[0].isdigit():
        return wyckoff
    return f"{multiplicity}{wyckoff}" if multiplicity or wyckoff else ""


def _split_formula_xyz(formula: object) -> list[str] | None:
    text = str(formula or "").strip()
    if text.startswith("(") and text.endswith(")"):
        text = text[1:-1]
    parts = [part.strip() for part in text.split(",")]
    return parts if len(parts) == 3 else None


def _safe_eval_fraction_expr(expr: str, params: dict[str, Any]) -> float:
    expr = re.sub(r"(?<![A-Za-z0-9_])([+-]?(?:\d+(?:/\d+)?|\d*\.\d+))([xyz])\b", r"\1*\2", expr)
    tree = ast.parse(expr, mode="eval")

    def visit(node: ast.AST) -> Fraction:
        if isinstance(node, ast.Expression):
            return visit(node.body)
        if isinstance(node, ast.Constant) and isinstance(node.value, (int, float)):
            return Fraction(str(node.value))
        if isinstance(node, ast.Name):
            if node.id not in params:
                raise ValueError(f"unknown Wyckoff parameter {node.id!r}")
            return Fraction(str(params[node.id]))
        if isinstance(node, ast.UnaryOp) and isinstance(node.op, ast.USub):
            return -visit(node.operand)
        if isinstance(node, ast.UnaryOp) and isinstance(node.op, ast.UAdd):
            return visit(node.operand)
        if isinstance(node, ast.BinOp):
            left = visit(node.left)
            right = visit(node.right)
            if isinstance(node.op, ast.Add):
                return left + right
            if isinstance(node.op, ast.Sub):
                return left - right
            if isinstance(node.op, ast.Mult):
                return left * right
            if isinstance(node.op, ast.Div):
                return left / right
        raise ValueError(f"unsupported Wyckoff expression {expr!r}")

    return float(visit(tree))


def parent_display_xyz(site: dict[str, Any]) -> list[float]:
    """Return the Source-formula representative for a mapped CIF site.

    The CIF may contain any symmetry-equivalent representative, while
    the local display uses the representative implied by the matched Source
    Wyckoff formula. Use the stored formula and fitted parameters when
    available, then fall back to the raw CIF fraction.
    """

    parts = _split_formula_xyz(site.get("wyckoff_formula"))
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


def _atoms_grouped_by_parent_site(
    atoms: list[dict[str, Any]],
    parent_sites: list[dict[str, Any]],
) -> list[dict[str, Any]]:
    """Group child rows by parent site and a proved presentation-grid order."""
    labels = [str(site.get("label") or "") for site in parent_sites]
    ranked = sorted(
        ((label, index) for index, label in enumerate(labels) if label),
        key=lambda item: len(item[0]),
        reverse=True,
    )

    def parent_rank(atom: dict[str, Any]) -> int:
        label = str(atom.get("label") or "")
        return next(
            (index for parent, index in ranked if label == parent or label.startswith(parent + "_")),
            len(labels),
        )

    def site_key(atom: dict[str, Any], index: int) -> tuple[int, int, str, int]:
        site = str(atom.get("site") or "").strip()
        match = re.match(r"^(\d+)\s*([A-Za-z].*)?$", site)
        if match is None:
            return 1, 10**9, site.lower(), index
        return 0, int(match.group(1)), str(match.group(2) or "").lower(), index

    presentation_order: dict[int, int] = {}
    for parent in range(len(labels)):
        members = [
            (index, atom)
            for index, atom in enumerate(atoms)
            if parent_rank(atom) == parent
        ]
        if not members:
            continue
        parsed = [
            _atom_presentation_grid_partition(atom)
            for _index, atom in members
        ]
        if any(item is None for item in parsed):
            continue
        partitions = [item for item in parsed if item is not None]
        sizes = {size for size, _indices in partitions}
        if len(sizes) != 1:
            continue
        size = sizes.pop()
        flat = [index for _size, indices in partitions for index in indices]
        if len(flat) != size or set(flat) != set(range(size)):
            continue
        for (source_index, _atom), (_size, indices) in zip(
            members,
            partitions,
            strict=True,
        ):
            presentation_order[source_index] = min(indices)

    ordered = [
        atom
        for _parent, _site, atom in sorted(
            (
                parent_rank(atom),
                (
                    (0, presentation_order[index], "", index)
                    if index in presentation_order
                    else (1, *site_key(atom, index)[1:])
                ),
                atom,
            )
            for index, atom in enumerate(atoms)
        )
    ]
    counts: dict[int, int] = {}
    out: list[dict[str, Any]] = []
    for atom in ordered:
        parent = parent_rank(atom)
        if parent < len(labels) and labels[parent]:
            counts[parent] = counts.get(parent, 0) + 1
            atom = {**atom, "label": f"{labels[parent]}_{counts[parent]}"}
        out.append(atom)
    return out


def _atom_presentation_grid_partition(
    atom: dict[str, Any],
) -> tuple[int, tuple[int, ...]] | None:
    """Return one atom orbit's complete-grid index witness."""

    raw_size = atom.get("_presentation_grid_size")
    raw_indices = atom.get("_presentation_grid_indices")
    if (
        isinstance(raw_size, bool)
        or not isinstance(raw_size, int)
        or raw_size <= 0
        or not isinstance(raw_indices, (list, tuple))
        or not raw_indices
    ):
        return None
    indices: list[int] = []
    for raw_index in raw_indices:
        if (
            isinstance(raw_index, bool)
            or not isinstance(raw_index, int)
            or not 0 <= raw_index < raw_size
        ):
            return None
        indices.append(raw_index)
    if len(set(indices)) != len(indices):
        return None
    site_match = re.fullmatch(
        r"\s*(\d+)\s*[A-Za-z].*",
        str(atom.get("site") or ""),
    )
    if site_match is None or len(indices) != int(site_match.group(1)):
        return None
    return raw_size, tuple(sorted(indices))


def _mode_rows_with_operation_labels(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    """Label the first supercell translation for each parent operation."""
    counts: dict[str, int] = {}
    seen: set[tuple[str, int]] = set()
    out: list[dict[str, Any]] = []
    for row in rows:
        atom = str(row.get("atom") or "")
        prefix = atom.rsplit("_", 1)[0] if "_" in atom else atom
        record = row.get("_operation_record")
        try:
            operation = int(record[4])
        except (TypeError, ValueError, IndexError):
            out.append(row)
            continue
        key = (prefix, operation)
        displayed = ""
        if prefix and key not in seen:
            seen.add(key)
            counts[prefix] = counts.get(prefix, 0) + 1
            displayed = f"{prefix}_{counts[prefix]}"
        out.append({**row, "atom": displayed or None})
    return out


def _mode_rows_with_complete_label_runs(
    rows: list[dict[str, Any]],
    atoms: list[dict[str, Any]],
    *,
    atom_prefix: str | None,
) -> list[dict[str, Any]] | None:
    """Preserve a complete child-site label partition already carried by rows."""

    expected_prefix = str(atom_prefix or "")
    if not expected_prefix or not rows or not atoms:
        return None

    def label_prefix(value: object) -> str:
        return re.sub(r"_\d+$", "", str(value or ""))

    expected: list[tuple[str, int, int, tuple[int, ...]]] = []
    for atom in atoms:
        if not isinstance(atom, dict):
            return None
        label = atom.get("label")
        if not isinstance(label, str) or label_prefix(label) != expected_prefix:
            continue
        site_match = re.fullmatch(r"\s*(\d+)\s*[A-Za-z].*", str(atom.get("site") or ""))
        if site_match is None:
            return None
        multiplicity = int(site_match.group(1))
        if multiplicity <= 0:
            return None
        partition = _atom_presentation_grid_partition(atom)
        if partition is None:
            return None
        grid_size, grid_indices = partition
        expected.append((label, multiplicity, grid_size, grid_indices))
    if not expected or len({label for label, *_rest in expected}) != len(expected):
        return None
    grid_sizes = {grid_size for _label, _multiplicity, grid_size, _indices in expected}
    if len(grid_sizes) != 1:
        return None
    grid_size = grid_sizes.pop()
    expected_indices = [
        index
        for _label, _multiplicity, _size, indices in expected
        for index in indices
    ]
    if len(expected_indices) != grid_size or set(expected_indices) != set(
        range(grid_size)
    ):
        return None

    starts: list[tuple[int, str]] = []
    presentation_order: list[int] = []
    for index, row in enumerate(rows):
        if not isinstance(row, dict):
            return None
        presentation_index = row.get("_presentation_grid_index")
        if isinstance(presentation_index, bool) or not isinstance(
            presentation_index, int
        ):
            return None
        presentation_order.append(presentation_index)
        label = row.get("atom")
        if label is None or label == "":
            continue
        if not isinstance(label, str):
            return None
        starts.append((index, label))
    if len(rows) != grid_size or sorted(presentation_order) != list(range(grid_size)):
        return None
    if [label for _index, label in starts] != [label for label, *_rest in expected]:
        return None
    if not starts or starts[0][0] != 0:
        return None

    ends = [index for index, _label in starts[1:]] + [len(rows)]
    if any(
        end - start != multiplicity
        for ((start, _label), end, (_expected_label, multiplicity, _size, _indices)) in zip(
            starts, ends, expected, strict=True
        )
    ):
        return None
    for ((start, _label), end, (_expected_label, _multiplicity, _size, indices)) in zip(
        starts,
        ends,
        expected,
        strict=True,
    ):
        if {
            presentation_order[index]
            for index in range(start, end)
        } != set(indices):
            return None
    return [dict(row) for row in rows]


def _mode_rows_grouped_by_child_orbit(
    rows: list[dict[str, Any]],
    atoms: list[dict[str, Any]],
    child_sg: int | None,
    atom_prefix: str | None = None,
    tol: float = 2e-4,
) -> list[dict[str, Any]]:
    """Group complete-mode rows by the displayed child Wyckoff orbit."""
    labeled_rows = _mode_rows_with_complete_label_runs(
        rows,
        atoms,
        atom_prefix=atom_prefix,
    )
    if labeled_rows is not None:
        return labeled_rows
    representative_rows = _mode_rows_grouped_by_presentation_orbits(
        rows,
        atoms,
        atom_prefix=atom_prefix,
        child_sg=child_sg,
        tol=tol,
    )
    if representative_rows is not None:
        return representative_rows
    if child_sg is None or not atoms:
        return _mode_rows_with_operation_labels(rows)
    try:
        operations = list(gemmi.find_spacegroup_by_number(int(child_sg)).operations())
    except Exception:
        return _mode_rows_with_operation_labels(rows)

    def same_point(left: list[float], right: list[float]) -> bool:
        return all(abs(((left[axis] - right[axis] + 0.5) % 1.0) - 0.5) <= tol for axis in range(3))

    def label_prefix(value: object) -> str:
        return re.sub(r"_\d+$", "", str(value or ""))

    def points_for_representative(representative: list[float]) -> list[list[float]]:
        points: list[list[float]] = []
        for operation in operations:
            point = [float(value) % 1.0 for value in operation.apply_to_xyz(representative)]
            if not any(same_point(point, existing) for existing in points):
                points.append(point)
        return points

    rows_with_prefix: list[tuple[dict[str, Any], str]] = []
    current_prefix = str(atom_prefix or "")
    for row in rows:
        explicit_prefix = label_prefix(row.get("atom"))
        if explicit_prefix and atom_prefix is None:
            current_prefix = explicit_prefix
        rows_with_prefix.append((row, current_prefix))

    orbit_atoms = list(atoms)
    orbit_points = [
        points_for_representative([float(value) for value in atom.get("xyz") or (0, 0, 0)])
        for atom in orbit_atoms
    ]
    known_labels = {str(atom.get("label") or "") for atom in orbit_atoms}
    known_prefixes = {
        label.rsplit("_", 1)[0]
        for label in known_labels
        if re.search(r"_\d+$", label)
    }
    for row, prefix in rows_with_prefix:
        label = str(row.get("atom") or "")
        match = re.match(r"^(.+)_\d+$", label)
        if not match or label in known_labels or match.group(1) not in known_prefixes:
            continue
        representative = [float(value) for value in row.get("xyz") or (0, 0, 0)]
        if any(
            label_prefix(atom.get("label")) == prefix
            and any(same_point(representative, point) for point in points)
            for atom, points in zip(orbit_atoms, orbit_points, strict=True)
        ):
            continue
        orbit_atoms.append({"label": label, "xyz": representative})
        orbit_points.append(points_for_representative(representative))
        known_labels.add(label)

    grouped: list[list[dict[str, Any]]] = [[] for _ in orbit_atoms]
    unmatched: list[dict[str, Any]] = []
    for row, prefix in rows_with_prefix:
        xyz = [float(value) for value in row.get("xyz") or (0, 0, 0)]
        orbit_index = next(
            (
                index
                for index, (atom, points) in enumerate(
                    zip(orbit_atoms, orbit_points, strict=True)
                )
                if (not prefix or label_prefix(atom.get("label")) == prefix)
                and any(same_point(xyz, point) for point in points)
            ),
            None,
        )
        if orbit_index is None:
            unmatched.append(row)
        else:
            grouped[orbit_index].append(row)

    out: list[dict[str, Any]] = []
    for atom, orbit_rows in zip(orbit_atoms, grouped):
        for index, row in enumerate(orbit_rows):
            out.append({**row, "atom": str(atom.get("label") or "") if index == 0 else None})
    out.extend(_mode_rows_with_operation_labels(unmatched))
    return out


def _mode_parent_atom_prefix(label: object) -> str | None:
    match = re.search(r"\[([^:\]]+):[^:\]]+:(?:dsp|mag)\]", str(label or ""))
    return None if match is None else match.group(1)


def _ordinary_child_space_group_number(
    mode_details: dict[str, Any],
    selected: dict[str, Any],
) -> int | None:
    """Resolve the ordinary child group used to orbit positional rows."""

    def ordinary_number(*candidates: Any) -> int | None:
        for candidate in candidates:
            if isinstance(candidate, bool):
                continue
            try:
                number = int(candidate)
            except (TypeError, ValueError):
                continue
            if 1 <= number <= 230:
                return number
        return None

    subgroup_details = mode_details.get("subgroup_details")
    if isinstance(subgroup_details, dict):
        candidates: list[Any] = [
            subgroup_details.get("number"),
            subgroup_details.get("subgroup"),
        ]
        magnetic = subgroup_details.get("magnetic_subgroup")
        if isinstance(magnetic, dict):
            candidates.append(magnetic.get("ordinary_number"))
        if (number := ordinary_number(*candidates)) is not None:
            return number

    subgroup_state = selected.get("subgroup_state")
    if isinstance(subgroup_state, dict):
        state_subgroup = subgroup_state.get("subgroup")
        candidates = [subgroup_state.get("number")]
        if isinstance(state_subgroup, dict):
            candidates.extend(
                [state_subgroup.get("ordinary_number"), state_subgroup.get("number")]
            )
        if (number := ordinary_number(*candidates)) is not None:
            return number

    orderparam = selected.get("orderparam")
    isotropy = (
        orderparam.get("isotropy")
        if isinstance(orderparam, dict) and isinstance(orderparam.get("isotropy"), dict)
        else orderparam
    )
    subgroup = isotropy.get("subgroup") if isinstance(isotropy, dict) else None
    if isinstance(subgroup, dict):
        return ordinary_number(subgroup.get("ordinary_number"), subgroup.get("number"))
    return None


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
    rows = [[sum(float(row[k]) * P[k][axis] for k in range(3)) for axis in range(3)] for row in basis]
    def norm(vector: list[float]) -> float:
        return math.sqrt(sum(value * value for value in vector))

    la, lb, lc = norm(rows[0]), norm(rows[1]), norm(rows[2])

    def angle(u: list[float], w: list[float], lu: float, lw: float) -> float:
        if not lu or not lw:
            return 0.0
        d = sum(u[k] * w[k] for k in range(3)) / (lu * lw)
        return math.degrees(math.acos(max(-1.0, min(1.0, d))))

    return {
        "a": la, "b": lb, "c": lc,
        "alpha": angle(rows[1], rows[2], lb, lc),
        "beta": angle(rows[0], rows[2], la, lc),
        "gamma": angle(rows[0], rows[1], la, lb),
    }


def render_mode_detail_text(
    backend_state: dict[str, Any], *, precision: int = 5
) -> str:
    """Render only the local comparison-relevant mode-detail sections.

    This deliberately does not try to reproduce the full public-web HTML/stdout.
    In particular, the ISODISTORT amplitude table is omitted because it is
    interactive amplitude state rather than the mode-basis kernel.
    """

    input_info = backend_state.get("input") or {}
    parent = backend_state.get("space_group") or input_info.get("parent") or {}
    selected = backend_state.get("selected") or {}
    mode_details = selected.get("mode_details") or {}
    parent_lattice = input_info.get("lattice") or {}
    _opd = selected.get("orderparam")
    _opd_iso = _opd.get("isotropy") if isinstance(_opd, dict) and isinstance(_opd.get("isotropy"), dict) else (_opd or {})
    supercell_lattice = _supercell_lattice(parent_lattice, _opd_iso.get("basis"))
    parent_sites = [site for site in input_info.get("atom_sites", []) if isinstance(site, dict)]
    atoms = _atoms_grouped_by_parent_site(
        list(mode_details.get("undistorted_atoms") or mode_details.get("distorted_atoms") or []),
        parent_sites,
    )
    definitions = mode_details.get("displacive_definitions") or []
    magnetic_definitions = mode_details.get("magnetic_definitions") or []
    strain_definitions = mode_details.get("strain_definitions") or []
    child_sg = _ordinary_child_space_group_number(mode_details, selected)

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
            f"{str(atom.get('label') or ''):<6} {str(site):<5} "
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
            child_sg,
            atom_prefix=_mode_parent_atom_prefix(mode.get("label")),
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
                child_sg,
                atom_prefix=_mode_parent_atom_prefix(mode.get("label")),
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
