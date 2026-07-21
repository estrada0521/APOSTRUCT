"""Undistorted-structure construction helpers for modes.

Extracted mechanically from the former monolithic runtime.
"""

from __future__ import annotations

import ast
from fractions import Fraction
from functools import lru_cache
import math
import re
from typing import Any
import gemmi
from ISODISTORT.Assembled.Backend.exactmath import (
    fraction_matrix_inverse3,
    fraction_matrix_multiply3,
    integer_determinant3,
)
from ISODISTORT.Assembled.Backend.source.tables import SourceTables
from ISODISTORT.Assembled.Backend.modes.structure.wyckoff_split import undistorted_rows_from_wyckoff_split
from ISODISTORT.Assembled.Backend.modes.structure.ordinary_presentation import (
    assign_formula_source_occurrences_to_presentation_grid,
    match_formula_assignments_to_geometric_grid_orbits,
)
from ISODISTORT.Assembled.Backend.modes.structure.magnetic_wyckoff import (
    group_ordinary_orbits_magnetic,
    identify_magnetic_wyckoff,
    identify_magnetic_wyckoff_branch,
    magnetic_group_setting,
    magnetic_orbit_points,
)
from ISODISTORT.Assembled.Backend.modes.structure.magnetic_presentation import (
    presentation_branch_labels_from_correspondence,
    presentation_grid_rows_for_magnetic_groups,
    selected_magnetic_correspondence_standard_points,
)
from ISODISTORT.Assembled.Backend.modes.engine.decoder import ModeDataDecoder

from ISODISTORT.Assembled.Backend.modes.common import (
    _assembled_data,
    _float_matrix_inverse_3,
    _fold01,
    _fold_fractional_xyz,
    _frac_close,
    _fraction_matmul,
    _fraction_matrix_inverse_3,
    _fraction_mod01,
    _fraction_row_multiply,
    _fraction_vecadd,
    _fraction_vecsub,
    _integer_basis_tuple,
    _isotropy_from_opd_row,
    _mat4_multiply_fraction,
    _matrix_from_basis_tuple,
    _mode_decoder,
    _normalize_setting_matrix,
    _origin_record_from_any,
    _origin_record_vector,
    _origin_vector,
    _row_multiply,
)
from ISODISTORT.Assembled.Backend.modes.site_transport import (
    _parent_point_from_default,
    _parent_point_to_default,
    _parent_setting_bridge,
    _source_default_wyckoff_params,
)


def _mode_rows_grouped_by_presentation_orbits(
    rows: list[dict[str, Any]],
    atoms: list[dict[str, Any]],
    *,
    atom_prefix: str | None,
    child_sg: int | None = None,
    tol: float = 2e-4,
) -> list[dict[str, Any]] | None:
    """Order complete mode rows by an exact Source orbit partition."""

    def prefix(value: object) -> str:
        return re.sub(r"_\d+$", "", str(value or ""))

    def finite_xyz(value: object) -> tuple[float, float, float] | None:
        if not isinstance(value, (list, tuple)) or len(value) != 3:
            return None
        try:
            xyz = tuple(float(component) for component in value)
        except (TypeError, ValueError):
            return None
        return xyz if all(math.isfinite(component) for component in xyz) else None

    def same_point(
        left: tuple[float, float, float],
        right: tuple[float, float, float],
    ) -> bool:
        return all(
            abs(((left[axis] - right[axis] + 0.5) % 1.0) - 0.5) <= tol
            for axis in range(3)
        )

    expected_prefix = str(atom_prefix or "")
    orbits: list[tuple[str, tuple[tuple[float, float, float], ...]]] = []
    for atom in atoms:
        if not isinstance(atom, dict):
            return None
        label = str(atom.get("label") or "")
        if not label or (expected_prefix and prefix(label) != expected_prefix):
            continue
        raw_orbit = atom.get("_mode_row_orbit_points")
        if raw_orbit is None:
            raw_orbit = atom.get("_presentation_orbit_points")
        if isinstance(raw_orbit, (list, tuple)) and raw_orbit:
            parsed_orbit = tuple(finite_xyz(point) for point in raw_orbit)
            if any(point is None for point in parsed_orbit):
                return None
        else:
            representative = finite_xyz(atom.get("xyz"))
            site_match = re.fullmatch(r"\s*(\d+)\s*[A-Za-z].*", str(atom.get("site") or ""))
            if representative is None or child_sg is None or site_match is None:
                return None
            try:
                operations = _child_symmetry_ops(int(child_sg))
            except (KeyError, RuntimeError, TypeError, ValueError):
                return None
            parsed: list[tuple[float, float, float]] = []
            for operation in operations:
                point = finite_xyz(_apply_child_op(operation, list(representative)))
                if point is None:
                    return None
                if not any(same_point(point, seen) for seen in parsed):
                    parsed.append(point)
            if len(parsed) != int(site_match.group(1)):
                return None
            parsed_orbit = tuple(parsed)
        orbits.append((label, parsed_orbit))  # type: ignore[arg-type]
    if not orbits or len({label for label, _orbit in orbits}) != len(orbits):
        return None

    row_points: list[tuple[float, float, float]] = []
    for row in rows:
        if not isinstance(row, dict):
            return None
        xyz = finite_xyz(row.get("xyz"))
        if xyz is None:
            return None
        row_points.append(xyz)

    def row_order_for(field: str) -> list[int]:
        order: list[int] = []
        for row in rows:
            value = row.get(field)
            if isinstance(value, bool) or not isinstance(value, int) or value < 0:
                return []
            order.append(value)
        return order

    presentation_order = row_order_for("_presentation_grid_index")
    source_order = row_order_for("_source_raw_index")
    if (
        len(presentation_order) == len(rows)
        and sorted(presentation_order) == list(range(len(rows)))
    ):
        row_order = presentation_order
    elif len(source_order) == len(rows):
        row_order = source_order
    else:
        return None

    unused = set(range(len(rows)))
    bucket_width: float | None = None
    if isinstance(tol, (int, float)) and not isinstance(tol, bool):
        candidate_width = 2.0 * float(tol)
        if (
            math.isfinite(candidate_width)
            and candidate_width > 0.0
            and math.isfinite(1.0 / candidate_width)
        ):
            bucket_width = candidate_width
    bucket_count = (
        max(1, int(math.floor(1.0 / bucket_width)))
        if bucket_width is not None
        else None
    )
    row_buckets: dict[tuple[int, int, int], list[int]] = {}
    if bucket_count is not None and bucket_width is not None:
        for index, point in enumerate(row_points):
            key = tuple(
                int(math.floor((point[axis] % 1.0) / bucket_width))
                % bucket_count
                for axis in range(3)
            )
            row_buckets.setdefault(key, []).append(index)  # type: ignore[arg-type]

    def indexed_matches(point: tuple[float, float, float]) -> list[int]:
        if bucket_count is None or bucket_width is None:
            return [
                index
                for index in sorted(unused)
                if same_point(point, row_points[index])
            ]
        center = tuple(
            int(math.floor((point[axis] % 1.0) / bucket_width)) % bucket_count
            for axis in range(3)
        )
        matches: list[int] = []
        visited_buckets: set[tuple[int, int, int]] = set()
        for first in (-1, 0, 1):
            for second in (-1, 0, 1):
                for third in (-1, 0, 1):
                    key = (
                        (center[0] + first) % bucket_count,
                        (center[1] + second) % bucket_count,
                        (center[2] + third) % bucket_count,
                    )
                    if key in visited_buckets:
                        continue
                    visited_buckets.add(key)
                    for index in row_buckets.get(key, ()):
                        if index in unused and same_point(
                            point, row_points[index]
                        ):
                            matches.append(index)
        return matches

    grouped: list[tuple[str, int, list[int]]] = []
    for label, orbit in orbits:
        indices: list[int] = []
        for point in orbit:
            matches = indexed_matches(point)
            if len(matches) != 1:
                return None
            unused.remove(matches[0])
            indices.append(matches[0])
        grouped.append((label, min(row_order[index] for index in indices), sorted(indices)))
    if unused:
        return None
    if len({raw_index for _label, raw_index, _indices in grouped}) != len(grouped):
        return None

    out: list[dict[str, Any]] = []
    for ordinal, (_label, _raw_index, indices) in enumerate(
        sorted(grouped, key=lambda item: item[1]),
        start=1,
    ):
        label = f"{expected_prefix}_{ordinal}"
        for orbit_index, row_index in enumerate(indices):
            out.append(
                {
                    **rows[row_index],
                    "atom": label if orbit_index == 0 else None,
                }
            )
    return out


def _split_formula_xyz(formula: Any) -> tuple[str, str, str] | None:
    if not isinstance(formula, str):
        return None
    text = formula.strip()
    if text.startswith("(") and text.endswith(")"):
        text = text[1:-1]
    parts = tuple(part.strip() for part in text.split(","))
    return parts if len(parts) == 3 else None



def _eval_fraction_formula(expr: str, params: dict[str, Any]) -> float:
    text = str(expr).strip().replace("−", "-")
    text = re.sub(
        r"(?<![A-Za-z0-9_])([+-]?(?:\d+(?:/\d+)?|\d*\.\d+))([xyz])\b",
        r"\1*\2",
        text,
    )
    tree = ast.parse(text, mode="eval")

    def eval_node(node: ast.AST) -> float:
        if isinstance(node, ast.Expression):
            return eval_node(node.body)
        if isinstance(node, ast.Constant) and isinstance(node.value, (int, float)):
            return float(node.value)
        if isinstance(node, ast.Name):
            return float(params.get(node.id, 0.0))
        if isinstance(node, ast.UnaryOp) and isinstance(node.op, ast.USub):
            return -eval_node(node.operand)
        if isinstance(node, ast.UnaryOp) and isinstance(node.op, ast.UAdd):
            return eval_node(node.operand)
        if isinstance(node, ast.BinOp):
            left = eval_node(node.left)
            right = eval_node(node.right)
            if isinstance(node.op, ast.Add):
                return left + right
            if isinstance(node.op, ast.Sub):
                return left - right
            if isinstance(node.op, ast.Mult):
                return left * right
            if isinstance(node.op, ast.Div):
                return left / right
        raise ValueError(f"unsupported Wyckoff formula expression: {expr!r}")

    return eval_node(tree) % 1.0



def _xyz_from_wyckoff_match(match: dict[str, Any], fallback_xyz: tuple[float, float, float]) -> list[float]:
    parts = _split_formula_xyz(match.get("formula"))
    params = match.get("params") or {}
    if parts is None or not isinstance(params, dict):
        return list(fallback_xyz)
    try:
        return [_eval_fraction_formula(part, params) for part in parts]
    except Exception:
        return list(fallback_xyz)

@lru_cache(maxsize=4096)
def _matched_child_site_display(
    child_sg: int,
    fract: tuple[float, float, float],
    *,
    multiplicity: int | None = None,
) -> tuple[str, tuple[float, float, float]] | None:
    query: dict[str, Any] = {"fract": [str(value) for value in fract]}
    if multiplicity is not None and int(multiplicity) > 0:
        query["multiplicity"] = str(int(multiplicity))
    data = _assembled_data()
    try:
        setting = int(data.default_inter_setting_record(int(child_sg))["id"])
        match = data.match_inter_wyckoff_site(
            int(child_sg),
            query,
            setting,
            tol=2e-5,
        )
    except Exception:
        match = None
    if not match:
        try:
            match = data.match_wyckoff_site(int(child_sg), query, tol=2e-5)
        except Exception:
            match = None
    if not match:
        return None
    label = str(match.get("label") or "")
    multiplicity = str(match.get("multiplicity") or "")
    site = f"{multiplicity}{label}" if multiplicity or label else ""
    if not site:
        return None
    display = _xyz_from_wyckoff_match(match, fract)
    return site, tuple(float(value) for value in display)


def _child_site_display(child_sg: int | None, xyz: Any, fallback: str) -> tuple[str, list[float]]:
    fract = _fold_fractional_xyz(xyz)
    if child_sg is None or fract is None:
        return fallback, list(fract) if fract is not None else list(xyz or [])
    match = _matched_child_site_display(int(child_sg), fract)
    return (match[0], list(match[1])) if match is not None else (fallback, list(fract))



def _display_basis_origin_matrix(
    basis: list[list[float]],
    origin: tuple[int, int, int, int] | Any,
) -> list[list[Fraction]]:
    origin_values = _origin_vector(origin)
    return [
        [Fraction(str(basis[row][col])).limit_denominator(1000000) for col in range(3)] + [Fraction(0)]
        for row in range(3)
    ] + [[Fraction(str(origin_values[col])).limit_denominator(1000000) for col in range(3)] + [Fraction(1)]]



def _internal_split_basis_origin(
    data: SourceTables,
    *,
    parent_sg: int,
    child_sg: int,
    display_basis: list[list[float]],
    display_origin: Any,
) -> tuple[list[list[float]], tuple[int, int, int, int]] | None:
    try:
        parent_choice = int(data.space["ispace_inter_choice"][int(parent_sg) - 1])
        child_choice = int(data.space["ispace_inter_choice"][int(child_sg) - 1])
        left = _mat4_multiply_fraction(
            _normalize_setting_matrix(data._inter_matrix_by_id(child_choice, 1)),  # noqa: SLF001
            _normalize_setting_matrix(data._space_cinter_base_matrix(int(child_sg), 0)),  # noqa: SLF001
        )
        right = _mat4_multiply_fraction(
            _normalize_setting_matrix(data._space_cinter_base_matrix(int(parent_sg), 1)),  # noqa: SLF001
            _normalize_setting_matrix(data._inter_matrix_by_id(parent_choice, 0)),  # noqa: SLF001
        )
        internal = _mat4_multiply_fraction(
            _mat4_multiply_fraction(data._mat4_inverse_fraction(left), _display_basis_origin_matrix(display_basis, display_origin)),  # noqa: SLF001
            data._mat4_inverse_fraction(right),  # noqa: SLF001
        )
    except Exception:
        return None
    basis = [[float(internal[row][col]) for col in range(3)] for row in range(3)]
    origin_values = [internal[3][col] for col in range(3)]
    den = 1
    for value in origin_values:
        den = math.lcm(den, value.denominator)
    origin = tuple(int(value * den) for value in origin_values) + (den,)
    return basis, origin  # type: ignore[return-value]



def _child_symmetry_ops(child_sg: int) -> tuple[gemmi.Op, ...]:
    return tuple(gemmi.find_spacegroup_by_number(int(child_sg)).operations())



def _apply_child_op(op: gemmi.Op, xyz: list[float]) -> list[float]:
    return [_fold01(value) for value in op.apply_to_xyz(list(xyz))]


def _periodic_position_matches(
    point: list[float],
    positions: list[list[float]],
    *,
    tolerance: float = 2e-5,
) -> list[int]:
    matches: list[int] = []
    for index, candidate in enumerate(positions):
        if all(
            abs(((float(point[axis]) - float(candidate[axis]) + 0.5) % 1.0) - 0.5)
            <= tolerance
            for axis in range(3)
        ):
            matches.append(index)
    return matches


def _child_orbit_display_layout(
    positions: list[list[float]],
    child_sg: int | None,
) -> tuple[list[int], set[int]] | None:
    """Return a stable orbit-major row order and its representative rows."""

    if child_sg is None or not positions:
        return None
    try:
        ops = _child_symmetry_ops(int(child_sg))
    except Exception:
        return None
    unused = set(range(len(positions)))
    orbit_rows: list[list[int]] = []
    while unused:
        seed = min(unused)
        orbit: set[int] = set()
        for operation in ops:
            matches = _periodic_position_matches(
                _apply_child_op(operation, positions[seed]),
                positions,
            )
            if len(matches) != 1:
                return None
            orbit.add(matches[0])
        if seed not in orbit or not orbit <= unused:
            return None
        ordered_orbit = sorted(orbit)
        orbit_rows.append(ordered_orbit)
        unused.difference_update(orbit)
    order = [index for orbit in orbit_rows for index in orbit]
    if len(order) != len(positions):
        return None
    return order, {orbit[0] for orbit in orbit_rows}


def _sort_child_orbit_rows_for_display(
    rows: list[dict[str, Any]],
    *,
    label_prefix: str,
) -> list[dict[str, Any]]:
    """Sort split sites by multiplicity/letter and renumber their display labels."""

    def key(item: tuple[int, dict[str, Any]]) -> tuple[int, str, int]:
        index, row = item
        site = str(row.get("site") or "")
        match = re.fullmatch(r"\s*(\d+)\s*([A-Za-z].*)?", site)
        multiplicity = int(match.group(1)) if match else 10**9
        letter = str(match.group(2) or "") if match else site
        return multiplicity, letter, index

    ordered = [dict(row) for _, row in sorted(enumerate(rows), key=key)]
    for ordinal, row in enumerate(ordered, start=1):
        row["label"] = f"{label_prefix}_{ordinal}"
    return ordered



def _child_orbit_representative_rows(
    *,
    label_prefix: str,
    positions: list[list[float]],
    child_sg: int | None,
    fallback_site: str,
    require_closed: bool = False,
) -> list[dict[str, Any]]:
    if child_sg is None or not positions:
        return []
    try:
        ops = _child_symmetry_ops(int(child_sg))
    except Exception:
        return []
    if require_closed and any(
        len(_periodic_position_matches(position, positions)) != 1
        for position in positions
    ):
        return []
    unused = set(range(len(positions)))
    rows: list[dict[str, Any]] = []
    while unused:
        seed = min(unused)
        orbit_points = [_apply_child_op(op, positions[seed]) for op in ops]
        orbit_indices: set[int] = set()
        for point in orbit_points:
            matches = set(_periodic_position_matches(point, positions))
            if require_closed and len(matches) != 1:
                return []
            orbit_indices.update(matches)
        if require_closed and (seed not in orbit_indices or not orbit_indices <= unused):
            return []
        orbit_indices.intersection_update(unused)
        if not orbit_indices:
            if require_closed:
                return []
            orbit_indices = {seed}
        if require_closed:
            representative = min(orbit_indices)
            fast_display = _matched_child_site_display(
                int(child_sg),
                tuple(float(value) for value in positions[representative]),
                multiplicity=len(orbit_indices),
            )
            if fast_display is not None:
                site_label, display_xyz = fast_display
                match = re.match(r"\s*(\d+)", str(site_label or ""))
                if match is not None and int(match.group(1)) == len(orbit_indices):
                    rows.append(
                        {
                            "label": f"{label_prefix}_{len(rows) + 1}",
                            "site": site_label,
                            "xyz": [_fold01(float(value)) for value in display_xyz],
                            "_source_index": representative,
                        }
                    )
                    unused.difference_update(orbit_indices)
                    continue
        display_by_index = {
            index: _child_site_display(child_sg, positions[index], fallback_site)
            for index in orbit_indices
        }

        def representative_key(index: int) -> tuple[int, int]:
            site_label = display_by_index[index][0]
            match = re.match(r"\s*(\d+)", str(site_label or ""))
            multiplicity = int(match.group(1)) if match else 10**9
            return multiplicity, index

        representative = min(orbit_indices, key=representative_key)
        site_label, display_xyz = display_by_index[representative]
        display_xyz = [_fold01(float(value)) for value in display_xyz]
        rows.append(
            {
                "label": f"{label_prefix}_{len(rows) + 1}",
                "site": site_label,
                "xyz": display_xyz,
                "_source_index": representative,
            }
        )
        unused.difference_update(orbit_indices)
    return rows



def _is_identity_basis_origin(
    basis: list[list[float]] | None,
    origin: Any,
) -> bool:
    if basis is None:
        return False
    identity = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
    return (
        all(abs(float(basis[row][col]) - identity[row][col]) <= 1e-8 for row in range(3) for col in range(3))
        and all(abs(_fold01(value)) <= 1e-8 for value in _origin_vector(origin))
    )



def _translation_child_shifts(
    basis: list[list[float]],
    *,
    max_size: int = 4096,
) -> list[tuple[Fraction, Fraction, Fraction]]:
    inverse = _fraction_matrix_inverse_3(basis)
    if inverse is None:
        return []
    generators = [
        tuple(_fraction_mod01(value) for value in _fraction_row_multiply(unit, inverse))
        for unit in (
            (Fraction(1), Fraction(0), Fraction(0)),
            (Fraction(0), Fraction(1), Fraction(0)),
            (Fraction(0), Fraction(0), Fraction(1)),
        )
    ]
    seen: set[tuple[Fraction, Fraction, Fraction]] = {(Fraction(0), Fraction(0), Fraction(0))}
    queue = [(Fraction(0), Fraction(0), Fraction(0))]
    while queue:
        current = queue.pop(0)
        for generator in generators:
            nxt = tuple(_fraction_mod01(current[axis] + generator[axis]) for axis in range(3))
            if nxt in seen:
                continue
            seen.add(nxt)
            queue.append(nxt)
            if len(seen) > max_size:
                return [(Fraction(0), Fraction(0), Fraction(0))]
    return sorted(seen)



def _parent_orbit_points(
    parent_sg: int,
    xyz: tuple[float, float, float],
    operations: list[str] | None = None,
    *,
    preserve_lattice_translation: bool = False,
) -> list[list[float]]:
    if operations:
        ops = []
        for triplet in operations:
            try:
                ops.append(gemmi.Op(str(triplet)))
            except (RuntimeError, TypeError, ValueError):
                continue
    else:
        try:
            ops = list(gemmi.find_spacegroup_by_number(int(parent_sg)).operations())
        except Exception:
            return [list(xyz)]
    points: list[list[float]] = []
    for op in ops:
        transformed = [float(value) for value in op.apply_to_xyz(list(xyz))]
        point = (
            transformed
            if preserve_lattice_translation
            else [_fold01(value) for value in transformed]
        )
        if not any(_frac_close(point, seen, tol=1e-7) for seen in points):
            points.append(point)
    return points



def _complete_child_atom_layout(
    *,
    parent_sg: int,
    child_sg: int | None,
    parent_xyz: tuple[float, float, float] | None,
    basis: list[list[float]] | None,
    origin: Any,
    label_prefix: str,
    fallback_site: str,
    parent_operations: list[str] | None = None,
) -> tuple[list[list[float]], dict[int, str]] | None:
    """Return the full selected-cell atom grid and child-orbit labels."""

    if child_sg is None or parent_xyz is None or basis is None:
        return None
    inverse = _float_matrix_inverse_3(basis)
    if inverse is None:
        return None
    origin_vector = _origin_vector(origin)
    positions: list[list[float]] = []
    for parent_point in _parent_orbit_points(
        parent_sg,
        parent_xyz,
        parent_operations,
        preserve_lattice_translation=True,
    ):
        child_base = _row_multiply(
            [parent_point[index] - origin_vector[index] for index in range(3)],
            inverse,
        )
        for shift in _translation_child_shifts(basis):
            child_point = [_fold01(child_base[index] + float(shift[index])) for index in range(3)]
            if not any(_frac_close(child_point, seen, tol=1e-7) for seen in positions):
                positions.append(child_point)
    representatives = _child_orbit_representative_rows(
        label_prefix=label_prefix,
        positions=positions,
        child_sg=child_sg,
        fallback_site=fallback_site,
    )
    labels = {
        int(row["_source_index"]): str(row["label"])
        for row in representatives
        if row.get("_source_index") is not None
    }
    return positions, labels


def _canonical_site_fractional_xyz(
    site: dict[str, Any],
) -> tuple[float, float, float] | None:
    """Keep the parent representative's lattice-translation provenance."""

    matched = site.get("representative_matched_fract")
    if isinstance(matched, (list, tuple)) and len(matched) == 3:
        try:
            parsed = tuple(float(value) for value in matched)
        except (TypeError, ValueError):
            parsed = ()
        if len(parsed) == 3 and all(math.isfinite(value) for value in parsed):
            return parsed  # type: ignore[return-value]
    return _fold_fractional_xyz(site.get("fract"))



def _mode_rows_on_complete_atom_layout(
    rows: list[dict[str, Any]],
    layout: tuple[list[list[float]], dict[int, str]] | None,
) -> list[dict[str, Any]]:
    """Join projected displacement support onto the complete selected-cell grid."""

    if layout is None:
        return rows
    positions, labels = layout
    tolerance = 1e-6
    # A 2*tol cell keeps every tol-close pair in the same or an adjacent bin.
    bucket_width = 2 * tolerance
    bucket_count = 500_000
    buckets: dict[tuple[int, int, int], list[int]] = {}
    row_xyz: dict[int, list[float]] = {}
    for row_index, row in enumerate(rows):
        xyz = row.get("xyz")
        if not isinstance(xyz, list):
            continue
        folded = tuple(float(xyz[index]) % 1.0 for index in range(3))
        if not all(math.isfinite(value) for value in folded):
            continue
        key = tuple(
            int(math.floor(value / bucket_width)) % bucket_count
            for value in folded
        )
        buckets.setdefault(key, []).append(row_index)  # type: ignore[arg-type]
        row_xyz[row_index] = xyz
    aligned: list[dict[str, Any]] = []
    for index, position in enumerate(positions):
        folded = tuple(float(position[axis]) % 1.0 for axis in range(3))
        candidate_indices: set[int] = set()
        if len(folded) == 3 and all(math.isfinite(value) for value in folded):
            center = tuple(
                int(math.floor(value / bucket_width)) % bucket_count
                for value in folded
            )
            for first in (-1, 0, 1):
                for second in (-1, 0, 1):
                    for third in (-1, 0, 1):
                        candidate_indices.update(
                            buckets.get(
                                (
                                    (center[0] + first) % bucket_count,
                                    (center[1] + second) % bucket_count,
                                    (center[2] + third) % bucket_count,
                                ),
                                (),
                            )
                        )
        matches = [
            row_index
            for row_index in sorted(candidate_indices)
            if _frac_close(position, row_xyz[row_index], tol=tolerance)
        ]
        source_index = max(
            matches,
            key=lambda row_index: sum(
                float(value) ** 2 for value in (rows[row_index].get("dxyz") or ())
            ),
            default=None,
        )
        source = None if source_index is None else rows[source_index]
        aligned.append(
            {
                "atom": labels.get(index),
                "xyz": list(position),
                "dxyz": (
                    [0.0, 0.0, 0.0]
                    if source is None
                    else [float(value) for value in source.get("dxyz") or (0.0, 0.0, 0.0)]
                ),
                "_presentation_grid_index": index,
                **(
                    {"_operation_record": source["_operation_record"]}
                    if source is not None and source.get("_operation_record") is not None
                    else {}
                ),
                **(
                    {"_source_raw_index": source["_source_raw_index"]}
                    if source is not None and source.get("_source_raw_index") is not None
                    else {}
                ),
            }
        )
    return aligned



def _web_magnetic_occurrence_gauge(
    decoder: ModeDataDecoder,
    sg: int,
    rows: list[dict[str, Any]],
) -> list[dict[str, Any]]:
    """Bridge faithful axial rows to Web's occurrence-frame convention.

    Each operation record is the Source-canonical-site to occurrence witness.
    Its determinant is invariant under the parent-to-child setting conjugation,
    so no point-operation ordinal is part of this presentation rule.
    """

    out: list[dict[str, Any]] = []
    for row in rows:
        record = row.get("_operation_record")
        if not isinstance(record, (list, tuple)) or len(record) != 5:
            out.append(row)
            continue
        lattice = int(decoder.iso.space["ispace_lattice"][int(sg) - 1])
        code = int(
            decoder.iso.space["ipoint_op_psettings"][(lattice - 1) * 72 + int(record[4]) - 1]
        )
        matrix = []
        for _ in range(9):
            matrix.append(code % 3 - 1)
            code //= 3
        determinant = integer_determinant3(matrix)
        if determinant >= 0:
            out.append(row)
        else:
            out.append(
                {
                    **row,
                    "dxyz": [-float(value) for value in row.get("dxyz") or (0.0, 0.0, 0.0)],
                }
            )
    return out



def _presentation_child_points(
    *,
    parent_sg: int,
    parent_xyz: tuple[float, float, float] | None,
    basis: list[list[float]] | None,
    origin: Any,
) -> list[list[float]]:
    """Expand one parent representative into the selected presentation cell."""

    if basis is None or parent_xyz is None:
        return []
    inverse = _float_matrix_inverse_3(basis)
    if inverse is None:
        return []
    origin_vector = _origin_vector(origin)
    shifts = _translation_child_shifts(basis)
    positions: list[list[float]] = []
    for parent_point in _parent_orbit_points(parent_sg, parent_xyz):
        child_base = _row_multiply(
            [parent_point[index] - origin_vector[index] for index in range(3)],
            inverse,
        )
        for shift in shifts:
            child_point = [_fold01(child_base[index] + float(shift[index])) for index in range(3)]
            if not any(_frac_close(child_point, seen, tol=1e-7) for seen in positions):
                positions.append(child_point)
    return positions



def _ordinary_rows_from_complete_child_grid(
    *,
    parent_sg: int,
    child_sg: int,
    label_prefix: str,
    parent_xyz: tuple[float, float, float] | None,
    fallback_site: str,
    source_basis: list[list[float]] | None,
    source_origin: Any,
    presentation_basis: list[list[float]] | None,
    presentation_origin: Any,
) -> list[dict[str, Any]]:
    """Partition the complete presentation grid in the standard child frame."""

    source_basis_values = _integer_basis_tuple(source_basis)
    source_inverse = _fraction_matrix_inverse_3(source_basis or [])
    source_origin_record = _origin_record_from_any(source_origin)
    if (
        source_basis_values is None
        or source_inverse is None
        or source_origin_record is None
        or presentation_basis is None
    ):
        return []
    presentation_points = _presentation_child_points(
        parent_sg=parent_sg,
        parent_xyz=parent_xyz,
        basis=presentation_basis,
        origin=presentation_origin,
    )
    if not presentation_points:
        return []
    decoder = _mode_decoder("Source")
    presentation_matrix = tuple(
        tuple(Fraction(str(value)) for value in row)
        for row in presentation_basis
    )
    presentation_origin_vector = tuple(
        Fraction(str(value)) for value in _origin_vector(presentation_origin)
    )
    source_origin_vector = _origin_record_vector(source_origin_record)

    standard_points: list[tuple[Fraction, Fraction, Fraction]] = []
    for point in presentation_points:
        parent_cinter = _fraction_vecadd(
            _fraction_row_multiply(
                tuple(Fraction(str(value)) for value in point), presentation_matrix
            ),
            presentation_origin_vector,
        )
        parent_pml = decoder.xyz_change_setting_point(
            int(parent_sg), "cinter", "pml", parent_cinter
        )
        child_pml = _fraction_row_multiply(
            _fraction_vecsub(parent_pml, source_origin_vector), source_inverse
        )
        standard_points.append(
            tuple(
                _fraction_mod01(value)
                for value in decoder.xyz_change_setting_point(
                    int(child_sg), "pml", "cinter", child_pml
                )
            )  # type: ignore[arg-type]
        )

    def periodic_match(left: Any, right: Any) -> bool:
        return all(
            abs(
                float(Fraction(left[axis]) - Fraction(right[axis]))
                - round(float(Fraction(left[axis]) - Fraction(right[axis])))
            )
            <= 1e-7
            for axis in range(3)
        )

    try:
        operations = _child_symmetry_ops(int(child_sg))
    except (KeyError, RuntimeError, TypeError, ValueError):
        return []
    unused = set(range(len(standard_points)))
    rows: list[dict[str, Any]] = []
    while unused:
        seed = min(unused)
        component: set[int] = set()
        for operation in operations:
            image = operation.apply_to_xyz(
                [float(value) for value in standard_points[seed]]
            )
            matches = [
                index
                for index, candidate in enumerate(standard_points)
                if periodic_match(image, candidate)
            ]
            if len(matches) != 1:
                return []
            component.add(matches[0])
        if seed not in component or not component <= unused:
            return []
        site_label, _display_xyz = _child_site_display(
            int(child_sg),
            [float(value) for value in standard_points[seed]],
            fallback_site,
        )
        rows.append(
            {
                "label": f"{label_prefix}_{len(rows) + 1}",
                "site": site_label,
                "xyz": list(presentation_points[seed]),
                "_source_index": seed,
            }
        )
        unused.difference_update(component)
    return rows



def _undistorted_rows_from_parent_wyckoff_split(
    *,
    parent_sg: int,
    child_sg: int | None,
    label_prefix: str,
    parent_xyz: tuple[float, float, float] | None,
    parent_site: str,
    basis: list[list[float]] | None,
    origin: Any,
) -> list[dict[str, Any]]:
    if child_sg is None or basis is None or parent_xyz is None:
        return []
    if int(parent_sg) == int(child_sg) and _is_identity_basis_origin(basis, origin):
        return [{"label": f"{label_prefix}_1", "site": parent_site, "xyz": list(parent_xyz)}]
    positions = _presentation_child_points(
        parent_sg=parent_sg,
        parent_xyz=parent_xyz,
        basis=basis,
        origin=origin,
    )
    return _child_orbit_representative_rows(
        label_prefix=label_prefix,
        positions=positions,
        child_sg=child_sg,
        fallback_site=parent_site,
        require_closed=True,
    )



def _undistorted_rows_for_site(
    *,
    sg: int,
    child_sg: int | None,
    site: dict[str, Any],
    label_prefix: str,
    fallback_site: str,
    parent_xyz: tuple[float, float, float] | None,
    split_basis: list[list[float]] | None,
    split_origin: Any,
    presentation_basis: list[list[float]] | None,
    presentation_origin: Any,
    parent_setting_id: int | None = None,
    symmetry_operations: list[str] | None = None,
) -> list[dict[str, Any]]:
    """Source-only undistorted-site rows for one parent atom-site group.

    Formula15 fixes the split topology and site order.  Coordinates are chosen
    from the same parent orbit in the public child setting; this avoids making
    an arbitrary Source-internal Wyckoff representative part of the display
    contract.  The geometric path is accepted only when it has exactly the
    formula15 site multiset, so an incompletely normalized setting cannot
    silently over-expand the structure.
    """

    parent_params = _source_default_wyckoff_params(
        int(sg),
        site,
        parent_setting_id,
        list(symmetry_operations or []),
    )
    formula_rows = undistorted_rows_from_wyckoff_split(
        _assembled_data(),
        parent_sg=sg,
        child_sg=child_sg,
        parent_wyckoff=str(site.get("wyckoff") or ""),
        label_prefix=label_prefix,
        parent_params=parent_params,
        subgroup_basis=split_basis,
        subgroup_origin=split_origin,
    )
    geometric_rows = _undistorted_rows_from_parent_wyckoff_split(
        parent_sg=sg,
        child_sg=child_sg,
        label_prefix=label_prefix,
        parent_xyz=parent_xyz,
        parent_site=fallback_site,
        basis=presentation_basis,
        origin=presentation_origin,
    )
    if formula_rows and parent_xyz is not None and presentation_basis is not None:
        try:
            data = _assembled_data()
            setting_id = (
                int(parent_setting_id)
                if parent_setting_id is not None
                else int(data.default_inter_setting_record(int(sg))["id"])
            )
            presentation_points = _presentation_child_points(
                parent_sg=int(sg),
                parent_xyz=parent_xyz,
                basis=presentation_basis,
                origin=presentation_origin,
            )
            assignments = assign_formula_source_occurrences_to_presentation_grid(
                child_sg=child_sg,
                formula_rows=formula_rows,
                presentation_grid_points=presentation_points,
                parent_pml_to_cinter_matrix=data.pml_to_cinter_matrix(int(sg), setting_id),
                parent_pml_to_cinter_origin=data.cml_to_cinter_origin(int(sg), setting_id),
                presentation_basis=presentation_basis,
                presentation_origin=presentation_origin,
            )
        except (IndexError, KeyError, TypeError, ValueError, ZeroDivisionError):
            assignments = None
        if assignments is not None:
            assigned_rows = [
                {
                    "label": f"{label_prefix}_{index}",
                    "site": assignment.formula_site,
                    "xyz": list(presentation_points[assignment.representative_grid_index]),
                    "_source_index": int(assignment.representative_grid_index),
                    "_wyckoff_formula15": assignment.formula15,
                    "_presentation_orbit_points": [
                        list(presentation_points[point_index])
                        for point_index in assignment.grid_indices
                    ],
                }
                for index, assignment in enumerate(assignments, start=1)
            ]
            geometric_order = (
                match_formula_assignments_to_geometric_grid_orbits(
                    child_sg=int(child_sg),
                    formula_assignments=assignments,
                    presentation_grid_points=presentation_points,
                    geometric_rows=geometric_rows,
                )
                if geometric_rows and child_sg is not None
                else None
            )
            if geometric_order is not None:
                return [
                    {
                        **geometric_rows[geometric_index],
                        "label": f"{label_prefix}_{index}",
                        "_wyckoff_formula15": assignment.formula15,
                        "_presentation_orbit_points": [
                            list(presentation_points[point_index])
                            for point_index in assignment.grid_indices
                        ],
                    }
                    for index, (assignment, geometric_index) in enumerate(
                        zip(assignments, geometric_order),
                        start=1,
                    )
                ]
            # B/M: Formula15 owns the complete split topology and Source row
            # order.  A geometric representative may replace it only through
            # the exact component bijection above; otherwise retain Formula15.
            if sorted(str(row.get("site") or "") for row in assigned_rows) != sorted(
                str(row.get("site") or "") for row in geometric_rows
            ):
                assigned_rows = [
                    {
                        **row,
                        "_source_formula_site": str(row.get("site") or ""),
                    }
                    for row in assigned_rows
                ]
            return assigned_rows
    if formula_rows and geometric_rows:
        unused = list(geometric_rows)
        ordered: list[dict[str, Any]] = []
        for index, formula_row in enumerate(formula_rows, start=1):
            match = next((row for row in unused if row.get("site") == formula_row.get("site")), None)
            if match is None:
                ordered = []
                break
            unused.remove(match)
            ordered.append(
                {
                    **match,
                    "label": f"{label_prefix}_{index}",
                    "_wyckoff_formula15": formula_row.get("_wyckoff_formula15"),
                }
            )
        if ordered and not unused:
            return ordered
        if len(formula_rows) == len(geometric_rows):
            if len(formula_rows) == 1:
                formula_site = str(formula_rows[0].get("site") or "")
                geometric_site = str(geometric_rows[0].get("site") or "")
                if formula_site and formula_site != geometric_site:
                    # M: retain Formula15 topology as an independent witness;
                    # the magnetic correspondence must still rederive it.
                    return [
                        {
                            **geometric_rows[0],
                            "_source_formula_site": formula_site,
                            "_wyckoff_formula15": formula_rows[0].get(
                                "_wyckoff_formula15"
                            ),
                        }
                    ]
            return geometric_rows
    if formula_rows:
        complete_rows = _ordinary_rows_from_complete_child_grid(
            parent_sg=sg,
            child_sg=int(child_sg) if child_sg is not None else 0,
            label_prefix=label_prefix,
            parent_xyz=parent_xyz,
            fallback_site=fallback_site,
            source_basis=split_basis,
            source_origin=split_origin,
            presentation_basis=presentation_basis,
            presentation_origin=presentation_origin,
        )
        if complete_rows and sorted(row["site"] for row in complete_rows) == sorted(
            row["site"] for row in formula_rows
        ):
            return complete_rows
        return formula_rows
    return geometric_rows



def _basis_from_opd_row(selected_opd: dict[str, Any] | None) -> list[list[float]] | None:
    if not isinstance(selected_opd, dict):
        return None
    iso = selected_opd.get("isotropy") or selected_opd
    basis = iso.get("basis") if isinstance(iso, dict) else None
    if not (
        isinstance(basis, list)
        and len(basis) == 3
        and all(isinstance(row, list) and len(row) == 3 for row in basis)
    ):
        return None
    try:
        return [[float(value) for value in row] for row in basis]
    except (TypeError, ValueError):
        return None



def _source_split_basis_from_opd_row(selected_opd: dict[str, Any] | None) -> list[list[float]] | None:
    if not isinstance(selected_opd, dict):
        return None
    iso = selected_opd.get("isotropy") or selected_opd
    basis = iso.get("source_basis_values") if isinstance(iso, dict) else None
    if not isinstance(basis, list) or len(basis) != 9:
        return None
    try:
        return [[float(basis[row * 3 + col]) for col in range(3)] for row in range(3)]
    except (TypeError, ValueError):
        return None



def _source_split_origin_from_opd_row(selected_opd: dict[str, Any] | None) -> tuple[int, int, int, int] | None:
    if not isinstance(selected_opd, dict):
        return None
    iso = selected_opd.get("isotropy") or selected_opd
    origin = iso.get("source_origin_values") if isinstance(iso, dict) else None
    if not isinstance(origin, list) or len(origin) != 4:
        return None
    try:
        return tuple(int(value) for value in origin)  # type: ignore[return-value]
    except (TypeError, ValueError):
        return None



def _point_rotation_matrix(
    data: SourceTables,
    child_sg: int,
    point_op: int,
) -> tuple[tuple[Fraction, Fraction, Fraction], ...]:
    units = (
        (Fraction(1), Fraction(0), Fraction(0)),
        (Fraction(0), Fraction(1), Fraction(0)),
        (Fraction(0), Fraction(0), Fraction(1)),
    )
    return tuple(tuple(data.vrot_fraction(int(child_sg), int(point_op), unit)) for unit in units)  # type: ignore[return-value]



def _fraction_vector_record(
    values: tuple[Fraction, Fraction, Fraction],
    point_op: int,
) -> tuple[int, int, int, int, int]:
    denominator = math.lcm(*(value.denominator for value in values), 1)
    return tuple(int(value * denominator) for value in values) + (denominator, int(point_op))  # type: ignore[return-value]



def _subgroup_parent_operation_records(
    data: Any,
    parent_sg: int,
    child_sg: int,
    basis: tuple[int, ...],
    origin: tuple[int, int, int, int],
) -> tuple[tuple[int, int, int, int, int], ...]:
    """Map child operations into the selected raw parent embedding."""

    basis_matrix = _matrix_from_basis_tuple(basis)
    basis_inverse = _fraction_matrix_inverse_3([[float(value) for value in row] for row in basis_matrix])
    if basis_inverse is None:
        return ()
    origin_vector = _origin_record_vector(origin)
    parent_point_ops: dict[tuple[tuple[Fraction, Fraction, Fraction], ...], int] = {}
    for record in data.generate_space_group_records(int(parent_sg)):
        parent_point_ops[_point_rotation_matrix(data, int(parent_sg), int(record[4]))] = int(record[4])

    out: list[tuple[int, int, int, int, int]] = []
    for x, y, z, den, child_point_op in data.generate_space_group_records(int(child_sg)):
        rotation = _point_rotation_matrix(data, int(child_sg), int(child_point_op))
        translation = (Fraction(int(x), int(den)), Fraction(int(y), int(den)), Fraction(int(z), int(den)))
        parent_rotation = _fraction_matmul(basis_inverse, _fraction_matmul(rotation, basis_matrix))
        point_op = parent_point_ops.get(parent_rotation)
        if point_op is None:
            raise ValueError(
                f"child SG{child_sg} operation {child_point_op} does not map to a point operation in SG{parent_sg}"
            )
        parent_translation = _fraction_vecadd(
            _fraction_vecsub(origin_vector, _fraction_row_multiply(origin_vector, parent_rotation)),
            _fraction_row_multiply(translation, basis_matrix),
        )
        out.append(_fraction_vector_record(parent_translation, point_op))
    return tuple(out)



def _split_basis_origin_for_wyckoff(
    *,
    parent_sg: int,
    child_sg: int | None,
    selected_opd: dict[str, Any] | None,
    presentation_basis: list[list[float]] | None,
    presentation_origin: Any,
) -> tuple[list[list[float]] | None, Any]:
    """Return the transform consumed by the Source-only get_new_wyckoff_ port.

    Complete-mode Wyckoff splitting follows ISO's raw ``data_isotropy``
    subgroup basis/origin, not the public/VALUE CELL presentation basis.  The
    two transforms can describe the same subgroup physically while only the raw
    transform closes the child Wyckoff orbit buffer used by ``get_new_wyckoff_``.
    """

    source_basis = _source_split_basis_from_opd_row(selected_opd)
    source_origin = _source_split_origin_from_opd_row(selected_opd)
    if source_basis is not None and source_origin is not None:
        return source_basis, source_origin
    if child_sg is not None and presentation_basis is not None:
        internal = _internal_split_basis_origin(
            _assembled_data(),
            parent_sg=int(parent_sg),
            child_sg=int(child_sg),
            display_basis=presentation_basis,
            display_origin=presentation_origin,
        )
        if internal is not None:
            return internal
    return presentation_basis, presentation_origin



def _selected_subgroup_number(selected_opd: dict[str, Any] | None) -> int | None:
    if not isinstance(selected_opd, dict):
        return None
    iso = selected_opd.get("isotropy") or selected_opd
    if not isinstance(iso, dict):
        return None
    subgroup = iso.get("subgroup")
    if isinstance(subgroup, dict):
        subgroup = subgroup.get("ordinary_number") or subgroup.get("number")
    try:
        return int(subgroup)
    except (TypeError, ValueError):
        return None



def _selected_magnetic_group_number(
    selected_opd: dict[str, Any] | None,
    child_sg: int | None,
) -> int | None:
    iso = _isotropy_from_opd_row(selected_opd)
    if not iso:
        return None
    subgroup = iso.get("subgroup")
    if not isinstance(subgroup, dict) or subgroup.get("ordinary_number") is None:
        return None
    try:
        group = int(subgroup["number"])
        setting = magnetic_group_setting(group)
    except (IndexError, KeyError, TypeError, ValueError):
        return None
    if child_sg is None or int(setting.ordinary_space_group) != int(child_sg):
        return None
    return group


def _magnetic_groups_in_presentation_order(
    groups: tuple[Any, ...],
    ordinary_rows: list[dict[str, Any]],
    presentation_points: list[list[float]],
    *,
    tolerance: float = 1e-7,
) -> tuple[tuple[Any, tuple[int, ...]], ...] | None:
    """Join magnetic orbit unions to a complete selected-cell grid partition."""

    if not groups or not ordinary_rows or not presentation_points:
        return None
    row_grid_indices: list[tuple[int, ...]] = []
    used_grid_indices: set[int] = set()
    for row in ordinary_rows:
        raw_orbit = row.get("_presentation_orbit_points")
        if not isinstance(raw_orbit, (list, tuple)) or not raw_orbit:
            return None
        orbit_indices: list[int] = []
        for raw_point in raw_orbit:
            if not isinstance(raw_point, (list, tuple)) or len(raw_point) != 3:
                return None
            try:
                matches = [
                    point_index
                    for point_index, point in enumerate(presentation_points)
                    if _frac_close(list(raw_point), point, tol=tolerance)
                ]
            except (OverflowError, TypeError, ValueError):
                return None
            if len(matches) != 1 or matches[0] in used_grid_indices:
                return None
            used_grid_indices.add(matches[0])
            orbit_indices.append(matches[0])
        if len(set(orbit_indices)) != len(orbit_indices):
            return None
        row_grid_indices.append(tuple(sorted(orbit_indices)))
    if used_grid_indices != set(range(len(presentation_points))):
        return None

    grouped: list[tuple[Any, tuple[int, ...]]] = []
    used_rows: set[int] = set()
    used_group_grid_indices: set[int] = set()
    try:
        for group in groups:
            raw_indices = group.ordinary_orbit_indices
            if not isinstance(raw_indices, (list, tuple)) or not raw_indices:
                return None
            ordinary_indices: list[int] = []
            for raw_index in raw_indices:
                if isinstance(raw_index, bool) or not isinstance(raw_index, int):
                    return None
                if not 0 <= raw_index < len(ordinary_rows) or raw_index in used_rows:
                    return None
                used_rows.add(raw_index)
                ordinary_indices.append(raw_index)
            grid_indices = tuple(
                sorted(
                    index
                    for ordinary_index in ordinary_indices
                    for index in row_grid_indices[ordinary_index]
                )
            )
            multiplicity = group.multiplicity
            if (
                not grid_indices
                or len(set(grid_indices)) != len(grid_indices)
                or any(index in used_group_grid_indices for index in grid_indices)
                or isinstance(multiplicity, bool)
                or not isinstance(multiplicity, int)
                or multiplicity != len(grid_indices)
            ):
                return None
            used_group_grid_indices.update(grid_indices)
            grouped.append((group, grid_indices))
    except (AttributeError, IndexError, TypeError, ValueError):
        return None
    if (
        used_rows != set(range(len(ordinary_rows)))
        or used_group_grid_indices != set(range(len(presentation_points)))
    ):
        return None
    return tuple(sorted(grouped, key=lambda item: min(item[1])))


def _rows_with_presentation_grid_witness(
    rows: list[dict[str, Any]],
    presentation_points: list[list[float]],
    expected_partitions: tuple[tuple[int, ...], ...],
    *,
    tolerance: float = 1e-7,
) -> list[dict[str, Any]] | None:
    """Attach grid partitions derived from each row's complete orbit points."""

    if not rows or not presentation_points or len(rows) != len(expected_partitions):
        return None
    witnessed_rows: list[dict[str, Any]] = []
    witnessed_partitions: list[tuple[int, ...]] = []
    used_indices: set[int] = set()
    for row in rows:
        raw_orbit = row.get("_mode_row_orbit_points")
        if not isinstance(raw_orbit, (list, tuple)) or not raw_orbit:
            return None
        row_indices: list[int] = []
        for raw_point in raw_orbit:
            if not isinstance(raw_point, (list, tuple)) or len(raw_point) != 3:
                return None
            try:
                matches = [
                    point_index
                    for point_index, point in enumerate(presentation_points)
                    if _frac_close(list(raw_point), point, tol=tolerance)
                ]
            except (OverflowError, TypeError, ValueError):
                return None
            if len(matches) != 1 or matches[0] in used_indices:
                return None
            used_indices.add(matches[0])
            row_indices.append(matches[0])
        partition = tuple(sorted(row_indices))
        witnessed_partitions.append(partition)
        witnessed_rows.append(
            {
                **row,
                "_presentation_grid_indices": list(partition),
                "_presentation_grid_size": len(presentation_points),
            }
        )
    if (
        used_indices != set(range(len(presentation_points)))
        or sorted(witnessed_partitions) != sorted(expected_partitions)
    ):
        return None
    return witnessed_rows


def _rows_with_optional_presentation_branch_labels(
    rows: list[dict[str, Any]],
    presentation_labels: tuple[str, ...],
) -> list[dict[str, Any]] | None:
    """Transport Wyckoff letters without discarding a proved orbit partition."""

    if not presentation_labels:
        return rows
    if len(presentation_labels) != len(rows):
        return None
    for row, label in zip(rows, presentation_labels, strict=True):
        match = re.match(r"^(\d+)", str(row.get("site") or ""))
        if match is None or not isinstance(label, str) or not label:
            return None
        row["site"] = f"{match.group(1)}{label}"
    return rows



def _magnetic_undistorted_rows_for_site(
    *,
    decoder: ModeDataDecoder,
    magnetic_group: int,
    parent_sg: int,
    child_sg: int,
    label_prefix: str,
    parent_xyz: tuple[float, float, float] | None,
    ordinary_rows: list[dict[str, Any]],
    source_basis: list[list[float]] | None,
    source_origin: Any,
    presentation_basis: list[list[float]] | None,
    presentation_origin: Any,
    magnetic_subgroup_selection: dict[str, Any] | None = None,
    parent_setting_id: int | None = None,
    presentation_grid_points: list[list[float]] | None = None,
) -> list[dict[str, Any]]:
    """Merge ordinary child orbits using the selected magnetic group."""

    source_basis_values = _integer_basis_tuple(source_basis)
    source_origin_record = _origin_record_from_any(source_origin)
    source_inverse = _fraction_matrix_inverse_3(source_basis or [])
    presentation_inverse = _fraction_matrix_inverse_3(presentation_basis or [])
    if (
        not ordinary_rows
        or source_basis_values is None
        or source_origin_record is None
        or source_inverse is None
        or presentation_basis is None
        or presentation_inverse is None
    ):
        return []
    source_matrix = _matrix_from_basis_tuple(source_basis_values)
    source_origin_vector = _origin_record_vector(source_origin_record)
    presentation_matrix = tuple(
        tuple(Fraction(str(value)) for value in row)
        for row in presentation_basis
    )
    presentation_origin_vector = tuple(
        Fraction(str(value)) for value in _origin_vector(presentation_origin)
    )
    try:
        parent_setting_bridge = _parent_setting_bridge(
            int(parent_sg), parent_setting_id
        )
    except (IndexError, KeyError, TypeError, ValueError, ZeroDivisionError):
        return []

    def selected_presentation_points() -> list[list[float]]:
        if presentation_grid_points is not None:
            try:
                points = [
                    [float(point[axis]) for axis in range(3)]
                    for point in presentation_grid_points
                    if len(point) == 3
                ]
            except (IndexError, OverflowError, TypeError, ValueError):
                return []
            if len(points) != len(presentation_grid_points) or not all(
                math.isfinite(value) for point in points for value in point
            ):
                return []
            return points
        return _presentation_child_points(
            parent_sg=int(parent_sg),
            parent_xyz=parent_xyz,
            basis=presentation_basis,
            origin=presentation_origin,
        )

    def to_standard(point: Any) -> tuple[Fraction, Fraction, Fraction]:
        parent_cinter = _fraction_vecadd(
            _fraction_row_multiply(
                tuple(Fraction(str(value)) for value in point),
                presentation_matrix,
            ),
            presentation_origin_vector,
        )
        parent_pml = decoder.xyz_change_setting_point(
            int(parent_sg),
            "cinter",
            "pml",
            _parent_point_to_default(parent_cinter, parent_setting_bridge),
        )
        child_pml = _fraction_row_multiply(
            _fraction_vecsub(parent_pml, source_origin_vector), source_inverse
        )
        return tuple(
            _fraction_mod01(value)
            for value in decoder.xyz_change_setting_point(
                int(child_sg), "pml", "cinter", child_pml
            )
        )  # type: ignore[return-value]

    def to_presentation(point: Any) -> list[float]:
        child_pml = decoder.xyz_change_setting_point(
            int(child_sg), "cinter", "pml", point
        )
        parent_pml = _fraction_vecadd(
            _fraction_row_multiply(child_pml, source_matrix), source_origin_vector
        )
        parent_default_cinter = decoder.xyz_change_setting_point(
            int(parent_sg), "pml", "cinter", parent_pml
        )
        parent_cinter = _parent_point_from_default(
            parent_default_cinter, parent_setting_bridge
        )
        displayed = _fraction_row_multiply(
            _fraction_vecsub(parent_cinter, presentation_origin_vector),
            presentation_inverse,
        )
        return [float(_fraction_mod01(value)) for value in displayed]

    def selected_correspondence_rows() -> list[dict[str, Any]]:
        selection = magnetic_subgroup_selection or {}
        correspondence = selection.get("full_operation_correspondence")
        if not isinstance(correspondence, (list, tuple)) or not correspondence:
            return []
        presentation_points = selected_presentation_points()
        if not presentation_points:
            return []
        standard_points = selected_magnetic_correspondence_standard_points(
            parent_sg=int(parent_sg),
            child_sg=int(child_sg),
            presentation_points=presentation_points,
            presentation_basis=presentation_basis,
            presentation_origin=presentation_origin,
            selected_basis=source_basis,
            selected_origin=source_origin,
            full_operation_correspondence=correspondence,
            parent_setting_id=parent_setting_id,
        )
        if len(standard_points) != len(presentation_points):
            return []

        expected_multiplicity = 0
        for row in ordinary_rows:
            source_formula_site = row.get("_source_formula_site")
            match = re.match(
                r"^(\d+)",
                str(source_formula_site or row.get("site") or ""),
            )
            if match is None:
                return []
            expected_multiplicity += int(match.group(1))
        if expected_multiplicity != len(presentation_points):
            return []

        def periodic_match(left: Any, right: Any) -> bool:
            return all(
                abs(
                    float(Fraction(left[axis]) - Fraction(right[axis]))
                    - round(float(Fraction(left[axis]) - Fraction(right[axis])))
                )
                <= 1e-7
                for axis in range(3)
            )

        unused = set(range(len(standard_points)))
        rows: list[dict[str, Any]] = []
        branch_witnesses: list[tuple[str, Any, int]] = []
        while unused:
            seed = min(unused)
            try:
                identification = identify_magnetic_wyckoff_branch(
                    int(magnetic_group),
                    standard_points[seed],
                    setting="cinter",
                )
                orbit = magnetic_orbit_points(
                    int(magnetic_group),
                    identification.representative,
                    record_setting="cinter",
                )
            except (IndexError, KeyError, TypeError, ValueError, ZeroDivisionError):
                return []
            component: set[int] = set()
            for point in orbit:
                matches = [
                    index
                    for index, candidate in enumerate(standard_points)
                    if periodic_match(point, candidate)
                ]
                if len(matches) != 1:
                    return []
                component.add(matches[0])
            if seed not in component or not component <= unused:
                return []
            representative_matches = [
                index
                for index in component
                if periodic_match(
                    standard_points[index], identification.representative
                )
            ]
            if len(representative_matches) != 1:
                return []
            representative_index = representative_matches[0]
            branch_witnesses.append(
                (
                    str(identification.row.label),
                    presentation_points[representative_index],
                    len(component),
                )
            )
            rows.append(
                {
                    "label": f"{label_prefix}_{len(rows) + 1}",
                    "site": f"{len(component)}{identification.row.label}",
                    "xyz": list(presentation_points[representative_index]),
                    "_mode_row_orbit_points": [
                        list(presentation_points[index])
                        for index in sorted(component)
                    ],
                }
            )
            unused.difference_update(component)
        if len(rows) == len(ordinary_rows) == 1:
            source_formula_site = str(
                ordinary_rows[0].get("_source_formula_site") or ""
            )
            if (
                source_formula_site
                and str(rows[0].get("site") or "") == source_formula_site
            ):
                return rows
        presentation_labels = presentation_branch_labels_from_correspondence(
            magnetic_group=int(magnetic_group),
            branches=branch_witnesses,
        )
        return (
            _rows_with_optional_presentation_branch_labels(
                rows, presentation_labels
            )
            or []
        )

    selected_parent_frame = parent_setting_bridge is not None

    if selected_parent_frame:
        selected_rows = selected_correspondence_rows()
        if selected_rows:
            return selected_rows

    def reference_frame_grid_rows() -> list[dict[str, Any]]:
        """Use a type-IV reference embedding only when it is the display frame."""

        try:
            setting = magnetic_group_setting(int(magnetic_group))
        except (IndexError, KeyError, TypeError, ValueError):
            return []
        reference_basis = tuple(
            tuple(Fraction(int(setting.reference_basis[row * 3 + col])) for col in range(3))
            for row in range(3)
        )
        if (
            int(setting.magnetic_type) != 4
            or int(setting.reference_space_group) != int(parent_sg)
            or presentation_matrix != reference_basis
        ):
            return []
        presentation_points = selected_presentation_points()
        if not presentation_points:
            return []

        def periodic_match(left: Any, right: Any) -> bool:
            for axis in range(3):
                delta = float(Fraction(left[axis]) - Fraction(str(right[axis])))
                delta -= round(delta)
                if abs(delta) > 1e-7:
                    return False
            return True

        unused = set(range(len(presentation_points)))
        rows: list[dict[str, Any]] = []
        while unused:
            seed = min(unused)
            try:
                orbit = magnetic_orbit_points(
                    int(magnetic_group),
                    presentation_points[seed],
                    record_setting="cinter",
                )
                row, _representative = identify_magnetic_wyckoff(
                    int(magnetic_group), presentation_points[seed]
                )
            except (IndexError, KeyError, TypeError, ValueError, ZeroDivisionError):
                return []
            component: set[int] = set()
            for point in orbit:
                matches = [
                    index
                    for index, candidate in enumerate(presentation_points)
                    if periodic_match(point, candidate)
                ]
                if len(matches) != 1:
                    return []
                component.add(matches[0])
            if len(component) != len(orbit) or seed not in component or not component <= unused:
                return []
            rows.append(
                {
                    "label": f"{label_prefix}_{len(rows) + 1}",
                    "site": f"{len(orbit)}{row.label}",
                    "xyz": list(presentation_points[seed]),
                    "_mode_row_orbit_points": [
                        list(presentation_points[index])
                        for index in sorted(component)
                    ],
                }
            )
            unused.difference_update(component)
        return rows

    try:
        groups = group_ordinary_orbits_magnetic(
            int(magnetic_group),
            [to_standard(row["xyz"]) for row in ordinary_rows],
            presentation_tolerance=1e-8,
        )
    except (IndexError, KeyError, TypeError, ValueError, ZeroDivisionError):
        return selected_correspondence_rows() or reference_frame_grid_rows()
    coverage = sorted(index for group in groups for index in group.ordinary_orbit_indices)
    if coverage != list(range(len(ordinary_rows))):
        return selected_correspondence_rows()

    presentation_points = selected_presentation_points()
    if not presentation_points:
        return selected_correspondence_rows()
    ordered_groups = _magnetic_groups_in_presentation_order(
        tuple(groups),
        ordinary_rows,
        presentation_points,
    )
    if ordered_groups is None:
        return selected_correspondence_rows()
    groups = tuple(group for group, _grid_indices in ordered_groups)
    group_grid_indices = tuple(
        grid_indices for _group, grid_indices in ordered_groups
    )

    grouped_rows: list[dict[str, Any]] = []
    for index, (group, grid_indices) in enumerate(
        zip(groups, group_grid_indices, strict=True),
        start=1,
    ):
        row = {
            "label": f"{label_prefix}_{index}",
            "site": f"{int(group.multiplicity)}{group.wyckoff_label}",
            "xyz": to_presentation(group.standard_representative),
            "_presentation_grid_indices": list(grid_indices),
            "_presentation_grid_size": len(presentation_points),
            "_mode_row_orbit_points": [
                list(presentation_points[grid_index])
                for grid_index in grid_indices
            ],
        }
        grouped_rows.append(row)
    presentation_rows = presentation_grid_rows_for_magnetic_groups(
        magnetic_group=int(magnetic_group),
        label_prefix=label_prefix,
        presentation_points=presentation_points,
        groups=groups,
        to_standard=to_standard,
    )
    if presentation_rows:
        witnessed_rows = _rows_with_presentation_grid_witness(
            presentation_rows,
            presentation_points,
            group_grid_indices,
        )
        if witnessed_rows is None:
            return selected_correspondence_rows()
        return witnessed_rows
    if not selected_parent_frame:
        selected_rows = selected_correspondence_rows()
        if selected_rows:
            return selected_rows
    return grouped_rows



def _presentation_basis_candidate(
    basis: list[list[float]] | None,
) -> tuple[list[list[float]] | None, str]:
    if basis is None:
        return None, "none"
    identity = ((1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, 1.0))
    is_identity = all(abs(float(basis[row][col]) - identity[row][col]) <= 1e-12 for row in range(3) for col in range(3))
    return basis, "identity" if is_identity else "child_basis"



def _basis_cinter_to_pml(
    decoder: ModeDataDecoder,
    sg: int,
    basis_cinter: list[list[float]] | None,
    selected_setting_id: int | None = None,
) -> tuple[int, ...] | None:
    if basis_cinter is None:
        return None
    rows_list: list[tuple[Fraction, Fraction, Fraction]] = []
    for row in basis_cinter:
        values: list[Fraction] = []
        for value in row:
            raw = Fraction(str(value))
            snapped = raw.limit_denominator(384)
            values.append(snapped if abs(float(raw - snapped)) <= 5e-9 else raw)
        rows_list.append(tuple(values))
    rows = tuple(rows_list)
    setting_bridge = _parent_setting_bridge(int(sg), selected_setting_id)
    if setting_bridge is not None:
        selected_to_default = fraction_matrix_inverse3(setting_bridge[0])
        rows = fraction_matrix_multiply3(rows, selected_to_default)
    inverse = fraction_matrix_inverse3(decoder.pml_to_cinter_matrix(sg))
    basis_pml = fraction_matrix_multiply3(rows, inverse)
    if any(value.denominator != 1 for row in basis_pml for value in row):
        return None
    return tuple(int(value) for row in basis_pml for value in row)



def _to_child_fractionals(
    positions: list[Any],
    basis: list[list[float]] | None,
    origin: Any,
) -> list[list[float]] | None:
    if basis is None:
        return None
    inverse = _float_matrix_inverse_3(basis)
    if inverse is None:
        return None
    shift = _origin_vector(origin)
    out: list[list[float]] = []
    for position in positions:
        try:
            parent = [float(Fraction(str(value))) for value in position]
        except (TypeError, ValueError):
            return None
        child = _row_multiply([parent[index] - shift[index] for index in range(3)], inverse)
        out.append([_fold01(value) for value in child])
    return out



def _presentation_atom_layout(
    positions: list[Any],
    basis: list[list[float]] | None,
    origin: Any,
    rule: str,
    child_sg: int | None = None,
) -> tuple[list[list[float]], list[int], set[int]] | None:
    child_positions = _to_child_fractionals(positions, basis, origin)
    if child_positions is None:
        return None
    orbit_layout = _child_orbit_display_layout(child_positions, child_sg)
    if orbit_layout is None:
        order = list(range(len(child_positions)))
        representatives = set(order)
    else:
        order, representatives = orbit_layout
    return child_positions, order, representatives



def _presentation_mode_vectors(
    modes: list[Any],
    basis: list[list[float]] | None,
) -> list[list[list[float]]]:
    inverse = _float_matrix_inverse_3(basis) if basis is not None else None
    return [
        [
            (
                [float(value) for value in vector[:3]]
                if inverse is None
                else _row_multiply([float(value) for value in vector[:3]], inverse)
            )
            for vector in mode
        ]
        for mode in modes
        if isinstance(mode, list)
    ]
