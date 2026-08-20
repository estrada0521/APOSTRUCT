"""Compare displacive, magnetic, and strain mode vectors."""

from __future__ import annotations

import math
import re
from collections import OrderedDict
from typing import Any, Sequence

import numpy as np

from Verification.comparison.basis_mode import (
    BasisEquivalentModeFrameTransport,
    basis_equivalent_mode_frame_transport,
)
from Verification.comparison.structure import (
    _exact_integer_value,
    _frac_vec_close,
    _validated_basis_equivalent_atom_result,
    _validated_basis_equivalent_lattice_result,
    _vec_close,
)
from Verification.result import physical_k_label
from Verification.parsers.complete_mode_text import (
    DisplaciveModeDefinition,
    StrainModeDefinition,
)

_VISIBLE_VECTOR_DECIMALS = 4
_VISIBLE_VECTOR_HALF_QUANTUM = math.nextafter(
    0.5 * 10.0**-_VISIBLE_VECTOR_DECIMALS,
    math.inf,
)
_MODE_MATCH_CACHE_MAX_ENTRIES = 65_536


def _visible_vector_tolerance(tol: float) -> float:
    """Respect the precision of vectors printed on the Web surface."""

    return max(tol, _VISIBLE_VECTOR_HALF_QUANTUM)


def _row_matching(
    web_rows: list[Any],
    local_rows: list[dict[str, Any]],
    predicate,
    candidate_indices: Sequence[Sequence[int]] | None = None,
) -> list[int] | None:
    """Return a complete bipartite assignment from Web rows to local rows."""

    if len(web_rows) != len(local_rows):
        return None
    if candidate_indices is not None and len(candidate_indices) != len(web_rows):
        return None
    local_to_web = [-1] * len(local_rows)

    def assign(web_index: int, seen: set[int]) -> bool:
        candidates = (
            range(len(local_rows))
            if candidate_indices is None
            else candidate_indices[web_index]
        )
        for local_index in candidates:
            local_row = local_rows[local_index]
            if local_index in seen or not predicate(web_rows[web_index], local_row):
                continue
            seen.add(local_index)
            previous = local_to_web[local_index]
            if previous < 0 or assign(previous, seen):
                local_to_web[local_index] = web_index
                return True
        return False

    if not all(assign(web_index, set()) for web_index in range(len(web_rows))):
        return None
    web_to_local = [-1] * len(web_rows)
    for local_index, web_index in enumerate(local_to_web):
        web_to_local[web_index] = local_index
    return web_to_local


def _periodic_xyz_candidate_indices(
    web_rows: list[Any],
    local_rows: list[dict[str, Any]],
    tol: float,
) -> list[list[int]] | None:
    """Index tol-close periodic positions without changing assignment order."""

    if not math.isfinite(tol) or tol <= 0.0:
        return None
    bucket_width = 2.0 * tol
    if bucket_width <= 0.0:
        return None
    inverse_width = 1.0 / bucket_width
    if not math.isfinite(inverse_width):
        return None
    bucket_count = max(1, int(math.floor(inverse_width)))
    buckets: dict[tuple[int, int, int], list[int]] = {}
    for local_index, row in enumerate(local_rows):
        xyz = row.get("xyz") or []
        if len(xyz) != 3:
            continue
        try:
            folded = tuple(float(value) % 1.0 for value in xyz)
        except (TypeError, ValueError):
            continue
        if not all(math.isfinite(value) for value in folded):
            continue
        key = tuple(
            min(bucket_count - 1, int(math.floor(value / bucket_width)))
            for value in folded
        )
        buckets.setdefault(key, []).append(local_index)  # type: ignore[arg-type]

    result: list[list[int]] = []
    for web_row in web_rows:
        xyz = getattr(web_row, "xyz", ())
        if len(xyz) != 3:
            return None
        try:
            folded = tuple(float(value) % 1.0 for value in xyz)
        except (TypeError, ValueError):
            return None
        if not all(math.isfinite(value) for value in folded):
            return None
        center = tuple(
            min(bucket_count - 1, int(math.floor(value / bucket_width)))
            for value in folded
        )
        candidates: set[int] = set()
        for first in (-1, 0, 1):
            for second in (-1, 0, 1):
                for third in (-1, 0, 1):
                    candidates.update(
                        buckets.get(
                            (
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
                if _frac_vec_close(xyz, local_rows[local_index].get("xyz") or [], tol)
            ]
        )
    return result


def _folded_position_key(
    rows: Sequence[Any],
    *,
    web: bool,
) -> tuple[tuple[float, float, float], ...] | None:
    positions: list[tuple[float, float, float]] = []
    for row in rows:
        try:
            xyz = getattr(row, "xyz", ()) if web else row.get("xyz") or []
            if len(xyz) != 3:
                return None
            position = tuple(float(value) % 1.0 for value in xyz)
        except (AttributeError, TypeError, ValueError):
            return None
        if len(position) != 3 or not all(math.isfinite(value) for value in position):
            return None
        positions.append(position)  # type: ignore[arg-type]
    return tuple(positions)


class _PeriodicCandidateCache:
    """Reuse immutable position candidates within one vector comparison."""

    def __init__(
        self,
        tol: float,
        *,
        max_mode_match_entries: int = _MODE_MATCH_CACHE_MAX_ENTRIES,
    ) -> None:
        self._tol = tol
        self._max_mode_match_entries = max(0, max_mode_match_entries)
        self._position_keys: dict[
            tuple[bool, int],
            tuple[Sequence[Any], tuple[tuple[float, float, float], ...] | None],
        ] = {}
        self._interned_position_keys: dict[
            tuple[tuple[float, float, float], ...],
            tuple[tuple[float, float, float], ...],
        ] = {}
        self._values: dict[
            tuple[
                tuple[tuple[float, float, float], ...],
                tuple[tuple[float, float, float], ...],
            ],
            tuple[tuple[int, ...], ...] | None,
        ] = {}
        self._mode_match_facts: OrderedDict[
            tuple[int, int],
            tuple[
                DisplaciveModeDefinition,
                dict[str, Any],
                tuple[bool, bool, bool],
            ],
        ] = OrderedDict()

    def _position_key(
        self,
        rows: Sequence[Any],
        *,
        web: bool,
    ) -> tuple[tuple[float, float, float], ...] | None:
        identity = (web, id(rows))
        cached = self._position_keys.get(identity)
        if cached is not None and cached[0] is rows:
            return cached[1]
        key = _folded_position_key(rows, web=web)
        if key is not None:
            key = self._interned_position_keys.setdefault(key, key)
        self._position_keys[identity] = (rows, key)
        return key

    def get(
        self,
        web_rows: Sequence[Any],
        local_rows: Sequence[dict[str, Any]],
    ) -> tuple[tuple[int, ...], ...] | None:
        web_key = self._position_key(web_rows, web=True)
        local_key = self._position_key(local_rows, web=False)
        if web_key is None or local_key is None:
            result = _periodic_xyz_candidate_indices(
                list(web_rows), list(local_rows), self._tol
            )
            return None if result is None else tuple(tuple(row) for row in result)
        key = (web_key, local_key)
        if key not in self._values:
            result = _periodic_xyz_candidate_indices(
                list(web_rows), list(local_rows), self._tol
            )
            self._values[key] = (
                None if result is None else tuple(tuple(row) for row in result)
            )
        return self._values[key]

    def get_mode_match_facts(
        self,
        web_mode: DisplaciveModeDefinition,
        local_mode: dict[str, Any],
    ) -> tuple[bool, bool, bool] | None:
        key = (id(web_mode), id(local_mode))
        cached = self._mode_match_facts.get(key)
        if cached is None or cached[0] is not web_mode or cached[1] is not local_mode:
            return None
        self._mode_match_facts.move_to_end(key)
        return cached[2]

    def put_mode_match_facts(
        self,
        web_mode: DisplaciveModeDefinition,
        local_mode: dict[str, Any],
        value: tuple[bool, bool, bool],
    ) -> None:
        if self._max_mode_match_entries == 0:
            return
        key = (id(web_mode), id(local_mode))
        self._mode_match_facts[key] = (
            web_mode,
            local_mode,
            value,
        )
        self._mode_match_facts.move_to_end(key)
        while len(self._mode_match_facts) > self._max_mode_match_entries:
            self._mode_match_facts.popitem(last=False)


def _optional_float_close(left: float | None, right: object, tol: float) -> bool:
    if left is None or right is None:
        return left is None and right is None
    try:
        return math.isclose(float(left), float(right), rel_tol=0.0, abs_tol=tol)
    except (TypeError, ValueError):
        return False


def _forward_filled_atom_labels(rows: list[Any], *, web: bool) -> list[str | None]:
    labels: list[str | None] = []
    raw_labels = [
        getattr(row, "atom_label", None) if web else row.get("atom")
        for row in rows
    ]
    current = next(
        (value for value in raw_labels if isinstance(value, str) and value),
        None,
    )
    for raw in raw_labels:
        if isinstance(raw, str) and raw:
            current = raw
        labels.append(current)
    return labels


def _mode_row_assignments(
    web_mode: DisplaciveModeDefinition,
    local_mode: dict[str, Any],
    tol: float,
    candidate_cache: _PeriodicCandidateCache | None = None,
) -> tuple[
    list[Any],
    list[dict[str, Any]],
    list[int] | None,
    list[int] | None,
    list[int] | None,
]:
    """Return the three physical row assignments without materializing rows."""

    web_source_rows = web_mode.rows
    local_source_rows = local_mode.get("rows") or []
    web_rows = list(web_source_rows)
    local_rows = list(local_source_rows)
    vector_tol = _visible_vector_tolerance(tol)

    def xyz(web_row: Any, local_row: dict[str, Any]) -> bool:
        return _frac_vec_close(web_row.xyz, local_row.get("xyz") or [], tol)

    position_candidates = (
        _periodic_xyz_candidate_indices(web_rows, local_rows, tol)
        if candidate_cache is None
        else candidate_cache.get(web_source_rows, local_source_rows)
    )
    position_is_indexed = position_candidates is not None

    def candidate_xyz(web_row: Any, local_row: dict[str, Any]) -> bool:
        return position_is_indexed or xyz(web_row, local_row)

    exact = _row_matching(
        web_rows,
        local_rows,
        lambda web_row, local_row: candidate_xyz(web_row, local_row)
        and _vec_close(web_row.dxyz, local_row.get("dxyz") or [], vector_tol),
        position_candidates,
    )
    sign_flip = _row_matching(
        web_rows,
        local_rows,
        lambda web_row, local_row: candidate_xyz(web_row, local_row)
        and _vec_close(
            web_row.dxyz,
            [-(float(x)) for x in (local_row.get("dxyz") or [])],
            vector_tol,
        ),
        position_candidates,
    )
    positions = _row_matching(
        web_rows, local_rows, candidate_xyz, position_candidates
    )
    if candidate_cache is not None:
        candidate_cache.put_mode_match_facts(
            web_mode,
            local_mode,
            (exact is not None, sign_flip is not None, positions is not None),
        )
    return web_rows, local_rows, exact, sign_flip, positions


def _mode_row_match_facts(
    web_mode: DisplaciveModeDefinition,
    local_mode: dict[str, Any],
    tol: float,
    candidate_cache: _PeriodicCandidateCache | None = None,
) -> tuple[bool, bool, bool]:
    if candidate_cache is not None:
        cached = candidate_cache.get_mode_match_facts(web_mode, local_mode)
        if cached is not None:
            return cached
    _, _, exact, sign_flip, positions = _mode_row_assignments(
        web_mode,
        local_mode,
        tol,
        candidate_cache,
    )
    return exact is not None, sign_flip is not None, positions is not None


def _mode_row_payload(
    web_mode: DisplaciveModeDefinition,
    local_mode: dict[str, Any],
    tol: float,
    candidate_cache: _PeriodicCandidateCache | None = None,
) -> dict[str, Any]:
    web_rows, local_rows, exact, sign_flip, positions = _mode_row_assignments(
        web_mode,
        local_mode,
        tol,
        candidate_cache,
    )
    summary = _mode_row_summary(
        web_mode,
        local_mode,
        web_rows=web_rows,
        local_rows=local_rows,
        exact=exact,
        sign_flip=sign_flip,
        positions=positions,
        tol=tol,
    )
    web_atoms = _forward_filled_atom_labels(web_rows, web=True)
    local_atoms = _forward_filled_atom_labels(local_rows, web=False)
    vector_tol = _visible_vector_tolerance(tol)
    best_assignment = exact or sign_flip or positions
    rows = []
    for web_index, web_row in enumerate(web_rows):
        local_index = (
            best_assignment[web_index]
            if best_assignment is not None
            else web_index
        )
        local_row = local_rows[local_index] if local_index < len(local_rows) else {}
        xyz_match = _frac_vec_close(web_row.xyz, local_row.get("xyz") or [], tol)
        dxyz_match = _vec_close(
            web_row.dxyz,
            local_row.get("dxyz") or [],
            vector_tol,
        )
        dxyz_sign_match = _vec_close(
            web_row.dxyz,
            [-(float(x)) for x in (local_row.get("dxyz") or [])],
            vector_tol,
        )
        rows.append(
            {
                "status": "ok" if xyz_match and dxyz_match else ("sign" if xyz_match and dxyz_sign_match else "fail"),
                "local_row_index": local_index,
                "web": {
                    "atom": web_atoms[web_index],
                    "atom_raw": web_row.atom_label,
                    "xyz": web_row.xyz,
                    "dxyz": web_row.dxyz,
                },
                "local": {
                    "atom": local_atoms[local_index] if local_index < len(local_atoms) else None,
                    "atom_raw": local_row.get("atom"),
                    "xyz": local_row.get("xyz"),
                    "dxyz": local_row.get("dxyz"),
                },
            }
        )
    return {
        **summary,
        "rows": rows,
    }


def _mode_row_summary(
    web_mode: DisplaciveModeDefinition,
    local_mode: dict[str, Any],
    *,
    web_rows: list[Any],
    local_rows: list[dict[str, Any]],
    exact: list[int] | None,
    sign_flip: list[int] | None,
    positions: list[int] | None,
    tol: float,
) -> dict[str, Any]:
    """Return every per-mode field except the discarded row diagnostics."""

    count_match = len(web_rows) == len(local_rows)
    if not count_match:
        status = "count"
    elif exact is not None:
        status = "ok"
    elif sign_flip is not None:
        status = "sign"
    elif positions is not None:
        status = "vector"
    else:
        status = "position"
    label_match = web_mode.label == local_mode.get("label")
    normfactor_match = _optional_float_close(
        web_mode.normfactor,
        local_mode.get("normfactor"),
        tol,
    )
    if status != "ok":
        strict_status = status
    elif not label_match:
        strict_status = "label"
    elif not normfactor_match:
        strict_status = "normfactor"
    else:
        strict_status = "ok"
    return {
        "status": status,
        "strict_status": strict_status,
        "count_match": count_match,
        "web_label": web_mode.label,
        "local_label": local_mode.get("label"),
        "label_match": label_match,
        "web_normfactor": web_mode.normfactor,
        "local_normfactor": local_mode.get("normfactor"),
        "normfactor_match": normfactor_match,
        "row_count": {"web": len(web_rows), "local": len(local_rows)},
    }


def _mode_row_probe(
    web_mode: DisplaciveModeDefinition,
    local_mode: dict[str, Any],
    tol: float,
    candidate_cache: _PeriodicCandidateCache | None = None,
) -> dict[str, Any]:
    web_rows = list(web_mode.rows)
    local_rows = list(local_mode.get("rows") or [])
    exact, sign_flip, positions = _mode_row_match_facts(
        web_mode,
        local_mode,
        tol,
        candidate_cache,
    )
    return _mode_row_summary(
        web_mode,
        local_mode,
        web_rows=web_rows,
        local_rows=local_rows,
        exact=[] if exact else None,
        sign_flip=[] if sign_flip else None,
        positions=[] if positions else None,
        tol=tol,
    )


def _mode_physical_block_key(
    label: str, *, expected_kind: str | None = None
) -> tuple[str, str] | None:
    identity = _mode_physical_block_identity(label, expected_kind=expected_kind)
    return None if identity is None else identity[0]


def _mode_physical_block_identity(
    label: str, *, expected_kind: str | None = None
) -> tuple[tuple[str, str], tuple[str, str], str] | None:
    """Return exact and direction-neutral identities for one mode block."""

    match = re.fullmatch(
        r"(?P<formal>.+)(?P<direction>\([^()]*\))"
        r"(?P<site>\[[^\[\]]+:(?P<kind>dsp|mag)\])(?P<tail>.+)",
        str(label),
    )
    if match is None:
        return None
    if expected_kind is not None and match.group("kind") != expected_kind:
        return None
    tail = re.fullmatch(r"(?P<base>.*?)(?:_\d+)?\([^()]*\)", match.group("tail"))
    if tail is None:
        return None
    physical_formal = physical_k_label(match.group("formal"))
    exact_prefix = _canonical_mode_block_decimals(
        physical_formal + match.group("direction") + match.group("site")
    )
    physical_prefix = _canonical_mode_block_decimals(
        physical_formal + match.group("site")
    )
    return (
        (exact_prefix, tail.group("base")),
        (physical_prefix, tail.group("base")),
        match.group("direction"),
    )


_MODE_BLOCK_DECIMAL = re.compile(
    r"(?<![A-Za-z0-9_.])(?P<whole>\d*)\.(?P<fraction>\d+)(?![\d.])"
)


def _canonical_mode_block_decimals(value: str) -> str:
    """Ignore only redundant decimal zeroes in visible formal labels."""

    def canonical(match: re.Match[str]) -> str:
        whole = (match.group("whole") or "0").lstrip("0") or "0"
        fraction = match.group("fraction").rstrip("0")
        return whole if not fraction else f"{whole}.{fraction}"

    return _MODE_BLOCK_DECIMAL.sub(canonical, value)


def _signed_permutation_modes_match(
    web_modes: list[DisplaciveModeDefinition] | tuple[DisplaciveModeDefinition, ...],
    local_modes: list[dict[str, Any]],
    tol: float,
) -> bool:
    """Match only exact per-column signs inside one semantic mode block."""

    return _signed_permutation_mode_assignment(web_modes, local_modes, tol) is not None


def _signed_permutation_mode_assignment(
    web_modes: list[DisplaciveModeDefinition] | tuple[DisplaciveModeDefinition, ...],
    local_modes: list[dict[str, Any]],
    tol: float,
    *,
    require_unambiguous_positions: bool = False,
    candidate_cache: _PeriodicCandidateCache | None = None,
) -> tuple[int, ...] | None:
    """Return a complete within-block Web-to-local column assignment."""

    if len(web_modes) != len(local_modes) or not web_modes:
        return None
    web_groups: dict[tuple[str, str], list[tuple[int, DisplaciveModeDefinition]]] = {}
    local_groups: dict[tuple[str, str], list[tuple[int, dict[str, Any]]]] = {}
    for index, mode in enumerate(web_modes):
        key = _mode_physical_block_key(mode.label)
        if key is None:
            return None
        web_groups.setdefault(key, []).append((index, mode))
    for index, mode in enumerate(local_modes):
        key = _mode_physical_block_key(str(mode.get("label") or ""))
        if key is None:
            return None
        local_groups.setdefault(key, []).append((index, mode))
    if set(web_groups) != set(local_groups):
        return None
    assignment = [-1] * len(web_modes)
    for key, web_group in web_groups.items():
        local_group = local_groups[key]
        if len(web_group) != len(local_group):
            return None
        edges = []
        for _, web_mode in web_group:
            candidates = []
            for local_slot, (_, local_mode) in enumerate(local_group):
                exact, sign_flip, _ = _mode_row_match_facts(
                    web_mode,
                    local_mode,
                    tol,
                    candidate_cache,
                )
                if not exact and not sign_flip:
                    continue
                if require_unambiguous_positions and not _mode_position_vectors_are_unambiguous(
                    web_mode,
                    local_mode,
                    tol,
                    candidate_cache,
                ):
                    continue
                candidates.append(local_slot)
            edges.append(candidates)
        matched_web_by_local = [-1] * len(local_group)

        def assign(web_index: int, seen: set[int]) -> bool:
            for local_slot in edges[web_index]:
                if local_slot in seen:
                    continue
                seen.add(local_slot)
                prior = matched_web_by_local[local_slot]
                if prior < 0 or assign(prior, seen):
                    matched_web_by_local[local_slot] = web_index
                    return True
            return False

        if not all(assign(web_index, set()) for web_index in range(len(web_group))):
            return None
        for local_slot, web_slot in enumerate(matched_web_by_local):
            web_index = web_group[web_slot][0]
            assignment[web_index] = local_group[local_slot][0]
    if any(index < 0 for index in assignment):
        return None
    return tuple(assignment)


def _rows_in_reference_order(
    reference_rows: list[Any],
    rows: list[dict[str, Any]],
    tol: float,
    candidate_cache: _PeriodicCandidateCache | None = None,
) -> list[dict[str, Any]] | None:
    candidates = (
        _periodic_xyz_candidate_indices(reference_rows, rows, tol)
        if candidate_cache is None
        else candidate_cache.get(reference_rows, rows)
    )
    if (
        candidates is None
        or len(candidates) != len(reference_rows)
        or any(len(indices) != 1 for indices in candidates)
    ):
        return None
    assignment = [indices[0] for indices in candidates]
    if len(set(assignment)) != len(rows):
        return None
    return [rows[index] for index in assignment]


def _componentwise_minimax_quotient(
    local: np.ndarray,
    web: np.ndarray,
    initial: np.ndarray,
    bound: float,
) -> np.ndarray | None:
    """Find a basis quotient inside every visible component's rounding box."""

    if (
        local.ndim != 2
        or web.ndim != 2
        or local.shape != web.shape
        or initial.shape != (local.shape[1], web.shape[1])
        or local.shape[0] == 0
        or local.shape[1] == 0
        or not np.isfinite(local).all()
        or not np.isfinite(web).all()
        or not np.isfinite(initial).all()
        or not math.isfinite(bound)
        or bound <= 0.0
    ):
        return None
    initial_residual = local @ initial - web
    if not np.isfinite(initial_residual).all():
        return None
    unresolved = tuple(
        column
        for column in range(web.shape[1])
        if float(np.max(np.abs(initial_residual[:, column]), initial=0.0)) > bound
    )
    if not unresolved:
        return initial.copy()
    try:
        from scipy.optimize import linprog
    except ImportError:
        return None
    row_count, dimension = local.shape
    quotient = initial.copy()
    ones = np.ones(row_count, dtype=float)
    positive = np.column_stack((local, -ones))
    negative = np.column_stack((-local, -ones))
    constraints = np.vstack((positive, negative))
    objective = np.concatenate((np.zeros(dimension, dtype=float), [1.0]))
    bounds = [(None, None)] * dimension + [(0.0, None)]
    for column in unresolved:
        target = web[:, column]
        limits = np.concatenate((target, -target))
        try:
            solved = linprog(
                objective,
                A_ub=constraints,
                b_ub=limits,
                bounds=bounds,
                method="highs",
            )
            success = solved.success
            if not isinstance(success, (bool, np.bool_)) or not bool(success):
                return None
            solution = np.asarray(solved.x, dtype=float)
        except (ArithmeticError, AttributeError, TypeError, ValueError):
            return None
        if (
            solution.shape != (dimension + 1,)
            or not np.isfinite(solution).all()
            or solution[-1] < 0.0
        ):
            return None
        quotient[:, column] = solution[:dimension]
    return quotient if np.isfinite(quotient).all() else None


def _family_span_proof(
    web_modes: list[DisplaciveModeDefinition]
    | tuple[DisplaciveModeDefinition, ...],
    local_modes: list[dict[str, Any]],
    tol: float,
    *,
    vector_kind: str,
    candidate_cache: _PeriodicCandidateCache | None = None,
    allow_parent_direction_gauge: bool = False,
    allow_site_irrep_basis: bool = False,
) -> dict[str, Any] | None:
    """Prove equal visible mode spans inside each selected physical block."""

    if (
        vector_kind not in {"dsp", "mag"}
        or len(web_modes) != len(local_modes)
        or not web_modes
        or not math.isfinite(tol)
        or tol <= 0.0
    ):
        return None
    web_groups: dict[tuple[str, ...], list[DisplaciveModeDefinition]] = {}
    local_groups: dict[tuple[str, ...], list[dict[str, Any]]] = {}
    web_directions: dict[tuple[str, ...], set[str]] = {}
    local_directions: dict[tuple[str, ...], set[str]] = {}
    for mode in web_modes:
        identity = _mode_physical_block_identity(
            mode.label,
            expected_kind=vector_kind,
        )
        if identity is None:
            return None
        exact_key, physical_key, direction = identity
        selected_key = physical_key if allow_parent_direction_gauge else exact_key
        key = selected_key[:1] if allow_site_irrep_basis else selected_key
        web_groups.setdefault(key, []).append(mode)
        web_directions.setdefault(key, set()).add(direction)
    for mode in local_modes:
        identity = _mode_physical_block_identity(
            str(mode.get("label") or ""),
            expected_kind=vector_kind,
        )
        if identity is None:
            return None
        exact_key, physical_key, direction = identity
        selected_key = physical_key if allow_parent_direction_gauge else exact_key
        key = selected_key[:1] if allow_site_irrep_basis else selected_key
        local_groups.setdefault(key, []).append(mode)
        local_directions.setdefault(key, set()).add(direction)
    if set(web_groups) != set(local_groups):
        return None
    direction_gauge_keys: set[tuple[str, ...]] = set()
    if allow_parent_direction_gauge:
        for key in web_groups:
            if len(web_directions[key]) != 1 or len(local_directions[key]) != 1:
                return None
            if web_directions[key] != local_directions[key]:
                direction_gauge_keys.add(key)
        if not direction_gauge_keys:
            return None

    block_count = 0
    signed_permutation_block_count = 0
    family_basis_block_count = 0
    componentwise_minimax_block_count = 0
    max_condition = 0.0
    max_residual = 0.0
    max_rounding_ratio = 0.0
    for key, web_group in web_groups.items():
        local_group = local_groups[key]
        dimension = len(web_group)
        if dimension == 0 or len(local_group) != dimension:
            return None
        if _signed_permutation_mode_assignment(
            web_group,
            local_group,
            tol,
            candidate_cache=candidate_cache,
        ) is not None:
            block_count += 1
            signed_permutation_block_count += 1
            continue
        reference_rows = list(web_group[0].rows)
        if not reference_rows:
            return None
        web_columns: list[list[float]] = []
        local_columns: list[list[float]] = []
        for web_mode in web_group:
            rows = [
                {"xyz": row.xyz, "dxyz": row.dxyz}
                for row in web_mode.rows
            ]
            ordered = _rows_in_reference_order(
                reference_rows, rows, tol, candidate_cache
            )
            if ordered is None:
                return None
            try:
                vectors = [tuple(float(value) for value in row["dxyz"]) for row in ordered]
            except (KeyError, TypeError, ValueError):
                return None
            if any(len(vector) != 3 for vector in vectors):
                return None
            web_columns.append([value for vector in vectors for value in vector])
        for local_mode in local_group:
            rows = list(local_mode.get("rows") or [])
            if any(not isinstance(row, dict) for row in rows):
                return None
            ordered = _rows_in_reference_order(
                reference_rows, rows, tol, candidate_cache
            )
            if ordered is None:
                return None
            try:
                vectors = [
                    tuple(float(value) for value in (row.get("dxyz") or []))
                    for row in ordered
                ]
            except (TypeError, ValueError):
                return None
            if any(len(vector) != 3 for vector in vectors):
                return None
            local_columns.append([value for vector in vectors for value in vector])
        try:
            web = np.asarray(web_columns, dtype=float).T
            local = np.asarray(local_columns, dtype=float).T
        except (TypeError, ValueError):
            return None
        if (
            web.shape != local.shape
            or web.shape[1] != dimension
            or web.shape[0] != 3 * len(reference_rows)
            or not np.isfinite(web).all()
            or not np.isfinite(local).all()
        ):
            return None
        try:
            quotient, _, local_rank, singular_values = np.linalg.lstsq(
                local, web, rcond=None
            )
        except np.linalg.LinAlgError:
            return None
        if (
            int(local_rank) != dimension
            or quotient.shape != (dimension, dimension)
            or not np.isfinite(quotient).all()
            or len(singular_values) != dimension
            or singular_values[-1] <= 0.0
        ):
            return None
        local_condition = float(singular_values[0] / singular_values[-1])
        if not math.isfinite(local_condition) or local_condition > 1.0 / tol:
            return None
        residual_matrix = local @ quotient - web
        if not np.isfinite(residual_matrix).all():
            return None
        residual = float(np.max(np.abs(residual_matrix), initial=0.0))
        visible_rounding_bound = _visible_vector_tolerance(tol)
        if residual > visible_rounding_bound:
            maximum_box_l2 = math.nextafter(
                math.sqrt(local.shape[0]) * visible_rounding_bound,
                math.inf,
            )
            if np.any(
                np.linalg.norm(residual_matrix, axis=0) > maximum_box_l2
            ):
                return None
            minimax_quotient = _componentwise_minimax_quotient(
                local,
                web,
                quotient,
                visible_rounding_bound,
            )
            if minimax_quotient is None:
                return None
            quotient = minimax_quotient
            residual_matrix = local @ quotient - web
            if not np.isfinite(residual_matrix).all():
                return None
            residual = float(np.max(np.abs(residual_matrix), initial=0.0))
            if residual > visible_rounding_bound:
                return None
            componentwise_minimax_block_count += 1
        try:
            quotient_condition = float(np.linalg.cond(quotient))
        except np.linalg.LinAlgError:
            return None
        condition = max(quotient_condition, local_condition)
        if not math.isfinite(condition) or condition > 1.0 / tol:
            return None
        block_count += 1
        family_basis_block_count += 1
        max_condition = max(max_condition, condition)
        max_residual = max(max_residual, residual)
        max_rounding_ratio = max(
            max_rounding_ratio,
            residual / visible_rounding_bound,
        )
    return {
        "block_scope": (
            "parent_irrep_direction_child_site"
            if allow_site_irrep_basis
            else "formal_family"
        ),
        "block_count": block_count,
        "signed_permutation_block_count": signed_permutation_block_count,
        "family_basis_block_count": family_basis_block_count,
        "componentwise_minimax_block_count": componentwise_minimax_block_count,
        "parent_direction_gauge_block_count": len(direction_gauge_keys),
        "max_condition": max_condition,
        "max_residual": max_residual,
        "visible_component_rounding_bound": _visible_vector_tolerance(tol),
        "max_rounding_ratio": max_rounding_ratio,
    }


def _site_irrep_basis_scope_is_broader(
    modes: Sequence[DisplaciveModeDefinition] | Sequence[dict[str, Any]],
    *,
    vector_kind: str,
) -> bool:
    """Return whether one child-site space spans distinct print families.

    The tail labels name a site-symmetry-adapted presentation basis.  They do
    not split the parent-IR distortion space carried by that child site.
    """

    grouped: dict[str, set[tuple[str, str]]] = {}
    for mode in modes:
        label = (
            mode.label
            if isinstance(mode, DisplaciveModeDefinition)
            else mode.get("label")
        )
        identity = _mode_physical_block_identity(
            str(label or ""), expected_kind=vector_kind
        )
        if identity is None:
            return False
        exact_key = identity[0]
        grouped.setdefault(exact_key[0], set()).add(exact_key)
    return any(len(keys) > 1 for keys in grouped.values())


def _vector_modes_payload(
    web_modes: list[DisplaciveModeDefinition] | tuple[DisplaciveModeDefinition, ...],
    local_modes: list[dict[str, Any]],
    tol: float,
    *,
    vector_kind: str,
) -> dict[str, Any]:
    payload, _ = _vector_modes_payload_with_assignment(
        web_modes,
        local_modes,
        tol,
        vector_kind=vector_kind,
    )
    return payload


def _vector_modes_payload_with_assignment(
    web_modes: list[DisplaciveModeDefinition] | tuple[DisplaciveModeDefinition, ...],
    local_modes: list[dict[str, Any]],
    tol: float,
    *,
    vector_kind: str,
    require_unambiguous_signed_permutation: bool = False,
    materialize_rows: bool = True,
) -> tuple[dict[str, Any], tuple[int, ...] | None]:
    count_match = len(web_modes) == len(local_modes)
    candidate_cache = _PeriodicCandidateCache(tol)
    mode_builder = _mode_row_payload if materialize_rows else _mode_row_probe
    modes = [
        mode_builder(web_mode, local_mode, tol, candidate_cache)
        for web_mode, local_mode in zip(web_modes, local_modes, strict=False)
    ]
    hard_statuses = {mode["status"] for mode in modes}
    strict_statuses = {mode["strict_status"] for mode in modes}
    signed_permutation_assignment = (
        _signed_permutation_mode_assignment(
            web_modes,
            local_modes,
            tol,
            require_unambiguous_positions=require_unambiguous_signed_permutation,
            candidate_cache=candidate_cache,
        )
        if count_match
        else None
    )
    signed_permutation = signed_permutation_assignment is not None
    family_span = None
    if count_match and signed_permutation_assignment is None:
        site_scope = (
            _site_irrep_basis_scope_is_broader(web_modes, vector_kind=vector_kind)
            and _site_irrep_basis_scope_is_broader(
                local_modes, vector_kind=vector_kind
            )
        )
        for direction_gauge, site_irrep_basis in (
            (False, False),
            (False, True),
            (True, False),
            (True, True),
        ):
            if site_irrep_basis and not site_scope:
                continue
            family_span = _family_span_proof(
                web_modes,
                local_modes,
                tol,
                vector_kind=vector_kind,
                candidate_cache=candidate_cache,
                allow_parent_direction_gauge=direction_gauge,
                allow_site_irrep_basis=site_irrep_basis,
            )
            if family_span is not None:
                break
    if not count_match:
        status = "count"
    elif hard_statuses <= {"ok"}:
        status = "ok"
    elif hard_statuses <= {"ok", "sign"}:
        status = "sign"
    elif signed_permutation:
        status = "signed_permutation"
    elif family_span is not None:
        status = "family_basis"
    elif "position" in hard_statuses:
        status = "position"
    else:
        status = "vector"
    if not count_match:
        strict_status = "count"
    elif signed_permutation and status == "signed_permutation":
        strict_status = "signed_permutation"
    elif family_span is not None and status == "family_basis":
        strict_status = "family_basis"
    elif strict_statuses <= {"ok"}:
        strict_status = "ok"
    else:
        strict_status = next(
            (candidate for candidate in ("position", "vector", "sign", "label", "normfactor") if candidate in strict_statuses),
            "fail",
        )
    payload = {
        "status": status,
        "strict_status": strict_status,
        "count_match": count_match,
        "web_count": len(web_modes),
        "local_count": len(local_modes),
        "unpaired_web_tail_labels": [
            mode.label for mode in web_modes[len(local_modes):]
        ],
        "unpaired_local_tail_labels": [
            str(mode.get("label") or "") for mode in local_modes[len(web_modes):]
        ],
        "modes": modes,
    }
    if family_span is not None:
        payload["physical_family_basis"] = family_span
    return payload, signed_permutation_assignment


def _vector_modes_probe(
    web_modes: list[DisplaciveModeDefinition]
    | tuple[DisplaciveModeDefinition, ...],
    local_modes: list[dict[str, Any]],
    tol: float,
    *,
    vector_kind: str,
) -> dict[str, Any]:
    """Compute aggregate adoption status without retaining row diagnostics."""

    payload, _ = _vector_modes_payload_with_assignment(
        web_modes,
        local_modes,
        tol,
        vector_kind=vector_kind,
        materialize_rows=False,
    )
    return payload


def _transported_mode_dicts(
    local_modes: list[dict[str, Any]],
    transport: BasisEquivalentModeFrameTransport,
) -> list[dict[str, Any]] | None:
    if len(local_modes) != len(transport.modes):
        return None
    out: list[dict[str, Any]] = []
    for mode_index, (source, mode) in enumerate(zip(local_modes, transport.modes, strict=True)):
        source_rows = source.get("rows")
        if (
            mode.source_mode_index != mode_index
            or not isinstance(source_rows, (list, tuple))
            or len(source_rows) != len(mode.rows)
        ):
            return None
        rows = []
        for row_index, (source_row, row) in enumerate(zip(source_rows, mode.rows, strict=True)):
            if not isinstance(source_row, dict) or row.source_row_index != row_index:
                return None
            rows.append({**source_row, "xyz": row.xyz, "dxyz": row.dxyz})
        out.append({**source, "label": mode.label, "normfactor": mode.normfactor, "rows": rows})
    return out


def _mode_position_vectors_are_unambiguous(
    web_mode: DisplaciveModeDefinition,
    local_mode: dict[str, Any],
    tol: float,
    candidate_cache: _PeriodicCandidateCache | None = None,
) -> bool:
    local_source_rows = local_mode.get("rows") or []
    local_rows = list(local_source_rows)
    candidates = (
        _periodic_xyz_candidate_indices(list(web_mode.rows), local_rows, tol)
        if candidate_cache is None
        else candidate_cache.get(web_mode.rows, local_source_rows)
    )
    if candidates is None or len(candidates) != len(web_mode.rows):
        return False
    for indices in candidates:
        if not indices:
            return False
        vectors = [local_rows[index].get("dxyz") or [] for index in indices]
        if any(not _vec_close(vectors[0], vector, tol) for vector in vectors[1:]):
            return False
    return True


def _basis_equivalent_mode_payload(
    web_modes: list[DisplaciveModeDefinition] | tuple[DisplaciveModeDefinition, ...],
    local_modes: list[dict[str, Any]],
    *,
    vector_kind: str,
    selected_state: object,
    subgroup: dict[str, Any],
    lattice: dict[str, Any],
    atoms: dict[str, Any],
    tol: float,
) -> dict[str, Any] | None:
    """Adopt mode rows only after the complete equivalent-frame proof closes."""

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
        or len(web_modes) != len(local_modes)
        or not web_modes
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
    transport = basis_equivalent_mode_frame_transport(
        local_modes,
        basis_change=lattice.get("basis_change"),
        request_origin=request[3],
        local_origin=local_selected[3],
        request_child_space_group=request_child,
        local_child_space_group=local_child,
        vector_kind=vector_kind,
        request_basis=request[2],
    )
    if transport is None:
        return None
    transformed_modes = _transported_mode_dicts(local_modes, transport)
    if transformed_modes is None:
        return None
    payload, signed_permutation_assignment = _vector_modes_payload_with_assignment(
        web_modes,
        transformed_modes,
        tol,
        vector_kind=vector_kind,
        require_unambiguous_signed_permutation=True,
    )
    mode_payloads = payload.get("modes")
    common_invalid = (
        payload.get("count_match") is not True
        or not isinstance(mode_payloads, list)
        or len(mode_payloads) != len(web_modes)
    )
    if common_invalid:
        return None
    if payload.get("status") in {"ok", "sign"}:
        if (
            any(
                not _mode_position_vectors_are_unambiguous(web_mode, local_mode, tol)
                for web_mode, local_mode in zip(web_modes, transformed_modes, strict=True)
            )
            or any(
                mode.get("status") not in {"ok", "sign"}
                or mode.get("count_match") is not True
                for mode in mode_payloads
            )
        ):
            return None
    elif payload.get("status") == "signed_permutation":
        if signed_permutation_assignment is None:
            return None
    elif payload.get("status") == "family_basis":
        if not isinstance(payload.get("physical_family_basis"), dict):
            return None
    else:
        return None
    return {
        **payload,
        "validated": True,
        "provenance": transport.provenance,
        "basis_change": transport.basis_change,
        "condition_number": transport.condition_number,
        "vector_kind": transport.vector_kind,
        **(
            {"physical_signed_permutation": True}
            if payload.get("status") == "signed_permutation"
            else {}
        ),
    }


def _strain_mode_payload(web_mode: StrainModeDefinition, local_mode: dict[str, Any], tol: float) -> dict[str, Any]:
    local_components = tuple(local_mode.get("components") or local_mode.get("tensor") or ())
    vector_tol = _visible_vector_tolerance(tol)
    components_match = _vec_close(web_mode.components, local_components, vector_tol)
    sign_match = _vec_close(
        web_mode.components,
        tuple(-float(value) for value in local_components),
        vector_tol,
    )
    label_match = web_mode.label == local_mode.get("label")
    normfactor_match = _optional_float_close(web_mode.normfactor, local_mode.get("normfactor"), tol)
    physical_status = "ok" if components_match else ("sign" if sign_match else "tensor")
    if physical_status != "ok":
        strict_status = physical_status
    elif not label_match:
        strict_status = "label"
    elif not normfactor_match:
        strict_status = "normfactor"
    else:
        strict_status = "ok"
    return {
        "status": physical_status,
        "strict_status": strict_status,
        "label_match": label_match,
        "normfactor_match": normfactor_match,
        "web": {"label": web_mode.label, "normfactor": web_mode.normfactor, "components": web_mode.components},
        "local": {
            "label": local_mode.get("label"),
            "normfactor": local_mode.get("normfactor"),
            "components": local_components,
        },
    }


def _strain_modes_payload(web_modes: tuple[StrainModeDefinition, ...], local_modes: list[dict[str, Any]], tol: float) -> dict[str, Any]:
    count_match = len(web_modes) == len(local_modes)
    modes = [
        _strain_mode_payload(web_mode, local_mode, tol)
        for web_mode, local_mode in zip(web_modes, local_modes, strict=False)
    ]
    statuses = {mode["status"] for mode in modes}
    strict_statuses = {mode["strict_status"] for mode in modes}
    if not count_match:
        status = strict_status = "count"
    else:
        status = "ok" if statuses <= {"ok"} else ("sign" if statuses <= {"ok", "sign"} else "tensor")
        strict_status = "ok" if strict_statuses <= {"ok"} else next(
            (candidate for candidate in ("tensor", "sign", "label", "normfactor") if candidate in strict_statuses),
            "fail",
        )
    return {
        "status": status,
        "strict_status": strict_status,
        "count_match": count_match,
        "web_count": len(web_modes),
        "local_count": len(local_modes),
        "modes": modes,
    }
