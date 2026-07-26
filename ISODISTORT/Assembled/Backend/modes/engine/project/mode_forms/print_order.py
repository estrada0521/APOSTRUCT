"""mode-kernel print-order and real-basis shaping helpers."""

from __future__ import annotations

import math

from ISODISTORT.Assembled.Backend.modes.engine.project.mode_forms.normalization import (
    _active_mode_indices,
    _dominant_mode_component,
    _normalize_mode_vectors,
    _same_mode,
)

def _regroup_repeated_component_modes(modes: list[list[list[float]]]) -> list[list[list[float]]]:
    """Match mode-kernel print order when a component block is split by another.

    In the guarded type-3 real-form path, paired modes are printed by vector
    component.  The projection bridge can emit ``a,a,b,c,b,c`` when two site
    vector irreps share the same source row.  Re-group only that exact shape so
    alternating legitimate patterns such as ``a,b,a,b`` are left untouched.
    """

    if len(modes) != 6:
        return modes
    components = [_dominant_mode_component(mode) for mode in modes]
    if (
        components[0] == components[1]
        and components[2] == components[4]
        and components[3] == components[5]
        and len({components[0], components[2], components[3]}) == 3
    ):
        return [modes[index] for index in (0, 1, 2, 4, 3, 5)]
    return modes
def _type3_real_print_tail(mode: list[list[float]]) -> list[list[float]]:
    active = [
        index
        for index, vector in enumerate(mode)
        if any(abs(float(component)) > 1e-12 for component in vector)
    ]
    if len(active) < 3 or len(active) % 3 != 0:
        return mode
    drop_count = len(active) // 3
    drop = set(active[:drop_count])
    negative = set(active[drop_count:2 * drop_count])
    out = [[0.0, 0.0, 0.0] for _ in mode]
    for index in active[drop_count:]:
        if index not in drop:
            scale = -2.0 if index in negative else 2.0
            out[index] = [scale * float(component) for component in mode[index]]
    return _normalize_mode_vectors(out)
def _type3_real_print_modes(modes: list[list[list[float]]]) -> list[list[list[float]]]:
    """Split repeated type-3 real-form pairs into mode-kernel printed partners."""

    if len(modes) % 4 != 0:
        return modes
    out: list[list[list[float]]] = []
    changed = False
    for offset in range(0, len(modes), 4):
        chunk = modes[offset:offset + 4]
        active_counts = [
            sum(
                1
                for vector in mode
                if any(abs(float(component)) > 1e-12 for component in vector)
            )
            for mode in chunk
        ]
        if (
            len(chunk) == 4
            and active_counts == [12, 12, 12, 12]
            and _same_mode(chunk[0], chunk[2])
            and _same_mode(chunk[1], chunk[3])
        ):
            out.extend([
                chunk[0],
                _type3_real_print_tail(chunk[0]),
                chunk[1],
                _type3_real_print_tail(chunk[1]),
            ])
            changed = True
        elif (
            len(chunk) == 4
            and active_counts == [12, 12, 12, 12]
            and _same_mode(chunk[0], chunk[1])
            and _same_mode(chunk[2], chunk[3])
        ):
            out.extend([
                chunk[0],
                _type3_real_print_tail(chunk[0]),
                chunk[2],
                _type3_real_print_tail(chunk[2]),
            ])
            changed = True
        else:
            out.extend(chunk)
    return out if changed else modes
def _type3_parametric_component_print_order(modes: list[list[list[float]]]) -> list[list[list[float]]]:
    """Group alternating type-3 parametric print partners by vector component."""

    if len(modes) % 4 != 0:
        return modes
    out: list[list[list[float]]] = []
    changed = False
    for offset in range(0, len(modes), 4):
        chunk = modes[offset:offset + 4]
        components = [_dominant_mode_component(mode) for mode in chunk]
        if (
            len(chunk) == 4
            and components[0] == components[2]
            and components[1] == components[3]
            and components[0] != components[1]
        ):
            out.extend([chunk[0], chunk[2], chunk[1], chunk[3]])
            changed = True
        else:
            out.extend(chunk)
    return out if changed else modes
def _type3_parametric_scalar_plane_print_order(modes: list[list[list[float]]]) -> list[list[list[float]]]:
    if len(modes) != 6:
        return modes
    counts = [len(_active_mode_indices(mode)) for mode in modes]
    if counts == [12, 12, 8, 4, 8, 4]:
        return [modes[0], modes[1], modes[2], modes[4], modes[3], modes[5]]
    return modes
def _type3_parametric_second_arm_phase(modes: list[list[list[float]]]) -> list[list[list[float]]]:
    """Apply the second-arm sign convention for paired type-3 modes."""

    if len(modes) != 6:
        return modes
    counts = [len(_active_mode_indices(mode)) for mode in modes]
    if counts != [12, 12, 8, 8, 8, 8]:
        return modes
    half_index = len(modes[0]) // 2 if modes and modes[0] else 0
    out: list[list[list[float]]] = []
    for mode_index, mode in enumerate(modes):
        if mode_index % 2 == 0:
            out.append(mode)
            continue
        phased = []
        for atom_index, vector in enumerate(mode):
            sign = -1.0 if atom_index >= half_index else 1.0
            phased.append([sign * float(component) for component in vector])
        out.append(phased)
    return out
def _combine_modes(
    modes: list[list[list[float]]],
    terms: tuple[tuple[int, float], ...],
) -> list[list[float]]:
    atom_count = len(modes[0]) if modes else 0
    out = [[0.0, 0.0, 0.0] for _ in range(atom_count)]
    for mode_index, coefficient in terms:
        mode = modes[mode_index]
        for atom_index, vector in enumerate(mode):
            for axis in range(3):
                out[atom_index][axis] += coefficient * float(vector[axis])
    return _normalize_mode_vectors(out)
def _type3_parametric_kdim2_print_basis(modes: list[list[list[float]]]) -> list[list[list[float]]]:
    """Apply mode-kernel real-basis print convention for two-parameter k pairs.

    For hexagonal/trigonal two-parameter k labels, the printed ``vmode`` basis
    is rotated after project-vector expansion: the first scalar/plane partners are kept,
    while the following plane partners are mixed by ``2/sqrt(3)``.  This is the
    same final-basis transform used by C/D-family paired type-3 blocks; it is
    guarded by the component patterns rather than space-group labels.
    """

    components = [_dominant_mode_component(mode) for mode in modes]
    factor = 2.0 / math.sqrt(3.0)
    if len(modes) == 6 and components == [2, 2, 0, 0, 0, 0]:
        return [
            modes[0],
            modes[1],
            modes[2],
            modes[3],
            _combine_modes(modes, ((2, factor), (5, 1.0))),
            _combine_modes(modes, ((3, factor), (4, -1.0))),
        ]
    if (
        len(modes) == 12
        and components[:8] == [2, 2, 2, 2, 0, 0, 1, 1]
        and components[8:] in ([1, 1, 0, 0], [0, 0, 1, 1])
    ):
        active_counts = [len(_active_mode_indices(mode)) for mode in modes]
        if max(active_counts, default=0) < 4 and len(set(active_counts)) != 1:
            return modes
        return [
            modes[0],
            modes[1],
            modes[2],
            modes[3],
            modes[4],
            modes[5],
            _combine_modes(modes, ((4, -factor), (6, -1.0))),
            _combine_modes(modes, ((5, -factor), (7, -1.0))),
            modes[8],
            modes[9],
            modes[10],
            modes[11],
        ]
    return modes
def _type1_parametric_scalar_plane_print_order(
    modes: list[list[list[float]]],
    *,
    allow_support_pairing: bool = True,
) -> list[list[list[float]]]:
    if len(modes) != 6:
        return modes
    components = [_dominant_mode_component(mode) for mode in modes]
    counts = [len(_active_mode_indices(mode)) for mode in modes]
    if (
        components[0] == components[1]
        and components[2] == components[3]
        and components[4] == components[5]
        and len({components[0], components[2], components[4]}) == 3
    ):
        active = [_active_mode_indices(mode) for mode in modes]
        cross_component_families = (
            (active[2] <= active[4] or active[4] <= active[2])
            and (active[3] <= active[5] or active[5] <= active[3])
        )
        if cross_component_families and (
            allow_support_pairing or counts == [12, 12, 8, 4, 8, 4]
        ):
            return [modes[0], modes[1], modes[2], modes[4], modes[3], modes[5]]
    return modes
def _has_scalar_plane_vector_split(
    preps: list[dict[str, object]],
    vector_slots: tuple[int, int, int],
) -> bool:
    """Return true when site-vector irreps split one axis from a plane."""

    if len(set(vector_slots)) != 2:
        return False
    counts = {slot: vector_slots.count(slot) for slot in set(vector_slots)}
    if sorted(counts.values()) != [1, 2]:
        return False
    prep_pg_irreps = {
        int(prep.get("pg_irrep", 0))
        for prep in preps
        if int(prep.get("project_count", 0)) == 1 and int(prep.get("vector_dim", 0)) == 1
    }
    return set(counts) <= prep_pg_irreps
