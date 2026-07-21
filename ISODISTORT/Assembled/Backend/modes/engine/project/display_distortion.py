"""Source DISPLAY DISTORTION projection on an explicit occurrence cell."""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from numbers import Integral
from collections.abc import Mapping, Sequence
from typing import Any, Callable, TypeVar

import numpy as np

from ISODISTORT.Assembled.Backend.exactmath import integer_determinant3

_LayoutT = TypeVar("_LayoutT")


@dataclass(frozen=True)
class OccurrenceSublattice:
    """Exact finite-index relation between Source and display cell bases."""

    index: int
    transition: tuple[int, ...]


@dataclass(frozen=True)
class DisplayProjectionFamily:
    """One Source project family resolved over every display occurrence."""

    family: int
    modes: tuple[tuple[tuple[float, float, float], ...], ...]


@dataclass(frozen=True)
class DisplayProjectionBlock:
    """Complete projected output for one Source project preparation."""

    prep_index: int
    pg_irrep: int
    project_count: int
    vector_dim: int
    bridge_project_counts: tuple[int, ...]
    families: tuple[DisplayProjectionFamily, ...]


@dataclass(frozen=True)
class _DisplayProjectionPrep:
    prep_index: int
    pg_irrep: int
    project_count: int
    site_pg: int
    vector_basis_id: int
    vector_dim: int
    bridge_project_counts: tuple[int, ...]
    active_values: tuple[tuple[object, object], ...]


def _exact_basis_determinant(values: Sequence[object] | None) -> int | None:
    if (
        values is None
        or isinstance(values, (str, bytes, bytearray))
        or len(values) != 9
        or any(
            isinstance(value, bool) or not isinstance(value, Integral)
            for value in values
        )
    ):
        return None
    return integer_determinant3(values)  # type: ignore[arg-type]


def occurrence_sublattice(
    source_basis: Sequence[object] | None,
    occurrence_basis: Sequence[object] | None,
) -> OccurrenceSublattice | None:
    """Prove ``occurrence_basis = Q @ source_basis`` for integral index ``Q``."""

    source_det = _exact_basis_determinant(source_basis)
    occurrence_det = _exact_basis_determinant(occurrence_basis)
    if (
        source_basis is None
        or occurrence_basis is None
        or source_det in (None, 0)
        or occurrence_det in (None, 0)
    ):
        return None
    source = tuple(
        tuple(Fraction(int(source_basis[row * 3 + col])) for col in range(3))
        for row in range(3)
    )
    occurrence = tuple(
        tuple(Fraction(int(occurrence_basis[row * 3 + col])) for col in range(3))
        for row in range(3)
    )
    a, b, c = source[0]
    d, e, f = source[1]
    g, h, i = source[2]
    inverse = (
        (
            Fraction(e * i - f * h, source_det),
            Fraction(c * h - b * i, source_det),
            Fraction(b * f - c * e, source_det),
        ),
        (
            Fraction(f * g - d * i, source_det),
            Fraction(a * i - c * g, source_det),
            Fraction(c * d - a * f, source_det),
        ),
        (
            Fraction(d * h - e * g, source_det),
            Fraction(b * g - a * h, source_det),
            Fraction(a * e - b * d, source_det),
        ),
    )
    transition = tuple(
        tuple(
            sum(occurrence[row][axis] * inverse[axis][col] for axis in range(3))
            for col in range(3)
        )
        for row in range(3)
    )
    if any(value.denominator != 1 for row in transition for value in row):
        return None
    transition_values = tuple(int(value) for row in transition for value in row)
    transition_det = _exact_basis_determinant(transition_values)
    if transition_det is None:
        return None
    index = abs(transition_det)
    if index <= 1 or abs(occurrence_det) != abs(source_det) * index:
        return None
    return OccurrenceSublattice(index=index, transition=transition_values)


def source_display_projection_families(
    decoder: Any,
    case: Any,
    *,
    gid: object,
    full_dim: object,
    site_pg: object,
    pg_irrep: object,
    vector_basis_id: object,
    vector_dim: object,
    vector_setting: object,
    project_count: object,
    project_active_values: Sequence[Sequence[object]],
    direction_rows: Sequence[Sequence[float]],
    occurrence_records: Sequence[Sequence[object]],
) -> tuple[DisplayProjectionFamily, ...] | None:
    """Apply Source ``mode_b0`` and ``project_vector_`` on every occurrence."""

    integer_args = (
        gid,
        full_dim,
        site_pg,
        pg_irrep,
        vector_basis_id,
        vector_dim,
        vector_setting,
        project_count,
    )
    if any(
        isinstance(value, bool) or not isinstance(value, Integral)
        for value in integer_args
    ):
        return None
    gid = int(gid)
    full_dim = int(full_dim)
    site_pg = int(site_pg)
    pg_irrep = int(pg_irrep)
    vector_basis_id = int(vector_basis_id)
    vector_dim = int(vector_dim)
    vector_setting = int(vector_setting)
    project_count = int(project_count)
    if (
        gid <= 0
        or full_dim <= 0
        or full_dim > 48
        or site_pg <= 0
        or pg_irrep <= 0
        or vector_basis_id <= 0
        or vector_dim not in (1, 2, 3)
        or vector_setting not in (1, 2)
        or project_count <= 0
        or not occurrence_records
    ):
        return None
    try:
        direction = np.asarray(direction_rows, dtype=float)
    except (OverflowError, TypeError, ValueError):
        return None
    if (
        direction.ndim != 2
        or direction.shape[0] <= 0
        or direction.shape[0] > 48
        or direction.shape[1] != full_dim
        or not np.isfinite(direction).all()
    ):
        return None

    family_stride = 144
    basis_length = max(4608, project_count * family_stride)
    basis = np.zeros(basis_length, dtype=float)
    seen_indices: set[int] = set()
    nonzero_families: set[int] = set()
    try:
        for active in project_active_values:
            if isinstance(active, (str, bytes, bytearray)) or len(active) != 2:
                return None
            index, value = active
            if (
                isinstance(index, bool)
                or not isinstance(index, Integral)
                or isinstance(value, bool)
            ):
                return None
            index = int(index)
            numeric = float(value)
            if (
                not 0 <= index < basis_length
                or index // family_stride >= project_count
                or index % family_stride >= full_dim * 3
                or not np.isfinite(numeric)
                or index in seen_indices
            ):
                return None
            seen_indices.add(index)
            basis[index] = numeric
            if numeric != 0.0:
                nonzero_families.add(index // family_stride)
    except (OverflowError, TypeError, ValueError):
        return None
    if not seen_indices or nonzero_families != set(range(project_count)):
        return None

    family_modes = [
        [
            [[0.0, 0.0, 0.0] for _record in occurrence_records]
            for _row in range(direction.shape[0])
        ]
        for _family in range(project_count)
    ]
    for occurrence_index, values in enumerate(occurrence_records):
        try:
            valid_record = (
                not isinstance(values, (str, bytes, bytearray))
                and len(values) == 5
                and all(
                    not isinstance(value, bool) and isinstance(value, Integral)
                    for value in values
                )
            )
        except TypeError:
            valid_record = False
        if not valid_record:
            return None
        record = tuple(int(value) for value in values)
        if record[3] <= 0 or record[4] <= 0:
            return None
        try:
            irrep = np.asarray(
                decoder.little_phase_matrix_by_gid_record_for_case(
                    gid,
                    record,
                    case,
                ),
                dtype=complex,
            )
        except (IndexError, KeyError, OverflowError, TypeError, ValueError):
            return None
        if irrep.shape != (full_dim, full_dim) or not np.isfinite(irrep).all():
            return None

        mode_b0 = np.zeros(basis_length, dtype=complex)
        for family in range(project_count):
            family_base = family * family_stride
            for component in range(3):
                source_vector = np.asarray(
                    [
                        basis[family_base + row * 3 + component]
                        for row in range(full_dim)
                    ],
                    dtype=complex,
                )
                selected = direction @ (irrep @ source_vector)
                for row, value in enumerate(selected):
                    mode_b0[family_base + row * 3 + component] = value
        if not np.isfinite(mode_b0).all() or np.max(np.abs(mode_b0.imag)) > 1e-9:
            return None
        try:
            emitted, output = decoder.project_vector_from_boundary(
                site_pg=site_pg,
                vector_basis_id=vector_basis_id,
                target_vector_rep=pg_irrep,
                atom_count=int(direction.shape[0]),
                vector_dim=vector_dim,
                vector_setting=vector_setting,
                point_op=int(record[4]),
                project_count=project_count,
                basis_function=mode_b0.real,
                output_length=basis_length,
            )
        except (IndexError, KeyError, OverflowError, TypeError, ValueError):
            return None
        try:
            valid_output = (
                not isinstance(emitted, bool)
                and isinstance(emitted, Integral)
                and int(emitted) == project_count
                and len(output) >= project_count * family_stride
            )
        except TypeError:
            valid_output = False
        if not valid_output:
            return None
        for family in range(project_count):
            for row in range(direction.shape[0]):
                offset = family * family_stride + row * 3
                try:
                    vector = tuple(float(output[offset + axis]) for axis in range(3))
                except (IndexError, OverflowError, TypeError, ValueError):
                    return None
                if not all(np.isfinite(value) for value in vector):
                    return None
                family_modes[family][row][occurrence_index] = list(vector)

    for modes in family_modes:
        tensor = np.asarray(modes, dtype=float)
        if tensor.size == 0 or float(np.max(np.abs(tensor))) <= 1e-12:
            return None

    return tuple(
        DisplayProjectionFamily(
            family=family,
            modes=tuple(
                tuple(tuple(float(value) for value in vector) for vector in mode)
                for mode in modes
            ),
        )
        for family, modes in enumerate(family_modes)
    )


def _exact_integral(value: object) -> bool:
    return not isinstance(value, bool) and isinstance(value, Integral)


def _freeze_projection_preps(
    *,
    gid: int,
    raw_preps: Sequence[object],
    project_basis_by_branch: Mapping[tuple[int, int], object],
    bridge_sources_by_branch: Mapping[tuple[int, int], Sequence[object]],
) -> tuple[_DisplayProjectionPrep, ...] | None:
    if not raw_preps:
        return None
    frozen: list[_DisplayProjectionPrep] = []
    for prep_index, prep in enumerate(raw_preps):
        if not isinstance(prep, Mapping):
            return None
        fields = (
            "gid",
            "pg_irrep",
            "project_count",
            "site_pg",
            "vector_basis_id",
            "vector_dim",
        )
        if any(not _exact_integral(prep.get(field)) for field in fields):
            return None
        if int(prep["gid"]) != gid:
            return None
        pg_irrep = int(prep["pg_irrep"])
        project_count = int(prep["project_count"])
        site_pg = int(prep["site_pg"])
        vector_basis_id = int(prep["vector_basis_id"])
        vector_dim = int(prep["vector_dim"])
        if (
            pg_irrep <= 0
            or project_count <= 0
            or site_pg <= 0
            or vector_basis_id <= 0
            or vector_dim not in (1, 2, 3)
        ):
            return None
        project_basis = project_basis_by_branch.get((gid, pg_irrep))
        sources = bridge_sources_by_branch.get((gid, pg_irrep))
        if not isinstance(project_basis, Mapping) or not isinstance(sources, Sequence):
            return None
        if (
            not _exact_integral(project_basis.get("count"))
            or int(project_basis["count"]) != project_count
            or not _exact_integral(project_basis.get("family_stride"))
            or int(project_basis["family_stride"]) != 144
        ):
            return None
        active_values = project_basis.get("active_values")
        if not isinstance(active_values, Sequence) or isinstance(
            active_values, (str, bytes, bytearray)
        ):
            return None
        frozen_active: list[tuple[object, object]] = []
        for active in active_values:
            if (
                not isinstance(active, Sequence)
                or isinstance(active, (str, bytes, bytearray))
                or len(active) != 2
            ):
                return None
            frozen_active.append((active[0], active[1]))
        if len(sources) != project_count:
            return None
        bridge_counts: list[int] = []
        for family, source in enumerate(sources):
            if (
                not isinstance(source, Mapping)
                or not _exact_integral(source.get("gid"))
                or int(source["gid"]) != gid
                or not _exact_integral(source.get("pg_irrep"))
                or int(source["pg_irrep"]) != pg_irrep
                or not _exact_integral(source.get("family"))
                or int(source["family"]) != family
                or not _exact_integral(source.get("bridge_project_count"))
                or int(source["bridge_project_count"]) < 0
            ):
                return None
            bridge_counts.append(int(source["bridge_project_count"]))
        frozen.append(
            _DisplayProjectionPrep(
                prep_index=prep_index,
                pg_irrep=pg_irrep,
                project_count=project_count,
                site_pg=site_pg,
                vector_basis_id=vector_basis_id,
                vector_dim=vector_dim,
                bridge_project_counts=tuple(bridge_counts),
                active_values=tuple(frozen_active),
            )
        )
    return tuple(frozen)


def source_display_projection_plan(
    decoder: Any,
    case: Any,
    *,
    gid: int,
    full_dim: int,
    vector_setting: int,
    direction_rows: Sequence[Sequence[float]],
    occurrence_records: Sequence[Sequence[object]],
    raw_preps: Sequence[object],
    project_basis_by_branch: Mapping[tuple[int, int], object],
    bridge_sources_by_branch: Mapping[tuple[int, int], Sequence[object]],
) -> tuple[DisplayProjectionBlock, ...] | None:
    """Build every Source display block or decline the optional route."""

    preps = _freeze_projection_preps(
        gid=gid,
        raw_preps=raw_preps,
        project_basis_by_branch=project_basis_by_branch,
        bridge_sources_by_branch=bridge_sources_by_branch,
    )
    if preps is None:
        return None
    blocks: list[DisplayProjectionBlock] = []
    for prep in preps:
        families = source_display_projection_families(
            decoder,
            case,
            gid=gid,
            full_dim=full_dim,
            site_pg=prep.site_pg,
            pg_irrep=prep.pg_irrep,
            vector_basis_id=prep.vector_basis_id,
            vector_dim=prep.vector_dim,
            vector_setting=vector_setting,
            project_count=prep.project_count,
            project_active_values=prep.active_values,
            direction_rows=direction_rows,
            occurrence_records=occurrence_records,
        )
        if families is None or tuple(family.family for family in families) != tuple(
            range(prep.project_count)
        ):
            return None
        blocks.append(
            DisplayProjectionBlock(
                prep_index=prep.prep_index,
                pg_irrep=prep.pg_irrep,
                project_count=prep.project_count,
                vector_dim=prep.vector_dim,
                bridge_project_counts=prep.bridge_project_counts,
                families=families,
            )
        )
    return tuple(blocks)


def resolve_display_projection(
    candidate_layout: _LayoutT | None,
    build_candidate: Callable[[_LayoutT], tuple[DisplayProjectionBlock, ...] | None],
    build_baseline: Callable[[], _LayoutT | None],
) -> tuple[_LayoutT | None, tuple[DisplayProjectionBlock, ...], bool]:
    """Commit a complete candidate layout and projection together, or neither."""

    if candidate_layout is not None:
        try:
            projected = build_candidate(candidate_layout)
        except (
            IndexError,
            KeyError,
            OverflowError,
            TypeError,
            ValueError,
            ZeroDivisionError,
        ):
            projected = None
        if projected:
            return candidate_layout, projected, True
    return build_baseline(), (), False
