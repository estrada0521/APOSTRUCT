"""Source-only dynamic-k subgroup-to-order-parameter selection.

This module computes the dynamic invariant-direction surface for ISO
``DISPLAY DIRECTION``.  It consumes a selected child subgroup embedding and
returns invariant directions in Source k-family/irrep order.  Web presentation
and FINAL values are deliberately outside this module.
"""

from __future__ import annotations

from collections import Counter
from collections.abc import Callable, Hashable, Iterable, Sequence
from dataclasses import dataclass, replace
from fractions import Fraction
from functools import lru_cache
from itertools import product
import math
from typing import Any

import numpy as np

from APOSTRUCT.Backend.exactmath import (
    fraction_matrix_inverse3,
)
from APOSTRUCT.Backend.lattice_quotient import (
    reciprocal_mesh_source_order as _ordered_reciprocal_mesh,
)
from APOSTRUCT.Backend.modes.engine.input import Case


OperationRecord = tuple[int, int, int, int, int]


def _strict_integral_values(values: Sequence[Any], *, size: int = 4) -> tuple[int, ...] | None:
    """Normalize proof k-parameters without accepting lossy ``int()`` casts."""
    try:
        if isinstance(values, (str, bytes, bytearray)) or len(values) != size:
            return None
    except (TypeError, ValueError):
        return None
    normalized: list[int] = []
    for value in values:
        if isinstance(value, bool):
            return None
        try:
            fraction = Fraction(str(value))
        except (TypeError, ValueError, ZeroDivisionError):
            return None
        if fraction.denominator != 1:
            return None
        normalized.append(int(fraction))
    return tuple(normalized)


@dataclass(frozen=True)
class DynamicKOccurrence:
    """One exact reciprocal-mesh occurrence before presentation selection."""

    mesh_order: int
    k_label: str
    source_kparam: tuple[int, int, int, int]
    reciprocal_vector_pml: tuple[Fraction, Fraction, Fraction]
    inversion_bucket: tuple[tuple[int, int], ...]


@dataclass(frozen=True)
class DynamicSubductionRow:
    gid: int
    irrep_label: str
    k_label: str
    source_kparam: tuple[int, int, int, int]
    reciprocal_vector_pml: tuple[Fraction, Fraction, Fraction]
    direction_matrix: tuple[tuple[float, ...], ...]
    source_occurrences: tuple[DynamicKOccurrence, ...] = ()
    carrier_source_kparam: tuple[int, int, int, int] | None = None


@dataclass(frozen=True)
class OccurrenceAliasCandidate:
    """One exact stored-occurrence class that can be replayed as a heading."""

    representative_gid: int
    representative_source_kparam: tuple[int, int, int, int]
    candidate_source_kparam: tuple[int, int, int, int]
    heading_identity: Hashable
    candidate_row: DynamicSubductionRow


def _reciprocal_translation_key(
    values: Sequence[Fraction],
) -> tuple[Fraction, Fraction, Fraction]:
    """Identify one PML reciprocal vector modulo reciprocal translations."""

    if len(values) != 3:
        raise ValueError("Source occurrence reciprocal vector must have length 3")
    return tuple(Fraction(value) % 1 for value in values)  # type: ignore[return-value]


def _source_parameter_translation_key(
    values: tuple[int, int, int, int],
) -> tuple[Fraction, Fraction, Fraction]:
    if values[3] == 0:
        raise ValueError("Source occurrence k-parameter denominator must be nonzero")
    return tuple(
        Fraction(values[index], values[3]) % 1 for index in range(3)
    )  # type: ignore[return-value]


def stored_occurrence_alias_candidates(
    decoder: Any,
    rows: Sequence[DynamicSubductionRow],
    *,
    aliases_for_class: Callable[
        [DynamicSubductionRow, tuple[int, int, int, int]],
        Sequence[DynamicSubductionRow],
    ],
    heading_identity_for_row: Callable[[DynamicSubductionRow], Hashable | None],
    represented_occurrence_for_row: Callable[
        [DynamicSubductionRow],
        tuple[
            tuple[int, int, int, int],
            tuple[Fraction, Fraction, Fraction],
            bool,
        ]
        | None,
    ],
    represented_heading_identities: Iterable[Hashable] = (),
    allowed_irrep_types: Iterable[int] | None = None,
) -> tuple[OccurrenceAliasCandidate, ...]:
    """Return unambiguous aliases named by stored Source occurrence classes.

    Dynamic family discovery can merge several exact ``source_kparam``
    classes into one row. The caller supplies the occurrence currently used
    to represent each raw family; only the other stored classes are eligible.
    Aliases discovered by re-solving arbitrary reciprocal vectors are not.
    Occupied-site emission remains a separate downstream admission gate.
    """

    seen_headings: set[Hashable] = set()
    for heading in represented_heading_identities:
        try:
            hash(heading)
            seen_headings.add(heading)
        except TypeError:
            continue

    allowed_types = (
        None
        if allowed_irrep_types is None
        else frozenset(int(value) for value in allowed_irrep_types)
    )
    candidates: list[OccurrenceAliasCandidate] = []
    for row in rows:
        try:
            gid_values = _strict_integral_values((row.gid,), size=1)
            if gid_values is None or gid_values[0] <= 0:
                continue
            little = decoder.little_record_by_gid(gid_values[0])
            if int(little.old_id) != 0 or (
                allowed_types is not None and int(little.irrep_type) not in allowed_types
            ):
                continue
        except (IndexError, KeyError, TypeError, ValueError):
            continue
        represented = represented_occurrence_for_row(row)
        if represented is None:
            raise ValueError("Source occurrence family has no represented row")
        represented_values = _strict_integral_values(represented[0])
        if represented_values is None:
            raise ValueError("represented Source occurrence has invalid k parameters")
        represented_source_kparam = represented_values
        represented_pml = _reciprocal_translation_key(represented[1])
        promoted = bool(represented[2])
        represented_source_key = _source_parameter_translation_key(
            represented_source_kparam
        )

        classes: list[
            tuple[
                tuple[int, int, int, int],
                list[DynamicKOccurrence],
            ]
        ] = []
        class_indices: dict[tuple[int, int, int, int], int] = {}
        for occurrence in row.source_occurrences:
            try:
                source_kparam_values = _strict_integral_values(occurrence.source_kparam)
                if source_kparam_values is None:
                    continue
                source_kparam = (
                    source_kparam_values[0],
                    source_kparam_values[1],
                    source_kparam_values[2],
                    source_kparam_values[3],
                )
            except (AttributeError, TypeError, ValueError):
                continue
            index = class_indices.get(source_kparam)
            if index is None:
                class_indices[source_kparam] = len(classes)
                classes.append((source_kparam, [occurrence]))
            else:
                classes[index][1].append(occurrence)

        if len(classes) <= 1:
            continue
        source_classes = {
            source_kparam
            for source_kparam, _occurrences in classes
            if _source_parameter_translation_key(source_kparam)
            == represented_source_key
        }
        # A promoted heading may use a standardized k absent from the raw
        # Source classes; its exact mesh occurrence remains authoritative.
        if len(source_classes) > 1 or (not promoted and not source_classes):
            raise ValueError(
                "represented k parameters must identify one Source occurrence class"
            )
        represented_classes = set(source_classes)
        if promoted:
            exact_pml_classes = {
                source_kparam
                for source_kparam, occurrences in classes
                if any(
                    tuple(occurrence.reciprocal_vector_pml) == tuple(represented[1])
                    for occurrence in occurrences
                )
            }
            pml_classes = exact_pml_classes or {
                source_kparam
                for source_kparam, occurrences in classes
                if any(
                    _reciprocal_translation_key(occurrence.reciprocal_vector_pml)
                    == represented_pml
                    for occurrence in occurrences
                )
            }
            if len(pml_classes) != 1:
                raise ValueError(
                    "represented reciprocal vector must identify one Source mesh class"
                )
            represented_classes.update(pml_classes)
        for source_kparam, occurrences in classes:
            if source_kparam in represented_classes:
                continue
            try:
                exact_aliases = [
                    alias
                    for alias in aliases_for_class(row, source_kparam)
                    if int(alias.gid) == int(row.gid)
                    and _strict_integral_values(alias.source_kparam) == source_kparam
                ]
            except (IndexError, KeyError, TypeError, ValueError):
                continue
            if len(exact_aliases) != 1:
                continue
            alias = replace(
                exact_aliases[0],
                source_kparam=source_kparam,
                reciprocal_vector_pml=tuple(occurrences[0].reciprocal_vector_pml),
                source_occurrences=tuple(occurrences),
            )
            try:
                direction = np.asarray(alias.direction_matrix, dtype=float)
                if (
                    direction.ndim != 2
                    or direction.size == 0
                    or not np.all(np.isfinite(direction))
                    or int(np.linalg.matrix_rank(direction, tol=1e-10)) <= 0
                ):
                    continue
                heading = heading_identity_for_row(alias)
                if heading is None:
                    continue
                hash(heading)
            except (TypeError, ValueError, np.linalg.LinAlgError):
                continue
            if heading in seen_headings:
                continue
            candidates.append(
                OccurrenceAliasCandidate(
                    representative_gid=int(row.gid),
                    representative_source_kparam=represented_source_kparam,
                    candidate_source_kparam=source_kparam,
                    heading_identity=heading,
                    candidate_row=alias,
                )
            )
    heading_counts = Counter(candidate.heading_identity for candidate in candidates)
    return tuple(
        candidate
        for candidate in candidates
        if heading_counts[candidate.heading_identity] == 1
    )


@dataclass(frozen=True)
class _SelectedFamily:
    k_label: str
    source_kparam: tuple[int, int, int, int]
    reciprocal_vector_pml: tuple[Fraction, Fraction, Fraction]
    source_occurrences: tuple[DynamicKOccurrence, ...]


@dataclass(frozen=True)
class _KFamilyMatch:
    label: str
    kparam: tuple[int, int, int, int]
    kdim: int
    source_order: int
    representative_gid: int
    representative_vector: tuple[Fraction, Fraction, Fraction]


@dataclass(frozen=True)
class KvecStandardEmbedding:
    subgroup_number: int
    basis_values: tuple[int, ...]
    origin_values: tuple[int, int, int, int]


@dataclass(frozen=True)
class KvecStandardProvenance:
    """Source-order provenance retained from the k-vector scan."""

    sg: int
    representative_gid: int
    source_kparam: tuple[int, int, int, int]
    standard_kparam: tuple[int, int, int, int]
    selected_source_order: int
    selected_point_op: int
    selected_point_matrix: tuple[int, ...]
    selected_operation_record: OperationRecord | None
    source_star: tuple[tuple[Fraction, Fraction, Fraction], ...]
    standard_star: tuple[tuple[Fraction, Fraction, Fraction], ...]
    arm_permutation: tuple[int, ...]
    arm_conjugated: tuple[bool, ...]
    reciprocal_shifts: tuple[tuple[int, int, int], ...]
    source_embedding: KvecStandardEmbedding | None = None
    standard_embedding: KvecStandardEmbedding | None = None


def _key(vector: Sequence[Fraction]) -> tuple[tuple[int, int], ...]:
    return tuple(((Fraction(value) % 1).numerator, (Fraction(value) % 1).denominator) for value in vector[:3])


def _reciprocal_point_group_orbit_key(
    decoder: Any,
    sg: int,
    vector: Sequence[Fraction],
) -> tuple[tuple[int, int], ...]:
    """Return the canonical PML reciprocal orbit under the parent point group."""

    source = tuple(Fraction(value) for value in vector[:3])
    keys: list[tuple[tuple[int, int], ...]] = []
    seen_ops: set[int] = set()
    units = (
        (Fraction(1), Fraction(0), Fraction(0)),
        (Fraction(0), Fraction(1), Fraction(0)),
        (Fraction(0), Fraction(0), Fraction(1)),
    )
    for record in decoder.generate_space_group_records(int(sg)):
        point_op = int(record[4])
        if point_op in seen_ops:
            continue
        seen_ops.add(point_op)
        columns = [decoder.vrot_fraction(int(sg), point_op, unit) for unit in units]
        direct = tuple(
            tuple(columns[col][row] for col in range(3))
            for row in range(3)
        )
        inverse = fraction_matrix_inverse3(direct)
        reciprocal = tuple(
            sum(inverse[col][row] * source[col] for col in range(3)) % 1
            for row in range(3)
        )
        keys.append(_key(reciprocal))
    return min(keys) if keys else _key(source)


def _solve(rows: Sequence[Sequence[Fraction]], rhs: Sequence[Fraction], width: int) -> tuple[Fraction, ...] | None:
    matrix = [list(row[:width]) + [Fraction(value)] for row, value in zip(rows, rhs, strict=True)]
    pivot_row = 0
    pivots: list[int] = []
    for col in range(width):
        pivot = next((row for row in range(pivot_row, len(matrix)) if matrix[row][col]), None)
        if pivot is None:
            continue
        matrix[pivot_row], matrix[pivot] = matrix[pivot], matrix[pivot_row]
        divisor = matrix[pivot_row][col]
        matrix[pivot_row] = [value / divisor for value in matrix[pivot_row]]
        for row, values in enumerate(matrix):
            if row == pivot_row or not values[col]:
                continue
            factor = values[col]
            matrix[row] = [values[index] - factor * matrix[pivot_row][index] for index in range(width + 1)]
        pivots.append(col)
        pivot_row += 1
    if any(all(not row[col] for col in range(width)) and row[width] for row in matrix):
        return None
    if len(pivots) < width:
        return None
    solution = [Fraction(0)] * width
    for row, col in enumerate(pivots):
        solution[col] = matrix[row][width]
    return tuple(solution)


def _k_dimension(decoder: Any, sg: int, gid: int) -> int:
    lattice = int(decoder.iso.space["ispace_lattice"][sg - 1])
    slot = (lattice - 1) * 27 + int(decoder.iso.little["little_irr_k"][gid - 1]) - 1
    return int(decoder.iso.little["little_k_dim"][slot])


def _source_kparam_for_vector(
    decoder: Any,
    sg: int,
    gid: int,
    vector: Sequence[Fraction],
) -> tuple[int, int, int, int] | None:
    return _source_kparam_for_vector_values(
        decoder,
        int(sg),
        int(gid),
        tuple(Fraction(value) for value in vector[:3]),
    )


def _source_kparam_for_vector_values(
    decoder: Any,
    sg: int,
    gid: int,
    vector: tuple[Fraction, Fraction, Fraction],
) -> tuple[int, int, int, int] | None:
    dim = _k_dimension(decoder, sg, gid)
    if dim <= 0:
        vectors = decoder.little_k_vectors_by_gid(gid).vectors
        return (0, 0, 0, 1) if any(_key(item) == _key(vector) for item in vectors) else None
    kslot = int(decoder.iso.little["little_irr_k"][gid - 1])
    slot = (sg - 1) * 27 + kslot - 1
    count = int(decoder.iso.little["little_k_star_count"][slot])
    pointer = int(decoder.iso.little["little_k_star_ml_pointer"][slot])
    raw = decoder.iso.little["little_k_star"]
    target = tuple(Fraction(value) for value in vector[:3])
    for arm in range(count):
        record = tuple(int(value) for value in raw[16 * (pointer - 1 + arm):16 * (pointer + arm)])
        pieces = []
        for offset in range(0, 16, 4):
            den = record[offset + 3]
            pieces.append(
                (Fraction(0), Fraction(0), Fraction(0))
                if den == 0
                else tuple(Fraction(record[offset + axis], den) for axis in range(3))
            )
        rows = [[pieces[index + 1][axis] for index in range(dim)] for axis in range(3)]
        for sx in (-1, 0, 1):
            for sy in (-1, 0, 1):
                for sz in (-1, 0, 1):
                    shifted = (target[0] + sx, target[1] + sy, target[2] + sz)
                    solved = _solve(rows, [shifted[axis] - pieces[0][axis] for axis in range(3)], dim)
                    if solved is None:
                        continue
                    den = math.lcm(*(value.denominator for value in solved), 1)
                    padded = tuple(solved) + (Fraction(0),) * (3 - dim)
                    return tuple(int(value * den) for value in padded[:3]) + (den,)  # type: ignore[return-value]
    return None


def _little_k_vectors_by_gid_kparam(
    decoder: Any,
    gid: int,
    kparam: Sequence[Fraction | int],
) -> tuple[tuple[Fraction, Fraction, Fraction], ...]:
    values = tuple(kparam)
    return _little_k_vectors_by_gid_kparam_values(decoder, int(gid), values)


def _little_k_vectors_by_gid_kparam_values(
    decoder: Any,
    gid: int,
    values: tuple[Fraction | int, ...],
) -> tuple[tuple[Fraction, Fraction, Fraction], ...]:
    if len(values) == 4 and int(values[3]) != 0:
        params = tuple(Fraction(int(values[index]), int(values[3])) for index in range(3))
    else:
        params = tuple(Fraction(value) for value in values[:3])
    params += (Fraction(0),) * (3 - len(params))
    sg = int(decoder.iso.little["little_irr_space_group"][gid - 1])
    kslot = int(decoder.iso.little["little_irr_k"][gid - 1])
    slot = (sg - 1) * 27 + kslot - 1
    count = int(decoder.iso.little["little_k_star_count"][slot])
    pointer = int(decoder.iso.little["little_k_star_ml_pointer"][slot])
    raw = decoder.iso.little["little_k_star"]
    vectors: list[tuple[Fraction, Fraction, Fraction]] = []
    for arm in range(count):
        record = tuple(int(value) for value in raw[16 * (pointer - 1 + arm):16 * (pointer + arm)])
        pieces: list[tuple[Fraction, Fraction, Fraction]] = []
        for offset in range(0, 16, 4):
            denominator = int(record[offset + 3])
            pieces.append(
                (Fraction(0), Fraction(0), Fraction(0))
                if denominator == 0
                else tuple(Fraction(record[offset + axis], denominator) for axis in range(3))
            )
        vectors.append(
            tuple(
                pieces[0][axis]
                + sum(params[index] * pieces[index + 1][axis] for index in range(3))
                for axis in range(3)
            )
        )
    return tuple(vectors)


@lru_cache(maxsize=4096)
def _first_little_k_affine_by_gid_cached(
    decoder: Any,
    gid: int,
) -> tuple[tuple[Fraction, Fraction, Fraction], ...] | None:
    """Return the exact affine pieces of the first Source little-k star arm."""

    sg = int(decoder.iso.little["little_irr_space_group"][gid - 1])
    kslot = int(decoder.iso.little["little_irr_k"][gid - 1])
    slot = (sg - 1) * 27 + kslot - 1
    count = int(decoder.iso.little["little_k_star_count"][slot])
    if count <= 0:
        return None
    pointer = int(decoder.iso.little["little_k_star_ml_pointer"][slot])
    raw = decoder.iso.little["little_k_star"]
    record = tuple(int(value) for value in raw[16 * (pointer - 1):16 * pointer])
    pieces: list[tuple[Fraction, Fraction, Fraction]] = []
    for offset in range(0, 16, 4):
        denominator = int(record[offset + 3])
        pieces.append(
            (Fraction(0), Fraction(0), Fraction(0))
            if denominator == 0
            else tuple(Fraction(record[offset + axis], denominator) for axis in range(3))
        )
    return tuple(pieces)


def _evaluate_first_little_k_affine(
    pieces: tuple[tuple[Fraction, Fraction, Fraction], ...],
    params: Sequence[Fraction],
) -> tuple[Fraction, Fraction, Fraction]:
    return tuple(
        pieces[0][axis]
        + sum(params[index] * pieces[index + 1][axis] for index in range(len(params)))
        for axis in range(3)
    )  # type: ignore[return-value]


def _first_little_k_vector_by_gid_kparam(
    decoder: Any,
    gid: int,
    values: tuple[Fraction | int, ...],
) -> tuple[Fraction, Fraction, Fraction] | None:
    """B: evaluate only the first Source star arm used by kvec_standard_."""

    pieces = _first_little_k_affine_by_gid_cached(decoder, int(gid))
    if pieces is None:
        return None
    if len(values) == 4 and int(values[3]) != 0:
        params = tuple(Fraction(int(values[index]), int(values[3])) for index in range(3))
    else:
        params = tuple(value if isinstance(value, Fraction) else Fraction(value) for value in values[:3])
    params += (Fraction(0),) * (3 - len(params))
    return _evaluate_first_little_k_affine(pieces, params)


def _rotate_kvector_by_point_matrix(
    decoder: Any,
    sg: int,
    vector: Sequence[Fraction],
    point: Sequence[int],
) -> tuple[Fraction, Fraction, Fraction]:
    return tuple(
        value % 1
        for value in _rotate_kvector_by_point_matrix_raw(decoder, sg, vector, point)
    )  # type: ignore[return-value]


def _rotate_kvector_by_point_matrix_raw(
    decoder: Any,
    sg: int,
    vector: Sequence[Fraction],
    point: Sequence[int],
) -> tuple[Fraction, Fraction, Fraction]:
    lattice = int(decoder.iso.space["ispace_lattice"][int(sg) - 1])
    record = tuple(
        int(value)
        for value in decoder.iso.space["lattice_ml"][(lattice - 1) * 36:lattice * 36]
    )
    denominators = tuple(
        int(value)
        for value in decoder.iso.space["lattice_ml_denom"][(lattice - 1) * 2:lattice * 2]
    )
    denominator = int(denominators[1])
    if denominator == 0:
        raise ValueError(f"zero reciprocal lattice_ml denominator for lattice {lattice}")
    rotated = decoder._vmlt_fraction(  # noqa: SLF001
        record[18:27],
        decoder._vmlt_fraction(  # noqa: SLF001
            tuple(int(value) for value in point),
            decoder._vmlt_fraction(record[27:36], tuple(Fraction(value) for value in vector[:3])),  # noqa: SLF001
        ),
    )
    return tuple(value / denominator for value in rotated)  # type: ignore[return-value]


@lru_cache(maxsize=256)
def _lattice_point_operations_for_sg(
    decoder: Any,
    sg: int,
) -> tuple[tuple[int, tuple[int, ...]], ...]:
    lattice = int(decoder.iso.space["ispace_lattice"][int(sg) - 1])
    representatives = [
        candidate
        for candidate, candidate_lattice in enumerate(decoder.iso.space["ispace_lattice"], start=1)
        if int(candidate_lattice) == lattice
    ]
    if not representatives:
        raise ValueError(f"no lattice representative for SG{sg} lattice={lattice}")
    representative_sg = representatives[-1]
    return tuple(
        (
            int(record[4]),
            tuple(
                int(value)
                for value in decoder.iso.space["ipoint_op"][9 * (int(record[4]) - 1):9 * int(record[4])]
            ),
        )
        for record in decoder.generate_space_group_records(representative_sg)
    )


@lru_cache(maxsize=256)
def _kvec_standard_point_operations_for_sg(
    decoder: Any,
    sg: int,
) -> tuple[tuple[int, tuple[int, ...]], ...]:
    """Return Source parent operations in standardization-adapter scan order."""

    operations: list[int] = []
    for record in decoder.generate_space_group_records(int(sg)):
        point_op = int(record[4])
        if point_op not in operations:
            operations.append(point_op)
    if not operations:
        raise ValueError(f"SG{sg} has no generated point operations")
    inversion = 25 if int(operations[0]) == 1 else 61
    if inversion not in operations:
        original = tuple(operations)
        products = decoder.iso.space["ipoint_op_mlt"]
        for point_op in original:
            mate = int(products[(inversion - 1) * 72 + point_op - 1])
            if mate not in operations:
                operations.append(mate)
    return tuple(
        (
            point_op,
            tuple(
                int(value)
                for value in decoder.iso.space["ipoint_op"][(point_op - 1) * 9:point_op * 9]
            ),
        )
        for point_op in operations
    )


@lru_cache(maxsize=256)
def _lattice_point_matrices_for_sg(decoder: Any, sg: int) -> tuple[tuple[int, ...], ...]:
    return tuple(matrix for _point_op, matrix in _lattice_point_operations_for_sg(decoder, sg))


def _kvector_pml_to_cinter(
    decoder: Any,
    sg: int,
    vector: Sequence[Fraction],
) -> tuple[Fraction, Fraction, Fraction]:
    raw = tuple(vector[:3])
    values = (
        raw
        if all(isinstance(value, Fraction) for value in raw)
        else tuple(Fraction(value) for value in raw)
    )
    return _kvector_pml_to_cinter_cached(decoder, int(sg), values)


@lru_cache(maxsize=4096)
def _kvector_pml_to_cinter_cached(
    decoder: Any,
    sg: int,
    values: tuple[Fraction, ...],
) -> tuple[Fraction, Fraction, Fraction]:
    """B: reuse the exact Source frame conversion for repeated k vectors."""

    matrix = _cinter_to_pml_matrix(decoder, int(sg))
    return tuple(
        sum(matrix[row][col] * values[col] for col in range(3))
        for row in range(3)
    )  # type: ignore[return-value]


@lru_cache(maxsize=4096)
def _first_little_k_cinter_affine_cached(
    decoder: Any,
    sg: int,
    gid: int,
) -> tuple[tuple[Fraction, Fraction, Fraction], ...] | None:
    """Transport the first Source arm's affine pieces once into cinter."""

    pieces = _first_little_k_affine_by_gid_cached(decoder, int(gid))
    if pieces is None:
        return None
    return tuple(
        _kvector_pml_to_cinter_cached(decoder, int(sg), tuple(piece))
        for piece in pieces
    )


@lru_cache(maxsize=256)
def _cinter_to_pml_matrix(decoder: Any, sg: int) -> tuple[tuple[Fraction, ...], ...]:
    return tuple(
        tuple(Fraction(value) for value in row)
        for row in decoder.cinter_to_pml_matrix(int(sg))
    )


def _kvec_parameter_periods(decoder: Any, representative_gid: int) -> tuple[int, ...]:
    gid = int(representative_gid)
    sg = int(decoder.iso.little["little_irr_space_group"][gid - 1])
    lattice = int(decoder.iso.space["ispace_lattice"][sg - 1])
    kslot = int(decoder.iso.little["little_irr_k"][gid - 1])
    kdim = _k_dimension(decoder, sg, gid)
    slot = (lattice - 1) * 27 + kslot - 1
    record = tuple(int(value) for value in decoder.iso.little["little_k"][slot * 16:(slot + 1) * 16])
    periods: list[int] = []
    for index in range(kdim):
        start = 4 * (index + 1)
        denominator = abs(int(record[start + 3]))
        if denominator == 0:
            periods.append(1)
            continue
        period = 1
        for numerator in record[start:start + 3]:
            period = math.lcm(period, denominator // math.gcd(denominator, abs(int(numerator))))
        periods.append(max(1, int(period)))
    return tuple(periods)


@lru_cache(maxsize=8192)
def _kvec_parameter_solutions(
    decoder: Any,
    representative_gid: int,
    vector: tuple[Fraction, Fraction, Fraction],
) -> tuple[tuple[Fraction, ...], ...]:
    gid = int(representative_gid)
    sg = int(decoder.iso.little["little_irr_space_group"][gid - 1])
    kdim = _k_dimension(decoder, sg, gid)
    if kdim <= 0:
        return (tuple(),)
    target = tuple(Fraction(value) for value in vector[:3])
    kslot = int(decoder.iso.little["little_irr_k"][gid - 1])
    slot = (sg - 1) * 27 + kslot - 1
    count = int(decoder.iso.little["little_k_star_count"][slot])
    pointer = int(decoder.iso.little["little_k_star_ml_pointer"][slot])
    raw = decoder.iso.little["little_k_star"]
    out: list[tuple[Fraction, ...]] = []
    seen: set[tuple[Fraction, ...]] = set()
    for arm in range(count):
        record = tuple(int(value) for value in raw[16 * (pointer - 1 + arm):16 * (pointer + arm)])
        pieces: list[tuple[Fraction, Fraction, Fraction]] = []
        for offset in range(0, 16, 4):
            denominator = int(record[offset + 3])
            pieces.append(
                (Fraction(0), Fraction(0), Fraction(0))
                if denominator == 0
                else tuple(Fraction(record[offset + axis], denominator) for axis in range(3))
            )
        rows = [[pieces[index + 1][axis] for index in range(kdim)] for axis in range(3)]
        for reciprocal_shift in product((-1, 0, 1), repeat=3):
            shifted = tuple(target[axis] + reciprocal_shift[axis] for axis in range(3))
            solution = _solve(
                rows,
                [shifted[axis] - pieces[0][axis] for axis in range(3)],
                kdim,
            )
            if solution is None:
                continue
            rebuilt = tuple(
                pieces[0][axis]
                + sum(solution[index] * pieces[index + 1][axis] for index in range(kdim))
                for axis in range(3)
            )
            if _key(rebuilt) != _key(target):
                continue
            canonical = tuple(Fraction(value) % 1 for value in solution)
            if canonical not in seen:
                seen.add(canonical)
                out.append(canonical)
    return tuple(out)


def _kparam_from_fractions(vector: Sequence[Fraction]) -> tuple[int, int, int, int]:
    values = tuple(Fraction(value) for value in vector[:3])
    denominator = math.lcm(*(value.denominator for value in values), 1)
    return tuple(int(value * denominator) for value in values) + (denominator,)  # type: ignore[return-value]


@lru_cache(maxsize=8192)
def _kvec_standard_selection(
    decoder: Any,
    sg: int,
    representative_gid: int,
    kparam: tuple[int, int, int, int],
) -> tuple[tuple[int, int, int, int], int, int, tuple[int, ...]]:
    """Return the faithful parameter result and selected lattice operation.

    Request-scoped: ``decoder`` identity is part of the cache key (see
    ``ModeDataDecoder``). Callers such as ``_kvec_standard_kparam`` previously
    re-entered the full holohedry search on every site/family touch; memoize
    the pure ``(sg, gid, kparam)`` selection within that request.
    """

    point_operations = _kvec_standard_point_operations_for_sg(decoder, int(sg))
    identity = (1, 0, 0, 0, 1, 0, 0, 0, 1)
    try:
        selected_source_order, (selected_point_op, selected_point_matrix) = next(
            (index, item)
            for index, item in enumerate(point_operations)
            if item[1] == identity
        )
    except StopIteration as exc:
        raise ValueError(f"SG{sg} lattice holohedry has no identity operation") from exc
    original_vector = _first_little_k_vector_by_gid_kparam(
        decoder,
        int(representative_gid),
        tuple(kparam),
    )
    if original_vector is None:
        return (
            tuple(int(value) for value in kparam),
            int(selected_source_order),
            int(selected_point_op),
            selected_point_matrix,
        )
    current_param = tuple(int(value) for value in kparam)
    current_key = tuple(
        component
        for value in _kvector_pml_to_cinter(decoder, int(sg), original_vector)
        for component in (abs(value), int(value < 0))
    )
    kdim = _k_dimension(decoder, int(sg), int(representative_gid))
    periods = _kvec_parameter_periods(decoder, int(representative_gid))
    cinter_affine = _first_little_k_cinter_affine_cached(
        decoder,
        int(sg),
        int(representative_gid),
    )
    if cinter_affine is None:
        return current_param, selected_source_order, selected_point_op, selected_point_matrix
    shifted_deltas = tuple(
        (
            parameter_shift,
            tuple(
                sum(
                    int(parameter_shift[index])
                    * int(periods[index])
                    * cinter_affine[index + 1][axis]
                    for index in range(kdim)
                )
                for axis in range(3)
            ),
        )
        for parameter_shift in product((-1, 0, 1), repeat=kdim)
    )
    for source_order, (point_op, point_matrix) in enumerate(point_operations):
        rotated = _rotate_kvector_by_point_matrix(decoder, int(sg), original_vector, point_matrix)
        for solved in _kvec_parameter_solutions(decoder, int(representative_gid), rotated):
            solved_cinter = _evaluate_first_little_k_affine(
                cinter_affine,
                tuple(Fraction(solved[index]) for index in range(kdim)),
            )
            for parameter_shift, shifted_delta in shifted_deltas:
                candidate_cinter = tuple(
                    solved_cinter[axis] + shifted_delta[axis]
                    for axis in range(3)
                )
                candidate_key = tuple(
                    component
                    for value in candidate_cinter
                    for component in (abs(value), int(value < 0))
                )
                if candidate_key < current_key:
                    signed = tuple(
                        Fraction(solved[index])
                        + int(parameter_shift[index]) * int(periods[index])
                        for index in range(kdim)
                    )
                    current_param = _kparam_from_fractions(
                        tuple(signed) + (Fraction(0),) * (3 - kdim)
                    )
                    current_key = candidate_key
                    selected_source_order = int(source_order)
                    selected_point_op = int(point_op)
                    selected_point_matrix = tuple(int(value) for value in point_matrix)
    return current_param, selected_source_order, selected_point_op, selected_point_matrix


def _kvec_standard_arm_map(
    decoder: Any,
    sg: int,
    source_star: tuple[tuple[Fraction, Fraction, Fraction], ...],
    standard_star: tuple[tuple[Fraction, Fraction, Fraction], ...],
    point_matrix: tuple[int, ...],
) -> tuple[tuple[int, ...], tuple[bool, ...], tuple[tuple[int, int, int], ...]]:
    if len(source_star) != len(standard_star):
        raise ValueError(
            f"kvec_standard star size mismatch: source={len(source_star)}, "
            f"standard={len(standard_star)}"
        )
    candidates: list[list[tuple[int, bool, tuple[int, int, int]]]] = []
    for source_vector in source_star:
        rotated = _rotate_kvector_by_point_matrix_raw(decoder, int(sg), source_vector, point_matrix)
        source_candidates: list[tuple[int, bool, tuple[int, int, int]]] = []
        for is_conjugated in (False, True):
            for standard_index, standard_vector in enumerate(standard_star):
                sign = -1 if is_conjugated else 1
                delta = tuple(
                    rotated[axis] - sign * standard_vector[axis]
                    for axis in range(3)
                )
                if all(value.denominator == 1 for value in delta):
                    source_candidates.append(
                        (
                            standard_index,
                            is_conjugated,
                            tuple(int(value) for value in delta),  # type: ignore[arg-type]
                        )
                    )
        if not source_candidates:
            raise ValueError(
                f"kvec standard arm {len(candidates)} has no match: "
                f"SG{sg} source={source_vector} rotated={rotated}"
            )
        candidates.append(source_candidates)

    dead: set[tuple[int, int]] = set()

    def assign(
        source_index: int,
        used_mask: int,
    ) -> tuple[tuple[int, bool, tuple[int, int, int]], ...] | None:
        if source_index == len(candidates):
            return ()
        state = (source_index, used_mask)
        if state in dead:
            return None
        for candidate in candidates[source_index]:
            standard_index = int(candidate[0])
            bit = 1 << standard_index
            if used_mask & bit:
                continue
            suffix = assign(source_index + 1, used_mask | bit)
            if suffix is not None:
                return (candidate,) + suffix
        dead.add(state)
        return None

    assignment = assign(0, 0)
    if assignment is None:
        raise ValueError(
            f"kvec_standard has no bijective arm assignment: SG{sg} "
            f"source_arms={len(source_star)} standard_arms={len(standard_star)}"
        )
    return (
        tuple(int(candidate[0]) for candidate in assignment),
        tuple(bool(candidate[1]) for candidate in assignment),
        tuple(candidate[2] for candidate in assignment),
    )


@lru_cache(maxsize=256)
def kvec_standard_provenance(
    decoder: Any,
    sg: int,
    representative_gid: int,
    kparam: tuple[int, int, int, int],
    source_embedding: KvecStandardEmbedding | None = None,
    standard_embedding: KvecStandardEmbedding | None = None,
) -> KvecStandardProvenance:
    """Expose Source-only standardization provenance to APOSTRUCT callers."""

    values = tuple(int(value) for value in kparam)
    standard_param, source_order, point_op, point_matrix = _kvec_standard_selection(
        decoder, int(sg), int(representative_gid), values
    )
    source_star = _little_k_vectors_by_gid_kparam(decoder, int(representative_gid), values)
    standard_star = _little_k_vectors_by_gid_kparam(decoder, int(representative_gid), standard_param)
    arm_permutation, arm_conjugated, reciprocal_shifts = _kvec_standard_arm_map(
        decoder, int(sg), source_star, standard_star, point_matrix
    )
    selected_operation_record = next(
        (
            tuple(int(value) for value in record)
            for record in decoder.generate_space_group_records(int(sg))
            if int(record[4]) == int(point_op)
        ),
        None,
    )
    return KvecStandardProvenance(
        sg=int(sg),
        representative_gid=int(representative_gid),
        source_kparam=values,
        standard_kparam=standard_param,
        selected_source_order=int(source_order),
        selected_point_op=int(point_op),
        selected_point_matrix=point_matrix,
        selected_operation_record=selected_operation_record,  # type: ignore[arg-type]
        source_star=source_star,
        standard_star=standard_star,
        arm_permutation=arm_permutation,
        arm_conjugated=arm_conjugated,
        reciprocal_shifts=reciprocal_shifts,
        source_embedding=source_embedding,
        standard_embedding=standard_embedding,
    )




def _kvec_standard_kparam(
    decoder: Any,
    sg: int,
    representative_gid: int,
    kparam: tuple[int, int, int, int],
) -> tuple[int, int, int, int]:
    return _kvec_standard_selection(decoder, sg, representative_gid, kparam)[0]


@lru_cache(maxsize=256)
def _family_rows(decoder: Any, sg: int) -> tuple[tuple[str, int, int], ...]:
    out: list[tuple[str, int, int]] = []
    seen: set[str] = set()
    for gid, row_sg in enumerate(decoder.iso.little["little_irr_space_group"], start=1):
        if int(row_sg) != sg:
            continue
        kslot = int(decoder.iso.little["little_irr_k"][gid - 1])
        lattice = int(decoder.iso.space["ispace_lattice"][sg - 1])
        label = str(decoder.iso.little["little_k_label"][(lattice - 1) * 27 + kslot - 1]).strip()
        if label in seen:
            continue
        seen.add(label)
        out.append((label, gid, _k_dimension(decoder, sg, gid)))
    return tuple(out)


@lru_cache(maxsize=8192)
def _k_family_candidates_for_standard_vector(
    decoder: Any,
    sg: int,
    target: tuple[Fraction, Fraction, Fraction],
) -> tuple[_KFamilyMatch, ...]:
    matches: list[_KFamilyMatch] = []
    for source_order, (label, representative_gid, kdim) in enumerate(
        _family_rows(decoder, int(sg)),
        start=1,
    ):
        solved = _source_kparam_for_vector(
            decoder,
            int(sg),
            int(representative_gid),
            target,
        )
        if solved is None:
            continue
        kparam = (
            (0, 0, 0, 1)
            if int(kdim) <= 0
            else _kvec_standard_kparam(
                decoder,
                int(sg),
                int(representative_gid),
                tuple(int(value) for value in solved),
            )
        )
        arms = (
            tuple(decoder.little_k_vectors_by_gid(int(representative_gid)).vectors)
            if int(kdim) <= 0
            else _little_k_vectors_by_gid_kparam(
                decoder,
                int(representative_gid),
                kparam,
            )
        )
        if not any(_key(arm) == _key(target) for arm in arms):
            continue
        matches.append(
            _KFamilyMatch(
                label=str(label),
                kparam=tuple(int(value) for value in kparam),
                kdim=int(kdim),
                source_order=int(source_order),
                representative_gid=int(representative_gid),
                representative_vector=target,
            )
        )
    return tuple(matches)


def _k_family_matches_for_vector(
    decoder: Any,
    sg: int,
    vector: Sequence[Fraction],
    preferred_kparams: dict[str, Sequence[int] | Sequence[Sequence[int]]] | None = None,
) -> tuple[_KFamilyMatch, ...]:
    matches: dict[tuple[str, int], _KFamilyMatch] = {}
    seen: set[tuple[tuple[int, int], ...]] = set()
    rotated_vectors = [tuple(Fraction(value) for value in vector[:3])]
    rotated_vectors.extend(
        _rotate_kvector_by_point_matrix(decoder, int(sg), vector, point_matrix)
        for point_matrix in _lattice_point_matrices_for_sg(decoder, int(sg))
    )
    for rotated in rotated_vectors:
        key = _key(rotated)
        if key in seen:
            continue
        seen.add(key)
        for match in _k_family_candidates_for_standard_vector(
            decoder, int(sg), tuple(Fraction(value) % 1 for value in rotated[:3])
        ):
            matches.setdefault((str(match.label), int(match.representative_gid)), match)
    if not matches:
        return ()
    selected_class = min(
        (
            not _matches_preferred_kparam(
                preferred_kparams,
                match.label,
                match.kparam,
            ),
            int(match.kdim),
        )
        for match in matches.values()
    )
    target = _key(vector)
    selected_matches: list[_KFamilyMatch] = []
    for selected in sorted(matches.values(), key=lambda item: item.source_order):
        if (
            not _matches_preferred_kparam(
                preferred_kparams,
                selected.label,
                selected.kparam,
            ),
            int(selected.kdim),
        ) != selected_class:
            continue
        if int(selected.kdim) > 0:
            occurrence_kparam = _source_kparam_for_vector(
                decoder,
                int(sg),
                int(selected.representative_gid),
                tuple(Fraction(value) for value in vector[:3]),
            )
            if occurrence_kparam is not None:
                standardized = _kvec_standard_kparam(
                    decoder,
                    int(sg),
                    int(selected.representative_gid),
                    tuple(int(value) for value in occurrence_kparam),
                )
                standard_star = _little_k_vectors_by_gid_kparam(
                    decoder,
                    int(selected.representative_gid),
                    standardized,
                )
                if any(
                    _key(arm) == target
                    or _key(tuple(-value for value in arm)) == target
                    for arm in standard_star
                ):
                    selected = replace(selected, kparam=standardized)
        selected_matches.append(selected)
    return tuple(selected_matches)


def _lattice_k_orbit_key(
    decoder: Any,
    sg: int,
    vector: Sequence[Fraction],
) -> tuple[tuple[int, int], ...]:
    keys = [
        _key(
            _rotate_kvector_by_point_matrix(
                decoder,
                int(sg),
                vector,
                point_matrix,
            )
        )
        for point_matrix in _lattice_point_matrices_for_sg(decoder, int(sg))
    ]
    return min(keys) if keys else _key(vector)


def _selected_family_records(
    decoder: Any,
    sg: int,
    basis: Sequence[int],
    preferred_kparams: dict[str, Sequence[int] | Sequence[Sequence[int]]] | None = None,
) -> tuple[_SelectedFamily, ...]:
    family_orbit = tuple[str, int, tuple[tuple[int, int], ...]]
    selected: dict[family_orbit, tuple[int, _KFamilyMatch]] = {}
    occurrences: dict[family_orbit, list[DynamicKOccurrence]] = {}

    for mesh_order, vector in enumerate(_ordered_reciprocal_mesh(basis)):
        families = _k_family_matches_for_vector(
            decoder,
            int(sg),
            vector,
            preferred_kparams,
        )
        if not families:
            continue
        orbit_key = _lattice_k_orbit_key(decoder, int(sg), vector)
        for family in families:
            key = (str(family.label), int(family.representative_gid), orbit_key)
            occurrences.setdefault(key, []).append(
                DynamicKOccurrence(
                    mesh_order=mesh_order,
                    k_label=str(family.label),
                    source_kparam=tuple(int(value) for value in family.kparam),
                    reciprocal_vector_pml=tuple(Fraction(value) for value in vector),
                    inversion_bucket=orbit_key,
                )
            )
            selected.setdefault(key, (mesh_order, family))
    families: list[_SelectedFamily] = []
    for key, (_mesh_order, family) in sorted(
        selected.items(),
        key=lambda item: (item[1][1].source_order, item[1][0]),
    ):
        by_kparam: dict[tuple[int, int, int, int], list[DynamicKOccurrence]] = {}
        for occurrence in occurrences[key]:
            by_kparam.setdefault(occurrence.source_kparam, []).append(occurrence)
        selected_kparam = tuple(int(value) for value in family.kparam)
        ordered_kparams = [selected_kparam] if selected_kparam in by_kparam else []
        ordered_kparams.extend(
            kparam for kparam in by_kparam if kparam != selected_kparam
        )
        for source_kparam in ordered_kparams:
            source_occurrences = tuple(by_kparam[source_kparam])
            families.append(
                _SelectedFamily(
                    k_label=str(family.label),
                    source_kparam=source_kparam,
                    reciprocal_vector_pml=(
                        tuple(family.representative_vector)
                        if source_kparam == selected_kparam
                        else tuple(source_occurrences[0].reciprocal_vector_pml)
                    ),
                    source_occurrences=source_occurrences,
                )
            )
    return tuple(families)


def _display_direction_family_records(
    decoder: Any,
    sg: int,
    basis: Sequence[int],
) -> tuple[_SelectedFamily, ...]:
    """Return the one-representative-per-Bravais-star DISPLAY DIRECTION mesh."""

    selected: dict[
        tuple[tuple[int, int], ...],
        tuple[int, _KFamilyMatch, DynamicKOccurrence],
    ] = {}
    for mesh_order, vector in enumerate(_ordered_reciprocal_mesh(basis)):
        matches = _k_family_matches_for_vector(decoder, int(sg), vector)
        if not matches:
            continue
        family = min(matches, key=lambda item: item.source_order)
        orbit_key = _lattice_k_orbit_key(decoder, int(sg), vector)
        occurrence = DynamicKOccurrence(
            mesh_order=mesh_order,
            k_label=str(family.label),
            source_kparam=tuple(int(value) for value in family.kparam),
            reciprocal_vector_pml=tuple(Fraction(value) for value in vector),
            inversion_bucket=orbit_key,
        )
        selected.setdefault(orbit_key, (mesh_order, family, occurrence))

    return tuple(
        _SelectedFamily(
            k_label=str(family.label),
            source_kparam=tuple(int(value) for value in family.kparam),
            reciprocal_vector_pml=(
                tuple(decoder.little_k_vectors_by_gid(int(family.representative_gid)).vectors[0])
                if int(family.kdim) <= 0
                else _first_little_k_vector_by_gid_kparam(
                    decoder,
                    int(family.representative_gid),
                    tuple(int(value) for value in family.kparam),
                )
                or occurrence.reciprocal_vector_pml
            ),
            source_occurrences=(occurrence,),
        )
        for _mesh_order, family, occurrence in sorted(
            selected.values(),
            key=lambda item: (item[1].source_order, item[0]),
        )
    )


def _matches_preferred_kparam(
    preferred_kparams: dict[str, Sequence[int] | Sequence[Sequence[int]]] | None,
    label: str,
    kparam: tuple[int, int, int, int],
) -> bool:
    raw = (preferred_kparams or {}).get(label)
    if raw is None:
        return False
    values = list(raw)
    targets = [values] if len(values) == 4 and all(not isinstance(value, Sequence) for value in values) else values
    candidate = tuple(Fraction(kparam[index], kparam[3]) % 1 for index in range(3))
    return any(
        isinstance(target, Sequence)
        and len(target) >= 4
        and candidate
        == tuple(Fraction(int(target[index]), int(target[3])) % 1 for index in range(3))
        for target in targets
    )


def _subduction_rows_for_families(
    decoder: Any,
    *,
    sg: int,
    basis_values: tuple[int, ...],
    operations: Sequence[OperationRecord],
    families: Sequence[_SelectedFamily],
    preferred_kparams: dict[str, Sequence[int] | Sequence[Sequence[int]]] | None = None,
    irrep_source: Any | None = None,
) -> tuple[DynamicSubductionRow, ...]:
    identity_op = int(decoder.generate_space_group_records(sg)[0][4])
    constrained_operations = tuple(operations) + tuple(
        (basis_values[3 * row], basis_values[3 * row + 1], basis_values[3 * row + 2], 1, identity_op)
        for row in range(3)
    )
    rows: list[DynamicSubductionRow] = []
    seen_rows: dict[tuple[int, tuple[tuple[int, int], ...]], int] = {}
    row_preferred: list[bool] = []
    for family in families:
        k_label = family.k_label
        discovered_kparam = family.source_kparam
        reciprocal_vector_pml = family.reciprocal_vector_pml
        source_kparam = discovered_kparam
        dim = _k_dimension(
            decoder,
            sg,
            next(
                gid
                for gid, row_sg in enumerate(decoder.iso.little["little_irr_space_group"], start=1)
                if int(row_sg) == sg
                and str(
                    decoder.iso.little["little_k_label"][
                        (int(decoder.iso.space["ispace_lattice"][sg - 1]) - 1) * 27
                        + int(decoder.iso.little["little_irr_k"][gid - 1])
                        - 1
                    ]
                ).strip()
                == k_label
            ),
        )
        case_params = tuple(Fraction(source_kparam[index], source_kparam[3]) for index in range(dim))
        for gid, row_sg in enumerate(decoder.iso.little["little_irr_space_group"], start=1):
            if int(row_sg) != sg:
                continue
            kslot = int(decoder.iso.little["little_irr_k"][gid - 1])
            lattice = int(decoder.iso.space["ispace_lattice"][sg - 1])
            if str(decoder.iso.little["little_k_label"][(lattice - 1) * 27 + kslot - 1]).strip() != k_label:
                continue
            if irrep_source is None:
                case = Case(
                    sg=sg,
                    wyckoff="a",
                    k_label=k_label,
                    params=(1, 1, 1, 90, 90, 90),
                    k_params=case_params,
                )
                matrices = [decoder.get_irreps_matrix_for_case(gid, record, case) for record in constrained_operations]
            else:
                matrices = [
                    irrep_source.get_irreps_matrix_for_record(gid, record, source_kparam)
                    for record in constrained_operations
                ]
            reduced = decoder.xrowop2_like(decoder.subgroup_constraint_matrix_from_irrep_matrices(matrices))
            count, raw_direction = decoder.subgroup_direction_from_reduced_matrix(reduced)
            if count <= 0:
                continue
            full_dim = int(decoder.little_record_by_gid(gid).full_dim)
            direction = tuple(
                tuple(float(raw_direction[free, coordinate]) for free in range(count))
                for coordinate in range(full_dim)
            )
            if direction is None:
                continue
            row_key = (
                int(gid),
                _reciprocal_point_group_orbit_key(decoder, sg, reciprocal_vector_pml),
            )
            candidate = DynamicSubductionRow(
                gid=gid,
                irrep_label=str(decoder.iso.little["little_irr_full_label"][gid - 1]).strip(),
                k_label=k_label,
                source_kparam=source_kparam,
                reciprocal_vector_pml=reciprocal_vector_pml,
                direction_matrix=direction,
                source_occurrences=family.source_occurrences,
            )
            candidate_preferred = _matches_preferred_kparam(preferred_kparams, k_label, source_kparam)
            if row_key in seen_rows:
                index = seen_rows[row_key]
                if candidate_preferred and not row_preferred[index]:
                    rows[index] = replace(
                        candidate,
                        source_occurrences=rows[index].source_occurrences + candidate.source_occurrences,
                    )
                    row_preferred[index] = True
                else:
                    rows[index] = replace(
                        rows[index],
                        source_occurrences=rows[index].source_occurrences + candidate.source_occurrences,
                    )
                continue
            seen_rows[row_key] = len(rows)
            rows.append(candidate)
            row_preferred.append(candidate_preferred)
    return tuple(rows)


def dynamic_subduction_rows(
    decoder: Any,
    *,
    sg: int,
    basis: Sequence[int],
    operations: Sequence[OperationRecord],
    preferred_kparams: dict[str, Sequence[int] | Sequence[Sequence[int]]] | None = None,
    irrep_source: Any | None = None,
) -> tuple[DynamicSubductionRow, ...]:
    """Return ordinary dynamic invariant directions in Source print order."""

    basis_values = tuple(int(value) for value in basis)
    return _subduction_rows_for_families(
        decoder,
        sg=sg,
        basis_values=basis_values,
        operations=operations,
        families=_selected_family_records(
            decoder,
            sg,
            basis_values,
            preferred_kparams,
        ),
        preferred_kparams=preferred_kparams,
        irrep_source=irrep_source,
    )


def display_direction_subduction_rows(
    decoder: Any,
    *,
    sg: int,
    basis: Sequence[int],
    operations: Sequence[OperationRecord],
    irrep_source: Any | None = None,
) -> tuple[DynamicSubductionRow, ...]:
    """Return binary-faithful DISPLAY DIRECTION rows in Source print order."""

    basis_values = tuple(int(value) for value in basis)
    return _subduction_rows_for_families(
        decoder,
        sg=sg,
        basis_values=basis_values,
        operations=operations,
        families=_display_direction_family_records(decoder, sg, basis_values),
        irrep_source=irrep_source,
    )


def dynamic_subduction_occurrence_alias_rows(
    decoder: Any,
    *,
    sg: int,
    basis: Sequence[int],
    operations: Sequence[OperationRecord],
    row: DynamicSubductionRow,
    irrep_source: Any,
    requested_source_kparam: Sequence[int] | None = None,
) -> tuple[DynamicSubductionRow, ...]:
    """Recover invariant raw Cases hidden by DISPLAY DIRECTION orbit dedupe."""

    basis_values = tuple(int(value) for value in basis)
    identity_op = int(decoder.generate_space_group_records(int(sg))[0][4])
    constrained_operations = tuple(operations) + tuple(
        (
            basis_values[3 * axis],
            basis_values[3 * axis + 1],
            basis_values[3 * axis + 2],
            1,
            identity_op,
        )
        for axis in range(3)
    )

    def parameter_key(kparam: Sequence[int]) -> tuple[Fraction, Fraction, Fraction]:
        """Keep raw Source occurrences distinct even when congruent modulo 1."""

        return tuple(
            Fraction(int(kparam[index]), int(kparam[3]))
            for index in range(3)
        )  # type: ignore[return-value]

    seen = {parameter_key(row.source_kparam)}
    aliases: list[DynamicSubductionRow] = []

    def append_alias(
        source_kparam: Sequence[int],
        occurrence: DynamicKOccurrence,
    ) -> None:
        source_kparam = tuple(int(value) for value in source_kparam)
        key = parameter_key(source_kparam)
        if key in seen:
            return
        seen.add(key)
        matrices = [
            irrep_source.get_irreps_matrix_for_record(
                int(row.gid),
                record,
                source_kparam,
            )
            for record in constrained_operations
        ]
        reduced = decoder.xrowop2_like(
            decoder.subgroup_constraint_matrix_from_irrep_matrices(matrices)
        )
        count, direction = decoder.subgroup_direction_from_reduced_matrix(reduced)
        if count <= 0:
            return
        full_dim = int(decoder.little_record_by_gid(int(row.gid)).full_dim)
        aliases.append(
            DynamicSubductionRow(
                gid=int(row.gid),
                irrep_label=str(row.irrep_label),
                k_label=str(row.k_label),
                source_kparam=source_kparam,
                reciprocal_vector_pml=tuple(occurrence.reciprocal_vector_pml),
                direction_matrix=tuple(
                    tuple(float(direction[free, coordinate]) for free in range(count))
                    for coordinate in range(full_dim)
                ),
                source_occurrences=(occurrence,),
            )
        )

    for occurrence in row.source_occurrences:
        source_kparam = _source_kparam_for_vector(
            decoder,
            int(sg),
            int(row.gid),
            occurrence.reciprocal_vector_pml,
        )
        if source_kparam is None:
            continue
        append_alias(source_kparam, occurrence)

    if requested_source_kparam is not None:
        requested = tuple(int(value) for value in requested_source_kparam)
        matching_occurrence = next(
            (
                occurrence
                for occurrence in row.source_occurrences
                if tuple(int(value) for value in occurrence.source_kparam)
                == requested
            ),
            None,
        )
        if matching_occurrence is not None:
            append_alias(requested, matching_occurrence)
    return tuple(aliases)


def dynamic_magnetic_subduction_occurrence_alias_rows(
    decoder: Any,
    magnetic_source: Any,
    *,
    sg: int,
    basis: Sequence[int],
    operations: Sequence[OperationRecord],
    row: DynamicSubductionRow,
    requested_source_kparam: Sequence[int] | None = None,
) -> tuple[DynamicSubductionRow, ...]:
    """Recover invariant magnetic Cases hidden by physical-orbit dedupe."""

    basis_values = tuple(int(value) for value in basis)
    table = magnetic_source._magnetic_table() if hasattr(magnetic_source, "_magnetic_table") else None
    if table is None:
        from APOSTRUCT.Backend.source.magnetic import data as magnetic_data

        table = magnetic_data().table
    identity_point_op = int(decoder.generate_space_group_records(int(sg))[0][4])
    identity_magnetic_op = int(table["mag_point_op_nonmag2mag"][identity_point_op - 1])
    constrained_operations = tuple(operations) + tuple(
        (
            basis_values[3 * axis],
            basis_values[3 * axis + 1],
            basis_values[3 * axis + 2],
            1,
            identity_magnetic_op,
        )
        for axis in range(3)
    )

    def parameter_key(kparam: Sequence[int]) -> tuple[Fraction, Fraction, Fraction]:
        return tuple(
            Fraction(int(kparam[index]), int(kparam[3]))
            for index in range(3)
        )  # type: ignore[return-value]

    def periodic_key(kparam: Sequence[int]) -> tuple[Fraction, Fraction, Fraction]:
        return tuple(value % 1 for value in parameter_key(kparam))  # type: ignore[return-value]

    seen = {parameter_key(row.source_kparam)}
    aliases: list[DynamicSubductionRow] = []
    requested = (
        tuple(int(value) for value in requested_source_kparam)
        if requested_source_kparam is not None
        else None
    )
    requested_key = parameter_key(requested) if requested is not None else None
    exact_requested_occurrence = (
        next(
            (
                occurrence
                for occurrence in row.source_occurrences
                if tuple(int(value) for value in occurrence.source_kparam)
                == requested
            ),
            None,
        )
        if requested is not None
        else None
    )

    def append_alias(
        source_kparam: Sequence[int],
        occurrence: DynamicKOccurrence,
    ) -> None:
        source_kparam = tuple(int(value) for value in source_kparam)
        key = parameter_key(source_kparam)
        if key in seen:
            return
        seen.add(key)
        matrices = [
            magnetic_source.get_irreps_matrix_for_magnetic_record(
                int(row.gid),
                record,
                source_kparam,
            )
            for record in constrained_operations
        ]
        reduced = decoder.xrowop2_like(
            decoder.subgroup_constraint_matrix_from_irrep_matrices(matrices)
        )
        count, direction = decoder.subgroup_direction_from_reduced_matrix(reduced)
        if count <= 0:
            return
        full_dim = int(decoder.little_record_by_gid(int(row.gid)).full_dim)
        aliases.append(
            DynamicSubductionRow(
                gid=int(row.gid),
                irrep_label=str(row.irrep_label),
                k_label=str(row.k_label),
                source_kparam=tuple(int(value) for value in source_kparam),
                reciprocal_vector_pml=tuple(occurrence.reciprocal_vector_pml),
                direction_matrix=tuple(
                    tuple(float(direction[free, coordinate]) for free in range(count))
                    for coordinate in range(full_dim)
                ),
                source_occurrences=(occurrence,),
            )
        )

    occurrence_kparams: list[tuple[DynamicKOccurrence, tuple[int, int, int, int]]] = []
    for occurrence in row.source_occurrences:
        source_kparam = _source_kparam_for_vector(
            decoder,
            int(sg),
            int(row.gid),
            occurrence.reciprocal_vector_pml,
        )
        if source_kparam is None:
            continue
        occurrence_kparams.append((occurrence, source_kparam))
        if (
            exact_requested_occurrence is not None
            and parameter_key(source_kparam) == requested_key
        ):
            append_alias(requested, exact_requested_occurrence)
        else:
            append_alias(source_kparam, occurrence)

    if requested is not None and requested_key not in seen:
        matching_occurrence = exact_requested_occurrence
        if matching_occurrence is None:
            requested_periodic = periodic_key(requested)
            matching_occurrence = next(
                (
                    occurrence
                    for occurrence, source_kparam in occurrence_kparams
                    if periodic_key(source_kparam) == requested_periodic
                ),
                None,
            )
        if matching_occurrence is not None:
            append_alias(requested, matching_occurrence)
    return tuple(aliases)


def _magnetic_subduction_rows_for_families(
    decoder: Any,
    magnetic_source: Any,
    *,
    sg: int,
    basis: Sequence[int],
    operations: Sequence[OperationRecord],
    families: Sequence[_SelectedFamily],
    preferred_kparams: dict[str, Sequence[int] | Sequence[Sequence[int]]] | None = None,
) -> tuple[DynamicSubductionRow, ...]:
    """Return magnetic fixed directions for an ordered set of k families."""

    basis_values = tuple(int(value) for value in basis)
    table = magnetic_source._magnetic_table() if hasattr(magnetic_source, "_magnetic_table") else None
    if table is None:
        from APOSTRUCT.Backend.source.magnetic import data as magnetic_data

        table = magnetic_data().table
    identity_point_op = int(decoder.generate_space_group_records(int(sg))[0][4])
    identity_magnetic_op = int(table["mag_point_op_nonmag2mag"][identity_point_op - 1])
    constrained_operations = tuple(operations) + tuple(
        (basis_values[3 * row], basis_values[3 * row + 1], basis_values[3 * row + 2], 1, identity_magnetic_op)
        for row in range(3)
    )
    rows: list[DynamicSubductionRow] = []
    seen_rows: dict[tuple[int, tuple[tuple[int, int], ...]], int] = {}
    row_preferred: list[bool] = []
    for family in families:
        k_label = family.k_label
        discovered_kparam = family.source_kparam
        reciprocal_vector_pml = family.reciprocal_vector_pml
        source_kparam = discovered_kparam
        for gid, row_sg in enumerate(decoder.iso.little["little_irr_space_group"], start=1):
            if int(row_sg) != sg:
                continue
            kslot = int(decoder.iso.little["little_irr_k"][gid - 1])
            lattice = int(decoder.iso.space["ispace_lattice"][sg - 1])
            if str(decoder.iso.little["little_k_label"][(lattice - 1) * 27 + kslot - 1]).strip() != k_label:
                continue
            matrices = [
                magnetic_source.get_irreps_matrix_for_magnetic_record(gid, record, source_kparam)
                for record in constrained_operations
            ]
            reduced = decoder.xrowop2_like(decoder.subgroup_constraint_matrix_from_irrep_matrices(matrices))
            count, direction = decoder.subgroup_direction_from_reduced_matrix(reduced)
            if count <= 0:
                continue
            full_dim = int(decoder.little_record_by_gid(gid).full_dim)
            row_key = (
                int(gid),
                _reciprocal_point_group_orbit_key(decoder, sg, reciprocal_vector_pml),
            )
            candidate = DynamicSubductionRow(
                gid=gid,
                irrep_label=str(decoder.iso.little["little_irr_full_label"][gid - 1]).strip(),
                k_label=k_label,
                source_kparam=source_kparam,
                reciprocal_vector_pml=reciprocal_vector_pml,
                direction_matrix=tuple(
                    tuple(float(direction[free, coordinate]) for free in range(count))
                    for coordinate in range(full_dim)
                ),
                source_occurrences=family.source_occurrences,
            )
            candidate_preferred = _matches_preferred_kparam(preferred_kparams, k_label, source_kparam)
            if row_key in seen_rows:
                index = seen_rows[row_key]
                if candidate_preferred and not row_preferred[index]:
                    rows[index] = replace(
                        candidate,
                        source_occurrences=rows[index].source_occurrences + candidate.source_occurrences,
                    )
                    row_preferred[index] = True
                else:
                    rows[index] = replace(
                        rows[index],
                        source_occurrences=rows[index].source_occurrences + candidate.source_occurrences,
                    )
                continue
            seen_rows[row_key] = len(rows)
            rows.append(candidate)
            row_preferred.append(candidate_preferred)
    return tuple(rows)


def dynamic_magnetic_subduction_rows(
    decoder: Any,
    magnetic_source: Any,
    *,
    sg: int,
    basis: Sequence[int],
    operations: Sequence[OperationRecord],
    preferred_kparams: dict[str, Sequence[int] | Sequence[Sequence[int]]] | None = None,
) -> tuple[DynamicSubductionRow, ...]:
    """Return magnetic invariant directions using time-reversal matrices."""

    basis_values = tuple(int(value) for value in basis)
    return _magnetic_subduction_rows_for_families(
        decoder,
        magnetic_source,
        sg=int(sg),
        basis=basis_values,
        operations=operations,
        families=_selected_family_records(
            decoder,
            int(sg),
            basis_values,
            preferred_kparams,
        ),
        preferred_kparams=preferred_kparams,
    )


def display_direction_magnetic_subduction_rows(
    decoder: Any,
    magnetic_source: Any,
    *,
    sg: int,
    basis: Sequence[int],
    operations: Sequence[OperationRecord],
) -> tuple[DynamicSubductionRow, ...]:
    """Return magnetic reverse-direction rows in binary family order."""

    basis_values = tuple(int(value) for value in basis)
    return _magnetic_subduction_rows_for_families(
        decoder,
        magnetic_source,
        sg=int(sg),
        basis=basis_values,
        operations=operations,
        families=_display_direction_family_records(decoder, int(sg), basis_values),
    )
