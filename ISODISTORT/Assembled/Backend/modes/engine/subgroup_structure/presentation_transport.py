"""Source-only child-orbit to parent-site coset factorization.

The local projection adapter uses parent Wyckoff coset representatives from
Source records.  Complete-mode presentation organizes those occurrences under
the selected child subgroup.  A child operation generally does not equal the
stored parent occurrence record: it factors through a parent-site stabilizer
and a selected-child lattice translation.  This module exposes that metadata
without choosing a Web print basis.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
import math
from typing import Any, Iterable, Sequence

import numpy as np

from ISODISTORT.Assembled.Backend.exactmath import (
    fraction_matrix_inverse3 as _matrix_inverse,
    fraction_matrix_multiply3 as _matrix_multiply,
    fraction_row_multiply3 as _row_multiply,
)


OperationRecord = tuple[int, int, int, int, int]
FractionVector = tuple[Fraction, Fraction, Fraction]
FractionMatrix = tuple[FractionVector, FractionVector, FractionVector]


@dataclass(frozen=True)
class PresentationTransportFactor:
    """One exact setting-aware site-coset factorization row."""

    orbit_representative_index: int
    child_operation_index: int
    child_point_op: int
    matched_occurrence_index: int
    child_operation_record: OperationRecord
    presentation_record: OperationRecord
    stabilizer_record: OperationRecord
    lattice_translation: tuple[int, int, int]
    lattice_coordinates: tuple[int, int, int]
    affine_rotation_equal: bool
    affine_translation_equal: bool


@dataclass(frozen=True)
class ChildSiteIrrepProjector:
    """One Source child-site irrep contained in a restricted parent irrep."""

    child_site_pg: int
    child_pg_irrep: int
    child_irrep_old_id: int
    child_irrep_label: str
    child_irrep_dimension: int
    rank: int
    projector: np.ndarray
    basis: np.ndarray


@dataclass(frozen=True)
class SiteVectorPrintColumn:
    """One canonical Source vector column for a site-PG irrep occurrence."""

    site_pg: int
    pg_irrep: int
    vector_setting: int
    vector_basis_id: int
    component_start: int
    basis_column: int
    vector: tuple[float, float, float]


@dataclass(frozen=True)
class SourcePrintIntertwiner:
    """A unique Source-column change of basis for one printed family."""

    matrix: np.ndarray
    source_rank: int
    target_rank: int
    residual: float
    condition_number: float


@dataclass(frozen=True)
class CanonicalRankOnePrintScalar:
    """One relative rank-one print gauge shared by canonical child orbits."""

    scalar: float
    contributing_orbit_indices: tuple[int, ...]
    orbit_scalars: tuple[float, ...]
    residual: float


@dataclass(frozen=True)
class ChildPrintSubductionBasis:
    """A projector basis oriented by Source child-site print columns."""

    basis: np.ndarray
    print_change: np.ndarray
    transported_columns: np.ndarray
    residual: float


@dataclass(frozen=True)
class ChildPrintSubductionEmbedding:
    """A projector image embedded in Source repeated child print columns."""

    basis: np.ndarray
    print_coefficients: np.ndarray
    transported_columns: np.ndarray
    child_print_columns: np.ndarray
    residual: float


@dataclass(frozen=True)
class MinimalPrintedProjectBlock:
    """The runtime print-topology identity needed by the Assembled kernel."""

    printed_block_index: int
    target_occurrence: int
    vector_component: int
    target_vector_rep: int
    project_family_index: int
    selected_candidate_index: int | None = None
    raw_star_arm: int | None = None


@dataclass(frozen=True)
class FactorCarrierPullback:
    """Factor carriers pulled back to their raw target-row gauges."""

    coefficients: np.ndarray
    pullback_matrices: np.ndarray
    target_records: tuple[OperationRecord, ...]


@dataclass(frozen=True)
class CanonicalChildPrintFactors:
    """Identity-child-operation factors defining Source print representatives."""

    factors: tuple[PresentationTransportFactor, ...]
    source_factor_indices: tuple[int, ...]
    orbit_representative_indices: tuple[int, ...]


@dataclass(frozen=True)
class PrintedProjectBlockSource:
    """Source project and site columns for one runtime print block."""

    printed_block_index: int
    project_family_index: int
    vector_component: int
    target_vector_rep: int
    project_matrix: np.ndarray
    site_vectors: np.ndarray
    selected_candidate_index: int | None = None
    raw_star_arm: int | None = None


def _point_matrix(decoder: Any, sg: int, point_op: int) -> FractionMatrix:
    units: FractionMatrix = (
        (Fraction(1), Fraction(0), Fraction(0)),
        (Fraction(0), Fraction(1), Fraction(0)),
        (Fraction(0), Fraction(0), Fraction(1)),
    )
    return tuple(
        tuple(decoder.vrot_fraction(int(sg), int(point_op), unit))
        for unit in units
    )  # type: ignore[return-value]


def _record_affine(decoder: Any, sg: int, record: OperationRecord) -> tuple[FractionMatrix, FractionVector]:
    denominator = int(record[3])
    if denominator == 0:
        raise ValueError(f"zero operation denominator: {record}")
    return (
        _point_matrix(decoder, int(sg), int(record[4])),
        tuple(Fraction(int(record[axis]), denominator) for axis in range(3)),  # type: ignore[arg-type]
    )


def _compose_affine(
    first: tuple[FractionMatrix, FractionVector],
    then: tuple[FractionMatrix, FractionVector],
) -> tuple[FractionMatrix, FractionVector]:
    first_matrix, first_translation = first
    then_matrix, then_translation = then
    rotated_translation = _row_multiply(first_translation, then_matrix)
    return (
        _matrix_multiply(first_matrix, then_matrix),
        tuple(
            rotated_translation[axis] + then_translation[axis]
            for axis in range(3)
        ),  # type: ignore[arg-type]
    )


def _fraction_record(values: FractionVector, point_op: int) -> OperationRecord:
    denominator = math.lcm(*(value.denominator for value in values), 1)
    numerators = [int(value * denominator) for value in values]
    divisor = math.gcd(denominator, math.gcd(abs(numerators[0]), math.gcd(abs(numerators[1]), abs(numerators[2]))))
    if divisor > 1:
        denominator //= divisor
        numerators = [value // divisor for value in numerators]
    return numerators[0], numerators[1], numerators[2], denominator, int(point_op)


def _record_from_affine(
    decoder: Any,
    parent_sg: int,
    affine: tuple[FractionMatrix, FractionVector],
) -> OperationRecord:
    matrix, translation = affine
    point_op_by_matrix = {
        _point_matrix(decoder, int(parent_sg), int(record[4])): int(record[4])
        for record in decoder.generate_space_group_records(int(parent_sg))
    }
    point_op = point_op_by_matrix.get(matrix)
    if point_op is None:
        raise KeyError(f"presentation rotation is not represented in parent SG{parent_sg}: {matrix}")
    return _fraction_record(translation, point_op)


def representative_operation_record_from_cinter(
    decoder: Any,
    sg: int,
    rotation: Sequence[Sequence[object]],
    translation_minus_lattice: Sequence[object],
) -> OperationRecord:
    """Convert a cinter representative operation to a full PML record.

    ``translation_minus_lattice`` is intentionally not reduced modulo the
    parent lattice: its integer part contributes to non-GM Bloch phases.
    Coordinates use the Source row convention ``x' = x R + t``.
    """

    if len(rotation) != 3 or any(len(row) != 3 for row in rotation):
        raise ValueError(f"expected a 3x3 representative rotation, got {rotation!r}")
    if len(translation_minus_lattice) != 3:
        raise ValueError(
            f"expected a representative translation triplet, got {translation_minus_lattice!r}"
        )
    cinter_rotation: FractionMatrix = tuple(
        tuple(Fraction(value) for value in row)
        for row in rotation
    )  # type: ignore[assignment]
    cinter_translation: FractionVector = tuple(
        Fraction(value) for value in translation_minus_lattice
    )  # type: ignore[assignment]
    pml_to_cinter: FractionMatrix = tuple(
        tuple(Fraction(value) for value in row)
        for row in decoder.pml_to_cinter_matrix(int(sg))
    )  # type: ignore[assignment]
    cinter_to_pml: FractionMatrix = tuple(
        tuple(Fraction(value) for value in row)
        for row in decoder.cinter_to_pml_matrix(int(sg))
    )  # type: ignore[assignment]
    origin: FractionVector = tuple(
        Fraction(value) for value in decoder._cml_to_cinter_origin(int(sg))  # noqa: SLF001
    )  # type: ignore[assignment]
    pml_rotation = _matrix_multiply(
        _matrix_multiply(pml_to_cinter, cinter_rotation),
        cinter_to_pml,
    )
    affine_translation = tuple(
        _row_multiply(origin, cinter_rotation)[axis]
        + cinter_translation[axis]
        - origin[axis]
        for axis in range(3)
    )
    pml_translation = _row_multiply(affine_translation, cinter_to_pml)  # type: ignore[arg-type]
    return _record_from_affine(
        decoder,
        int(sg),
        (pml_rotation, pml_translation),
    )


def compose_seitz_records(
    decoder: Any,
    sg: int,
    *,
    first: Sequence[int],
    then: Sequence[int],
) -> OperationRecord:
    """Compose two full PML records in Source row-vector order.

    Applying canonicalizer ``c`` before raw occurrence ``h`` is therefore
    ``compose_seitz_records(first=c, then=h)``.
    """

    first_record = tuple(int(value) for value in first)
    then_record = tuple(int(value) for value in then)
    if len(first_record) != 5 or len(then_record) != 5:
        raise ValueError(f"expected two five-field operation records: {first!r}, {then!r}")
    return _record_from_affine(
        decoder,
        int(sg),
        _compose_affine(
            _record_affine(decoder, int(sg), first_record),  # type: ignore[arg-type]
            _record_affine(decoder, int(sg), then_record),  # type: ignore[arg-type]
        ),
    )


def embedded_child_operation_records(
    decoder: Any,
    *,
    parent_sg: int,
    child_sg: int,
    subgroup_basis: Sequence[int],
    subgroup_origin: Sequence[int],
) -> tuple[OperationRecord, ...]:
    """Map Source-ordered child operations into the raw parent PML embedding."""

    basis_values = tuple(int(value) for value in subgroup_basis)
    if len(basis_values) != 9:
        raise ValueError(f"expected 9 subgroup basis values, got {len(basis_values)}")
    origin_values = tuple(int(value) for value in subgroup_origin)
    if len(origin_values) != 4 or origin_values[3] == 0:
        raise ValueError(f"invalid subgroup origin: {subgroup_origin}")
    basis: FractionMatrix = tuple(
        tuple(Fraction(basis_values[row * 3 + col]) for col in range(3))
        for row in range(3)
    )  # type: ignore[assignment]
    inverse = _matrix_inverse(basis)
    origin: FractionVector = tuple(Fraction(origin_values[axis], origin_values[3]) for axis in range(3))  # type: ignore[assignment]
    parent_point_ops = {
        _point_matrix(decoder, int(parent_sg), int(record[4])): int(record[4])
        for record in decoder.generate_space_group_records(int(parent_sg))
    }
    out: list[OperationRecord] = []
    for record in decoder.generate_space_group_records(int(child_sg)):
        child_matrix, child_translation = _record_affine(decoder, int(child_sg), record)
        parent_matrix = _matrix_multiply(inverse, _matrix_multiply(child_matrix, basis))
        point_op = parent_point_ops.get(parent_matrix)
        if point_op is None:
            raise KeyError(
                f"child SG{child_sg} operation {record[4]} does not embed in parent SG{parent_sg}"
            )
        parent_translation = tuple(
            origin[axis]
            - _row_multiply(origin, parent_matrix)[axis]
            + _row_multiply(child_translation, basis)[axis]
            for axis in range(3)
        )
        out.append(_fraction_record(parent_translation, point_op))
    return tuple(out)


def _lattice_coordinates(vector: FractionVector, basis_inverse: FractionMatrix) -> tuple[int, int, int] | None:
    coordinates = _row_multiply(vector, basis_inverse)
    if not all(value.denominator == 1 for value in coordinates):
        return None
    return tuple(int(value) for value in coordinates)  # type: ignore[return-value]


def _site_cosets_use_stabilizer_first(
    decoder: Any,
    parent_sg: int,
    parent_wyckoff_row: Any,
    stabilizers: Sequence[OperationRecord],
) -> bool:
    """Recover the affine side used by the Source Wyckoff coset table."""

    cosets = tuple(decoder.wyc_pg_cosets_records(int(parent_sg), stabilizers))
    generated = tuple(decoder.generate_space_group_records(int(parent_sg)))

    def modulo_parent_lattice(affine: tuple[FractionMatrix, FractionVector]) -> tuple[Any, ...]:
        matrix, translation = affine
        return matrix, tuple(value % 1 for value in translation)

    generated_set = {
        modulo_parent_lattice(_record_affine(decoder, int(parent_sg), record))
        for record in generated
    }
    coset_affines = tuple(_record_affine(decoder, int(parent_sg), record) for record in cosets)
    stabilizer_affines = tuple(
        _record_affine(decoder, int(parent_sg), record)
        for record in stabilizers
    )

    def covers_generated(stabilizer_first: bool) -> bool:
        products = {
            modulo_parent_lattice(
                _compose_affine(stabilizer, coset)
                if stabilizer_first
                else _compose_affine(coset, stabilizer)
            )
            for coset in coset_affines
            for stabilizer in stabilizer_affines
        }
        return products == generated_set

    stabilizer_first_valid = covers_generated(True)
    occurrence_first_valid = covers_generated(False)
    if stabilizer_first_valid != occurrence_first_valid:
        return stabilizer_first_valid

    # Some site subgroups are normal in the generated parent group, so both
    # sides cover the same operations.  Preserve the established alternate
    # R-setting convention in that genuinely ambiguous case.
    if stabilizer_first_valid:
        return int(generated[0][4]) != 1
    raise ValueError(
        f"Source Wyckoff cosets do not reconstruct SG{parent_sg} "
        f"for row {getattr(parent_wyckoff_row, 'label', '?')}"
    )


def factor_presentation_occurrences(
    decoder: Any,
    *,
    parent_sg: int,
    parent_wyckoff_row: Any,
    child_sg: int,
    subgroup_basis: Sequence[int],
    subgroup_origin: Sequence[int],
    occurrence_records: Iterable[Sequence[int]],
) -> tuple[PresentationTransportFactor, ...]:
    """Factor child actions through stored parent occurrence/site cosets.

    All ordering is Source ordering: occurrence records, child operations, and
    parent-site stabilizer records.  Every returned row verifies the exact
    affine identity modulo a selected-child lattice translation.  Multiple
    factors may target the same occurrence when the child subgroup stabilizes
    that child orbit; callers must retain them for site-irrep subduction rather
    than selecting one operation heuristically.
    """

    occurrences = tuple(tuple(int(value) for value in record) for record in occurrence_records)
    basis_values = tuple(int(value) for value in subgroup_basis)
    basis: FractionMatrix = tuple(
        tuple(Fraction(basis_values[row * 3 + col]) for col in range(3))
        for row in range(3)
    )  # type: ignore[assignment]
    basis_inverse = _matrix_inverse(basis)
    child_records = embedded_child_operation_records(
        decoder,
        parent_sg=int(parent_sg),
        child_sg=int(child_sg),
        subgroup_basis=basis_values,
        subgroup_origin=subgroup_origin,
    )
    raw_child_records = tuple(decoder.generate_space_group_records(int(child_sg)))
    if len(raw_child_records) != len(child_records):
        raise ValueError(f"child operation count changed during embedding: {len(raw_child_records)} != {len(child_records)}")
    stabilizers = tuple(decoder.wyc_pg_elements_records(int(parent_sg), parent_wyckoff_row))
    stabilizer_first = _site_cosets_use_stabilizer_first(
        decoder,
        int(parent_sg),
        parent_wyckoff_row,
        stabilizers,
    )
    occurrence_affines = tuple(_record_affine(decoder, int(parent_sg), record) for record in occurrences)
    stabilizer_affines = tuple(_record_affine(decoder, int(parent_sg), record) for record in stabilizers)
    child_affines = tuple(_record_affine(decoder, int(parent_sg), record) for record in child_records)

    # This is an exact runtime index over the Source-record candidate scan.
    # Candidate occurrence/site products do not depend on the orbit
    # representative or child operation.  Preserve their original nested
    # Source order while composing each exact affine only once.
    candidates_by_rotation: dict[
        FractionMatrix,
        list[tuple[int, int]],
    ] = {}
    for occurrence_index, occurrence_affine in enumerate(occurrence_affines):
        for stabilizer_index, stabilizer_affine in enumerate(stabilizer_affines):
            candidate_rotation = (
                _matrix_multiply(stabilizer_affine[0], occurrence_affine[0])
                if stabilizer_first
                else _matrix_multiply(occurrence_affine[0], stabilizer_affine[0])
            )
            candidates_by_rotation.setdefault(candidate_rotation, []).append(
                (occurrence_index, stabilizer_index)
            )

    factors: list[PresentationTransportFactor] = []
    unassigned = set(range(len(occurrences)))
    for representative_index in range(len(occurrences)):
        if representative_index not in unassigned:
            continue
        orbit_indices: set[int] = set()
        representative = occurrence_affines[representative_index]
        for child_index, (child_record, child_affine) in enumerate(zip(child_records, child_affines)):
            presentation_affine = _compose_affine(representative, child_affine)
            presentation_record = _record_from_affine(decoder, int(parent_sg), presentation_affine)
            match = None
            for occurrence_index, stabilizer_index in (
                candidates_by_rotation.get(presentation_affine[0], ())
            ):
                occurrence_affine = occurrence_affines[occurrence_index]
                stabilizer_affine = stabilizer_affines[stabilizer_index]
                if stabilizer_first:
                    rotated_translation = _row_multiply(
                        stabilizer_affine[1], occurrence_affine[0]
                    )
                    candidate_translation = tuple(
                        rotated_translation[axis] + occurrence_affine[1][axis]
                        for axis in range(3)
                    )
                else:
                    rotated_translation = _row_multiply(
                        occurrence_affine[1], stabilizer_affine[0]
                    )
                    candidate_translation = tuple(
                        rotated_translation[axis] + stabilizer_affine[1][axis]
                        for axis in range(3)
                    )
                lattice_translation = tuple(
                    presentation_affine[1][axis] - candidate_translation[axis]
                    for axis in range(3)
                )
                lattice_coordinates = _lattice_coordinates(lattice_translation, basis_inverse)
                if lattice_coordinates is None:
                    continue
                match = (
                    occurrence_index,
                    stabilizers[stabilizer_index],
                    lattice_translation,
                    lattice_coordinates,
                    candidate_translation,
                )
                break
            if match is None:
                raise KeyError(
                    f"no parent occurrence/site factor for child operation {child_record} "
                    f"from representative {representative_index}"
                )
            occurrence_index, stabilizer_record, lattice_translation, lattice_coordinates, candidate_translation = match
            orbit_indices.add(int(occurrence_index))
            translated_candidate = tuple(
                candidate_translation[axis] + lattice_translation[axis]
                for axis in range(3)
            )
            factors.append(
                PresentationTransportFactor(
                    orbit_representative_index=int(representative_index),
                    child_operation_index=int(child_index),
                    child_point_op=int(raw_child_records[child_index][4]),
                    matched_occurrence_index=int(occurrence_index),
                    child_operation_record=child_record,
                    presentation_record=presentation_record,
                    stabilizer_record=tuple(int(value) for value in stabilizer_record),
                    lattice_translation=tuple(int(value) for value in lattice_translation),
                    lattice_coordinates=lattice_coordinates,
                    affine_rotation_equal=True,
                    affine_translation_equal=translated_candidate == presentation_affine[1],
                )
            )
        unassigned.difference_update(orbit_indices)
    return tuple(factors)


def site_stabilizer_irrep_matrix(
    decoder: Any,
    *,
    site_pg: int,
    pg_irrep: int,
    parent_site_records: Iterable[Sequence[int]],
    stabilizer_record: OperationRecord,
) -> np.ndarray:
    """Return the faithful site-irrep matrix associated with one factor row."""

    parent_records = tuple(tuple(int(value) for value in record) for record in parent_site_records)
    mapped = decoder.map_parent_ops_to_site_project_records(
        int(site_pg),
        [int(record[4]) for record in parent_records],
    )
    if len(mapped) != len(parent_records):
        raise ValueError(f"site operation mapping lost records: {len(mapped)} != {len(parent_records)}")
    try:
        parent_index = next(
            index
            for index, record in enumerate(parent_records)
            if int(record[4]) == int(stabilizer_record[4])
        )
    except StopIteration as exc:
        raise KeyError(f"stabilizer {stabilizer_record} is not in the parent site setting") from exc
    matrices = decoder.site_little_sparse_matrices(int(site_pg), int(pg_irrep))
    return np.asarray(matrices[int(mapped[parent_index][4])], dtype=float)


def _projector_basis(projector: np.ndarray, *, tolerance: float) -> np.ndarray:
    """Choose a deterministic Source-column basis for a projector image."""

    accepted: list[np.ndarray] = []
    for column in range(projector.shape[1]):
        candidate = np.asarray(projector[:, column], dtype=complex)
        if np.linalg.norm(candidate) <= tolerance:
            continue
        trial = np.column_stack([*accepted, candidate]) if accepted else candidate.reshape(-1, 1)
        if np.linalg.matrix_rank(trial, tol=tolerance) <= len(accepted):
            continue
        pivot = next((value for value in candidate if abs(value) > tolerance), 1.0 + 0.0j)
        candidate = candidate / pivot
        candidate[np.abs(candidate) < tolerance] = 0.0
        accepted.append(candidate)
    if not accepted:
        return np.zeros((projector.shape[0], 0), dtype=complex)
    return np.column_stack(accepted)


def child_site_irrep_projectors(
    decoder: Any,
    *,
    parent_site_pg: int,
    parent_pg_irrep: int,
    parent_site_records: Iterable[Sequence[int]],
    child_site_pg: int,
    child_stabilizer_point_ops: Sequence[int],
    factors: Iterable[PresentationTransportFactor],
    tolerance: float = 1e-10,
) -> tuple[ChildSiteIrrepProjector, ...]:
    """Subduce one parent-site irrep onto Source child-site irrep projectors.

    ``factors`` must contain the child operations that stabilize one selected
    child orbit representative.  Their parent-site stabilizer records provide
    the restricted parent representation ``S(s_h)``.  Child operations are
    mapped as one complete setting before characters are read, preserving the
    Source setting correspondence that a single operation cannot identify.
    """

    factor_by_child_op = {int(item.child_point_op): item for item in factors}
    child_ops = tuple(int(value) for value in child_stabilizer_point_ops)
    if not child_ops:
        return ()
    missing = [point_op for point_op in child_ops if point_op not in factor_by_child_op]
    if missing:
        raise KeyError(f"missing presentation factors for child operations {missing}")
    mapped_child = decoder.map_parent_ops_to_site_project_records(int(child_site_pg), child_ops)
    if len(mapped_child) != len(child_ops):
        raise ValueError(f"child-site setting mapping lost records: {len(mapped_child)} != {len(child_ops)}")
    parent_matrices = tuple(
        site_stabilizer_irrep_matrix(
            decoder,
            site_pg=int(parent_site_pg),
            pg_irrep=int(parent_pg_irrep),
            parent_site_records=parent_site_records,
            stabilizer_record=factor_by_child_op[point_op].stabilizer_record,
        )
        for point_op in child_ops
    )
    count = int(decoder.iso.wyckoff["iwyckoff_pg_irrep_count"][int(child_site_pg) - 1])
    out: list[ChildSiteIrrepProjector] = []
    for child_pg_irrep in range(1, count + 1):
        child_matrices = decoder.site_little_sparse_matrices(int(child_site_pg), int(child_pg_irrep))
        canonical_ops = [int(record[4]) for record in mapped_child]
        representation = [np.asarray(child_matrices[point_op], dtype=complex) for point_op in canonical_ops]
        dimension = int(representation[0].shape[0])
        projector = sum(
            np.conjugate(np.trace(child_matrix)) * np.asarray(parent_matrix, dtype=complex)
            for child_matrix, parent_matrix in zip(representation, parent_matrices)
        ) * (float(dimension) / float(len(child_ops)))
        projector[np.abs(projector) < tolerance] = 0.0
        rank = int(np.linalg.matrix_rank(projector, tol=tolerance))
        if rank <= 0:
            continue
        basis = _projector_basis(projector, tolerance=tolerance)
        out.append(
            ChildSiteIrrepProjector(
                child_site_pg=int(child_site_pg),
                child_pg_irrep=int(child_pg_irrep),
                child_irrep_old_id=int(decoder.site_pg_irrep_old_id(int(child_site_pg), int(child_pg_irrep))),
                child_irrep_label=str(decoder.site_pg_irrep_label(int(child_site_pg), int(child_pg_irrep))),
                child_irrep_dimension=dimension,
                rank=rank,
                projector=projector,
                basis=basis,
            )
        )
    return tuple(out)


def site_vector_print_columns(
    decoder: Any,
    *,
    site_pg: int,
    pg_irrep: int,
    site_operation_records: Iterable[Sequence[int]],
    vector_setting: int,
) -> tuple[SiteVectorPrintColumn, ...]:
    """Return Source's canonical repeated-component print basis.

    Character projectors determine an invariant subspace but cannot determine
    the sign of a rank-one basis (or a basis inside repeated copies).  The
    missing convention is stored separately in
    ``iwyckoff_pg_vector_basis``.  This helper exposes those exact columns in
    repeated component order, using the full site-operation setting to select
    the correct vector-basis row.
    """

    vector_setting = int(vector_setting)
    if vector_setting not in (1, 2):
        raise ValueError(f"unsupported vector setting: {vector_setting}")
    records = tuple(tuple(int(value) for value in record) for record in site_operation_records)
    vector_basis_id = int(decoder.vector_basis_id_for_site_records(int(site_pg), records))
    reps = tuple(int(value) for value in decoder.site_vector_reps(int(site_pg)))
    selected_reps = reps[(vector_setting - 1) * 3:vector_setting * 3]
    old_id = int(decoder.site_pg_irrep_old_id(int(site_pg), int(pg_irrep)))
    dimension = int(decoder.image_record(old_id).dimension)
    basis_codes = decoder.iso.wyckoff["iwyckoff_pg_vector_basis"]
    out: list[SiteVectorPrintColumn] = []
    for component_start, component_irrep in enumerate(selected_reps, start=1):
        if int(component_irrep) != int(pg_irrep):
            continue
        for basis_column in range(dimension):
            code_base = (
                (vector_basis_id - 1) * 18
                + (vector_setting - 1) * 9
                + (component_start - 1) * 3
                + basis_column * 3
            )
            vector = tuple(
                float(decoder.iso.const[int(basis_codes[code_base + axis])])
                for axis in range(3)
            )
            out.append(
                SiteVectorPrintColumn(
                    site_pg=int(site_pg),
                    pg_irrep=int(pg_irrep),
                    vector_setting=vector_setting,
                    vector_basis_id=vector_basis_id,
                    component_start=int(component_start),
                    basis_column=int(basis_column),
                    vector=vector,  # type: ignore[arg-type]
                )
            )
    return tuple(out)


def solve_source_print_intertwiner(
    source_columns: Sequence[Sequence[complex | float]],
    target_columns: Sequence[Sequence[complex | float]],
    *,
    tolerance: float = 1e-10,
    condition_limit: float = 1e10,
) -> SourcePrintIntertwiner | None:
    """Solve ``source @ T = target`` only when the print gauge is unique.

    This helper applies no first-nonzero or label heuristic.  It returns
    ``None`` when source columns are dependent, source/target spans differ,
    column counts differ, or the solve is ill-conditioned.  Rank-one signs and
    higher-dimensional column remixes therefore share the same Source-only
    contract.
    """

    source = np.asarray(source_columns, dtype=complex)
    target = np.asarray(target_columns, dtype=complex)
    if source.ndim != 2 or target.ndim != 2:
        raise ValueError(f"print columns must be matrices, got {source.shape} and {target.shape}")
    if source.shape != target.shape or source.shape[1] == 0:
        return None
    source_rank = int(np.linalg.matrix_rank(source, tol=tolerance))
    target_rank = int(np.linalg.matrix_rank(target, tol=tolerance))
    column_count = int(source.shape[1])
    if source_rank != column_count or target_rank != column_count:
        return None
    combined_rank = int(np.linalg.matrix_rank(np.column_stack((source, target)), tol=tolerance))
    if combined_rank != column_count:
        return None
    singular_values = np.linalg.svd(source, compute_uv=False)
    if not len(singular_values) or float(singular_values[-1]) <= tolerance:
        return None
    condition_number = float(singular_values[0] / singular_values[-1])
    if not np.isfinite(condition_number) or condition_number > float(condition_limit):
        return None
    matrix, _residuals, _rank, _singular = np.linalg.lstsq(source, target, rcond=None)
    residual = float(np.max(np.abs(source @ matrix - target)))
    if residual > tolerance:
        return None
    matrix[np.abs(matrix) < tolerance] = 0.0
    return SourcePrintIntertwiner(
        matrix=matrix,
        source_rank=source_rank,
        target_rank=target_rank,
        residual=residual,
        condition_number=condition_number,
    )


def solve_canonical_rank1_print_scalar(
    current_orbit_columns: Sequence[Sequence[complex | float]],
    target_orbit_columns: Sequence[Sequence[complex | float]],
    *,
    orbit_representative_indices: Sequence[int],
    tolerance: float = 1e-10,
    condition_limit: float = 1e10,
) -> CanonicalRankOnePrintScalar | None:
    """Solve one relative scalar shared by all nonzero canonical orbits.

    Callers first select Source canonical child factors and construct the
    current emitted column ``E`` and aligned child-print target ``T`` in the
    same carrier and row gauge.  This helper then solves ``E_o c = T_o`` for
    every nonzero orbit independently and accepts the result only when one
    real scalar is shared by every orbit.  Zero/zero orbits are ignored;
    zero/nonzero mismatches, duplicate orbit identities, rank failures, and
    inconsistent scalars fail closed.

    The returned value is relative to the current emitted column.  It is not
    the absolute child-print alignment matrix, so applying it does not
    double-apply a sign or remix already present in the faithful column.
    """

    current = tuple(np.asarray(column, dtype=complex) for column in current_orbit_columns)
    target = tuple(np.asarray(column, dtype=complex) for column in target_orbit_columns)
    orbits = tuple(int(value) for value in orbit_representative_indices)
    if len(current) != len(target) or len(current) != len(orbits) or not current:
        return None
    if len(set(orbits)) != len(orbits):
        return None

    contributing: list[int] = []
    scalars: list[float] = []
    residual = 0.0
    for orbit, source_column, target_column in zip(orbits, current, target, strict=True):
        if source_column.shape != target_column.shape or source_column.size == 0:
            return None
        source_flat = source_column.reshape(-1, 1)
        target_flat = target_column.reshape(-1, 1)
        source_scale = float(np.max(np.abs(source_flat)))
        target_scale = float(np.max(np.abs(target_flat)))
        if source_scale <= float(tolerance) and target_scale <= float(tolerance):
            continue
        if source_scale <= float(tolerance) or target_scale <= float(tolerance):
            return None
        solved = solve_source_print_intertwiner(
            source_flat,
            target_flat,
            tolerance=float(tolerance),
            condition_limit=float(condition_limit),
        )
        if solved is None or solved.matrix.shape != (1, 1):
            return None
        value = complex(solved.matrix[0, 0])
        scale = max(1.0, abs(value.real))
        if abs(value.imag) > float(tolerance) * scale:
            return None
        contributing.append(int(orbit))
        scalars.append(float(value.real))
        residual = max(residual, float(solved.residual))

    if not scalars:
        return None
    scalar = scalars[0]
    for value in scalars[1:]:
        scale = max(1.0, abs(scalar), abs(value))
        if abs(value - scalar) > float(tolerance) * scale:
            return None
    return CanonicalRankOnePrintScalar(
        scalar=float(scalar),
        contributing_orbit_indices=tuple(contributing),
        orbit_scalars=tuple(scalars),
        residual=float(residual),
    )


def site_print_basis_intertwiner(
    decoder: Any,
    *,
    parent_sg: int,
    parent_wyckoff_row: Any,
    pg_irrep: int,
    vector_setting: int,
    representative_operation_record: Sequence[int],
    component_start: int | None = None,
    expected_canonical_columns: Sequence[Sequence[complex | float]] | None = None,
    tolerance: float = 1e-10,
    condition_limit: float = 1e10,
) -> SourcePrintIntertwiner | None:
    """Return the Source site-print gauge induced by the canonicalizer.

    The representative operation maps presentation/CIF coordinates to the
    canonical Source representative.  Its conjugation of the canonical site
    stabilizer selects another Source ``iwyckoff_pg_vector_basis`` setting.
    The returned matrix ``J`` is the unique nondegenerate solution of
    ``A_c U_presentation J = U_canonical``.  No determinant or positive-lead
    sign heuristic is used; unsupported or degenerate settings fail closed.
    """

    record = tuple(int(value) for value in representative_operation_record)
    if len(record) != 5 or int(record[3]) == 0:
        raise ValueError(
            f"invalid representative operation record: {representative_operation_record!r}"
        )
    if int(vector_setting) not in (1, 2):
        raise ValueError(f"unsupported vector setting: {vector_setting}")

    canonical_records = tuple(
        tuple(int(value) for value in item)
        for item in decoder.wyc_pg_elements_records(int(parent_sg), parent_wyckoff_row)
    )
    canonical_columns = site_vector_print_columns(
        decoder,
        site_pg=int(parent_wyckoff_row.site_pg),
        pg_irrep=int(pg_irrep),
        site_operation_records=canonical_records,
        vector_setting=int(vector_setting),
    )
    if component_start is not None:
        canonical_columns = tuple(
            item for item in canonical_columns
            if int(item.component_start) == int(component_start)
        )
    if not canonical_columns:
        return None

    canonicalizer = _record_affine(decoder, int(parent_sg), record)  # type: ignore[arg-type]
    inverse_matrix = _matrix_inverse(canonicalizer[0])
    inverse_translation = tuple(
        -value
        for value in _row_multiply(canonicalizer[1], inverse_matrix)
    )
    inverse = (inverse_matrix, inverse_translation)  # type: ignore[arg-type]
    presentation_records: list[OperationRecord] = []
    try:
        for site_record in canonical_records:
            conjugated = _compose_affine(
                _compose_affine(
                    canonicalizer,
                    _record_affine(decoder, int(parent_sg), site_record),
                ),
                inverse,
            )
            presentation_records.append(
                _record_from_affine(decoder, int(parent_sg), conjugated)
            )
        presentation_columns = site_vector_print_columns(
            decoder,
            site_pg=int(parent_wyckoff_row.site_pg),
            pg_irrep=int(pg_irrep),
            site_operation_records=presentation_records,
            vector_setting=int(vector_setting),
        )
        if component_start is not None:
            presentation_columns = tuple(
                item for item in presentation_columns
                if int(item.component_start) == int(component_start)
            )
    except (KeyError, ValueError):
        return None
    if len(presentation_columns) != len(canonical_columns):
        return None

    canonical = np.asarray([item.vector for item in canonical_columns], dtype=complex).T
    presentation = np.asarray([item.vector for item in presentation_columns], dtype=complex).T
    if expected_canonical_columns is not None:
        expected = np.asarray(expected_canonical_columns, dtype=complex)
        if expected.shape != canonical.shape or not np.allclose(
            expected,
            canonical,
            atol=float(tolerance),
            rtol=0.0,
        ):
            return None
    pml_to_cinter = tuple(
        tuple(Fraction(value) for value in row)
        for row in decoder.pml_to_cinter_matrix(int(parent_sg))
    )
    cinter_to_pml = tuple(
        tuple(Fraction(value) for value in row)
        for row in decoder.cinter_to_pml_matrix(int(parent_sg))
    )
    cinter_rotation = _matrix_multiply(
        _matrix_multiply(cinter_to_pml, canonicalizer[0]),
        pml_to_cinter,
    )
    # Source site-vector print columns use the internal-conventional component
    # frame, whereas operation records use PML direct coordinates.
    rotation = np.asarray(cinter_rotation, dtype=float).T
    if int(vector_setting) == 2:
        rotation = rotation * float(round(np.linalg.det(rotation)))
    transported = rotation @ presentation
    return solve_source_print_intertwiner(
        transported,
        canonical,
        tolerance=float(tolerance),
        condition_limit=float(condition_limit),
    )


def presentation_carrier_coefficients(
    decoder: Any,
    *,
    gid: int,
    presentation_records: Iterable[Sequence[int]],
    presentation_case: Any,
    direction_matrix: Sequence[Sequence[complex | float]],
    source_vectors: Sequence[Sequence[complex | float]],
) -> np.ndarray:
    """Evaluate ``W D_k(p_i) V`` for every presented occurrence.

    ``source_vectors`` is ``(full_dim, family_columns)`` and retains the exact
    project family/component identity.  ``direction_matrix`` is the already
    selected OPD ``W`` with shape ``(printed_rows, full_dim)``.  The output has
    shape ``(occurrences, printed_rows, family_columns)``.  In particular, the
    target Case is explicit, so equivalent ``+k``/``-k`` representatives do
    not silently share a real-carrier gauge.
    """

    vectors = np.asarray(source_vectors, dtype=complex)
    direction = np.asarray(direction_matrix, dtype=complex)
    if vectors.ndim != 2:
        raise ValueError(f"expected source vector matrix, got shape {vectors.shape}")
    if direction.ndim != 2:
        raise ValueError(f"expected direction matrix, got shape {direction.shape}")
    if direction.shape[1] != vectors.shape[0]:
        raise ValueError(
            f"carrier dimension mismatch: direction {direction.shape[1]} != source {vectors.shape[0]}"
        )
    out: list[np.ndarray] = []
    for values in presentation_records:
        record = tuple(int(value) for value in values)
        matrix = decoder._bridge_irrep_matrix_for_record(int(gid), record, presentation_case)
        if matrix.shape[1] != vectors.shape[0]:
            raise ValueError(f"irrep/source dimension mismatch: {matrix.shape} vs {vectors.shape}")
        out.append(direction @ np.asarray(matrix, dtype=complex) @ vectors)
    if not out:
        return np.zeros((0, direction.shape[0], vectors.shape[1]), dtype=complex)
    return np.stack(out, axis=0)


def factor_resolved_parent_site_vectors(
    decoder: Any,
    *,
    parent_site_pg: int,
    parent_pg_irrep: int,
    parent_site_records: Iterable[Sequence[int]],
    occurrence_records: Iterable[Sequence[int]],
    factors: Iterable[PresentationTransportFactor],
    parent_site_vectors: Sequence[Sequence[complex | float]],
    child_subduction_basis: Sequence[Sequence[complex | float]] | None = None,
    vector_setting: int,
) -> np.ndarray:
    """Expand a parent site carrier on every exact presentation factor.

    For a factor ``q0 h = qi s L`` this returns

    ``P(qi) U_parent S(s) J_child``.

    ``U_parent`` has shape ``(3, parent_site_irrep_dim)`` and ``J_child`` is
    the Source subduction basis in that parent-site carrier.  Omitting
    ``J_child`` uses the identity.  ``P`` is the same point-operation action
    as ``project_vector_``; pseudo vectors additionally receive its
    determinant.  The output shape is ``(factor, 3, child_columns)``.

    Factor rows are deliberately not collapsed by matched occurrence.  A
    child orbit stabilizer can produce several exact factorizations of the
    same occurrence, and their ``S(s)`` actions are part of the printed branch
    topology rather than duplicate numerical rows.
    """

    if int(vector_setting) not in (1, 2):
        raise ValueError(f"unsupported vector setting: {vector_setting}")
    occurrences = tuple(tuple(int(value) for value in record) for record in occurrence_records)
    factor_rows = tuple(factors)
    vectors = np.asarray(parent_site_vectors, dtype=complex)
    if vectors.ndim != 2 or vectors.shape[0] != 3:
        raise ValueError(f"expected parent site vectors (3,dim), got {vectors.shape}")
    parent_records = tuple(tuple(int(value) for value in record) for record in parent_site_records)
    if child_subduction_basis is None:
        child_basis = np.eye(vectors.shape[1], dtype=complex)
    else:
        child_basis = np.asarray(child_subduction_basis, dtype=complex)
    if child_basis.ndim != 2 or child_basis.shape[0] != vectors.shape[1]:
        raise ValueError(
            f"site/subduction dimension mismatch: {vectors.shape} vs {child_basis.shape}"
        )

    out: list[np.ndarray] = []
    for factor in factor_rows:
        occurrence_index = int(factor.matched_occurrence_index)
        if not 0 <= occurrence_index < len(occurrences):
            raise IndexError(f"factor occurrence index out of range: {occurrence_index}")
        occurrence = occurrences[occurrence_index]
        point = np.asarray(decoder.iso.point_ops[int(occurrence[4]) - 1], dtype=float)
        # project_vector_ treats stored point-operation columns as output
        # components, hence the column-vector action is point.T.
        physical = point.T
        if int(vector_setting) == 2:
            physical = physical * float(round(np.linalg.det(point)))
        site_action = site_stabilizer_irrep_matrix(
            decoder,
            site_pg=int(parent_site_pg),
            pg_irrep=int(parent_pg_irrep),
            parent_site_records=parent_records,
            stabilizer_record=factor.stabilizer_record,
        )
        if site_action.shape[0] != vectors.shape[1]:
            raise ValueError(
                f"site action dimension mismatch: {site_action.shape} vs {vectors.shape}"
            )
        out.append(physical @ vectors @ np.asarray(site_action, dtype=complex) @ child_basis)
    if not out:
        return np.zeros((0, 3, child_basis.shape[1]), dtype=complex)
    return np.stack(out, axis=0)


def factor_target_row_indices(
    decoder: Any,
    *,
    child_sg: int,
    factors: Iterable[PresentationTransportFactor],
    raw_child_points: Iterable[Sequence[Fraction | float | int]],
    occurrence_records: Iterable[Sequence[int]] | None = None,
    tolerance: float = 1e-8,
) -> tuple[int, ...]:
    """Map every factor's child action back to one raw child-layout row.

    ``raw_child_points`` must use the same raw order as the occurrence records
    passed to :func:`factor_presentation_occurrences`, but its coordinates are
    in the child internal-conventional (``cinter``) setting.  For each factor
    this helper applies the Source child operation ``child_operation_index``
    to its orbit representative ``q0``.  The operation is evaluated in child
    PML coordinates.  Candidate rows are converted to the same PML setting
    before the periodic ``mod 1`` match.  Folding in cinter is not equivalent
    for centered settings because its integer cell is larger than the Source
    primitive translation lattice.

    The returned tuple is factor ordered.  It is intentionally separate from
    ``matched_occurrence_index``: the latter identifies ``qi`` modulo a child
    lattice translation in ``q0 h = qi s L`` and can collapse distinct raw
    supercell rows, while this function resolves their displayed coordinates.
    Optional ``occurrence_records`` expose the symbolic Source branch already
    selected by the exact factorization ``q0 h = qi s L``.  A globally unique
    matching full record is authoritative even when substituted site
    parameters create accidental coordinate coincidences.  Coordinates are
    used only to distinguish several translation copies of that same record.
    ``lattice_coordinates`` are not applied again here.
    """

    points = tuple(
        tuple(Fraction(str(value)) if isinstance(value, float) else Fraction(value) for value in point)
        for point in raw_child_points
    )
    if any(len(point) != 3 for point in points):
        raise ValueError("raw child points must be coordinate triplets")
    records = (
        tuple(tuple(int(value) for value in record) for record in occurrence_records)
        if occurrence_records is not None
        else None
    )
    if records is not None:
        if len(records) != len(points):
            raise ValueError(
                f"point/occurrence count mismatch: {len(points)} != {len(records)}"
            )
        if any(len(record) != 5 or int(record[3]) == 0 for record in records):
            raise ValueError("occurrence records must be valid five-field records")
    child_records = tuple(decoder.generate_space_group_records(int(child_sg)))

    def fold(values: Sequence[Fraction]) -> tuple[Fraction, Fraction, Fraction]:
        return tuple(value % 1 for value in values)  # type: ignore[return-value]

    def periodic_close(left: Sequence[Fraction], right: Sequence[Fraction]) -> bool:
        for lhs, rhs in zip(left, right):
            delta = float((lhs - rhs) % 1)
            distance = min(delta, 1.0 - delta)
            if distance > float(tolerance):
                return False
        return True

    folded_points_pml = tuple(
        fold(
            decoder.xyz_change_setting_point(
                int(child_sg),
                "cinter",
                "pml",
                point,
            )
        )
        for point in points
    )
    folded_points_cinter = tuple(fold(point) for point in points)
    out: list[int] = []
    for factor in factors:
        representative_index = int(factor.orbit_representative_index)
        if not 0 <= representative_index < len(points):
            raise IndexError(f"factor representative index out of range: {representative_index}")
        child_index = int(factor.child_operation_index)
        if not 0 <= child_index < len(child_records):
            raise IndexError(f"factor child operation index out of range: {child_index}")
        record = tuple(int(value) for value in child_records[child_index])
        denominator = int(record[3])
        if denominator == 0:
            raise ValueError(f"zero child operation denominator: {record}")
        point_pml = decoder.xyz_change_setting_point(
            int(child_sg),
            "cinter",
            "pml",
            points[representative_index],
        )
        rotated = decoder.vrot_fraction(int(child_sg), int(record[4]), point_pml)
        transformed_pml = tuple(
            rotated[axis] + Fraction(int(record[axis]), denominator)
            for axis in range(3)
        )
        transformed = fold(transformed_pml)
        matches = [
            index
            for index, candidate in enumerate(folded_points_pml)
            if periodic_close(transformed, candidate)
        ]
        if len(matches) > 1:
            # PML periodicity correctly identifies primitive-lattice
            # equivalence, but it can merge distinct centering copies when
            # the supplied raw layout explicitly retains them.  Re-evaluate
            # only those equivalent candidates in child cinter coordinates.
            transformed_cinter = fold(
                decoder.xyz_change_setting_point(
                    int(child_sg),
                    "pml",
                    "cinter",
                    transformed_pml,
                )
            )
            cinter_matches = [
                index
                for index in matches
                if periodic_close(transformed_cinter, folded_points_cinter[index])
            ]
            if len(cinter_matches) == 1:
                matches = cinter_matches
        if records is not None:
            matched_occurrence_index = int(factor.matched_occurrence_index)
            if not 0 <= matched_occurrence_index < len(records):
                raise IndexError(
                    f"factor occurrence index out of range: {matched_occurrence_index}"
                )
            matched_record = records[matched_occurrence_index]
            record_matches = [
                index for index, record_value in enumerate(records)
                if record_value == matched_record
            ]
            if len(record_matches) == 1:
                matches = record_matches
            elif len(record_matches) > 1:
                coordinate_record_matches = [
                    index for index in matches if index in set(record_matches)
                ]
                if len(coordinate_record_matches) == 1:
                    matches = coordinate_record_matches
                elif matched_occurrence_index in coordinate_record_matches:
                    # Exact factorization retains the raw Source branch even
                    # when substituted parameters make several rows identical
                    # in both coordinate and full-operation record.
                    matches = [matched_occurrence_index]
        if len(matches) != 1:
            raise KeyError(
                "child operation did not identify one raw row: "
                f"representative={representative_index} child_op={child_index} "
                f"point={transformed} matches={matches}"
            )
        out.append(int(matches[0]))
    return tuple(out)


def child_print_aligned_subduction_basis(
    decoder: Any,
    *,
    parent_sg: int,
    representative_record: Sequence[int],
    child_basis_rows_in_parent_cinter: Sequence[Sequence[Fraction | float | int]],
    parent_site_vectors: Sequence[Sequence[complex | float]],
    projector_basis: Sequence[Sequence[complex | float]],
    child_print_vectors: Sequence[Sequence[complex | float]],
    vector_setting: int,
    parent_vector_setting: str = "cinter",
    tolerance: float = 1e-10,
) -> ChildPrintSubductionBasis | None:
    """Orient a child subduction basis by Source's child print columns.

    Character projectors determine an image but leave an arbitrary basis in
    that image, including an arbitrary sign in rank one.  This helper carries
    ``U_parent J_projector`` through the actual parent representative and the
    selected child basis, then solves

    ``U_parent_in_child C = U_child_print``

    and returns ``J_projector C``.  It returns ``None`` when that change is not
    unique or the Source column spans differ.  No positive-lead or Web label
    convention is used.

    ``parent_site_vectors`` and ``child_print_vectors`` are column matrices
    with physical axis first.  ``child_basis_rows_in_parent_cinter`` satisfies
    ``x_parent_cinter = x_child * basis``.  Parent site-vector tables normally
    use cinter axes; callers with a PML surface may set
    ``parent_vector_setting='pml'`` explicitly.
    """

    if int(vector_setting) not in (1, 2):
        raise ValueError(f"unsupported vector setting: {vector_setting}")
    setting = str(parent_vector_setting).strip().lower()
    if setting not in {"pml", "cinter"}:
        raise ValueError(f"unsupported parent vector setting: {parent_vector_setting!r}")
    record = tuple(int(value) for value in representative_record)
    if len(record) != 5:
        raise ValueError(f"expected 5-int representative record, got {record}")
    parent_vectors = np.asarray(parent_site_vectors, dtype=complex)
    projector = np.asarray(projector_basis, dtype=complex)
    child_vectors = np.asarray(child_print_vectors, dtype=complex)
    if parent_vectors.ndim != 2 or parent_vectors.shape[0] != 3:
        raise ValueError(f"expected parent vectors (3,dim), got {parent_vectors.shape}")
    if projector.ndim != 2 or projector.shape[0] != parent_vectors.shape[1]:
        raise ValueError(
            f"parent/projector dimension mismatch: {parent_vectors.shape} vs {projector.shape}"
        )
    if child_vectors.shape != (3, projector.shape[1]):
        raise ValueError(
            f"child print column mismatch: {child_vectors.shape} != {(3, projector.shape[1])}"
        )
    child_basis = np.asarray(
        [[float(Fraction(str(value)) if isinstance(value, float) else Fraction(value)) for value in row]
         for row in child_basis_rows_in_parent_cinter],
        dtype=float,
    )
    if child_basis.shape != (3, 3):
        raise ValueError(f"expected child basis 3x3, got {child_basis.shape}")
    if abs(float(np.linalg.det(child_basis))) <= tolerance:
        raise ValueError(f"singular child presentation basis: {child_basis}")

    point = np.asarray(decoder.iso.point_ops[int(record[4]) - 1], dtype=float)
    pml_to_cinter = np.asarray(decoder.pml_to_cinter_matrix(int(parent_sg)), dtype=float)
    if setting == "cinter":
        cinter_to_pml = np.linalg.inv(pml_to_cinter)
        setting_transport = np.linalg.inv(child_basis)
        row_transport = cinter_to_pml @ point @ pml_to_cinter @ setting_transport
    else:
        setting_transport = pml_to_cinter @ np.linalg.inv(child_basis)
        row_transport = point @ setting_transport
    handedness = 1.0
    if int(vector_setting) == 2:
        handedness = float(round(np.linalg.det(point)))
        setting_det = float(np.linalg.det(setting_transport))
        if abs(setting_det) <= tolerance:
            raise ValueError("singular parent-to-child vector transport")
        handedness *= 1.0 if setting_det > 0 else -1.0
    transported = handedness * row_transport.T @ parent_vectors @ projector
    change = solve_source_print_intertwiner(
        transported,
        child_vectors,
        tolerance=float(tolerance),
    )
    if change is None:
        return None
    aligned = projector @ change.matrix
    residual = float(np.max(np.abs(transported @ change.matrix - child_vectors)))
    return ChildPrintSubductionBasis(
        basis=aligned,
        print_change=change.matrix,
        transported_columns=transported,
        residual=residual,
    )


def child_print_subduction_embedding(
    decoder: Any,
    *,
    parent_sg: int,
    representative_record: Sequence[int],
    child_basis_rows_in_parent_cinter: Sequence[Sequence[Fraction | float | int]],
    parent_site_vectors: Sequence[Sequence[complex | float]],
    projector_basis: Sequence[Sequence[complex | float]],
    child_print_vectors: Sequence[Sequence[complex | float]],
    vector_setting: int,
    parent_vector_setting: str = "cinter",
    tolerance: float = 1e-10,
) -> ChildPrintSubductionEmbedding | None:
    """Embed a subduced image in all repeated Source child print copies.

    Unlike :func:`child_print_aligned_subduction_basis`, this helper does not
    require one child print column per projector column.  It solves

    ``U_child_repeated C = U_parent_in_child``

    using the full independent repeated-copy basis.  This is the faithful
    boundary for cases such as SG119 where a rank-one parent branch maps to a
    unique linear combination of three child PG1 copies rather than to any
    single x/y/z copy.  The projector basis itself is not re-signed; the signs
    and remix live in ``print_coefficients``.
    """

    if int(vector_setting) not in (1, 2):
        raise ValueError(f"unsupported vector setting: {vector_setting}")
    setting = str(parent_vector_setting).strip().lower()
    if setting not in {"pml", "cinter"}:
        raise ValueError(f"unsupported parent vector setting: {parent_vector_setting!r}")
    record = tuple(int(value) for value in representative_record)
    if len(record) != 5:
        raise ValueError(f"expected 5-int representative record, got {record}")
    parent_vectors = np.asarray(parent_site_vectors, dtype=complex)
    projector = np.asarray(projector_basis, dtype=complex)
    child_vectors = np.asarray(child_print_vectors, dtype=complex)
    if parent_vectors.ndim != 2 or parent_vectors.shape[0] != 3:
        raise ValueError(f"expected parent vectors (3,dim), got {parent_vectors.shape}")
    if projector.ndim != 2 or projector.shape[0] != parent_vectors.shape[1]:
        raise ValueError(
            f"parent/projector dimension mismatch: {parent_vectors.shape} vs {projector.shape}"
        )
    if child_vectors.ndim != 2 or child_vectors.shape[0] != 3 or child_vectors.shape[1] == 0:
        raise ValueError(f"expected nonempty child print matrix (3,copies), got {child_vectors.shape}")
    child_rank = int(np.linalg.matrix_rank(child_vectors, tol=float(tolerance)))
    if child_rank != child_vectors.shape[1]:
        return None
    child_basis = np.asarray(
        [[float(Fraction(str(value)) if isinstance(value, float) else Fraction(value)) for value in row]
         for row in child_basis_rows_in_parent_cinter],
        dtype=float,
    )
    if child_basis.shape != (3, 3):
        raise ValueError(f"expected child basis 3x3, got {child_basis.shape}")
    if abs(float(np.linalg.det(child_basis))) <= tolerance:
        raise ValueError(f"singular child presentation basis: {child_basis}")

    point = np.asarray(decoder.iso.point_ops[int(record[4]) - 1], dtype=float)
    pml_to_cinter = np.asarray(decoder.pml_to_cinter_matrix(int(parent_sg)), dtype=float)
    if setting == "cinter":
        cinter_to_pml = np.linalg.inv(pml_to_cinter)
        setting_transport = np.linalg.inv(child_basis)
        row_transport = cinter_to_pml @ point @ pml_to_cinter @ setting_transport
    else:
        setting_transport = pml_to_cinter @ np.linalg.inv(child_basis)
        row_transport = point @ setting_transport
    handedness = 1.0
    if int(vector_setting) == 2:
        handedness = float(round(np.linalg.det(point)))
        setting_det = float(np.linalg.det(setting_transport))
        if abs(setting_det) <= tolerance:
            raise ValueError("singular parent-to-child vector transport")
        handedness *= 1.0 if setting_det > 0 else -1.0
    transported = handedness * row_transport.T @ parent_vectors @ projector
    combined_rank = int(
        np.linalg.matrix_rank(np.column_stack((child_vectors, transported)), tol=float(tolerance))
    )
    if combined_rank != child_rank:
        return None
    coefficients, _residuals, _rank, _singular = np.linalg.lstsq(
        child_vectors,
        transported,
        rcond=None,
    )
    residual = float(np.max(np.abs(child_vectors @ coefficients - transported)))
    if residual > float(tolerance):
        return None
    coefficients[np.abs(coefficients) < tolerance] = 0.0
    return ChildPrintSubductionEmbedding(
        basis=projector,
        print_coefficients=coefficients,
        transported_columns=transported,
        child_print_columns=child_vectors,
        residual=residual,
    )


def printed_project_block_sources(
    decoder: Any,
    *,
    site_pg: int,
    pg_irrep: int,
    site_operation_records: Iterable[Sequence[int]],
    vector_setting: int,
    full_dim: int,
    project_basis_active_values: Iterable[Sequence[float | int]],
    printed_blocks: Iterable[Any],
    family_stride: int = 144,
) -> tuple[PrintedProjectBlockSource, ...]:
    """Decode the Source columns for every printed project block.

    The project return stores each family as ``carrier * 3 + site_column`` in
    a 144-value family stride.  Runtime print topology over Source records maps
    a project family to an emitted block and vector component.  Combining those
    two Source surfaces yields ``project_matrix(full_dim, site_columns)`` and
    ``site_vectors(3, site_columns)`` once.  Both direct tensor generation and
    factor-aligned child targets must consume this same decode.
    """

    if int(vector_setting) not in (1, 2):
        raise ValueError(f"unsupported vector setting: {vector_setting}")
    full_dim = int(full_dim)
    if full_dim <= 0:
        raise ValueError(f"invalid carrier dimension: {full_dim}")
    stride = int(family_stride)
    if stride <= 0:
        raise ValueError(f"invalid project family stride: {stride}")
    site_records = tuple(tuple(int(value) for value in record) for record in site_operation_records)
    source_columns = site_vector_print_columns(
        decoder,
        site_pg=int(site_pg),
        pg_irrep=int(pg_irrep),
        site_operation_records=site_records,
        vector_setting=int(vector_setting),
    )
    active = tuple(tuple(item) for item in project_basis_active_values)
    blocks = tuple(printed_blocks)

    def field(item: Any, name: str, default: int = 0) -> int:
        if isinstance(item, dict):
            return int(item.get(name, default))
        return int(getattr(item, name, default))

    def optional_field(item: Any, name: str) -> int | None:
        value = item.get(name) if isinstance(item, dict) else getattr(item, name, None)
        return None if value is None else int(value)

    out: list[PrintedProjectBlockSource] = []
    for block_index, block in enumerate(blocks):
        target = field(block, "target_vector_rep", int(pg_irrep))
        if target != int(pg_irrep):
            raise ValueError(f"printed target {target} does not match site irrep {pg_irrep}")
        family = field(block, "project_family_index")
        vector_component = field(block, "vector_component")
        component_columns = tuple(
            item for item in source_columns if int(item.component_start) == vector_component + 1
        )
        if not component_columns:
            raise KeyError(
                f"no Source site columns for PG{site_pg} irrep {pg_irrep} component {vector_component}"
            )
        site_vectors = np.asarray([item.vector for item in component_columns], dtype=complex).T
        project = np.zeros((full_dim, len(component_columns)), dtype=complex)
        for item in active:
            if len(item) < 2:
                continue
            index = int(item[0])
            if index // stride != family:
                continue
            local = index % stride
            carrier = local // 3
            site_column = local % 3
            if carrier < full_dim and site_column < len(component_columns):
                project[carrier, site_column] = complex(float(item[1]))
        if not np.any(np.abs(project) > 1e-12):
            raise KeyError(f"project family {family} has no active columns for printed block {block_index + 1}")
        out.append(
            PrintedProjectBlockSource(
                printed_block_index=field(block, "printed_block_index", block_index + 1),
                project_family_index=int(family),
                vector_component=int(vector_component),
                target_vector_rep=int(target),
                project_matrix=project,
                site_vectors=site_vectors,
                selected_candidate_index=optional_field(block, "selected_candidate_index"),
                raw_star_arm=optional_field(block, "raw_star_arm"),
            )
        )
    return tuple(out)


def printed_project_blocks_from_trace(
    decoder: Any,
    *,
    gid: int,
    site_pg: int,
    vector_setting: int,
    project_return_basis_type1: Iterable[dict[str, Any]],
    pg_irrep: int | None = None,
) -> tuple[MinimalPrintedProjectBlock, ...]:
    """Return the minimal DISPLAY DISTORTION print-block ordering.

    Target irreps follow ``project_``/subduction order in the supplied trace,
    not the order of the three site-vector slots.  Within each target, project
    family is the outer loop and repeated site-vector component is the inner
    loop.  Only families actually returned by ``project_`` are emitted.

    This minimal runtime topology exposes only the identities required to
    build Assembled Source tensors.
    """

    if int(vector_setting) not in (1, 2):
        raise ValueError(f"unsupported vector setting: {vector_setting}")
    reps = tuple(int(value) for value in decoder.site_vector_reps(int(site_pg)))
    start = (int(vector_setting) - 1) * 3
    vector_reps = reps[start:start + 3]
    little = decoder.little_record_by_gid(int(gid))
    full_dim = int(little.full_dim)
    star_count = max(1, int(decoder.little_transform_block_count(int(gid))))
    star_block_dim = full_dim // star_count if full_dim > 0 and full_dim % star_count == 0 else 0
    out: list[MinimalPrintedProjectBlock] = []
    seen_targets: set[int] = set()
    for item in project_return_basis_type1:
        if not isinstance(item, dict) or int(item.get("gid", 0)) != int(gid):
            continue
        target = int(item.get("pg_irrep", 0))
        if pg_irrep is not None and target != int(pg_irrep):
            continue
        if target in seen_targets:
            continue
        seen_targets.add(target)
        components = [index for index, value in enumerate(vector_reps) if int(value) == target]
        if not components:
            continue
        selected = tuple(int(value) for value in item.get("selected_indices", ()))
        if selected:
            family_indices = tuple(range(len(selected)))
        else:
            stride = int(item.get("family_stride", 144) or 144)
            family_indices = tuple(
                sorted(
                    {
                        int(active[0]) // stride
                        for active in item.get("active_values", ())
                        if isinstance(active, (list, tuple)) and len(active) >= 2
                    }
                )
            )
        for family in family_indices:
            selected_candidate = selected[family] if family < len(selected) else None
            raw_star_arm = (
                int(selected_candidate) // star_block_dim
                if selected_candidate is not None and star_block_dim > 0
                else None
            )
            for target_occurrence, vector_component in enumerate(components):
                out.append(
                    MinimalPrintedProjectBlock(
                        printed_block_index=len(out) + 1,
                        target_occurrence=int(target_occurrence),
                        vector_component=int(vector_component),
                        target_vector_rep=int(target),
                        project_family_index=int(family),
                        selected_candidate_index=selected_candidate,
                        raw_star_arm=raw_star_arm,
                    )
                )
    return tuple(out)


def factor_carrier_coefficients_on_target_rows(
    decoder: Any,
    *,
    gid: int,
    factors: Iterable[PresentationTransportFactor],
    target_row_indices: Sequence[int],
    target_occurrence_records: Iterable[Sequence[int]],
    presentation_case: Any,
    source_case: Any | None = None,
    source_factors: Iterable[PresentationTransportFactor] | None = None,
    direction_matrix: Sequence[Sequence[complex | float]],
    source_vectors: Sequence[Sequence[complex | float]],
) -> FactorCarrierPullback:
    """Pull ``W D(p_f) V`` back to each factor's raw target-row gauge.

    ``factor_target_row_indices`` folds ``q0 h`` onto a raw child-layout row.
    Distinct factors can reach the same row through different child-lattice
    translations.  Their physical site action may be identical while
    ``D(p_f)`` differs by a Bloch phase.  For each factor this helper computes

    ``Q_f = D_source(q_target) inverse(D_source(p_f))``

    and returns ``W Q_f D_presentation(p_f) V``.  With one Case this reduces
    to ``W D(q_target) V``.  ``Q_f`` is a full real-carrier
    matrix; no scalar sign, first-factor selection, or averaging is used.
    ``source_factors`` may provide the untransported ``p_f`` records when the
    presentation factors include an additional representative operation.  It
    must preserve the factor order and orbit identity; when omitted, the
    presentation factors are also used as source factors.
    """

    factor_rows = tuple(factors)
    source_factor_rows = factor_rows if source_factors is None else tuple(source_factors)
    row_indices = tuple(int(value) for value in target_row_indices)
    records = tuple(tuple(int(value) for value in record) for record in target_occurrence_records)
    if len(row_indices) != len(factor_rows):
        raise ValueError(f"factor/target-row count mismatch: {len(factor_rows)} != {len(row_indices)}")
    if len(source_factor_rows) != len(factor_rows):
        raise ValueError(
            f"source/presentation factor count mismatch: {len(source_factor_rows)} != {len(factor_rows)}"
        )
    identity_fields = tuple(
        name
        for name in PresentationTransportFactor.__dataclass_fields__
        if name != "presentation_record"
    )
    for factor_index, (factor, source_factor) in enumerate(
        zip(factor_rows, source_factor_rows, strict=True)
    ):
        if any(
            getattr(factor, name) != getattr(source_factor, name)
            for name in identity_fields
        ):
            raise ValueError(
                f"source/presentation factor identity mismatch at index {factor_index}"
            )
    direction = np.asarray(direction_matrix, dtype=complex)
    vectors = np.asarray(source_vectors, dtype=complex)
    if direction.ndim != 2 or vectors.ndim != 2 or direction.shape[1] != vectors.shape[0]:
        raise ValueError(f"carrier shape mismatch: direction={direction.shape} source={vectors.shape}")
    coefficients: list[np.ndarray] = []
    pullbacks: list[np.ndarray] = []
    target_records: list[OperationRecord] = []
    pullback_case = presentation_case if source_case is None else source_case

    for factor, source_factor, row_index in zip(
        factor_rows, source_factor_rows, row_indices, strict=True
    ):
        if not 0 <= row_index < len(records):
            raise IndexError(f"factor target row out of range: {row_index}")
        target_record = records[row_index]
        source_presentation = np.asarray(
            decoder._bridge_irrep_matrix_for_record(
                int(gid),
                source_factor.presentation_record,
                pullback_case,
            ),
            dtype=complex,
        )
        source_target = np.asarray(
            decoder._bridge_irrep_matrix_for_record(
                int(gid),
                target_record,
                pullback_case,
            ),
            dtype=complex,
        )
        presentation = np.asarray(
            decoder._bridge_irrep_matrix_for_record(
                int(gid),
                factor.presentation_record,
                presentation_case,
            ),
            dtype=complex,
        )
        if (
            source_presentation.shape != source_target.shape
            or source_presentation.shape != presentation.shape
            or source_presentation.shape[0] != direction.shape[1]
        ):
            raise ValueError(
                f"factor carrier shape mismatch: source_p={source_presentation.shape} "
                f"source_q={source_target.shape} presentation={presentation.shape} "
                f"direction={direction.shape}"
            )
        pullback = source_target @ np.linalg.inv(source_presentation)
        pullback[np.abs(pullback) < 1e-12] = 0.0
        coefficient = direction @ pullback @ presentation @ vectors
        coefficient[np.abs(coefficient) < 1e-12] = 0.0
        pullbacks.append(pullback)
        coefficients.append(coefficient)
        target_records.append(target_record)
    if not coefficients:
        shape = (0, direction.shape[0], vectors.shape[1])
        matrix_shape = (0, direction.shape[1], direction.shape[1])
        return FactorCarrierPullback(
            coefficients=np.zeros(shape, dtype=complex),
            pullback_matrices=np.zeros(matrix_shape, dtype=complex),
            target_records=(),
        )
    return FactorCarrierPullback(
        coefficients=np.stack(coefficients, axis=0),
        pullback_matrices=np.stack(pullbacks, axis=0),
        target_records=tuple(target_records),
    )


def canonical_child_print_factors(
    decoder: Any,
    *,
    child_sg: int,
    factors: Iterable[PresentationTransportFactor],
) -> CanonicalChildPrintFactors:
    """Return one Source canonical print factor for every child orbit.

    Complete physical fields require all child-operation factors.  Their
    carrier and site actions can cancel a print-column sign, however, because
    gauge-equivalent fields are physically identical.  Runtime print order over
    Source records fixes the column at the canonical child representative: the
    factor whose child operation is the Source identity operation.  This helper
    verifies that operation zero is affine identity and then selects exactly
    one such factor per ``orbit_representative_index`` in Source factor order.

    The result is only a print-gauge observation surface.  It must not replace
    all-factor validation of the emitted physical mode field.
    """

    factor_rows = tuple(factors)
    child_records = tuple(
        tuple(int(value) for value in record)
        for record in decoder.generate_space_group_records(int(child_sg))
    )
    if not child_records:
        raise ValueError(f"child SG{child_sg} has no generated operations")
    identity_record = child_records[0]
    if int(identity_record[3]) == 0 or any(int(identity_record[axis]) != 0 for axis in range(3)):
        raise ValueError(f"first child operation is not translation identity: {identity_record}")
    units = (
        (Fraction(1), Fraction(0), Fraction(0)),
        (Fraction(0), Fraction(1), Fraction(0)),
        (Fraction(0), Fraction(0), Fraction(1)),
    )
    rotated = tuple(
        tuple(decoder.vrot_fraction(int(child_sg), int(identity_record[4]), unit))
        for unit in units
    )
    if rotated != units:
        raise ValueError(f"first child operation is not rotational identity: {identity_record}")

    selected: list[PresentationTransportFactor] = []
    selected_indices: list[int] = []
    seen_orbits: set[int] = set()
    expected_orbits = tuple(dict.fromkeys(int(item.orbit_representative_index) for item in factor_rows))
    for factor_index, factor in enumerate(factor_rows):
        if int(factor.child_operation_index) != 0:
            continue
        orbit = int(factor.orbit_representative_index)
        if orbit in seen_orbits:
            raise ValueError(f"multiple canonical child factors for orbit {orbit}")
        seen_orbits.add(orbit)
        selected.append(factor)
        selected_indices.append(int(factor_index))
    if tuple(int(item.orbit_representative_index) for item in selected) != expected_orbits:
        raise ValueError(
            "canonical child factors do not cover Source orbit order: "
            f"{[item.orbit_representative_index for item in selected]} != {expected_orbits}"
        )
    return CanonicalChildPrintFactors(
        factors=tuple(selected),
        source_factor_indices=tuple(selected_indices),
        orbit_representative_indices=expected_orbits,
    )


def expand_occurrence_print_columns(
    carrier_coefficients: Sequence[Sequence[Sequence[complex | float]]],
    site_vectors: Sequence[Sequence[Sequence[complex | float]]] | Sequence[Sequence[complex | float]],
) -> np.ndarray:
    """Expand carrier coefficients into occurrence-major physical columns.

    Carrier shape is ``(occurrences, printed_rows, family_columns)``.  Site
    vectors are either ``(occurrences, family_columns, 3)`` or one
    ``(family_columns, 3)`` basis broadcast to every occurrence after the
    caller has applied the child-setting/witness physical transform.  Rows are
    ordered ``occurrence, x, y, z``; columns are ordered
    ``printed_row, family_column``.  This is the matrix consumed by
    :func:`solve_source_print_intertwiner`.
    """

    carrier = np.asarray(carrier_coefficients, dtype=complex)
    vectors = np.asarray(site_vectors, dtype=complex)
    if carrier.ndim != 3:
        raise ValueError(f"expected carrier tensor rank 3, got shape {carrier.shape}")
    occurrence_count, printed_rows, family_columns = carrier.shape
    if vectors.ndim == 2:
        vectors = np.broadcast_to(vectors, (occurrence_count,) + vectors.shape)
    if vectors.shape != (occurrence_count, family_columns, 3):
        raise ValueError(
            "site vector shape mismatch: "
            f"{vectors.shape} != {(occurrence_count, family_columns, 3)}"
        )
    expanded = np.einsum("orf,ofa->oarf", carrier, vectors)
    return expanded.transpose(0, 1, 2, 3).reshape(
        occurrence_count * 3,
        printed_rows * family_columns,
    )
