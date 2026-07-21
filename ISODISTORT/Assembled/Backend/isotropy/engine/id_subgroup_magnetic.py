"""Magnetic ``id_subgroup_magnetic_`` support helpers.

This module is the next faithful-port frontier for parametric magnetic OPD
rows.  It intentionally does not try to finish the whole subgroup
identification routine in one step; the helpers here isolate the magnetic
point-code, multiplication, and candidate-generator layers that differ from
the ordinary ``id_subgroup_`` port.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from fractions import Fraction
import math
from typing import Sequence

from ISODISTORT.Assembled.Backend.exactmath import fraction_row_multiply3
from ISODISTORT.Assembled.Backend.source import magnetic as magnetic_data
from ISODISTORT.Assembled.Backend.isotropy.engine.id_subgroup import (
    _matinv_numerator_denominator,
    _matrix3_rows,
    _matmlt_flat,
    _nonpivot_count,
    _origin_after_basis,
    generator_diff_records,
    id_subgroup_basis_candidate_from_rowop,
    id_subgroup_rowop_setup_matrix,
    input_operation_records_from_inverse,
    nice_lattice2_basis,
    origin_equation_system,
    rowop_reduce_matrix,
)
from ISODISTORT.Assembled.Backend.isotropy.engine.solve_eqs import solve_eqs_mod_int_first
from ISODISTORT.Assembled.Backend.isotropy.engine.source_data import SourceData


MagneticOperationRecord = tuple[int, int, int, int, int]


@dataclass(frozen=True)
class MagneticGroupCandidate:
    magnetic_group: int
    ordinary_space_group: int
    point_group: int
    point_group_code: int
    generators: tuple[MagneticOperationRecord, ...]


@dataclass(frozen=True)
class MagneticGeneratorMapping:
    candidate_magnetic_group: int
    candidate_generator_count: int
    per_generator_match_counts: tuple[int, ...]
    per_generator_match_slots: tuple[tuple[int, ...], ...]
    selected_input_generator_slots: tuple[int, ...]
    selected_candidate_generator_ops: tuple[int, ...]


@dataclass(frozen=True)
class MagneticPointOperationCorrespondence:
    """One selected-input to canonical-BNS magnetic point-operation match."""

    input_slot: int
    input_point_op: int
    canonical_point_op: int
    selected_generator: bool


@dataclass(frozen=True)
class MagneticFullOperationCorrespondence:
    """One lifted full-Seitz correspondence in Source operation order."""

    input_slot: int
    input_record: MagneticOperationRecord
    local_input_record: MagneticOperationRecord
    transformed_input_record: MagneticOperationRecord
    canonical_record: MagneticOperationRecord


@dataclass(frozen=True)
class MagneticIdSubgroupSelection:
    """Faithful first-success identity retained by ``id_subgroup_magnetic_``."""

    base_value_count: int
    candidate_scan_ordinal: int
    mapping_ordinal: int
    basis_candidate_ordinal: int
    candidate_basis: tuple[tuple[int, int, int], ...]
    candidate_origin: tuple[int, int, int, int]
    generator_mapping: MagneticGeneratorMapping
    point_operation_correspondence: tuple[MagneticPointOperationCorrespondence, ...]
    full_operation_correspondence: tuple[MagneticFullOperationCorrespondence, ...]


@dataclass(frozen=True)
class MagneticIdSubgroupResult:
    magnetic_group: int
    ordinary_space_group: int
    basis: tuple[int, ...]
    origin: tuple[int, int, int, int]
    # Excluded from repr/equality so the historical four-field result contract
    # remains stable while faithful callers can retain the selected search path.
    selection: MagneticIdSubgroupSelection | None = field(
        default=None,
        compare=False,
        repr=False,
    )


def magnetic_point_code(magnetic_point_op: int) -> int:
    """Return ``mag_point_op_code`` for one 1-based magnetic point operation."""

    table = magnetic_data.data().table
    return int(table["mag_point_op_code"][int(magnetic_point_op) - 1])


def magnetic_nonmag_point_op(magnetic_point_op: int) -> int:
    """Return the nonmagnetic point operation paired with a magnetic one."""

    table = magnetic_data.data().table
    return int(table["mag_point_op_mag2nonmag"][int(magnetic_point_op) - 1])


def magnetic_point_mul(left: int, right: int) -> int:
    """Return the binary ``mag_point_op_mlt`` product.

    The table is right-major in the ISO common block, matching the already
    verified ``generate_spacegroup_magnetic_`` closure port.
    """

    table = magnetic_data.data().table
    return int(table["mag_point_op_mlt"][(int(right) - 1) * 144 + (int(left) - 1)])


def magnetic_group_generators(magnetic_group: int) -> tuple[MagneticOperationRecord, ...]:
    """Return native ``mag_gen`` records for one 1-based magnetic group."""

    table = magnetic_data.data().table
    group = int(magnetic_group)
    count = int(table["mag_gen_count"][group - 1])
    raw = table["mag_gen"][(group - 1) * 20 : group * 20]
    return tuple(
        tuple(int(value) for value in raw[index * 5 : index * 5 + 5])  # type: ignore[misc]
        for index in range(count)
    )


def magnetic_operation_code_sum(operations: Sequence[Sequence[int]]) -> int:
    """Return the input operation-code sum used before candidate scanning."""

    return sum(magnetic_point_code(int(record[4])) for record in operations)


def magnetic_group_candidates_for_operations(
    operations: Sequence[Sequence[int]],
    *,
    ordinary_space_group: int | None = None,
) -> tuple[MagneticGroupCandidate, ...]:
    """Return magnetic groups that pass the first point-code filter.

    ``ordinary_space_group`` is optional because the binary discovers it during
    the full basis/generator search.  It is useful for diagnostics after a
    nonmagnetic bridge has identified the ordinary subgroup.
    """

    table = magnetic_data.data().table
    target_code = magnetic_operation_code_sum(operations)
    out: list[MagneticGroupCandidate] = []
    for group in range(1, len(table["mag_space_group"]) + 1):
        ordinary = int(table["mag_space_group"][group - 1])
        if ordinary_space_group is not None and ordinary != int(ordinary_space_group):
            continue
        point_group = int(table["mag_point_group"][group - 1])
        code = int(table["mag_point_group_code"][point_group - 1])
        if code != target_code:
            continue
        out.append(
            MagneticGroupCandidate(
                magnetic_group=group,
                ordinary_space_group=ordinary,
                point_group=point_group,
                point_group_code=code,
                generators=magnetic_group_generators(group),
            )
        )
    return tuple(out)


def _magnetic_closure_correspondence(
    input_ops: Sequence[Sequence[int]],
    selected_input_slots: tuple[int, ...],
    selected_candidate_ops: tuple[int, ...],
) -> tuple[MagneticPointOperationCorrespondence, ...] | None:
    """Close a generator match and retain its ordered point-operation map."""

    input_slot_by_point_op = {int(record[4]): index + 1 for index, record in enumerate(input_ops)}
    mapped_index_by_input_slot: dict[int, int] = {}
    mapped_index_by_candidate_op: dict[int, int] = {}
    mapped_input_slots: list[int] = []
    mapped_candidate_ops: list[int] = []

    for slot, candidate_op in zip(selected_input_slots, selected_candidate_ops):
        if slot in mapped_index_by_input_slot:
            return None
        index = len(mapped_input_slots) + 1
        mapped_index_by_input_slot[int(slot)] = index
        mapped_index_by_candidate_op[int(candidate_op)] = index
        mapped_input_slots.append(int(slot))
        mapped_candidate_ops.append(int(candidate_op))

    changed = True
    while changed:
        changed = False
        mapped_count = len(mapped_input_slots)
        for left_index in range(mapped_count):
            for right_index in range(mapped_count):
                candidate_product = magnetic_point_mul(
                    mapped_candidate_ops[left_index],
                    mapped_candidate_ops[right_index],
                )
                input_left_op = int(input_ops[mapped_input_slots[left_index] - 1][4])
                input_right_op = int(input_ops[mapped_input_slots[right_index] - 1][4])
                input_product_op = magnetic_point_mul(input_left_op, input_right_op)
                input_product_slot = input_slot_by_point_op.get(input_product_op)
                if input_product_slot is None:
                    return None

                existing_candidate_index = mapped_index_by_candidate_op.get(candidate_product)
                existing_input_index = mapped_index_by_input_slot.get(input_product_slot)
                if existing_candidate_index is None:
                    if existing_input_index is not None:
                        return None
                    index = len(mapped_input_slots) + 1
                    mapped_index_by_candidate_op[candidate_product] = index
                    mapped_index_by_input_slot[input_product_slot] = index
                    mapped_candidate_ops.append(candidate_product)
                    mapped_input_slots.append(input_product_slot)
                    changed = True
                    continue
                if mapped_input_slots[existing_candidate_index - 1] != input_product_slot:
                    return None
                if existing_input_index != existing_candidate_index:
                    return None
    if len(mapped_input_slots) < len(input_ops):
        return None
    selected = set(int(value) for value in selected_input_slots)
    canonical_by_slot = {
        int(slot): int(candidate_op)
        for slot, candidate_op in zip(mapped_input_slots, mapped_candidate_ops)
    }
    return tuple(
        MagneticPointOperationCorrespondence(
            input_slot=index,
            input_point_op=int(record[4]),
            canonical_point_op=canonical_by_slot[index],
            selected_generator=index in selected,
        )
        for index, record in enumerate(input_ops, start=1)
    )


def _magnetic_closure_accepts(
    input_ops: Sequence[Sequence[int]],
    selected_input_slots: tuple[int, ...],
    selected_candidate_ops: tuple[int, ...],
) -> bool:
    """Return whether a generator mapping closes under magnetic multiplication."""

    return (
        _magnetic_closure_correspondence(
            input_ops,
            selected_input_slots,
            selected_candidate_ops,
        )
        is not None
    )


def candidate_magnetic_generator_mappings(
    input_ops: Sequence[Sequence[int]],
    candidate: MagneticGroupCandidate,
) -> tuple[MagneticGeneratorMapping, ...]:
    """Enumerate accepted magnetic generator mappings for one candidate group."""

    match_slots: list[tuple[int, ...]] = []
    for generator in candidate.generators:
        code = magnetic_point_code(int(generator[4]))
        slots = tuple(
            index + 1
            for index, operation in enumerate(input_ops)
            if magnetic_point_code(int(operation[4])) == code
        )
        match_slots.append(slots)
    if any(not slots for slots in match_slots):
        return ()

    accepted: list[MagneticGeneratorMapping] = []
    candidate_ops = tuple(int(record[4]) for record in candidate.generators)
    total = 1
    for slots in match_slots:
        total *= len(slots)
    for ordinal in range(total):
        value = ordinal
        selected_list: list[int] = []
        for slots in match_slots:
            selected_list.append(int(slots[value % len(slots)]))
            value //= len(slots)
        selected = tuple(selected_list)
        if len(set(selected)) != len(selected):
            continue
        if not _magnetic_closure_accepts(input_ops, selected, candidate_ops):
            continue
        accepted.append(
            MagneticGeneratorMapping(
                candidate_magnetic_group=int(candidate.magnetic_group),
                candidate_generator_count=len(candidate.generators),
                per_generator_match_counts=tuple(len(slots) for slots in match_slots),
                per_generator_match_slots=tuple(match_slots),
                selected_input_generator_slots=selected,
                selected_candidate_generator_ops=candidate_ops,
            )
        )
    return tuple(accepted)


def id_subgroup_magnetic_input_generator_matrices(
    data: SourceData,
    parent_sg: int,
    basis: Sequence[int],
    input_ops: Sequence[Sequence[int]],
) -> tuple[tuple[tuple[int, int, int], ...], ...]:
    """Magnetic variant of the input generator matrix transform.

    Magnetic operation records keep magnetic point-operation labels.  The
    spatial matrix used by ``id_subgroup_magnetic_`` is the paired nonmagnetic
    point operation, while the point-op labels used for generator matching stay
    magnetic.
    """

    inverse, inverse_denominator = _matinv_numerator_denominator(basis)
    lattice = int(data.space["ispace_lattice"][int(parent_sg) - 1])
    setting = data.ml_lattice_setting_record(lattice)
    denom0, _denom1 = data.ml_lattice_setting_denominators(lattice)
    left = setting[:9]
    right = setting[9:18]
    out: list[tuple[tuple[int, int, int], ...]] = []
    for record in input_ops:
        nonmag_op = magnetic_nonmag_point_op(int(record[4]))
        raw = data._matmlt_iso_storage(
            data._matmlt_iso_storage(
                data._matmlt_iso_storage(inverse, left),
                data.point_operation_matrix(nonmag_op),
            ),
            right,
        )
        raw = data._matmlt_iso_storage(raw, tuple(int(value) for value in basis[:9]))
        reduced, denominator = data._reduce_integer_matrix_denominator(
            raw,
            int(denom0) * int(inverse_denominator),
        )
        if denominator != 1:
            raise ValueError(
                f"id_subgroup_magnetic input matrix did not reduce to integer: "
                f"parent={parent_sg} op={record[4]} denominator={denominator}"
            )
        out.append(_matrix3_rows(reduced))
    return tuple(out)


def _fraction_row_matrix(
    vector: Sequence[Fraction],
    matrix: Sequence[Sequence[Fraction | int]],
) -> tuple[Fraction, Fraction, Fraction]:
    return fraction_row_multiply3(vector, matrix)  # type: ignore[arg-type]


def _fraction_record(values: Sequence[Fraction]) -> tuple[int, int, int, int]:
    denominator = 1
    for value in values:
        denominator = denominator * value.denominator // math.gcd(
            denominator,
            value.denominator,
        )
    numerators = [int(value * denominator) for value in values]
    divisor = abs(denominator)
    for value in numerators:
        divisor = math.gcd(divisor, abs(value))
    if divisor > 1:
        numerators = [value // divisor for value in numerators]
        denominator //= divisor
    return numerators[0], numerators[1], numerators[2], denominator


def _record_translation_mod1(record: Sequence[int]) -> tuple[Fraction, Fraction, Fraction]:
    return tuple(
        Fraction(int(record[axis]), int(record[3])) % 1
        for axis in range(3)
    )  # type: ignore[return-value]


def _lift_full_operation_correspondence(
    data: SourceData,
    *,
    active_ops: Sequence[MagneticOperationRecord],
    local_input_ops: Sequence[MagneticOperationRecord],
    basis_candidate: Sequence[Sequence[int]],
    origin_solution: Sequence[int],
    origin_denominator: int,
    magnetic_group: int,
    ordinary_space_group: int,
    point_correspondence: Sequence[MagneticPointOperationCorrespondence],
) -> tuple[MagneticFullOperationCorrespondence, ...]:
    """Lift a point-operation map through the accepted affine embedding.

    The origin solver uses ``t_local B^-1 - t_canonical = o - o R``.  Applying
    that same Source convention to every closure record retains antiunitary
    labels and translation lifts, not merely the generator point-operation
    correspondence.  The transformed translation must match the canonical BNS
    Source record modulo its primitive integer lattice; otherwise the full
    lift is intentionally withheld while the point-operation surface remains
    available.
    """

    basis_flat = tuple(int(value) for row in basis_candidate for value in row)
    inverse, inverse_denominator = _matinv_numerator_denominator(basis_flat)
    inverse_matrix = tuple(
        tuple(
            Fraction(int(inverse[row * 3 + col]), int(inverse_denominator))
            for col in range(3)
        )
        for row in range(3)
    )
    origin = tuple(
        Fraction(int(origin_solution[axis]), int(origin_denominator))
        for axis in range(3)
    )
    canonical_by_op = {
        int(record[4]): tuple(int(value) for value in record[:5])
        for record in data.generate_magnetic_space_group_records(int(magnetic_group))
    }
    point_by_slot = {
        int(item.input_slot): int(item.canonical_point_op)
        for item in point_correspondence
    }
    out: list[MagneticFullOperationCorrespondence] = []
    for slot, (active, local) in enumerate(
        zip(active_ops, local_input_ops, strict=True),
        start=1,
    ):
        canonical_op = point_by_slot[int(slot)]
        canonical = canonical_by_op.get(int(canonical_op))
        if canonical is None:
            return ()
        candidate_matrix = data.id_subgroup_candidate_generator_matrix(
            int(ordinary_space_group),
            magnetic_nonmag_point_op(int(canonical_op)),
        )
        rotated_origin = _fraction_row_matrix(origin, candidate_matrix)
        local_translation = tuple(
            Fraction(int(local[axis]), int(local[3]))
            for axis in range(3)
        )
        transformed_translation = _fraction_row_matrix(
            local_translation,
            inverse_matrix,
        )
        canonical_translation = tuple(
            transformed_translation[axis] - origin[axis] + rotated_origin[axis]
            for axis in range(3)
        )
        if tuple(value % 1 for value in canonical_translation) != _record_translation_mod1(
            canonical
        ):
            return ()
        translation_record = _fraction_record(canonical_translation)
        out.append(
            MagneticFullOperationCorrespondence(
                input_slot=int(slot),
                input_record=tuple(int(value) for value in active[:5]),
                local_input_record=tuple(int(value) for value in local[:5]),
                transformed_input_record=(
                    int(translation_record[0]),
                    int(translation_record[1]),
                    int(translation_record[2]),
                    int(translation_record[3]),
                    int(canonical_op),
                ),
                canonical_record=tuple(int(value) for value in canonical[:5]),
            )
        )
    return tuple(out)


def id_subgroup_magnetic_identify_with_generator_block(
    data: SourceData,
    parent_sg: int,
    basis: Sequence[int],
    input_ops: Sequence[Sequence[int]],
    flag: int = 0,
    *,
    ordinary_space_group_hint: int | None = None,
) -> MagneticIdSubgroupResult | None:
    """Port the magnetic candidate/generator branch of ``id_subgroup_magnetic_``.

    This mirrors the ordinary ``id_subgroup_`` search, but candidate groups,
    generator matching, and closure use magnetic point-operation labels.
    Matrix equations still use the paired nonmagnetic point operation.
    """

    parent = int(parent_sg)
    basis_in = tuple(int(value) for value in basis[:9])
    lattice = int(data.space["ispace_lattice"][parent - 1])
    basis_flat = nice_lattice2_basis(data, lattice, basis_in) if int(flag) >= 0 else basis_in
    active_ops = tuple(
        tuple(int(value) for value in record[:5])
        for record in input_ops
        if len(record) >= 5 and int(record[4]) != 0
    )
    inverse, inverse_denominator = _matinv_numerator_denominator(basis_flat)
    local_input_ops = input_operation_records_from_inverse(inverse, inverse_denominator, active_ops)
    local_input_matrices = id_subgroup_magnetic_input_generator_matrices(
        data,
        parent,
        basis_flat,
        local_input_ops,
    )

    lattice_constraint: int | None = None
    candidates = magnetic_group_candidates_for_operations(
        local_input_ops,
        ordinary_space_group=ordinary_space_group_hint,
    )
    for base_value_count in (3, 5, 7, 11):
        for candidate_scan_ordinal, candidate in enumerate(candidates, start=1):
            candidate_lattice = int(data.space["ispace_lattice"][candidate.ordinary_space_group - 1])
            if lattice_constraint is not None and candidate_lattice != lattice_constraint:
                continue
            mappings = candidate_magnetic_generator_mappings(local_input_ops, candidate)
            if not mappings:
                continue
            candidate_matrices_by_op = {
                int(record[4]): data.id_subgroup_candidate_generator_matrix(
                    candidate.ordinary_space_group,
                    magnetic_nonmag_point_op(int(record[4])),
                )
                for record in candidate.generators
            }
            for mapping_ordinal, mapping in enumerate(mappings, start=1):
                selected_input_matrices = [
                    local_input_matrices[int(slot) - 1]
                    for slot in mapping.selected_input_generator_slots
                ]
                candidate_matrices = [
                    candidate_matrices_by_op[int(op)]
                    for op in mapping.selected_candidate_generator_ops
                ]
                rowop_input = id_subgroup_rowop_setup_matrix(selected_input_matrices, candidate_matrices)
                rowop_matrix = rowop_reduce_matrix(rowop_input, len(candidate_matrices) * 9, 9)
                free_count = _nonpivot_count(rowop_matrix, len(candidate_matrices) * 9, 9)
                limit = int(base_value_count) ** int(free_count)
                for candidate_index in range(1, limit + 1):
                    basis_candidate = id_subgroup_basis_candidate_from_rowop(
                        rowop_matrix,
                        len(candidate_matrices) * 9,
                        9,
                        base_value_count,
                        candidate_index,
                    )
                    if basis_candidate is None:
                        continue
                    diffs = generator_diff_records(
                        data,
                        candidate.ordinary_space_group,
                        mapping.selected_input_generator_slots,
                        local_input_ops,
                        tuple(value for row in basis_candidate for value in row),
                        candidate.generators,
                    )
                    equations = origin_equation_system(candidate_matrices, diffs)
                    solve = solve_eqs_mod_int_first(equations.matrix, equations.rhs, equations.denominator)
                    lattice_constraint = candidate_lattice
                    if not solve.success:
                        continue
                    basis_out = _matmlt_flat(data, basis_flat, basis_candidate)
                    origin_out = _origin_after_basis(basis_out, solve.solution, solve.denominator)
                    point_correspondence = _magnetic_closure_correspondence(
                        local_input_ops,
                        mapping.selected_input_generator_slots,
                        mapping.selected_candidate_generator_ops,
                    )
                    if point_correspondence is None:
                        continue
                    full_correspondence = _lift_full_operation_correspondence(
                        data,
                        active_ops=active_ops,
                        local_input_ops=local_input_ops,
                        basis_candidate=basis_candidate,
                        origin_solution=solve.solution,
                        origin_denominator=solve.denominator,
                        magnetic_group=int(candidate.magnetic_group),
                        ordinary_space_group=int(candidate.ordinary_space_group),
                        point_correspondence=point_correspondence,
                    )
                    return MagneticIdSubgroupResult(
                        int(candidate.magnetic_group),
                        int(candidate.ordinary_space_group),
                        basis_out,
                        origin_out,
                        selection=MagneticIdSubgroupSelection(
                            base_value_count=int(base_value_count),
                            candidate_scan_ordinal=int(candidate_scan_ordinal),
                            mapping_ordinal=int(mapping_ordinal),
                            basis_candidate_ordinal=int(candidate_index),
                            candidate_basis=tuple(
                                tuple(int(value) for value in row)
                                for row in basis_candidate
                            ),
                            candidate_origin=(
                                int(solve.solution[0]),
                                int(solve.solution[1]),
                                int(solve.solution[2]),
                                int(solve.denominator),
                            ),
                            generator_mapping=mapping,
                            point_operation_correspondence=point_correspondence,
                            full_operation_correspondence=full_correspondence,
                        ),
                    )
    return None
