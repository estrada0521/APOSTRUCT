"""Coupled multi-irrep OPD subgroup generation.

The Web Method-2 multi-slot path is not a Cartesian display of independently
computed OPD rows.  It enumerates relative domains, intersects the stabilizers
of every selected order-parameter subspace, and identifies the resulting
common subgroup.  This module implements that local path from Source data.
"""

from __future__ import annotations

from dataclasses import asdict
from fractions import Fraction
from functools import lru_cache
from itertools import product
import math
import re
from typing import Any, Sequence

import numpy as np

from ISODISTORT.Assembled.Backend.exactmath import (
    fraction_matrix_inverse3,
    integer_determinant3,
)

from ISODISTORT.Assembled.Backend.isotropy import catalog as isotropy_catalog
from ISODISTORT.Assembled.Backend.isotropy.engine.get_isotropy import (
    _equation_key,
    _full_orderparam,
    _transform_orderparam,
    _transform_orderparam_sparse,
)
from ISODISTORT.Assembled.Backend.isotropy.engine.id_subgroup import (
    id_subgroup_identify_with_generator_block,
)
from ISODISTORT.Assembled.Backend.isotropy.engine.id_subgroup_magnetic import (
    id_subgroup_magnetic_identify_with_generator_block,
)
from ISODISTORT.Assembled.Backend.isotropy.engine.lattice import get_new_fractionals
from ISODISTORT.Assembled.Backend.isotropy.engine.orderparam import orderparam_check
from ISODISTORT.Assembled.Backend.isotropy.engine.source_data import SourceData
from ISODISTORT.Assembled.Backend.isotropy.domains import (
    domain_operation_records,
    domain_operation_records_from_stabilizer,
    mapped_subgroup_records,
)
from ISODISTORT.Assembled.Backend.source import magnetic as magnetic_data
from ISODISTORT.Assembled.Backend.source.tables import SourceTables


def _number(value: Any) -> float:
    text = str(value).strip()
    try:
        return float(Fraction(text))
    except (ValueError, ZeroDivisionError):
        return float(text)


class _MagneticCoupledData:
    """Use magnetic parent operations with one slot's time-reversal action."""

    def __init__(
        self,
        base: SourceData,
        parent_sg: int,
        parent_magnetic_group: int,
        *,
        magnetic_irrep: bool,
    ):
        self._base = base
        self._parent_sg = int(parent_sg)
        self._parent_magnetic_group = int(parent_magnetic_group)
        self._magnetic_irrep = bool(magnetic_irrep)

    def __getattr__(self, name: str):
        return getattr(self._base, name)

    def generate_space_group_records(self, sg: int) -> tuple[tuple[int, int, int, int, int], ...]:
        if int(sg) != self._parent_sg:
            raise ValueError(f"magnetic coupled parent mismatch: {sg} != {self._parent_sg}")
        return self._base.generate_magnetic_space_group_records(self._parent_magnetic_group)

    def get_irreps_matrix_for_record(
        self,
        gid: int,
        record: tuple[int, int, int, int, int],
        kparam: Sequence[int],
    ):
        if self._magnetic_irrep:
            return self._base.get_irreps_matrix_for_magnetic_record(
                int(gid), record, tuple(int(value) for value in kparam)
            )
        table = magnetic_data.data().table
        nonmag_record = (
            int(record[0]),
            int(record[1]),
            int(record[2]),
            int(record[3]),
            int(table["mag_point_op_mag2nonmag"][int(record[4]) - 1]),
        )
        return self._base.get_irreps_matrix_for_record(
            int(gid), nonmag_record, tuple(int(value) for value in kparam)
        )

    def orderparam_to_subgroup(
        self,
        gid: int,
        kparam: Sequence[int],
        orderparam: Sequence[float],
        row_count: int,
    ) -> tuple[tuple[int, ...], tuple[tuple[int, int, int, int, int], ...]]:
        gid = int(gid)
        values = tuple(float(value) for value in orderparam)
        params = tuple(int(value) for value in kparam)
        dim = int(self._base.little["little_irr_full_dim"][gid - 1])
        candidates = list(self.generate_space_group_records(self._parent_sg))
        basis: tuple[int, ...] | None = None
        fractions: tuple[tuple[int, int, int, int], ...] = ()
        for row_index in range(int(row_count)):
            row = tuple(values[row_index * 48 + col] for col in range(dim))
            if row_index == 0:
                basis = self._base.newlat_order(gid, params, row)
                fractions = get_new_fractionals(basis)
            else:
                if basis is None:
                    raise ValueError("missing previous magnetic coupled basis")
                next_basis = self._base.newlat_order3(gid, params, row, basis)
                transform = self._base._basis_transform_matrix(basis, next_basis)  # noqa: SLF001
                fractions = tuple(
                    self._base._reduce_fraction_parent(basis, fraction)  # noqa: SLF001
                    for fraction in get_new_fractionals(transform)
                )
                basis = next_basis
            next_candidates: list[tuple[int, int, int, int, int]] = []
            for operation in candidates:
                for fraction in fractions:
                    shifted = self._base._vadd_fraction_operation(fraction, operation)  # noqa: SLF001
                    matrix = self.get_irreps_matrix_for_record(gid, shifted, params)
                    if orderparam_check(dim, 1, self._base._matrix_to_stride48(matrix), row):  # noqa: SLF001
                        next_candidates.append(shifted)
                        break
            candidates = next_candidates
        return basis or (1, 0, 0, 0, 1, 0, 0, 0, 1), tuple(candidates)


def _orderparam(row: dict[str, Any], dim: int) -> tuple[int, tuple[float, ...]]:
    direction = row.get("direction") or {}
    isotropy = row.get("isotropy") or {}
    # Static data_isotropy directions act in the stored irrep frame.  Their
    # source_display_rows may be a presentation transform such as (0,a) ->
    # (a,0.577a), which must not be fed back into subgroup stabilization.
    # Dynamic .iso rows have no static direction record; their generated source
    # rows are already the binary order-parameter frame.
    vectors = (
        isotropy.get("source_display_rows")
        if bool(isotropy.get("dynamic"))
        else isotropy.get("source_numeric_rows")
    ) or isotropy.get("source_display_rows") or direction.get("vectors") or []
    free = len(vectors)
    out = [0.0] * (max(1, free) * 48)
    for row_index, values in enumerate(vectors):
        for col, value in enumerate(values[:dim]):
            out[row_index * 48 + col] = _number(value)
    return free, tuple(out)


_SubgroupResult = tuple[
    tuple[int, ...],
    tuple[tuple[int, int, int, int, int], ...],
]


def _cached_orderparam_to_subgroup(
    irrep_data: SourceData | _MagneticCoupledData,
    gid: int,
    kparam: Sequence[int],
    values: Sequence[float],
    free: int,
    cache: dict[tuple[object, ...], _SubgroupResult] | None,
) -> _SubgroupResult:
    params = tuple(int(value) for value in kparam)
    orderparam = tuple(float(value) for value in values)
    key = (id(irrep_data), int(gid), params, int(free), orderparam)
    if cache is not None:
        cached = cache.get(key)
        if cached is not None:
            return cached
    result = irrep_data.orderparam_to_subgroup(
        int(gid), params, orderparam, int(free)
    )
    if cache is not None:
        cache[key] = result
    return result


def _cached_newlat_step(
    data: SourceData | _MagneticCoupledData,
    gid: int,
    kparam: Sequence[int],
    vector: Sequence[float],
    previous_basis: tuple[int, ...] | None,
    cache: dict[tuple[object, ...], tuple[int, ...]] | None,
) -> tuple[int, ...]:
    params = tuple(int(value) for value in kparam)
    values = tuple(float(value) for value in vector)
    key = (id(data), int(gid), params, values, previous_basis)
    if cache is not None:
        cached = cache.get(key)
        if cached is not None:
            return cached
    result = (
        data.newlat_order(int(gid), params, values)
        if previous_basis is None
        else data.newlat_order3(int(gid), params, values, previous_basis)
    )
    if cache is not None:
        cache[key] = result
    return result


def _combined_stabilizer(
    data: SourceData,
    parent_sg: int,
    slots: Sequence[dict[str, Any]],
    orderparams: Sequence[tuple[int, Sequence[float]]],
    operation_orders: Sequence[Sequence[tuple[int, int, int, int, int]] | None] | None = None,
    slot_data: Sequence[SourceData | _MagneticCoupledData] | None = None,
    subgroup_cache: dict[tuple[object, ...], _SubgroupResult] | None = None,
    basis_cache: dict[tuple[object, ...], tuple[int, ...]] | None = None,
    newlat_cache: dict[tuple[object, ...], tuple[int, ...]] | None = None,
    membership_cache: dict[tuple[object, ...], bool] | None = None,
) -> tuple[tuple[int, ...], tuple[tuple[int, int, int, int, int], ...]]:
    basis: tuple[int, ...] | None = None
    individual_subgroups: list[
        tuple[tuple[int, ...], tuple[tuple[int, int, int, int, int], ...]]
    ] = []
    action_data = tuple(slot_data) if slot_data is not None else (data,) * len(slots)
    for slot, (free, values), irrep_data in zip(slots, orderparams, action_data):
        gid = int(slot["irrep"]["gid"])
        kparam = tuple(int(value) for value in slot["source_kparam"])
        individual_subgroups.append(
            _cached_orderparam_to_subgroup(
                irrep_data, gid, kparam, values, int(free), subgroup_cache
            )
        )
    basis_key = (
        id(data),
        tuple(
            (
                int(slot["irrep"]["gid"]),
                tuple(int(value) for value in slot["source_kparam"]),
                int(free),
                tuple(float(value) for value in values),
            )
            for slot, (free, values) in zip(slots, orderparams)
        ),
    )
    if basis_cache is not None:
        basis = basis_cache.get(basis_key)
    if basis is None:
        for slot, (free, values) in zip(slots, orderparams):
            gid = int(slot["irrep"]["gid"])
            kparam = tuple(int(value) for value in slot["source_kparam"])
            dim = int(slot["irrep"]["full_dim"])
            for row_index in range(int(free)):
                vector = tuple(float(values[row_index * 48 + col]) for col in range(dim))
                basis = _cached_newlat_step(
                    data, gid, kparam, vector, basis, newlat_cache
                )
        if basis_cache is not None and basis is not None:
            basis_cache[basis_key] = basis
    final_basis = basis or (1, 0, 0, 0, 1, 0, 0, 0, 1)
    fractions = get_new_fractionals(final_basis)
    candidates: list[tuple[int, int, int, int, int]] = []
    for operation in data.generate_space_group_records(int(parent_sg)):
        for fraction in fractions:
            shifted = data._vadd_fraction_operation(fraction, operation)  # noqa: SLF001
            included = True
            for subgroup_basis, subgroup_operations in individual_subgroups:
                membership_key = (
                    shifted,
                    tuple(int(value) for value in subgroup_basis),
                    tuple(subgroup_operations),
                )
                if membership_cache is not None and membership_key in membership_cache:
                    member = membership_cache[membership_key]
                else:
                    member = _operation_in_subgroup(
                        shifted, subgroup_basis, subgroup_operations
                    )
                    if membership_cache is not None:
                        membership_cache[membership_key] = member
                if not member:
                    included = False
                    break
            if included:
                candidates.append(shifted)
                break
    if operation_orders:
        # coupled_subgroups2_ builds intersections in reverse slot order and
        # subgroups_intersect_ retains the left operand's operation order.
        # Preserve the corresponding stored child-operation order instead of
        # re-sorting the common set by parent point-operation number.
        # Magnetic intersections retain the final slot's transformed order.
        # Only its base-domain cinter order is available here; using an earlier
        # slot's order for a non-base final domain changes the chosen embedding.
        retained_order = (
            operation_orders[-1]
            if isinstance(data, _MagneticCoupledData)
            else next((order for order in reversed(operation_orders) if order), None)
        )
        if retained_order:
            point_order = {int(record[4]): index for index, record in enumerate(retained_order)}
            if isinstance(data, _MagneticCoupledData):
                mag_to_nonmag = magnetic_data.data().table["mag_point_op_mag2nonmag"]
                candidates.sort(
                    key=lambda record: point_order.get(
                        int(mag_to_nonmag[int(record[4]) - 1]),
                        len(point_order),
                    )
                )
            else:
                candidates.sort(key=lambda record: point_order.get(int(record[4]), len(point_order)))
    return final_basis, tuple(candidates)


def _operation_in_subgroup(
    candidate: tuple[int, int, int, int, int],
    basis: Sequence[int],
    operations: Sequence[tuple[int, int, int, int, int]],
) -> bool:
    representative = next((operation for operation in operations if int(operation[4]) == int(candidate[4])), None)
    if representative is None:
        return False
    delta = tuple(
        Fraction(int(candidate[axis]), int(candidate[3]))
        - Fraction(int(representative[axis]), int(representative[3]))
        for axis in range(3)
    )
    denominator = math.lcm(*(value.denominator for value in delta))
    record = tuple(int(value * denominator) for value in delta) + (denominator,)
    return _fraction_in_lattice(record, basis)


def _same_lattice(left: Sequence[int], right: Sequence[int]) -> bool:
    for source, target in ((left, right), (right, left)):
        for offset in range(0, 9, 3):
            if not _fraction_in_lattice(
                (int(source[offset]), int(source[offset + 1]), int(source[offset + 2]), 1),
                target,
            ):
                return False
    return True


def _same_subgroup(
    left_basis: Sequence[int],
    left_operations: Sequence[tuple[int, int, int, int, int]],
    right_basis: Sequence[int],
    right_operations: Sequence[tuple[int, int, int, int, int]],
) -> bool:
    """Match the subgroup-embedding equivalence used by ``same_subgroup_``."""

    if not _same_lattice(left_basis, right_basis):
        return False
    if {int(record[4]) for record in left_operations} != {
        int(record[4]) for record in right_operations
    }:
        return False
    return all(
        _operation_in_subgroup(record, right_basis, right_operations)
        for record in left_operations
    ) and all(
        _operation_in_subgroup(record, left_basis, left_operations)
        for record in right_operations
    )


def _operation_records(
    data: SourceData,
    parent_sg: int,
    kernel_basis: Sequence[int],
) -> tuple[tuple[int, int, int, int, int], ...]:
    generated = data.generate_space_group_records(int(parent_sg))
    fractions = get_new_fractionals(tuple(int(value) for value in kernel_basis))
    if isinstance(data, _MagneticCoupledData):
        return tuple(
            data._vadd_fraction_operation(fraction, operation)  # noqa: SLF001
            for operation in generated
            for fraction in fractions
        )
    lattice = int(data.space["ispace_lattice"][int(parent_sg) - 1])
    return tuple(
        data._vadd_fraction_operation(  # noqa: SLF001
            data._rotate_kernel_fraction_by_space_operation(  # noqa: SLF001
                fraction,
                operation,
                lattice=lattice,
            ),
            operation,
        )
        for operation in generated
        for fraction in fractions
    )


def _magnetic_domain_operation_records_from_stabilizer(
    data: _MagneticCoupledData,
    *,
    parent_sg: int,
    basis: tuple[int, ...],
    subgroup_operations: tuple[tuple[int, int, int, int, int], ...],
) -> tuple[tuple[int, int, int, int, int], ...]:
    """Port the coset-order core of ``generate_domains2_magnetic_``."""

    table = magnetic_data.data().table
    parent_records = tuple(data.generate_space_group_records(int(parent_sg)))
    subgroup_points = {int(record[4]) for record in subgroup_operations}
    identity = (1, 0, 0, 0, 1, 0, 0, 0, 1)
    identity_index = next(
        index
        for index, record in enumerate(parent_records, start=1)
        if data.point_operation_matrix(
            int(table["mag_point_op_mag2nonmag"][int(record[4]) - 1])
        ) == identity
        and not bool(table["mag_point_op_r"][int(record[4]) - 1])
    )
    point_to_index = {
        int(record[4]): index for index, record in enumerate(parent_records, start=1)
    }
    marked = [False] * (len(parent_records) + 1)
    cosets = [identity_index]
    marked[identity_index] = True
    for point in subgroup_points:
        marked[point_to_index[point]] = True

    multiplication = table["mag_point_op_mlt"]

    def point_product(left: int, right: int) -> int:
        return int(multiplication[(int(right) - 1) * 144 + int(left) - 1])

    target = len(parent_records) // len(subgroup_points)
    while len(cosets) < target:
        chosen = next(
            (index for index in range(1, len(parent_records) + 1) if not marked[index]),
            0,
        )
        if chosen == 0:
            raise RuntimeError("generate_domains2_magnetic point-coset scan did not close")
        cosets.append(chosen)
        chosen_point = int(parent_records[chosen - 1][4])
        for subgroup_point in subgroup_points:
            marked[point_to_index[point_product(chosen_point, subgroup_point)]] = True

    lattice = int(data.space["ispace_lattice"][int(parent_sg) - 1])
    out: list[tuple[int, int, int, int, int]] = []
    for fraction in get_new_fractionals(basis):
        for coset in cosets:
            parent_record = parent_records[coset - 1]
            nonmag_record = (
                int(parent_record[0]),
                int(parent_record[1]),
                int(parent_record[2]),
                int(parent_record[3]),
                int(table["mag_point_op_mag2nonmag"][int(parent_record[4]) - 1]),
            )
            rotated = data._rotate_kernel_fraction_by_space_operation(  # noqa: SLF001
                fraction,
                nonmag_record,
                lattice=lattice,
            )
            out.append(data._vadd_fraction_operation(rotated, parent_record))  # noqa: SLF001
    return tuple(out)


def _domain_orbit(
    data: SourceData | _MagneticCoupledData,
    slot: dict[str, Any],
    source: tuple[int, Sequence[float]],
    operations: Sequence[tuple[int, int, int, int, int]],
    numbered_operations: Sequence[tuple[int, int, int, int, int]] | None = None,
    use_numbered_order: bool = False,
    subgroup_cache: dict[tuple[object, ...], _SubgroupResult] | None = None,
) -> tuple[
    tuple[tuple[float, ...], ...],
    dict[tuple[int, tuple[float, ...]], int],
    tuple[int, ...],
]:
    free, values = source
    gid = int(slot["irrep"]["gid"])
    dim = int(slot["irrep"]["full_dim"])
    kparam = tuple(int(value) for value in slot["source_kparam"])
    domains: list[tuple[float, ...]] = []
    domain_subgroups: list[
        tuple[tuple[int, ...], tuple[tuple[int, int, int, int, int], ...]]
    ] = []
    domain_numbers: list[int] = []
    indices: dict[tuple[int, tuple[float, ...]], int] = {}
    numbered_keys: dict[tuple[int, tuple[float, ...]], int] = {}
    if numbered_operations is not None and not use_numbered_order:
        for domain_number, operation in enumerate(numbered_operations, start=1):
            matrix = data.get_irreps_matrix_for_record(gid, operation, kparam)
            transformed = _transform_orderparam(dim, free, values, data._matrix_to_stride48(matrix))  # noqa: SLF001
            numbered_keys.setdefault(_equation_key(dim, free, transformed), domain_number)
    # Dynamic rows retain subspaces in their freshly generated domain order.
    # Static rows keep the established parent-kernel orbit order and use the
    # stored domain operations only to recover display ordinals.
    domain_operations = numbered_operations if use_numbered_order and numbered_operations is not None else operations
    for domain_number, operation in enumerate(domain_operations, start=1):
        matrix = data.get_irreps_matrix_for_record(gid, operation, kparam)
        transformed = _transform_orderparam(dim, free, values, data._matrix_to_stride48(matrix))  # noqa: SLF001
        key = _equation_key(dim, free, transformed)
        if key in indices:
            continue
        subgroup = _cached_orderparam_to_subgroup(
            data,
            gid,
            kparam,
            transformed,
            int(free),
            subgroup_cache,
        )
        equivalent = next(
            (
                index
                for index, previous in enumerate(domain_subgroups)
                if _same_subgroup(subgroup[0], subgroup[1], previous[0], previous[1])
            ),
            None,
        )
        if equivalent is not None:
            # ``generate_domains2_`` deduplicates subgroup embeddings, not OP
            # subspaces.  Keep the first direction but map every equivalent
            # subspace back to that retained domain.
            indices[key] = equivalent
            continue
        indices[key] = len(domains)
        domains.append(transformed)
        domain_subgroups.append(subgroup)
        domain_numbers.append(numbered_keys.get(key, domain_number))
    return tuple(domains), indices, tuple(domain_numbers)


def _domain_actions(
    data: SourceData,
    slot: dict[str, Any],
    free: int,
    domains: Sequence[Sequence[float]],
    indices: dict[tuple[int, tuple[float, ...]], int],
    operations: Sequence[tuple[int, int, int, int, int]],
) -> tuple[tuple[int, ...], ...]:
    gid = int(slot["irrep"]["gid"])
    dim = int(slot["irrep"]["full_dim"])
    kparam = tuple(int(value) for value in slot["source_kparam"])
    out: list[tuple[int, ...]] = []
    action_by_matrix: dict[tuple[float, ...], tuple[int, ...]] = {}
    retained_values = 0
    retained_value_limit = 1_000_000
    for operation in operations:
        matrix = data.get_irreps_matrix_for_record(gid, operation, kparam)
        stride = data._matrix_to_stride48(matrix)  # noqa: SLF001
        matrix_key = tuple(
            float(stride[source_row * 48 + col])
            for source_row in range(dim)
            for col in range(dim)
        )
        action = action_by_matrix.get(matrix_key)
        if action is None:
            stride_nonzero_cache: dict[int, tuple[tuple[tuple[int, float], ...], ...]] = {}
            transformed_domains = [
                _transform_orderparam_sparse(
                    dim,
                    free,
                    domain,
                    stride,
                    stride_nonzero_cache,
                )
                for domain in domains
            ]
            action = tuple(
                _domain_index(dim, free, transformed, domains, indices)
                for transformed in transformed_domains
            )
            value_count = len(matrix_key) + len(action)
            if retained_values + value_count <= retained_value_limit:
                action_by_matrix[matrix_key] = action
                retained_values += value_count
        out.append(action)
    return tuple(out)


def _domain_index(
    dim: int,
    free: int,
    transformed: Sequence[float],
    domains: Sequence[Sequence[float]],
    indices: dict[tuple[int, tuple[float, ...]], int],
) -> int:
    subspace_tolerance = 1e-6
    key = _equation_key(dim, free, transformed)
    exact = indices.get(key)
    if exact is not None:
        return exact
    left = np.asarray(
        [[float(transformed[row * 48 + col]) for col in range(dim)] for row in range(free)],
        dtype=float,
    )
    left_rank = int(np.linalg.matrix_rank(left, tol=subspace_tolerance))
    for index, domain in enumerate(domains):
        right = np.asarray(
            [[float(domain[row * 48 + col]) for col in range(dim)] for row in range(free)],
            dtype=float,
        )
        if int(np.linalg.matrix_rank(right, tol=subspace_tolerance)) != left_rank:
            continue
        if int(np.linalg.matrix_rank(np.vstack((left, right)), tol=subspace_tolerance)) == left_rank:
            return index
    raise KeyError(key)


def _relative_domain_tuples(actions: Sequence[Sequence[Sequence[int]]]) -> tuple[tuple[int, ...], ...]:
    ranges = [range(len(slot_actions[0])) for slot_actions in actions]
    retained: list[tuple[int, ...]] = []
    for candidate in product(*ranges):
        orbit = [
            tuple(actions[slot][operation][domain] for slot, domain in enumerate(candidate))
            for operation in range(len(actions[0]))
        ]
        if candidate == min(orbit):
            retained.append(candidate)
    return tuple(retained)


def _format_coeff(value: float, variable: str) -> str:
    if abs(value) < 5e-7:
        return ""
    if abs(value - 1.0) < 5e-7:
        return variable
    if abs(value + 1.0) < 5e-7:
        return f"-{variable}"
    return f"{value:.3f}{variable}"


def _join_opd_columns(columns: Sequence[str], group_size: int) -> str:
    if group_size <= 0:
        group_size = len(columns) or 1
    parts: list[str] = []
    for index, column in enumerate(columns):
        if index:
            parts.append(";" if index % group_size == 0 else ",")
        parts.append(str(column))
    return "(" + "".join(parts) + ")"


def _group_source_opd(opd: str, group_size: int) -> str:
    columns = [part.strip() for part in re.split(r"[,;]", str(opd).strip().strip("()"))]
    return _join_opd_columns(columns, group_size)


def _symbolic_opd(
    dim: int,
    free: int,
    values: Sequence[float],
    variable_offset: int = 0,
    group_size: int = 0,
) -> str:
    variables = "abcdefghijklmnopqrstuvwxyz"
    columns: list[str] = []
    for col in range(int(dim)):
        terms = [
            _format_coeff(float(values[row * 48 + col]), variables[variable_offset + row])
            for row in range(int(free))
        ]
        expression = "+".join(term for term in terms if term).replace("+-", "-")
        columns.append(expression or "0")
    return _join_opd_columns(columns, group_size or int(dim))


def _rename_opd_variables(opd: str, variable_offset: int) -> str:
    variables = "abcdefghijklmnopqrstuvwxyz"
    return "".join(
        variables[variable_offset + ord(char) - ord("a")] if "a" <= char <= "z" else char
        for char in str(opd)
    )


def _subduction_frequency(
    data: SourceData,
    slot: dict[str, Any],
    kernel_basis: Sequence[int],
    subgroup_basis: Sequence[int],
    operations: Sequence[tuple[int, int, int, int, int]],
    cache: dict[tuple[object, ...], int] | None = None,
) -> int:
    gid = int(slot["irrep"]["gid"])
    kparam = tuple(int(value) for value in slot["source_kparam"])
    cache_key = (
        id(data),
        gid,
        kparam,
        tuple(int(value) for value in kernel_basis),
        tuple(int(value) for value in subgroup_basis),
        tuple(operations),
    )
    if cache is not None and cache_key in cache:
        return cache[cache_key]
    if not operations:
        return 0
    fractions = [
        fraction
        for fraction in get_new_fractionals(tuple(int(value) for value in kernel_basis))
        if _fraction_in_lattice(fraction, subgroup_basis)
    ]
    character_sum = 0.0
    count = 0
    for operation in operations:
        for fraction in fractions:
            shifted = data._vadd_fraction_operation(fraction, operation)  # noqa: SLF001
            matrix = data.get_irreps_matrix_for_record(gid, shifted, kparam)
            character_sum += float(matrix.trace())
            count += 1
    result = int(round(character_sum / count)) if count else 0
    if cache is not None:
        cache[cache_key] = result
    return result


def _fraction_in_lattice(fraction: Sequence[int], basis: Sequence[int]) -> bool:
    values = [Fraction(int(value), int(fraction[3] or 1)) for value in fraction[:3]]
    inverse = _basis_inverse(tuple(int(value) for value in basis[:9]))
    if inverse is None:
        return False
    coordinates = [sum(values[row] * inverse[row * 3 + col] for row in range(3)) for col in range(3)]
    return all(value.denominator == 1 for value in coordinates)


@lru_cache(maxsize=512)
def _basis_inverse(basis: tuple[int, ...]) -> tuple[Fraction, ...] | None:
    matrix = tuple(
        tuple(Fraction(int(basis[3 * row + column])) for column in range(3))
        for row in range(3)
    )
    try:
        inverse = fraction_matrix_inverse3(matrix)
    except ValueError:
        return None
    return tuple(value for row in inverse for value in row)


def _active_k_vectors(
    display_data: SourceTables,
    parent_sg: int,
    slot: dict[str, Any],
    free: int,
    values: Sequence[float],
) -> list[str]:
    dim = int(slot["irrep"]["full_dim"])
    rows = [
        [float(values[row * 48 + col]) for col in range(dim)]
        for row in range(int(free))
    ]
    return isotropy_catalog._active_k_vectors_for_gid(  # noqa: SLF001
        display_data,
        int(parent_sg),
        int(slot["irrep"]["gid"]),
        dim,
        rows,
        slot.get("display_k_params") or slot.get("k_params") or {},
    )


def _coupled_source_rows(rows: Sequence[dict[str, Any]]) -> tuple[dict[str, Any], ...]:
    """Match the static data-isotropy or generated ``s*.iso`` traversal."""
    if rows and all(not bool((row.get("isotropy") or {}).get("dynamic")) for row in rows):
        return tuple(sorted(
            rows,
            key=lambda row: int((row.get("isotropy") or {}).get("row_id") or 0),
        ))
    return tuple(rows)


def coupled_opd_rows(
    data: SourceData,
    *,
    display_data: SourceTables,
    parent_sg: int,
    slots: Sequence[dict[str, Any]],
    parent_setting_id: int,
    parent_cell: tuple[float, float, float, float, float, float] | None,
    row_limit: int | None = None,
) -> list[dict[str, Any]]:
    """Enumerate coupled OPD rows in ordered source/domain traversal."""

    if len(slots) < 2:
        raise ValueError("coupled_opd_rows requires at least two slots")
    if row_limit is not None and int(row_limit) < 1:
        raise ValueError("row_limit must be positive")
    magnetic_slots = tuple(bool(slot["irrep"].get("magnetic")) for slot in slots)
    magnetic_request = any(magnetic_slots)
    coupled_data: SourceData | _MagneticCoupledData = data
    slot_coupled_data: tuple[SourceData | _MagneticCoupledData, ...] = (data,) * len(slots)
    parent_magnetic_group: int | None = None
    if magnetic_request:
        parent_groups = {
            int(data.magnetic_orderparam_group_for_irrep(int(slot["irrep"]["gid"])))
            for slot in slots
            if bool(slot["irrep"].get("magnetic"))
        }
        if len(parent_groups) != 1:
            raise ValueError(f"magnetic coupled slots use different parent groups: {sorted(parent_groups)}")
        parent_magnetic_group = next(iter(parent_groups))
        slot_coupled_data = tuple(
            _MagneticCoupledData(
                data,
                int(parent_sg),
                parent_magnetic_group,
                magnetic_irrep=magnetic,
            )
            for magnetic in magnetic_slots
        )
        coupled_data = slot_coupled_data[0]

    rows: list[dict[str, Any]] = []
    limited_seen: set[tuple[Any, ...]] = set()
    row_id = 0
    subgroup_cache: dict[tuple[object, ...], _SubgroupResult] = {}
    basis_cache: dict[tuple[object, ...], tuple[int, ...]] = {}
    newlat_cache: dict[tuple[object, ...], tuple[int, ...]] = {}
    membership_cache: dict[tuple[object, ...], bool] = {}
    subduction_cache: dict[tuple[object, ...], int] = {}
    kernel_sources = tuple(
        (int(slot["irrep"]["full_dim"]), _full_orderparam(int(slot["irrep"]["full_dim"])))
        for slot in slots
    )
    slot_kernel_bases = tuple(
        _combined_stabilizer(
            coupled_data,
            parent_sg,
            (slot,),
            (source,),
            slot_data=(irrep_data,),
            subgroup_cache=subgroup_cache,
            basis_cache=basis_cache,
            newlat_cache=newlat_cache,
            membership_cache=membership_cache,
        )[0]
        for slot, source, irrep_data in zip(slots, kernel_sources, slot_coupled_data)
    )
    kernel_basis, _kernel_operations = _combined_stabilizer(
        coupled_data,
        parent_sg,
        slots,
        kernel_sources,
        slot_data=slot_coupled_data,
        subgroup_cache=subgroup_cache,
        basis_cache=basis_cache,
        newlat_cache=newlat_cache,
        membership_cache=membership_cache,
    )
    parent_kernel_operations = _operation_records(coupled_data, parent_sg, kernel_basis)
    mixed_dynamic_request = any(bool(slot.get("k_params")) for slot in slots)
    # Fortran's coupled row arrays vary the first irrep index fastest.
    ordered_slot_rows = tuple(_coupled_source_rows(slot["opd_rows"]) for slot in slots)
    domain_cache: dict[
        tuple[int, int],
        tuple[
            tuple[
                tuple[tuple[float, ...], ...],
                dict[tuple[int, tuple[float, ...]], int],
                tuple[int, ...],
            ],
            tuple[tuple[int, ...], ...],
        ],
    ] = {}
    relative_domain_cache: dict[
        tuple[tuple[tuple[int, ...], ...], ...],
        tuple[tuple[int, ...], ...],
    ] = {}
    for reversed_rows in product(*(rows for rows in reversed(ordered_slot_rows))):
        source_rows = tuple(reversed(reversed_rows))
        canonical = tuple(
            _orderparam(row, int(slot["irrep"]["full_dim"]))
            for slot, row in zip(slots, source_rows)
        )
        operations = parent_kernel_operations
        individual_operations: list[tuple[tuple[int, int, int, int, int], ...]] = []
        individual_source_orders: list[bool] = []
        stored_operation_orders: list[tuple[tuple[int, int, int, int, int], ...] | None] = []
        for slot, source, source_row, irrep_data in zip(
            slots, canonical, source_rows, slot_coupled_data
        ):
            subgroup_basis, subgroup_operations = _cached_orderparam_to_subgroup(
                irrep_data,
                int(slot["irrep"]["gid"]),
                tuple(int(value) for value in slot["source_kparam"]),
                source[1],
                int(source[0]),
                subgroup_cache,
            )
            # ``coupled_subgroups2_`` calls ``generate_domains2_`` with the
            # per-direction cinter basis stored in the isotropy row.  The raw
            # translation kernel returned by ``orderparam_to_subgroup_`` spans
            # the same lattice, but its fractional-coset order is different and
            # therefore assigns different domain numbers.
            source_iso = source_row.get("isotropy") or {}
            source_basis = tuple(
                int(value)
                for value in (source_iso.get("source_basis_values") or ())
            )
            source_origin = tuple(
                int(value)
                for value in (source_iso.get("source_origin_values") or ())
            )
            source_subgroup = source_iso.get("subgroup") or {}
            source_child = int(
                (
                    source_subgroup.get("ordinary_number")
                    if bool(slot["irrep"].get("magnetic"))
                    else source_subgroup.get("number")
                )
                or 0
            )
            # A mixed fixed/dynamic request traverses every slot through the
            # freshly assembled coupled domain table.  Fixed-only requests keep
            # the established stabilizer orbit.  Some locally reconstructed
            # cinter tables are not closed under the parent action; those are
            # detected below and fall back to the complete stabilizer orbit.
            use_source_domain_path = (
                magnetic_request or mixed_dynamic_request
            ) and len(source_basis) == 9
            numbering_basis = (
                source_basis
                if use_source_domain_path
                else tuple(int(value) for value in subgroup_basis)
            )
            binary_domain_path = False
            if (
                not magnetic_request
                and
                source_child > 0
                and len(source_basis) == 9
                and len(source_origin) == 4
            ):
                try:
                    domain_operations = domain_operation_records(
                        data,
                        parent_sg=int(parent_sg),
                        child_sg=source_child,
                        basis=source_basis,
                        origin=source_origin,  # type: ignore[arg-type]
                    )
                    numbering_basis = source_basis
                    binary_domain_path = True
                except (KeyError, RuntimeError, ValueError):
                    binary_domain_path = False
            if not binary_domain_path:
                if magnetic_request:
                    domain_operations = _magnetic_domain_operation_records_from_stabilizer(
                        coupled_data,
                        parent_sg=int(parent_sg),
                        basis=numbering_basis,
                        subgroup_operations=tuple(subgroup_operations),
                    )
                    binary_domain_path = True
                else:
                    domain_operations = domain_operation_records_from_stabilizer(
                        data,
                        parent_sg=int(parent_sg),
                        basis=numbering_basis,
                        subgroup_operations=tuple(subgroup_operations),
                    )
            individual_operations.append(domain_operations)
            individual_source_orders.append(binary_domain_path or use_source_domain_path)
            stored_order = None
            if (
                source_child > 0
                and len(source_basis) == 9
                and len(source_origin) == 4
            ):
                try:
                    stored_order = mapped_subgroup_records(
                        data,
                        parent_sg=int(parent_sg),
                        child_sg=source_child,
                        basis=source_basis,
                        origin=source_origin,  # type: ignore[arg-type]
                    )
                except (KeyError, RuntimeError, ValueError):
                    stored_order = None
            stored_operation_orders.append(stored_order)
        orbits: list[
            tuple[tuple[tuple[float, ...], ...], dict[tuple[int, tuple[float, ...]], int], tuple[int, ...]]
        ] = []
        actions: list[tuple[tuple[int, ...], ...]] = []
        for slot_index, (slot, source, source_row, slot_operations, source_domain_order, irrep_data) in enumerate(zip(
            slots,
            canonical,
            source_rows,
            individual_operations,
            individual_source_orders,
            slot_coupled_data,
        )):
            cache_key = (slot_index, id(source_row))
            cached_domain = domain_cache.get(cache_key)
            if cached_domain is None:
                orbit = _domain_orbit(
                    irrep_data,
                    slot,
                    source,
                    operations,
                    slot_operations,
                    source_domain_order,
                    subgroup_cache,
                )
                try:
                    action = _domain_actions(irrep_data, slot, source[0], orbit[0], orbit[1], operations)
                except KeyError:
                    if not source_domain_order:
                        raise
                    orbit = _domain_orbit(
                        irrep_data,
                        slot,
                        source,
                        operations,
                        slot_operations,
                        False,
                        subgroup_cache,
                    )
                    action = _domain_actions(irrep_data, slot, source[0], orbit[0], orbit[1], operations)
                domain_cache[cache_key] = (orbit, action)
            else:
                orbit, action = cached_domain
            orbits.append(orbit)
            actions.append(action)
        action_key = tuple(actions)
        relative_domains = relative_domain_cache.get(action_key)
        if relative_domains is None:
            relative_domains = _relative_domain_tuples(actions)
            relative_domain_cache[action_key] = relative_domains
        for domain_tuple in relative_domains:
            transformed = tuple(orbits[index][0][domain] for index, domain in enumerate(domain_tuple))
            displayed_domains = tuple(orbits[index][2][domain] for index, domain in enumerate(domain_tuple))
            basis, subgroup_operations = _combined_stabilizer(
                coupled_data,
                parent_sg,
                slots,
                tuple((canonical[index][0], transformed[index]) for index in range(len(slots))),
                tuple(
                    stored_operation_orders[index] if domain == 0 else None
                    for index, domain in enumerate(domain_tuple)
                ),
                slot_data=slot_coupled_data,
                subgroup_cache=subgroup_cache,
                basis_cache=basis_cache,
                newlat_cache=newlat_cache,
                membership_cache=membership_cache,
            )
            expected_subduction: list[int] = []
            for slot, source in zip(slots, canonical):
                expected_subduction.append(int(source[0]))
            # The standalone 9.6.1 binary evaluates each irrep with its own
            # translation kernel.  Method-2 coupled tables instead retain the
            # common coupled kernel at this gate.  Keep both values visible so
            # the binary-faithful and Web-presentation boundaries do not merge.
            binary_subduction = [
                _subduction_frequency(
                    irrep_data,
                    slot,
                    slot_kernel_basis,
                    basis,
                    subgroup_operations,
                    subduction_cache,
                )
                for slot, slot_kernel_basis, irrep_data in zip(
                    slots, slot_kernel_bases, slot_coupled_data
                )
            ]
            coupled_subduction = [
                _subduction_frequency(
                    irrep_data,
                    slot,
                    kernel_basis,
                    basis,
                    subgroup_operations,
                    subduction_cache,
                )
                for slot, irrep_data in zip(slots, slot_coupled_data)
            ]
            if coupled_subduction != expected_subduction:
                continue
            if magnetic_request:
                result = id_subgroup_magnetic_identify_with_generator_block(
                    data,
                    parent_sg,
                    basis,
                    subgroup_operations,
                    0,
                )
            else:
                result = id_subgroup_identify_with_generator_block(
                    data, parent_sg, basis, subgroup_operations, 0
                )
            if result is None:
                continue
            row_id += 1
            direction_parts: list[str] = []
            opd_parts: list[str] = []
            active_parts: list[tuple[str, ...]] = []
            variable_offset = 0
            for index, (slot, source_row, domain, displayed_domain) in enumerate(
                zip(slots, source_rows, domain_tuple, displayed_domains)
            ):
                iso = source_row.get("isotropy") or {}
                label = str(iso.get("opd_label") or (source_row.get("direction") or {}).get("label") or "")
                direction_parts.append(f"{label}({displayed_domain})")
                little_dim = int(slot["irrep"].get("little_dim") or slot["irrep"]["full_dim"])
                if domain == 0:
                    opd_parts.append(
                        _rename_opd_variables(
                            _group_source_opd(
                                str(iso.get("display_opd") or iso.get("source_opd") or ""),
                                little_dim,
                            ),
                            variable_offset,
                        ).strip("()")
                    )
                else:
                    opd_parts.append(
                        _symbolic_opd(
                            int(slot["irrep"]["full_dim"]),
                            canonical[index][0],
                            transformed[index],
                            variable_offset,
                            little_dim,
                        ).strip("()")
                    )
                active_parts.append(tuple(
                    _active_k_vectors(
                        display_data,
                        parent_sg,
                        slot,
                        canonical[index][0],
                        transformed[index],
                    )
                ))
                variable_offset += int(canonical[index][0])
            subgroup = int(
                result.magnetic_group if magnetic_request else result.subgroup
            )
            ordinary_subgroup = int(
                result.ordinary_space_group if magnetic_request else result.subgroup
            )
            det = abs(integer_determinant3(result.basis[:9]))
            magnetic_table = magnetic_data.data().table if magnetic_request else None
            symbol = (
                str(magnetic_table["mag_bns_label"][subgroup - 1]).strip()
                if magnetic_table is not None
                else display_data.default_setting_space_symbol(subgroup)
            )
            displayed = display_data.subgroup_change_setting_cinter(
                int(parent_sg),
                ordinary_subgroup,
                tuple(int(value) for value in result.basis),
                tuple(int(value) for value in result.origin),
                parent_setting_id=int(parent_setting_id),
                subgroup_setting_id=isotropy_catalog._opd_subgroup_setting_id(  # noqa: SLF001
                    int(parent_sg), ordinary_subgroup, int(parent_setting_id)
                ),
            )
            basis_rows = isotropy_catalog._basis_fraction_rows(  # noqa: SLF001
                displayed["basis"],
                int(displayed["basis_denominator"]),
            )
            basis_rows = isotropy_catalog.present_opd_basis_rows(
                ordinary_subgroup,
                basis_rows,
                parent_cell,
                data=isotropy_catalog._subgroup_core_data(),  # noqa: SLF001
                parametric=any(bool(slot.get("k_params")) for slot in slots),
                coupled=True,
            )
            direction_label = "".join(direction_parts)
            display_opd = "(" + "|".join(opd_parts) + ")"
            seen_active_vectors: set[str] = set()
            unique_active_parts: list[str] = []
            for vectors in active_parts:
                fresh = [vector for vector in vectors if vector not in seen_active_vectors]
                seen_active_vectors.update(fresh)
                if fresh:
                    unique_active_parts.append(",".join(fresh))
            iso_row = {
                "row_id": row_id,
                "source_ordinal": row_id,
                "coupled": True,
                "canonical": row_id == 1,
                "opd_label": direction_label,
                "direction_label": direction_label,
                "direction_domains": list(displayed_domains),
                "display_opd": display_opd,
                "source_opd": display_opd,
                "op_rows": sum(source[0] for source in canonical),
                "free": sum(source[0] for source in canonical),
                "s": det,
                "i": (
                    isotropy_catalog._magnetic_subgroup_index_for_groups(  # noqa: SLF001
                        int(parent_magnetic_group), subgroup, det
                    )
                    if magnetic_request and parent_magnetic_group is not None
                    else isotropy_catalog._subgroup_index(  # noqa: SLF001
                        display_data, int(parent_sg), subgroup, det
                    )
                ),
                "k_active_vectors": unique_active_parts,
                "k_active": ";".join(unique_active_parts),
                "subgroup": {
                    "number": subgroup,
                    "symbol": symbol,
                    **(
                        {
                            "ordinary_number": ordinary_subgroup,
                            "display_label": str(magnetic_table["mag_nlabel"][subgroup - 1]).strip(),
                        }
                        if magnetic_table is not None
                        else {}
                    ),
                },
                "basis": isotropy_catalog._basis_values(basis_rows),  # noqa: SLF001
                "basis_text": isotropy_catalog._basis_text(basis_rows),  # noqa: SLF001
                "basis_denominator": int(displayed["basis_denominator"]),
                "source_basis_values": list(result.basis),
                "source_origin_values": list(result.origin),
                **(
                    {"magnetic_parent_group": int(parent_magnetic_group)}
                    if magnetic_request and parent_magnetic_group is not None
                    else {}
                ),
                **(
                    {
                        "magnetic_subgroup_selection": asdict(result.selection),
                    }
                    if magnetic_request and result.selection is not None
                    else {}
                ),
                # Preserve the exact common stabilizer embedding for downstream
                # complete-mode subduction. Magnetic requests retain magnetic
                # point-operation ids; consumers must map them explicitly when
                # evaluating ordinary irreps.
                "source_operation_records": [list(record) for record in subgroup_operations],
                "subduction": {
                    "expected": expected_subduction,
                    "binary_slot_kernel": binary_subduction,
                    "coupled_common_kernel": coupled_subduction,
                },
                "det": det,
                "origin": isotropy_catalog._display_origin_text(  # noqa: SLF001
                    parent_sg,
                    ordinary_subgroup,
                    basis_rows,
                    displayed["origin"],
                ),
            }
            output_row = {
                "direction": {
                    "label": direction_label,
                    "opd": display_opd,
                    "components": [
                        {"slot": index + 1, "label": part, "domain": displayed_domains[index]}
                        for index, part in enumerate(direction_parts)
                    ],
                },
                "isotropy": iso_row,
                "_dedup_key": (
                    tuple(str((source_row.get("isotropy") or {}).get("opd_label") or "") for source_row in source_rows),
                    subgroup,
                    tuple(int(value) for value in result.basis),
                    tuple(int(value) % int(result.origin[3] or 1) for value in result.origin[:3]),
                    int(result.origin[3] or 1),
                ),
            }
            if row_limit is not None:
                dedup_key = output_row.pop("_dedup_key")
                if dedup_key in limited_seen:
                    continue
                limited_seen.add(dedup_key)
                iso_row["row_id"] = len(rows) + 1
                iso_row["canonical"] = not rows
                rows.append(output_row)
                if len(rows) >= int(row_limit):
                    return rows
            else:
                rows.append(output_row)
    unique: list[dict[str, Any]] = []
    seen: set[tuple[Any, ...]] = set()
    for row in rows:
        key = row.pop("_dedup_key")
        if key in seen:
            continue
        seen.add(key)
        unique.append(row)
    for index, row in enumerate(unique, start=1):
        row["isotropy"]["row_id"] = index
        row["isotropy"]["canonical"] = index == 1
    return unique
