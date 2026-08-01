"""Local nonmagnetic ``get_isotropy_`` dynamic-row generation.

This module implements the DISPLAY-ISOTROPY path that materializes dynamic
``*.iso`` rows for parametric-k irreps.  The public display layer sorts these
rows before printing, but the rows themselves are generated here from
``Source`` data only:

    OPD subspace -> orderparam_to_subgroup -> id_subgroup -> dynamic row

The current scope is the nonmagnetic ``param_4=0`` path used by
``DISPLAY ISOTROPY`` and ISODISTORT Method-2 OPD selection.
"""

from __future__ import annotations

from collections import deque
from collections.abc import Sequence
from fractions import Fraction
from functools import lru_cache

from distortropy.Backend.isotropy.engine.dynamic_isotropy_file import DynamicIsotropyRow
from distortropy.Backend.isotropy.engine.id_subgroup import id_subgroup_identify_with_generator_block
from distortropy.Backend.isotropy.engine.id_subgroup_magnetic import (
    id_subgroup_magnetic_identify_with_generator_block,
)
from distortropy.Backend.isotropy.engine.lattice import get_new_fractionals
from distortropy.Backend.isotropy.engine.orderparam import (
    eqs_to_orderparam,
    nice_orderparam,
    orderparam_add_eqs,
    orderparam_check,
    orderparam_label,
    orderparam_to_eqs,
    orthogonal_orderparam,
)
from distortropy.Backend.isotropy.engine.numerics import xrowop2
from distortropy.Backend.isotropy.engine.source_data import SourceData
from distortropy.Backend.source import magnetic as magnetic_data


Pair = tuple[int, int]
PairProductCache = dict[tuple[Pair, Pair], Pair]
OperationCache = dict[Pair, tuple[int, int, int, int, int]]
StrideMatrixCache = dict[tuple[int, int, int, int, int], tuple[float, ...]]
StrideNonzeroCache = dict[int, tuple[tuple[tuple[int, float], ...], ...]]
# Bound retained float objects and tuple slots to roughly 40 MiB per orbit node.
_ORBIT_ACTIVE_MEMO_SCALAR_BUDGET = 1_000_000


class _FactorizedMagneticPairProduct:
    """Multiply Source ``(space, kernel-coset)`` pairs without flattening the coset."""

    def __init__(
        self,
        base: SourceData,
        *,
        space_operations: Sequence[tuple[int, int, int, int, int]],
        kernel_basis: Sequence[int],
        kernel_fractions: Sequence[tuple[int, int, int, int]],
        lattice: int,
    ) -> None:
        self.base = base
        self.space_operations = space_operations
        self.kernel_basis = kernel_basis
        self.kernel_fractions = kernel_fractions
        self._kernel_fraction_values = tuple(
            self.base._fraction_values(fraction) for fraction in kernel_fractions
        )
        self.lattice = lattice
        self._magnetic_table = magnetic_data.data().table
        self._space_products: dict[tuple[int, int], Pair] = {}
        self._space_fraction_actions: dict[tuple[int, int], int] = {}

    def _fraction_index(
        self,
        values: tuple[Fraction, Fraction, Fraction],
    ) -> int:
        return self.base._kernel_fraction_index_for_delta(
            kernel_basis=self.kernel_basis,  # type: ignore[arg-type]
            kernel_fractions=self.kernel_fractions,
            delta=values,
        )

    def _space_product(self, left_space: int, right_space: int) -> Pair:
        key = (left_space, right_space)
        cached = self._space_products.get(key)
        if cached is not None:
            return cached

        left_record = self.space_operations[left_space - 1]
        right_record = self.space_operations[right_space - 1]
        nonmagnetic_point_op = int(
            self._magnetic_table["mag_point_op_mag2nonmag"][int(left_record[4]) - 1]
        )
        rotated = self.base._rotate_kernel_fraction_by_space_operation(
            tuple(int(value) for value in right_record[:4]),  # type: ignore[arg-type]
            (*left_record[:4], nonmagnetic_point_op),
            lattice=self.lattice,
        )
        left_values = self.base._fraction_values(
            tuple(int(value) for value in left_record[:4])  # type: ignore[arg-type]
        )
        rotated_values = self.base._fraction_values(rotated)
        fraction = self.base._fraction_record_from_values(
            tuple(left_values[axis] + rotated_values[axis] for axis in range(3))  # type: ignore[arg-type]
        )
        point_op = int(
            self._magnetic_table["mag_point_op_mlt"]
            [(int(right_record[4]) - 1) * 144 + (int(left_record[4]) - 1)]
        )
        result = self.base._operation_pair_for_record(
            space_operations=self.space_operations,
            kernel_basis=self.kernel_basis,  # type: ignore[arg-type]
            kernel_fractions=self.kernel_fractions,
            record=(*fraction, point_op),
        )
        self._space_products[key] = result
        return result

    def _space_fraction(self, space_index: int, fraction_index: int) -> int:
        key = (space_index, fraction_index)
        cached = self._space_fraction_actions.get(key)
        if cached is not None:
            return cached

        space_record = self.space_operations[space_index - 1]
        nonmagnetic_point_op = int(
            self._magnetic_table["mag_point_op_mag2nonmag"][int(space_record[4]) - 1]
        )
        rotated = self.base._rotate_kernel_fraction_by_space_operation(
            self.kernel_fractions[fraction_index - 1],
            (*space_record[:4], nonmagnetic_point_op),
            lattice=self.lattice,
        )
        result = self._fraction_index(self.base._fraction_values(rotated))
        self._space_fraction_actions[key] = result
        return result

    def multiply(self, left: Pair, right: Pair) -> Pair:
        space_index, correction_index = self._space_product(left[0], right[0])
        rotated_index = self._space_fraction(left[0], right[1])
        left_values = self._kernel_fraction_values[left[1] - 1]
        correction_values = self._kernel_fraction_values[correction_index - 1]
        rotated_values = self._kernel_fraction_values[rotated_index - 1]
        return space_index, self._fraction_index(
            tuple(
                left_values[axis]
                + correction_values[axis]
                + rotated_values[axis]
                for axis in range(3)
            )
        )


class _MagneticDynamicData:
    """SourceData adapter for magnetic ``get_isotropy_`` dynamic traversal."""

    def __init__(self, base: SourceData, gid: int):
        self._base = base
        self._gid = int(gid)
        self._magnetic_group = base.magnetic_orderparam_group_for_irrep(gid)
        self._ordinary_sg = int(base.little["little_irr_space_group"][int(gid) - 1])
        self.space = base.space
        self.little = base.little
        self.const = base.const
        self._factorized_pair_product: _FactorizedMagneticPairProduct | None = None

    def __getattr__(self, name: str):
        return getattr(self._base, name)

    def generate_space_group_records(self, sg: int) -> tuple[tuple[int, int, int, int, int], ...]:
        if int(sg) != self._ordinary_sg:
            raise ValueError(f"magnetic adapter parent mismatch: {sg} != {self._ordinary_sg}")
        return self._base.generate_magnetic_space_group_records(self._magnetic_group)

    def get_irreps_matrix_for_record(
        self,
        gid: int,
        record: tuple[int, int, int, int, int],
        kparam: Sequence[int],
    ):
        return self._base.get_irreps_matrix_for_magnetic_record(int(gid), record, tuple(int(value) for value in kparam))

    def orderparam_to_subgroup(
        self,
        gid: int,
        kparam: Sequence[int],
        orderparam: Sequence[float],
        row_count: int,
    ) -> tuple[tuple[int, ...], tuple[tuple[int, int, int, int, int], ...]]:
        return self._base.orderparam_to_subgroup_magnetic(
            int(gid),
            tuple(int(value) for value in kparam),
            tuple(float(value) for value in orderparam),
            int(row_count),
        )

    def configure_pair_product(
        self,
        *,
        space_operations: Sequence[tuple[int, int, int, int, int]],
        kernel_basis: Sequence[int],
        kernel_fractions: Sequence[tuple[int, int, int, int]],
        lattice: int,
    ) -> None:
        self._factorized_pair_product = _FactorizedMagneticPairProduct(
            self._base,
            space_operations=space_operations,
            kernel_basis=kernel_basis,
            kernel_fractions=kernel_fractions,
            lattice=lattice,
        )

    def operation_pair_multiply(
        self,
        *,
        space_operations: Sequence[tuple[int, int, int, int, int]],
        kernel_basis: tuple[int, ...] | list[int],
        kernel_fractions: Sequence[tuple[int, int, int, int]],
        lattice: int = 1,
        left: Pair,
        right: Pair,
    ) -> Pair:
        factorized = self._factorized_pair_product
        if factorized is None:
            raise RuntimeError("magnetic pair product is not configured")
        return factorized.multiply(left, right)


def _full_orderparam(dim: int) -> tuple[float, ...]:
    """Return the full OPD space in upstream stride-48 row layout."""

    out = [0.0] * (int(dim) * 48)
    for row in range(int(dim)):
        out[row * 48 + row] = 1.0
    return tuple(out)


def _is_full_orderparam(dim: int, free: int, orderparam: Sequence[float]) -> bool:
    if int(free) != int(dim):
        return False
    for row in range(int(dim)):
        for col in range(int(dim)):
            expected = 1.0 if row == col else 0.0
            if abs(float(orderparam[row * 48 + col]) - expected) > 1e-10:
                return False
    return True


@lru_cache(maxsize=128)
def _all_operation_pairs(space_count: int, fraction_count: int) -> tuple[Pair, ...]:
    return tuple(
        (space_index, fraction_index)
        for space_index in range(1, int(space_count) + 1)
        for fraction_index in range(1, int(fraction_count) + 1)
    )


def _operation_for_pair(
    space_operations: Sequence[tuple[int, int, int, int, int]],
    kernel_fractions: Sequence[tuple[int, int, int, int]],
    pair: Pair,
    cache: OperationCache | None = None,
) -> tuple[int, int, int, int, int]:
    if cache is not None:
        cached = cache.get(pair)
        if cached is not None:
            return cached
    operation = SourceData._vadd_fraction_operation(
        kernel_fractions[int(pair[1]) - 1],
        space_operations[int(pair[0]) - 1],
    )
    if cache is not None:
        cache[pair] = operation
    return operation


def _stride_matrix_for_pair(
    data: SourceData,
    *,
    gid: int,
    kparam: Sequence[int],
    generated_space: Sequence[tuple[int, int, int, int, int]],
    kernel_fractions: Sequence[tuple[int, int, int, int]],
    pair: Pair,
    operation_cache: OperationCache,
    stride_matrix_cache: StrideMatrixCache,
) -> tuple[float, ...]:
    operation = _operation_for_pair(generated_space, kernel_fractions, pair, cache=operation_cache)
    operation_key = tuple(int(value) for value in operation)
    cached = stride_matrix_cache.get(operation_key)
    if cached is not None:
        return cached
    matrix = data.get_irreps_matrix_for_record(int(gid), operation, tuple(int(value) for value in kparam))
    stride_matrix = data._matrix_to_stride48(matrix)
    stride_matrix_cache[operation_key] = stride_matrix
    return stride_matrix


def _pair_product(
    data: SourceData,
    *,
    generated_space: Sequence[tuple[int, int, int, int, int]],
    kernel_basis: Sequence[int],
    kernel_fractions: Sequence[tuple[int, int, int, int]],
    lattice: int,
    cache: PairProductCache,
    left: Pair,
    right: Pair,
) -> Pair:
    key = ((int(left[0]), int(left[1])), (int(right[0]), int(right[1])))
    cached = cache.get(key)
    if cached is not None:
        return cached
    product = data.operation_pair_multiply(
        space_operations=generated_space,
        kernel_basis=kernel_basis,
        kernel_fractions=kernel_fractions,
        lattice=lattice,
        left=left,
        right=right,
    )
    cache[key] = product
    return product


def _candidate_pairs_with_parent_closure(
    data: SourceData,
    *,
    gid: int,
    kparam: Sequence[int],
    dim: int,
    free: int,
    orderparam: Sequence[float],
    generated_space: Sequence[tuple[int, int, int, int, int]],
    kernel_basis: Sequence[int],
    kernel_fractions: Sequence[tuple[int, int, int, int]],
    lattice: int,
    pair_product_cache: PairProductCache,
    operation_cache: OperationCache,
    stride_matrix_cache: StrideMatrixCache,
    initial_flags: dict[Pair, int] | None = None,
) -> tuple[Pair, ...]:
    """Return pair scan order after the parent operation-pair closure pass.

    ``get_isotropy_`` maintains flags over
    ``(ispace_index, kernel_fraction_index)``. Matrix evaluation uses
    ``vadd(kernel_fraction, ispace_op)``, but candidate enumeration marks
    ``selected * stabilizer`` products as excluded before trying the next
    representative. This is the parent-node half of the closure table and does
    not alter OPD matrix semantics.
    """

    pairs = _all_operation_pairs(len(generated_space), len(kernel_fractions))
    if initial_flags is None:
        flags = _stabilizer_pair_flags(
            data,
            gid=gid,
            kparam=kparam,
            dim=dim,
            free=free,
            orderparam=orderparam,
            generated_space=generated_space,
            kernel_fractions=kernel_fractions,
            operation_cache=operation_cache,
            stride_matrix_cache=stride_matrix_cache,
        )
    else:
        flags = {pair: int(initial_flags.get(pair, 0)) for pair in pairs}

    out: list[Pair] = []
    marked_pairs = [pair for pair in pairs if flags.get(pair) == -1]
    for pair in pairs:
        if flags[pair] != 0:
            continue
        for marked in marked_pairs:
            product = _pair_product(
                data,
                generated_space=generated_space,
                kernel_basis=kernel_basis,
                kernel_fractions=kernel_fractions,
                lattice=lattice,
                cache=pair_product_cache,
                left=pair,
                right=marked,
            )
            flags[product] = 1
        out.append(pair)
    return tuple(out)


def _stabilizer_pair_flags(
    data: SourceData,
    *,
    gid: int,
    kparam: Sequence[int],
    dim: int,
    free: int,
    orderparam: Sequence[float],
    generated_space: Sequence[tuple[int, int, int, int, int]],
    kernel_fractions: Sequence[tuple[int, int, int, int]],
    operation_cache: OperationCache,
    stride_matrix_cache: StrideMatrixCache,
) -> dict[Pair, int]:
    pairs = _all_operation_pairs(len(generated_space), len(kernel_fractions))
    flags: dict[Pair, int] = {}
    for pair in pairs:
        stride_matrix = _stride_matrix_for_pair(
            data,
            gid=gid,
            kparam=kparam,
            generated_space=generated_space,
            kernel_fractions=kernel_fractions,
            pair=pair,
            operation_cache=operation_cache,
            stride_matrix_cache=stride_matrix_cache,
        )
        flags[pair] = -1 if orderparam_check(int(dim), int(free), stride_matrix, orderparam) else 0
    return flags


def _close_child_pair_flags(
    data: SourceData,
    *,
    gid: int,
    kparam: Sequence[int],
    dim: int,
    free: int,
    orderparam: Sequence[float],
    generated_space: Sequence[tuple[int, int, int, int, int]],
    kernel_basis: Sequence[int],
    kernel_fractions: Sequence[tuple[int, int, int, int]],
    lattice: int,
    parent_flags: dict[Pair, int],
    selected_pair: Pair,
    seen: set[tuple[int, tuple[float, ...]]],
    pair_product_cache: PairProductCache,
    operation_cache: OperationCache,
    stride_matrix_cache: StrideMatrixCache,
    stride_nonzero_cache: StrideNonzeroCache,
    orbit_key_by_stride_id: dict[int, tuple[int, tuple[float, ...]]] | None = None,
) -> dict[Pair, int]:
    """Propagate operation-pair closure flags to a newly created child node."""

    pairs = _all_operation_pairs(len(generated_space), len(kernel_fractions))
    flags = {pair: 0 for pair in pairs}
    minus_list: list[Pair] = []

    def mark_minus(pair: Pair) -> None:
        if flags[pair] == -1:
            return
        flags[pair] = -1
        minus_list.append(pair)

    check_cache: dict[int, bool] = {}
    transformed_key_cache: dict[int, tuple[int, tuple[float, ...]]] = {}

    for marked in pairs:
        if parent_flags.get(marked) != -1:
            continue
        mark_minus(marked)
        mark_minus(
            _pair_product(
                data,
                generated_space=generated_space,
                kernel_basis=kernel_basis,
                kernel_fractions=kernel_fractions,
                lattice=lattice,
                cache=pair_product_cache,
                left=selected_pair,
                right=marked,
            )
        )
        mark_minus(
            _pair_product(
                data,
                generated_space=generated_space,
                kernel_basis=kernel_basis,
                kernel_fractions=kernel_fractions,
                lattice=lattice,
                cache=pair_product_cache,
                left=marked,
                right=selected_pair,
            )
        )

    for pair in pairs:
        if flags[pair] != 0:
            continue
        stride_matrix = _stride_matrix_for_pair(
            data,
            gid=gid,
            kparam=kparam,
            generated_space=generated_space,
            kernel_fractions=kernel_fractions,
            pair=pair,
            operation_cache=operation_cache,
            stride_matrix_cache=stride_matrix_cache,
        )
        stride_id = id(stride_matrix)
        preserves = check_cache.get(stride_id)
        if preserves is None:
            preserves = bool(orderparam_check(int(dim), int(free), stride_matrix, orderparam))
            check_cache[stride_id] = preserves
        if preserves:
            for marked in tuple(minus_list):
                mark_minus(
                    _pair_product(
                        data,
                        generated_space=generated_space,
                        kernel_basis=kernel_basis,
                        kernel_fractions=kernel_fractions,
                        lattice=lattice,
                        cache=pair_product_cache,
                        left=pair,
                        right=marked,
                    )
                )
                mark_minus(
                    _pair_product(
                        data,
                        generated_space=generated_space,
                        kernel_basis=kernel_basis,
                        kernel_fractions=kernel_fractions,
                        lattice=lattice,
                        cache=pair_product_cache,
                        left=marked,
                        right=pair,
                    )
                )
            continue

        key = None
        if orbit_key_by_stride_id is not None:
            key = orbit_key_by_stride_id.get(stride_id)
        if key is None:
            key = transformed_key_cache.get(stride_id)
        if key is None:
            transformed = _transform_orderparam_sparse(
                int(dim),
                int(free),
                orderparam,
                stride_matrix,
                stride_nonzero_cache,
            )
            eq_count, equations, key = _orderparam_equation_data(
                int(dim),
                int(free),
                transformed,
            )
            transformed_key_cache[stride_id] = key
        if key not in seen:
            continue
        for marked in tuple(minus_list):
            flags[
                _pair_product(
                    data,
                    generated_space=generated_space,
                    kernel_basis=kernel_basis,
                    kernel_fractions=kernel_fractions,
                    lattice=lattice,
                    cache=pair_product_cache,
                    left=pair,
                    right=marked,
                )
            ] = 1
    return flags


@lru_cache(maxsize=200_000)
def _reduced_equation_key_cached(
    dim: int,
    eq_count: int,
    active_equations: tuple[float, ...],
) -> tuple[int, tuple[float, ...]]:
    required = max(0, (int(dim) - 1) * 50 + int(eq_count))
    equations = [0.0] * required
    offset = 0
    for col in range(int(dim)):
        for row in range(int(eq_count)):
            equations[row + col * 50] = float(active_equations[offset])
            offset += 1
    reduced = equations if _equations_are_reduced(int(dim), int(eq_count), equations) else xrowop2(equations, int(eq_count), int(dim))
    active: list[float] = []
    for col in range(int(dim)):
        for row in range(int(eq_count)):
            active.append(round(float(reduced[row + col * 50]), 10))
    return int(eq_count), tuple(active)


def _equations_are_reduced(dim: int, eq_count: int, equations: Sequence[float]) -> bool:
    pivot_cols: list[int] = []
    for row in range(int(eq_count)):
        pivot_col = None
        for col in range(int(dim)):
            value = float(equations[row + col * 50]) if row + col * 50 < len(equations) else 0.0
            if abs(value) > 1e-6:
                pivot_col = col
                if abs(value - 1.0) > 1e-6:
                    return False
                break
        if pivot_col is None:
            return False
        if pivot_cols and pivot_col <= pivot_cols[-1]:
            return False
        pivot_cols.append(pivot_col)
        for other in range(int(eq_count)):
            if other == row:
                continue
            value = float(equations[other + pivot_col * 50]) if other + pivot_col * 50 < len(equations) else 0.0
            if abs(value) > 1e-6:
                return False
    return True


def _reduced_equation_key(dim: int, eq_count: int, equations: Sequence[float]) -> tuple[int, tuple[float, ...]]:
    """Canonical key used by ``get_isotropy_``'s point_simeqs list."""

    return _reduced_equation_key_cached(
        int(dim),
        int(eq_count),
        tuple(
            round(
                float(equations[row + col * 50])
                if row + col * 50 < len(equations)
                else 0.0,
                12,
            )
            for col in range(int(dim))
            for row in range(int(eq_count))
        ),
    )


@lru_cache(maxsize=4096)
def _orderparam_equation_data_cached(
    dim: int,
    free: int,
    active_orderparam: tuple[float, ...],
) -> tuple[int, tuple[float, ...], tuple[int, tuple[float, ...]]]:
    """Memoize the exact OPD-to-equation boundary without changing its key."""

    eq_count, equations = orderparam_to_eqs(
        int(dim),
        int(free),
        active_orderparam,
        op_stride=int(dim),
    )
    return (
        eq_count,
        equations,
        _reduced_equation_key(int(dim), eq_count, equations),
    )


def _active_orderparam_values(
    dim: int,
    free: int,
    orderparam: Sequence[float],
) -> tuple[float, ...]:
    return tuple(
        float(orderparam[row * 48 + col])
        if row * 48 + col < len(orderparam)
        else 0.0
        for row in range(int(free))
        for col in range(int(dim))
    )


def _orderparam_equation_data(
    dim: int,
    free: int,
    orderparam: Sequence[float],
) -> tuple[int, tuple[float, ...], tuple[int, tuple[float, ...]]]:
    active = _active_orderparam_values(int(dim), int(free), orderparam)
    return _orderparam_equation_data_cached(int(dim), int(free), active)


def _equation_key(dim: int, free: int, orderparam: Sequence[float]) -> tuple[int, tuple[float, ...]]:
    """Canonical key for one OPD subspace."""

    return _orderparam_equation_data(int(dim), int(free), orderparam)[2]


def _orderparam_from_equations(
    dim: int,
    eq_count: int,
    equations: Sequence[float],
    *,
    reduce_first: bool,
) -> tuple[int, tuple[float, ...]]:
    """Return the OPD reconstructed from its equation-row representation."""

    source = xrowop2(equations, int(eq_count), int(dim)) if reduce_first else equations
    return eqs_to_orderparam(int(dim), int(eq_count), source)


def _transform_orderparam(
    dim: int,
    free: int,
    orderparam: Sequence[float],
    irrep_matrix: Sequence[float],
) -> tuple[float, ...]:
    """Apply one irrep matrix to an OPD row set in upstream stride-48 layout."""

    out = [0.0] * (int(free) * 48)
    for row in range(int(free)):
        for col in range(int(dim)):
            total = 0.0
            for source_row in range(int(dim)):
                matrix_index = source_row * 48 + col
                op_index = row * 48 + source_row
                matrix_value = float(irrep_matrix[matrix_index]) if matrix_index < len(irrep_matrix) else 0.0
                op_value = float(orderparam[op_index]) if op_index < len(orderparam) else 0.0
                total += matrix_value * op_value
            out[row * 48 + col] = total
    return tuple(out)


def _stride_nonzeros(
    dim: int,
    irrep_matrix: Sequence[float],
    cache: StrideNonzeroCache,
) -> tuple[tuple[tuple[int, float], ...], ...]:
    key = id(irrep_matrix)
    cached = cache.get(key)
    if cached is not None:
        return cached
    columns: list[tuple[tuple[int, float], ...]] = []
    for col in range(int(dim)):
        entries: list[tuple[int, float]] = []
        for source_row in range(int(dim)):
            matrix_index = source_row * 48 + col
            matrix_value = float(irrep_matrix[matrix_index]) if matrix_index < len(irrep_matrix) else 0.0
            if matrix_value != 0.0:
                entries.append((source_row, matrix_value))
        columns.append(tuple(entries))
    cached = tuple(columns)
    cache[key] = cached
    return cached


def _transform_orderparam_sparse(
    dim: int,
    free: int,
    orderparam: Sequence[float],
    irrep_matrix: Sequence[float],
    nonzero_cache: StrideNonzeroCache,
) -> tuple[float, ...]:
    out = [0.0] * (int(free) * 48)
    nonzeros = _stride_nonzeros(int(dim), irrep_matrix, nonzero_cache)
    for row in range(int(free)):
        row_offset = row * 48
        for col, entries in enumerate(nonzeros):
            total = 0.0
            for source_row, matrix_value in entries:
                op_index = row_offset + source_row
                op_value = float(orderparam[op_index]) if op_index < len(orderparam) else 0.0
                total += matrix_value * op_value
            out[row_offset + col] = total
    return tuple(out)


def _orbit_keys_and_representative(
    data: SourceData,
    *,
    gid: int,
    kparam: Sequence[int],
    dim: int,
    free: int,
    orderparam: Sequence[float],
    generated_space: Sequence[tuple[int, int, int, int, int]],
    kernel_fractions: Sequence[tuple[int, int, int, int]],
    operation_cache: OperationCache,
    stride_matrix_cache: StrideMatrixCache,
    stride_nonzero_cache: StrideNonzeroCache,
) -> tuple[set[tuple[int, tuple[float, ...]]], tuple[float, ...], dict[int, tuple[int, tuple[float, ...]]]]:
    """Return point-simeqs keys and the ``nice_orderparam_`` orbit representative."""

    base_active = _active_orderparam_values(int(dim), int(free), orderparam)
    base_equations = _orderparam_equation_data_cached(int(dim), int(free), base_active)
    keys = {base_equations[2]}
    key_by_stride_id: dict[int, tuple[int, tuple[float, ...]]] = {}
    representative = tuple(float(value) for value in orderparam)
    if _is_full_orderparam(int(dim), int(free), representative):
        return keys, representative, key_by_stride_id
    key_by_stride: dict[tuple[float, ...], tuple[int, tuple[float, ...]]] = {}
    key_by_active = {base_active: base_equations[2]}
    active_scalar_count = len(base_active)
    for pair in _all_operation_pairs(len(generated_space), len(kernel_fractions)):
        stride_matrix = _stride_matrix_for_pair(
            data,
            gid=gid,
            kparam=kparam,
            generated_space=generated_space,
            kernel_fractions=kernel_fractions,
            pair=pair,
            operation_cache=operation_cache,
            stride_matrix_cache=stride_matrix_cache,
        )
        stride_id = id(stride_matrix)
        key = key_by_stride.get(stride_matrix)
        if key is None:
            transformed = _transform_orderparam_sparse(
                int(dim),
                int(free),
                orderparam,
                stride_matrix,
                stride_nonzero_cache,
            )
            active = _active_orderparam_values(int(dim), int(free), transformed)
            key = key_by_active.get(active)
            if key is None:
                eq_count, equations, key = _orderparam_equation_data_cached(int(dim), int(free), active)
                if active_scalar_count + len(active) <= _ORBIT_ACTIVE_MEMO_SCALAR_BUDGET:
                    key_by_active[active] = key
                    active_scalar_count += len(active)
                if key not in keys:
                    keys.add(key)
                    next_free, canonical = _orderparam_from_equations(
                        int(dim),
                        eq_count,
                        equations,
                        reduce_first=True,
                    )
                    if next_free == int(free) and nice_orderparam(
                        int(dim),
                        int(free),
                        representative,
                        canonical,
                    ) == 2:
                        representative = tuple(canonical)
            key_by_stride[stride_matrix] = key
        key_by_stride_id[stride_id] = key
    return keys, representative, key_by_stride_id


def _fmt_matrix_value(value: float) -> str:
    value = float(value)
    if abs(value) < 5e-7:
        value = 0.0
    rounded = round(value)
    if abs(value - rounded) < 5e-7:
        return str(int(rounded))
    return f"{value:.6f}".rstrip("0").rstrip(".")


def _row_matrix(free: int, dim: int, orderparam: Sequence[float]) -> tuple[tuple[str, ...], ...]:
    """Return the orthogonalized matrix written to dynamic ``*.iso`` rows."""

    _, orthogonal = orthogonal_orderparam(int(free), int(dim), orderparam)
    return tuple(
        tuple(_fmt_matrix_value(orthogonal[row * 48 + col]) for col in range(int(dim)))
        for row in range(int(free))
    )


def _row_from_orderparam(
    data: SourceData,
    *,
    parent_sg: int,
    gid: int,
    kparam: Sequence[int],
    free: int,
    orderparam: Sequence[float],
    counts: dict[int, int],
) -> DynamicIsotropyRow | None:
    basis, operations = data.orderparam_to_subgroup(gid, tuple(int(value) for value in kparam), tuple(orderparam), free)
    result = id_subgroup_identify_with_generator_block(data, parent_sg, basis, operations, 0)
    if result is None:
        return None
    counts[int(free)] = counts.get(int(free), 0) + 1
    return DynamicIsotropyRow(
        subgroup_number=int(result.subgroup),
        free=int(free),
        basis_values=tuple(int(value) for value in result.basis),
        origin_values=tuple(int(value) for value in result.origin),
        direction=orderparam_label(int(free), counts[int(free)]),
        matrix=_row_matrix(int(free), int(data.little["little_irr_full_dim"][int(gid) - 1]), orderparam),
    )


def _row_from_orderparam_magnetic(
    data: SourceData,
    *,
    parent_sg: int,
    gid: int,
    kparam: Sequence[int],
    free: int,
    orderparam: Sequence[float],
    counts: dict[int, int],
) -> DynamicIsotropyRow | None:
    basis, operations = data.orderparam_to_subgroup_magnetic(
        int(gid),
        tuple(int(value) for value in kparam),
        tuple(float(value) for value in orderparam),
        int(free),
    )
    result = id_subgroup_magnetic_identify_with_generator_block(
        data,
        int(parent_sg),
        basis,
        operations,
        0,
    )
    if result is None:
        return None
    counts[int(free)] = counts.get(int(free), 0) + 1
    return DynamicIsotropyRow(
        subgroup_number=int(result.magnetic_group),
        free=int(free),
        basis_values=tuple(int(value) for value in result.basis),
        origin_values=tuple(int(value) for value in result.origin),
        direction=orderparam_label(int(free), counts[int(free)]),
        matrix=_row_matrix(int(free), int(data.little["little_irr_full_dim"][int(gid) - 1]), orderparam),
        magnetic_subgroup_selection=result.selection,
        magnetic_operation_records=tuple(tuple(int(value) for value in record) for record in operations),
    )


def generate_dynamic_isotropy_rows(
    data: SourceData,
    *,
    gid: int,
    kparam: Sequence[int],
    _point_occurrence: int | None = None,
) -> tuple[DynamicIsotropyRow, ...]:
    """Generate dynamic ``*.iso`` rows for one nonmagnetic irrep.

    Rows are returned in upstream discovery order.  Use
    ``sort_dynamic_rows_for_file`` when comparing with final dynamic file order.
    """

    gid = int(gid)
    if _point_occurrence is not None and int(_point_occurrence) < 1:
        raise ValueError("point occurrence must be positive")
    parent_sg = int(data.little["little_irr_space_group"][gid - 1])
    lattice = int(data.space["ispace_lattice"][parent_sg - 1])
    dim = int(data.little["little_irr_full_dim"][gid - 1])
    generated_space = data.generate_space_group_records(parent_sg)
    rows: list[DynamicIsotropyRow] = []
    counts: dict[int, int] = {}
    seen: set[tuple[int, tuple[float, ...]]] = set()
    queue: deque[tuple[int, tuple[float, ...], tuple[float, ...], dict[Pair, int]]] = deque()
    pair_product_cache: PairProductCache = {}
    operation_cache: OperationCache = {}
    stride_matrix_cache: StrideMatrixCache = {}
    stride_nonzero_cache: StrideNonzeroCache = {}

    first = _full_orderparam(dim)
    first_basis, _ = data.orderparam_to_subgroup(gid, tuple(int(value) for value in kparam), first, dim)
    kernel_fractions = get_new_fractionals(first_basis)
    first_keys, first_rep, _ = _orbit_keys_and_representative(
        data,
        gid=gid,
        kparam=kparam,
        dim=dim,
        free=dim,
        orderparam=first,
        generated_space=generated_space,
        kernel_fractions=kernel_fractions,
        operation_cache=operation_cache,
        stride_matrix_cache=stride_matrix_cache,
        stride_nonzero_cache=stride_nonzero_cache,
    )
    seen.update(first_keys)
    first_flags = _stabilizer_pair_flags(
        data,
        gid=gid,
        kparam=kparam,
        dim=dim,
        free=dim,
        orderparam=first_rep,
        generated_space=generated_space,
        kernel_fractions=kernel_fractions,
        operation_cache=operation_cache,
        stride_matrix_cache=stride_matrix_cache,
    )
    queue.append((dim, first, first_rep, first_flags))

    while queue:
        free, node_orderparam, row_orderparam, pair_flags = queue.popleft()
        row = (
            _row_from_orderparam(
                data,
                parent_sg=parent_sg,
                gid=gid,
                kparam=kparam,
                free=free,
                orderparam=row_orderparam,
                counts=counts,
            )
            # A selected P occurrence only depends on accepted free=1 rows.
            # Higher-free subgroup presentation does not affect BFS expansion.
            if _point_occurrence is None or int(free) == 1
            else None
        )
        if row is not None:
            rows.append(row)
            # Dynamic file labels count accepted free=1 rows in Source BFS
            # discovery order.  Once P<n> is accepted, later nodes cannot
            # change that selected row or its final file-order ordinal.
            if (
                _point_occurrence is not None
                and int(row.free) == 1
                and counts.get(1) == int(_point_occurrence)
            ):
                return tuple(rows)

        for pair in _candidate_pairs_with_parent_closure(
            data,
            gid=gid,
            kparam=kparam,
            dim=dim,
            free=free,
            orderparam=node_orderparam,
            generated_space=generated_space,
            kernel_basis=first_basis,
            kernel_fractions=kernel_fractions,
            lattice=lattice,
            pair_product_cache=pair_product_cache,
            operation_cache=operation_cache,
            stride_matrix_cache=stride_matrix_cache,
            initial_flags=pair_flags,
        ):
            stride_matrix = _stride_matrix_for_pair(
                data,
                gid=gid,
                kparam=kparam,
                generated_space=generated_space,
                kernel_fractions=kernel_fractions,
                pair=pair,
                operation_cache=operation_cache,
                stride_matrix_cache=stride_matrix_cache,
            )
            if orderparam_check(dim, free, stride_matrix, node_orderparam):
                continue
            eq_count, equations, _ = _orderparam_equation_data(
                dim,
                free,
                node_orderparam,
            )
            next_eq_count, next_equations = orderparam_add_eqs(dim, eq_count, equations, stride_matrix)
            if next_eq_count <= eq_count or next_eq_count >= dim:
                continue
            key = _reduced_equation_key(dim, next_eq_count, next_equations)
            if key in seen:
                continue
            next_free, next_orderparam = _orderparam_from_equations(
                dim,
                next_eq_count,
                next_equations,
                reduce_first=False,
            )
            orbit_keys, representative, orbit_key_by_stride_id = _orbit_keys_and_representative(
                data,
                gid=gid,
                kparam=kparam,
                dim=dim,
                free=next_free,
                orderparam=next_orderparam,
                generated_space=generated_space,
                kernel_fractions=kernel_fractions,
                operation_cache=operation_cache,
                stride_matrix_cache=stride_matrix_cache,
                stride_nonzero_cache=stride_nonzero_cache,
            )
            if orbit_keys & seen:
                continue
            child_flags = _close_child_pair_flags(
                data,
                gid=gid,
                kparam=kparam,
                dim=dim,
                free=next_free,
                orderparam=next_orderparam,
                generated_space=generated_space,
                kernel_basis=first_basis,
                kernel_fractions=kernel_fractions,
                lattice=lattice,
                parent_flags=pair_flags,
                selected_pair=pair,
                seen=seen,
                pair_product_cache=pair_product_cache,
                operation_cache=operation_cache,
                stride_matrix_cache=stride_matrix_cache,
                stride_nonzero_cache=stride_nonzero_cache,
                orbit_key_by_stride_id=orbit_key_by_stride_id,
            )
            seen.update(orbit_keys)
            queue.append((next_free, next_orderparam, representative, child_flags))

    return tuple(rows)


def generate_dynamic_point_isotropy_row(
    data: SourceData,
    *,
    gid: int,
    kparam: Sequence[int],
    occurrence: int,
) -> DynamicIsotropyRow | None:
    """Return exact Source-selected ``P<n>`` without completing later BFS nodes."""

    occurrence = int(occurrence)
    rows = generate_dynamic_isotropy_rows(
        data,
        gid=int(gid),
        kparam=kparam,
        _point_occurrence=occurrence,
    )
    label = orderparam_label(1, occurrence)
    return next(
        (row for row in rows if int(row.free) == 1 and row.direction == label),
        None,
    )

def generate_dynamic_isotropy_rows_magnetic(
    data: SourceData,
    *,
    gid: int,
    kparam: Sequence[int],
    _point_occurrence: int | None = None,
) -> tuple[DynamicIsotropyRow, ...]:
    """Generate dynamic magnetic ``DISPLAY ISOTROPY`` rows.

    This follows the same OPD-subspace traversal as the ordinary dynamic
    generator, but evaluates matrices, operation products, and subgroup
    identification through magnetic Source-table operations.
    """

    gid = int(gid)
    if _point_occurrence is not None and int(_point_occurrence) < 1:
        raise ValueError("point occurrence must be positive")
    magnetic_data = _MagneticDynamicData(data, gid)
    parent_sg = int(data.little["little_irr_space_group"][gid - 1])
    lattice = int(data.space["ispace_lattice"][parent_sg - 1])
    dim = int(data.little["little_irr_full_dim"][gid - 1])
    generated_space = magnetic_data.generate_space_group_records(parent_sg)
    rows: list[DynamicIsotropyRow] = []
    counts: dict[int, int] = {}
    seen: set[tuple[int, tuple[float, ...]]] = set()
    queue: deque[tuple[int, tuple[float, ...], tuple[float, ...], dict[Pair, int]]] = deque()
    pair_product_cache: PairProductCache = {}
    operation_cache: OperationCache = {}
    stride_matrix_cache: StrideMatrixCache = {}
    stride_nonzero_cache: StrideNonzeroCache = {}

    first = _full_orderparam(dim)
    first_basis, _ = magnetic_data.orderparam_to_subgroup(gid, tuple(int(value) for value in kparam), first, dim)
    kernel_fractions = get_new_fractionals(first_basis)
    magnetic_data.configure_pair_product(
        space_operations=generated_space,
        kernel_basis=first_basis,
        kernel_fractions=kernel_fractions,
        lattice=lattice,
    )
    first_keys, first_rep, _ = _orbit_keys_and_representative(
        magnetic_data,  # type: ignore[arg-type]
        gid=gid,
        kparam=kparam,
        dim=dim,
        free=dim,
        orderparam=first,
        generated_space=generated_space,
        kernel_fractions=kernel_fractions,
        operation_cache=operation_cache,
        stride_matrix_cache=stride_matrix_cache,
        stride_nonzero_cache=stride_nonzero_cache,
    )
    seen.update(first_keys)
    first_flags = _stabilizer_pair_flags(
        magnetic_data,  # type: ignore[arg-type]
        gid=gid,
        kparam=kparam,
        dim=dim,
        free=dim,
        orderparam=first_rep,
        generated_space=generated_space,
        kernel_fractions=kernel_fractions,
        operation_cache=operation_cache,
        stride_matrix_cache=stride_matrix_cache,
    )
    queue.append((dim, first, first_rep, first_flags))

    while queue:
        free, node_orderparam, row_orderparam, pair_flags = queue.popleft()
        row = (
            _row_from_orderparam_magnetic(
                data,
                parent_sg=parent_sg,
                gid=gid,
                kparam=kparam,
                free=free,
                orderparam=row_orderparam,
                counts=counts,
            )
            if _point_occurrence is None or int(free) == 1
            else None
        )
        if row is not None:
            rows.append(row)
            # Magnetic dynamic rows use the same Source P<n> discovery
            # count as ordinary rows, so the accepted occurrence is final.
            if (
                _point_occurrence is not None
                and int(row.free) == 1
                and counts.get(1) == int(_point_occurrence)
            ):
                return tuple(rows)

        for pair in _candidate_pairs_with_parent_closure(
            magnetic_data,  # type: ignore[arg-type]
            gid=gid,
            kparam=kparam,
            dim=dim,
            free=free,
            orderparam=node_orderparam,
            generated_space=generated_space,
            kernel_basis=first_basis,
            kernel_fractions=kernel_fractions,
            lattice=lattice,
            pair_product_cache=pair_product_cache,
            operation_cache=operation_cache,
            stride_matrix_cache=stride_matrix_cache,
            initial_flags=pair_flags,
        ):
            stride_matrix = _stride_matrix_for_pair(
                magnetic_data,  # type: ignore[arg-type]
                gid=gid,
                kparam=kparam,
                generated_space=generated_space,
                kernel_fractions=kernel_fractions,
                pair=pair,
                operation_cache=operation_cache,
                stride_matrix_cache=stride_matrix_cache,
            )
            if orderparam_check(dim, free, stride_matrix, node_orderparam):
                continue
            eq_count, equations, _ = _orderparam_equation_data(
                dim,
                free,
                node_orderparam,
            )
            next_eq_count, next_equations = orderparam_add_eqs(dim, eq_count, equations, stride_matrix)
            if next_eq_count <= eq_count or next_eq_count >= dim:
                continue
            key = _reduced_equation_key(dim, next_eq_count, next_equations)
            if key in seen:
                continue
            next_free, next_orderparam = _orderparam_from_equations(
                dim,
                next_eq_count,
                next_equations,
                reduce_first=False,
            )
            orbit_keys, representative, orbit_key_by_stride_id = _orbit_keys_and_representative(
                magnetic_data,  # type: ignore[arg-type]
                gid=gid,
                kparam=kparam,
                dim=dim,
                free=next_free,
                orderparam=next_orderparam,
                generated_space=generated_space,
                kernel_fractions=kernel_fractions,
                operation_cache=operation_cache,
                stride_matrix_cache=stride_matrix_cache,
                stride_nonzero_cache=stride_nonzero_cache,
            )
            if orbit_keys & seen:
                continue
            child_flags = _close_child_pair_flags(
                magnetic_data,  # type: ignore[arg-type]
                gid=gid,
                kparam=kparam,
                dim=dim,
                free=next_free,
                orderparam=next_orderparam,
                generated_space=generated_space,
                kernel_basis=first_basis,
                kernel_fractions=kernel_fractions,
                lattice=lattice,
                parent_flags=pair_flags,
                selected_pair=pair,
                seen=seen,
                pair_product_cache=pair_product_cache,
                operation_cache=operation_cache,
                stride_matrix_cache=stride_matrix_cache,
                stride_nonzero_cache=stride_nonzero_cache,
                orbit_key_by_stride_id=orbit_key_by_stride_id,
            )
            seen.update(orbit_keys)
            queue.append((next_free, next_orderparam, representative, child_flags))

    return tuple(rows)


def generate_dynamic_point_isotropy_row_magnetic(
    data: SourceData,
    *,
    gid: int,
    kparam: Sequence[int],
    occurrence: int,
) -> DynamicIsotropyRow | None:
    """Return exact magnetic ``P<n>`` without completing later BFS nodes."""

    occurrence = int(occurrence)
    rows = generate_dynamic_isotropy_rows_magnetic(
        data,
        gid=int(gid),
        kparam=kparam,
        _point_occurrence=occurrence,
    )
    label = orderparam_label(1, occurrence)
    return next(
        (row for row in rows if int(row.free) == 1 and row.direction == label),
        None,
    )
