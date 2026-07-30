"""Exact fixed-K OPD/domain restriction for ``DISPLAY INVARIANT``.

The binary does not replace the parent action by the setwise stabilizer of a
selected OPD.  For every synchronized parent-group matrix ``D``, it projects
``D`` onto the selected, domain-transported direction columns ``B``:

``(B^T B)^-1 B^T D B``.

The direction columns are required to be mutually orthogonal, matching the
explicit check in ``invar_subspace_``.  The projected matrices need not form a
group; they are the exact matrix list averaged by the downstream invariant
kernel.
"""

from __future__ import annotations

from fractions import Fraction
from typing import Sequence

from .domains import (
    domain_count_from_isotropy_row,
    domain_count_from_subgroup,
    domain_operation_record_from_isotropy_row,
    domain_operation_record_from_subgroup,
)
from .algebra import ExactMatrix, _coefficient
from .authority import InvariantSource, invariant_source
from .source import (
    _exact_record_matrix,
    _quadratic3_record_matrix,
    _rational_entry,
    _source_display_entry,
    _source_gid,
    coupled_fixed_irrep_matrices,
    coupled_parametric_irrep_matrices,
)


def _exact_direction_matrix(
    projection: InvariantSource,
    gid: int,
    direction: str,
) -> tuple[int, ExactMatrix]:
    row_id, raw = projection.direction_matrix_by_gid_label(int(gid), str(direction))
    if row_id is None:
        raise ValueError("fixed DISPLAY INVARIANT requires a Source OPD row")
    matrix = tuple(
        tuple(_rational_entry(complex(value)) for value in row)
        for row in raw
    )
    if not matrix or not matrix[0]:
        raise ValueError("OPD direction matrix is empty")
    width = len(matrix[0])
    if any(len(row) != width for row in matrix):
        raise ValueError("OPD direction matrix is ragged")
    return int(row_id), matrix


def _exact_dynamic_direction_matrix(
    projection: InvariantSource,
    gid: int,
    direction: str,
    kparam: Sequence[float | Fraction | int],
) -> ExactMatrix:
    row_id, raw = projection.direction_matrix_by_gid_label(
        int(gid),
        str(direction),
        kparam=kparam,
    )
    if row_id is not None:
        raise ValueError("parametric DISPLAY INVARIANT expected a dynamic OPD row")
    matrix = tuple(
        tuple(_source_display_entry(value) for value in row)
        for row in raw
    )
    if not matrix or not matrix[0] or any(len(row) != len(matrix[0]) for row in matrix):
        raise ValueError("dynamic OPD direction matrix is empty or ragged")
    return matrix


def _rectangular_block_diagonal(blocks: Sequence[ExactMatrix]) -> ExactMatrix:
    row_count = sum(len(block) for block in blocks)
    column_count = sum(len(block[0]) for block in blocks)
    rows = [[Fraction(0) for _column in range(column_count)] for _row in range(row_count)]
    row_offset = 0
    column_offset = 0
    for block in blocks:
        width = len(block[0])
        for row, values in enumerate(block):
            if len(values) != width:
                raise ValueError("OPD direction matrix is ragged")
            for column, value in enumerate(values):
                rows[row_offset + row][column_offset + column] = _coefficient(value)  # type: ignore[assignment]
        row_offset += len(block)
        column_offset += width
    return tuple(tuple(row) for row in rows)


def _matrix_product(left: ExactMatrix, right: ExactMatrix) -> ExactMatrix:
    if not left or not right or len(left[0]) != len(right):
        raise ValueError("matrix dimensions differ")
    return tuple(
        tuple(
            sum(left[row][inner] * right[inner][column] for inner in range(len(right)))
            for column in range(len(right[0]))
        )
        for row in range(len(left))
    )


def _project_action(
    full_action: Sequence[ExactMatrix],
    direction_blocks: Sequence[ExactMatrix],
) -> tuple[tuple[ExactMatrix, ...], tuple[int, ...]]:
    direction_matrix = _rectangular_block_diagonal(direction_blocks)
    row_count = len(direction_matrix)
    column_count = len(direction_matrix[0])
    norms = []
    for column in range(column_count):
        norm = sum(
            (
                direction_matrix[row][column] * direction_matrix[row][column]  # type: ignore[operator]
                for row in range(row_count)
            ),
            Fraction(0),
        )
        norms.append(_coefficient(norm))
    for left in range(column_count):
        for right in range(left):
            overlap = sum(
                (
                    direction_matrix[row][left] * direction_matrix[row][right]  # type: ignore[operator]
                    for row in range(row_count)
                ),
                Fraction(0),
            )
            if overlap:
                raise ValueError("combined OPD direction columns are not orthogonal")
    if any(norm == 0 for norm in norms):
        raise ValueError("combined OPD direction contains a zero column")

    projected = []
    for full_matrix in full_action:
        transformed = _matrix_product(full_matrix, direction_matrix)
        rows = []
        for left in range(column_count):
            values = []
            for right in range(column_count):
                coefficient = sum(
                    (
                        direction_matrix[row][left] * transformed[row][right]  # type: ignore[operator]
                        for row in range(row_count)
                    ),
                    Fraction(0),
                )
                values.append(coefficient / norms[left])  # type: ignore[operator]
            rows.append(tuple(values))
        projected.append(tuple(rows))
    return tuple(projected), tuple(len(block[0]) for block in direction_blocks)


def restricted_fixed_irrep_action(
    space_group: int,
    irreps: Sequence[str],
    directions: Sequence[str],
    domains: Sequence[int],
    *,
    projection_source: InvariantSource | None = None,
) -> tuple[tuple[ExactMatrix, ...], tuple[int, ...]]:
    """Return binary-projected matrices and free dimensions for fixed irreps."""

    labels = tuple(str(value) for value in irreps)
    opd_labels = tuple(str(value) for value in directions)
    domain_values = tuple(domains)
    if not labels or not (len(labels) == len(opd_labels) == len(domain_values)):
        raise ValueError("irreps, directions, and domains must have the same nonzero length")
    if any(type(value) is not int or value < 1 for value in domain_values):
        raise ValueError("domain values must be positive integers")

    projection = projection_source or invariant_source()
    data = projection.source_data
    gids = tuple(_source_gid(data, int(space_group), label) for label in labels)
    direction_blocks: list[ExactMatrix] = []
    block_dimensions: list[int] = []
    for gid, direction, domain in zip(gids, opd_labels, domain_values, strict=True):
        row_id, block = _exact_direction_matrix(projection, gid, direction)
        upper = domain_count_from_isotropy_row(projection, data, int(space_group), row_id)
        if domain > upper:
            raise ValueError(f"domain {domain} exceeds the Source upper bound {upper}")
        if domain > 1:
            record = domain_operation_record_from_isotropy_row(
                projection,
                data,
                sg=int(space_group),
                row_id=row_id,
                domain=domain,
            )
            if record is None:
                raise ValueError(f"Source domain {domain} has no operation record")
            block = _matrix_product(_exact_record_matrix(data, gid, record), block)
        direction_blocks.append(block)
        block_dimensions.append(len(block[0]))

    full_action = coupled_fixed_irrep_matrices(
        int(space_group),
        labels,
        source_data=data,
    )
    projected, projected_dimensions = _project_action(full_action, direction_blocks)
    if projected_dimensions != tuple(block_dimensions):
        raise AssertionError("projected OPD dimensions changed")
    return projected, projected_dimensions


def restricted_parametric_irrep_action(
    space_group: int,
    irreps: Sequence[str],
    k_parameters: Sequence[Sequence[float | Fraction | int]],
    directions: Sequence[str],
    domains: Sequence[int],
    *,
    projection_source: InvariantSource | None = None,
) -> tuple[tuple[ExactMatrix, ...], tuple[int, ...]]:
    """Return the binary OPD/domain projection when any factor is parametric."""

    labels = tuple(str(value) for value in irreps)
    parameter_rows = tuple(tuple(values) for values in k_parameters)
    opd_labels = tuple(str(value) for value in directions)
    domain_values = tuple(domains)
    if not labels or not (
        len(labels) == len(parameter_rows) == len(opd_labels) == len(domain_values)
    ):
        raise ValueError("parametric factor fields must have the same nonzero length")
    if any(type(value) is not int or value < 1 for value in domain_values):
        raise ValueError("domain values must be positive integers")

    projection = projection_source or invariant_source()
    data = projection.source_data
    gids = tuple(_source_gid(data, int(space_group), label) for label in labels)
    source_parameters = []
    for gid, values in zip(gids, parameter_rows, strict=True):
        little = projection.little_record_by_gid(gid)
        if int(little.old_id) > 0:
            if values:
                raise ValueError("fixed Source irrep does not accept explicit k parameters")
            source_parameters.append(None)
        else:
            if not values:
                raise ValueError("parametric Source irrep requires explicit k parameters")
            source_parameters.append(projection.source_kparam_for_gid(gid, values))
    direction_blocks = []
    for gid, kparam, direction, domain in zip(
        gids,
        source_parameters,
        opd_labels,
        domain_values,
        strict=True,
    ):
        little = projection.little_record_by_gid(gid)
        if int(little.old_id) > 0:
            row_id, block = _exact_direction_matrix(projection, gid, direction)
            upper = domain_count_from_isotropy_row(
                projection,
                data,
                int(space_group),
                row_id,
            )
            if domain > upper:
                raise ValueError(f"domain {domain} exceeds the Source upper bound {upper}")
            if domain > 1:
                record = domain_operation_record_from_isotropy_row(
                    projection,
                    data,
                    sg=int(space_group),
                    row_id=row_id,
                    domain=domain,
                )
                if record is None:
                    raise ValueError(f"Source domain {domain} has no operation record")
                block = _matrix_product(_exact_record_matrix(data, gid, record), block)
            direction_blocks.append(block)
            continue
        if kparam is None:
            raise AssertionError("parametric factor lost its exact k parameter record")
        block = _exact_dynamic_direction_matrix(projection, gid, direction, kparam)
        dynamic_row = projection._dynamic_row_for_gid_kparam_direction(gid, kparam, direction)
        if dynamic_row is None:
            raise ValueError(f"dynamic OPD row not found for gid={gid}: {direction}")
        basis = tuple(int(value) for value in dynamic_row.basis_values)
        origin = tuple(int(value) for value in dynamic_row.origin_values)
        upper = domain_count_from_subgroup(
            data,
            sg=int(space_group),
            child_sg=int(dynamic_row.subgroup_number),
            basis=basis,
        )
        if domain > upper:
            raise ValueError(f"domain {domain} exceeds the Source upper bound {upper}")
        if domain > 1:
            record = domain_operation_record_from_subgroup(
                data,
                sg=int(space_group),
                child_sg=int(dynamic_row.subgroup_number),
                basis=basis,
                origin=origin,  # type: ignore[arg-type]
                domain=domain,
            )
            if record is None:
                raise ValueError(f"Source domain {domain} has no operation record")
            block = _matrix_product(
                _quadratic3_record_matrix(projection, gid, record, kparam),
                block,
            )
        direction_blocks.append(block)

    full_action = coupled_parametric_irrep_matrices(
        int(space_group),
        labels,
        parameter_rows,
        projection_source=projection,
    )
    return _project_action(full_action, direction_blocks)
