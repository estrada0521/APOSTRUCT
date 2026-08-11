"""Prepare dense coefficient buffers for ``project_vector_``."""

from __future__ import annotations

from collections.abc import Sequence


def dense_source_e8(source: dict[str, object]) -> list[float]:
    dense = [0.0] * 768
    for active in source.get("source_active_values", []):  # type: ignore[union-attr]
        if isinstance(active, list) and len(active) >= 2:
            dense[int(active[0])] = float(active[1])
    return dense


def dense_orderparam_source(
    source_rows: list[list[float]],
    term_count: int,
    row_count: int,
) -> list[float]:
    dense = [0.0] * max(512, row_count * 48 + 16)
    for row_index, row_values in enumerate(source_rows[:row_count]):
        for term_index in range(min(term_count, len(row_values))):
            dense[row_index * 48 + (16 - term_count) + term_index] = float(row_values[term_index])
    return dense


def effective_orderparam_rows(source_rows: list[list[float]]) -> int:
    if not source_rows:
        return 0
    first = tuple(round(float(value), 12) for value in source_rows[0])
    if all(tuple(round(float(value), 12) for value in row_values) == first for row_values in source_rows):
        return 1
    return len(source_rows)


def direct_bridge_coefficients(
    project_basis_by_branch: dict[tuple[int, int], dict[str, object]],
    gid: int,
    pg_irrep: int,
    source: dict[str, object],
    weight_buffer: Sequence[float],
    vector_loop_count: int,
    project_count: int,
    output_length: int,
) -> list[float] | None:
    basis_item = project_basis_by_branch.get((gid, pg_irrep), {})
    if not isinstance(basis_item, dict):
        return None
    active_values = basis_item.get("active_values", [])
    if not isinstance(active_values, list):
        return None
    family_stride = int(basis_item.get("family_stride", 144) or 144)
    source_family = int(source.get("family", 0))
    row_count = int(project_count) + 1
    active = [0.0] * max(family_stride, row_count * 3)
    found = False
    for item in active_values:
        if not (isinstance(item, list) and len(item) >= 2):
            continue
        index = int(item[0])
        if index // family_stride != source_family:
            continue
        local_index = index % family_stride
        if local_index < len(active):
            active[local_index] = float(item[1])
            found = True
    if not found:
        return None
    out = [0.0] * int(output_length)
    for row_index in range(row_count):
        for vector_index in range(max(0, int(vector_loop_count))):
            value = 0.0
            for column_index in range(row_count):
                source_index = column_index * 3 + vector_index
                weight_index = row_index + column_index * 48
                source_value = active[source_index] if source_index < len(active) else 0.0
                weight_value = weight_buffer[weight_index] if weight_index < len(weight_buffer) else 0.0
                value += source_value * weight_value
            output_index = row_index * 3 + vector_index
            if output_index < len(out):
                out[output_index] = value
    return out


def direct_bridge_basis_function(
    source_rows: list[list[float]],
    vector_loop_count: int,
    atom_count: int,
    project_count: int,
    coefficients: list[float],
    output_length: int,
) -> list[float]:
    out = [0.0] * int(output_length)
    term_count = int(project_count) + 1
    for row_index in range(max(0, int(atom_count))):
        row_values = source_rows[row_index] if row_index < len(source_rows) else []
        for vector_index in range(max(0, int(vector_loop_count))):
            value = 0.0
            for term_index in range(term_count):
                source_value = (
                    float(row_values[term_index])
                    if term_index < len(row_values)
                    else 0.0
                )
                coefficient_index = term_index * 3 + vector_index
                coefficient_value = (
                    coefficients[coefficient_index]
                    if coefficient_index < len(coefficients)
                    else 0.0
                )
                value += source_value * coefficient_value
            output_index = row_index * 3 + vector_index
            if output_index < len(out):
                out[output_index] = value
    return out
