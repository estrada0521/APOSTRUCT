"""Mode presentation layout helpers."""

from __future__ import annotations

import re
from typing import Any
import numpy as np
from ISODISTORT.Assembled.Backend.modes.engine.decoder import ModeDataDecoder
from ISODISTORT.Assembled.Backend.modes.engine.subgroup_structure.presentation_transport import (
    site_vector_print_columns,
)


def _mode_site_irrep_labels(
    decoder: ModeDataDecoder,
    trace: dict[str, Any],
    block: dict[str, Any],
    mode_count: int,
    *,
    opd_rows_as_components: bool = False,
) -> list[str] | None:
    summaries = [
        item
        for item in block.get("mode_block_summaries") or []
        if isinstance(item, dict)
    ]
    try:
        site_pg = int((trace.get("wyckoff") or {})["site_pg"])
        bases = [
            decoder.site_pg_irrep_label(site_pg, int(item["pg_irrep"]))
            for item in summaries
        ]
        dimensions = [
            int(
                decoder.image_record(
                    decoder.site_pg_irrep_old_id(site_pg, int(item["pg_irrep"]))
                ).dimension
            )
            for item in summaries
        ]
    except Exception:
        return None
    totals: dict[str, int] = {}
    occurrences: list[int] = []
    display_dimensions: list[int] = []
    for summary, base, dimension in zip(summaries, bases, dimensions, strict=True):
        count = int(summary.get("flat_mode_count") or 0)
        if dimension <= 0 or count <= 0:
            return None
        display_dimension = dimension
        source_row_count = int(summary.get("source_row_count") or 0)
        opd_row_count = sum(
            len(group.get("rows") or [])
            for group in summary.get("opd_groups") or []
            if isinstance(group, dict)
        )
        if (
            opd_rows_as_components
            and (
                int(summary.get("little_type") or 0) == 3
                or len(summary.get("opd_groups") or []) > 1
            )
            and opd_row_count > 0
            and count % opd_row_count == 0
        ):
            display_dimension = opd_row_count
        elif (
            source_row_count > 1
            and count == source_row_count
            and int(summary.get("opd_group_count") or 0) == source_row_count
        ):
            display_dimension = 1
        elif (
            int(summary.get("little_type") or 0) == 3
            and source_row_count > 1
            and int(summary.get("opd_group_count") or 0) != source_row_count
            and count % source_row_count == 0
        ):
            display_dimension = source_row_count
        if count % display_dimension:
            if count == 1:
                # Complex site irreps can be emitted one real print column per
                # independent project family.
                display_dimension = 1
            else:
                return None
        occurrence_count = count // display_dimension
        display_dimensions.append(display_dimension)
        occurrences.append(occurrence_count)
        totals[base] = totals.get(base, 0) + occurrence_count
    seen: dict[str, int] = {}
    labels: list[str] = []
    for summary, base, dimension, occurrence_count in zip(
        summaries, bases, display_dimensions, occurrences, strict=True
    ):
        source_row_count = int(summary.get("source_row_count") or 0)
        flat_mode_count = int(summary.get("flat_mode_count") or 0)
        opd_row_count = sum(
            len(group.get("rows") or [])
            for group in summary.get("opd_groups") or []
            if isinstance(group, dict)
        )
        if (
            not (
                opd_rows_as_components
                and (
                    int(summary.get("little_type") or 0) == 3
                    or len(summary.get("opd_groups") or []) > 1
                )
                and opd_row_count > 0
            )
            and source_row_count > 1
            and flat_mode_count == source_row_count * dimension
            and (
                int(summary.get("little_type") or 0) == 3
                or int(summary.get("opd_group_count") or 0) == source_row_count
            )
            and (
                occurrence_count == 1
                or int(summary.get("opd_group_count") or 0) == source_row_count
            )
        ):
            labels.extend(
                f"{base}({chr(ord('a') + component) if component < 26 else component + 1})"
                for component in range(flat_mode_count)
            )
            continue
        for _occurrence in range(occurrence_count):
            seen[base] = seen.get(base, 0) + 1
            displayed = f"{base}_{seen[base]}" if totals[base] > 1 else base
            labels.extend(
                f"{displayed}({chr(ord('a') + component) if component < 26 else component + 1})"
                for component in range(dimension)
            )
    if len(labels) != mode_count:
        return None
    return labels


def _mode_source_layout(
    block: dict[str, Any], mode_count: int
) -> list[tuple[int, int, int]] | None:
    """Return ``(family, component, family_width)`` in local emitted order."""

    layout = [
        (
            int(summary.get("source_family") or 0),
            component,
            int(summary.get("flat_mode_count") or 0),
        )
        for summary in block.get("mode_block_summaries") or []
        if isinstance(summary, dict)
        for component in _summary_source_print_components(summary)
    ]
    return layout if len(layout) == mode_count else None


def _summary_source_print_components(summary: dict[str, Any]) -> list[int]:
    """Map projected-block summaries to their displayed component identities."""

    mode_count = int(summary.get("flat_mode_count") or 0)
    if int(summary.get("little_type") or 0) != 3:
        return list(range(mode_count))
    return _summary_projected_print_components(summary)


def _summary_projected_print_components(summary: dict[str, Any]) -> list[int]:
    """Transpose group-major emitted ordinals to block-major print ordinals."""

    mode_count = int(summary.get("flat_mode_count") or 0)
    groups = [
        list(group.get("rows") or [])
        for group in summary.get("opd_groups") or []
        if isinstance(group, dict)
    ]
    row_count = sum(len(rows) for rows in groups)
    if mode_count <= 0 or row_count <= 0 or mode_count % row_count:
        return list(range(mode_count))
    block_count = mode_count // row_count
    components: list[int] = []
    row_offset = 0
    for rows in groups:
        for block_index in range(block_count):
            components.extend(
                block_index * row_count + row_offset + row_index
                for row_index in range(len(rows))
            )
        row_offset += len(rows)
    return components if len(components) == mode_count else list(range(mode_count))


def _mode_source_metadata_layout(
    block: dict[str, Any],
    mode_count: int,
    *,
    site_pg: int | None = None,
) -> list[dict[str, Any]] | None:
    """Keep the numeric Source branch identity for every emitted mode."""

    layout = [
        {
            "gid": int(summary.get("gid") or block.get("gid") or 0),
            "site_pg": int(site_pg or (block.get("wyckoff") or {}).get("site_pg") or 0),
            "pg_irrep": int(summary.get("pg_irrep") or 0),
            "family": int(summary.get("source_family") or 0),
            "component": component,
            "print_component": print_component,
            "family_width": int(summary.get("flat_mode_count") or 0),
            "source_row_count": int(summary.get("source_row_count") or 0),
            "little_type": int(summary.get("little_type") or 0),
            "opd_groups": list(summary.get("opd_groups") or []),
        }
        for summary in block.get("mode_block_summaries") or []
        if isinstance(summary, dict)
        for component, print_component in enumerate(
            _summary_projected_print_components(summary)
        )
    ]
    return layout if len(layout) == mode_count else None


def _source_family_vector_matrix(
    trace: dict[str, Any],
    *,
    gid: int,
    pg_irrep: int,
    family: int,
    full_dim: int,
) -> np.ndarray | None:
    """Decode one faithful project family as carrier-by-site columns."""

    item = next(
        (
            candidate
            for candidate in trace.get("project_return_basis_type1") or []
            if isinstance(candidate, dict)
            and int(candidate.get("gid") or 0) == int(gid)
            and int(candidate.get("pg_irrep") or 0) == int(pg_irrep)
        ),
        None,
    )
    source = next(
        (
            candidate
            for candidate in trace.get("project_vector_bridge_source") or []
            if isinstance(candidate, dict)
            and int(candidate.get("gid") or 0) == int(gid)
            and int(candidate.get("pg_irrep") or 0) == int(pg_irrep)
            and int(candidate.get("family") or 0) == int(family)
        ),
        None,
    )
    if item is None or source is None:
        return None
    site_columns = int(source.get("vector_loop_count") or 0)
    if full_dim <= 0 or site_columns <= 0:
        return None
    stride = int(item.get("family_stride") or 144)
    matrix = np.zeros((int(full_dim), site_columns), dtype=complex)
    found = False
    for active in item.get("active_values") or []:
        if not (isinstance(active, list) and len(active) >= 2):
            continue
        index = int(active[0])
        if index // stride != int(family):
            continue
        local = index % stride
        carrier = local // 3
        site_column = local % 3
        if carrier < full_dim and site_column < site_columns:
            matrix[carrier, site_column] = complex(float(active[1]))
            found = found or abs(float(active[1])) > 1e-12
    return matrix if found else None


def _parallel_source_child_vector(
    decoder: ModeDataDecoder,
    *,
    child_sg: int,
    child_wyckoff_site: str,
    reference: np.ndarray,
    vector_setting: int,
) -> np.ndarray | None:
    """Find the unique Source child-site print column spanning reference."""

    match = re.fullmatch(r"\s*\d+\s*([A-Za-z]+)\s*", child_wyckoff_site)
    if match is None:
        return None
    try:
        row = next(
            item
            for item in decoder.wyckoff_rows(int(child_sg))
            if item.label == match.group(1)
        )
    except StopIteration:
        return None
    records = decoder.wyc_pg_elements_records(int(child_sg), row)
    count = int(decoder.iso.wyckoff["iwyckoff_pg_irrep_count"][int(row.site_pg) - 1])
    reference_norm = float(np.linalg.norm(reference))
    if reference_norm <= 1e-12:
        return None
    candidates: list[np.ndarray] = []
    for child_pg_irrep in range(1, count + 1):
        for column in site_vector_print_columns(
            decoder,
            site_pg=int(row.site_pg),
            pg_irrep=int(child_pg_irrep),
            site_operation_records=records,
            vector_setting=int(vector_setting),
        ):
            vector = np.asarray(column.vector, dtype=float)
            norm = float(np.linalg.norm(vector))
            if norm <= 1e-12:
                continue
            if (
                float(np.linalg.norm(np.cross(reference, vector)))
                <= 1e-10 * reference_norm * norm
            ):
                if not any(
                    np.allclose(vector, seen, atol=1e-12) for seen in candidates
                ):
                    candidates.append(vector)
    return candidates[0] if len(candidates) == 1 else None
