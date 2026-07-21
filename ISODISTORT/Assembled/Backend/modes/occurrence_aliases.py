"""Fail-closed admission helpers for replayed dynamic occurrences."""

from __future__ import annotations

from collections.abc import Hashable, Mapping, Sequence
from dataclasses import dataclass
import math


@dataclass(frozen=True)
class OccurrenceSiteEmission:
    """Source emission identity observed for one spec at one occupied site."""

    site_index: int
    row_id: int
    source_signature: Hashable
    mode_count: int
    all_modes_nonzero: bool
    atom_operation_records: Hashable = ()
    atom_fractionals: Hashable = ()


def admitted_occurrence_alias_spec_orders(
    candidate_to_anchor: Mapping[int, int],
    *,
    eligible_sites_by_spec: Mapping[int, Sequence[int]],
    emissions: Mapping[tuple[int, int], OccurrenceSiteEmission],
) -> frozenset[int]:
    """Admit aliases whose Source emission matches their representative row.

    Candidate generation proves that a stored reciprocal occurrence has one
    invariant alias. This downstream gate separately proves that the alias is
    printable at every selected occupied site. Missing blocks, changed Source
    topology, zero logical modes, or an empty total emission all fail closed.
    """

    admitted: set[int] = set()
    for candidate_order, anchor_order in candidate_to_anchor.items():
        candidate_order = int(candidate_order)
        anchor_order = int(anchor_order)
        if candidate_order == anchor_order:
            continue
        candidate_sites = tuple(
            int(value) for value in eligible_sites_by_spec.get(candidate_order, ())
        )
        anchor_sites = tuple(
            int(value) for value in eligible_sites_by_spec.get(anchor_order, ())
        )
        if not candidate_sites or candidate_sites != anchor_sites:
            continue
        emitted = 0
        valid = True
        for site_index in candidate_sites:
            candidate = emissions.get((candidate_order, site_index))
            anchor = emissions.get((anchor_order, site_index))
            if (
                candidate is None
                or anchor is None
                or int(candidate.site_index) != site_index
                or int(anchor.site_index) != site_index
                or int(candidate.row_id) != int(anchor.row_id)
                or candidate.source_signature != anchor.source_signature
                or candidate.atom_operation_records != anchor.atom_operation_records
                or candidate.atom_fractionals != anchor.atom_fractionals
                or int(candidate.mode_count) < 0
                or int(candidate.mode_count) != int(anchor.mode_count)
                or not bool(candidate.all_modes_nonzero)
            ):
                valid = False
                break
            emitted += int(candidate.mode_count)
        if valid and emitted > 0:
            admitted.add(candidate_order)
    return frozenset(admitted)


MagneticOccurrenceSiteEmission = OccurrenceSiteEmission
admitted_magnetic_occurrence_alias_spec_orders = admitted_occurrence_alias_spec_orders


def _exact_source_direction_matrix(
    value: object,
) -> tuple[tuple[float, ...], ...] | None:
    if (
        isinstance(value, (str, bytes, bytearray))
        or not isinstance(value, Sequence)
        or not value
    ):
        return None
    rows: list[tuple[float, ...]] = []
    width: int | None = None
    for raw_row in value:
        if (
            isinstance(raw_row, (str, bytes, bytearray))
            or not isinstance(raw_row, Sequence)
            or not raw_row
        ):
            return None
        row: list[float] = []
        for component in raw_row:
            if type(component) is not float or not math.isfinite(component):
                return None
            row.append(component)
        if width is None:
            width = len(row)
        elif len(row) != width:
            return None
        rows.append(tuple(row))
    return tuple(rows)


def direction_distinct_occurrence_alias_spec_orders(
    render_specs: Sequence[tuple[Mapping[str, object], object, object, object]],
    candidate_to_anchor: Mapping[int, int],
) -> frozenset[int]:
    """Select exact Source aliases whose direction is not the representative's."""

    selected: set[int] = set()
    for raw_candidate, raw_anchor in candidate_to_anchor.items():
        if (
            type(raw_candidate) is not int
            or type(raw_anchor) is not int
            or raw_candidate < 0
            or raw_anchor < 0
            or raw_candidate >= len(render_specs)
            or raw_anchor >= len(render_specs)
        ):
            continue
        candidate = _exact_source_direction_matrix(
            render_specs[raw_candidate][0].get("direction_matrix")
        )
        anchor = _exact_source_direction_matrix(
            render_specs[raw_anchor][0].get("direction_matrix")
        )
        if (
            candidate is None
            or anchor is None
            or len(candidate) != len(anchor)
            or any(len(left) != len(right) for left, right in zip(candidate, anchor))
            or any(
                type(left) is not type(right)
                for candidate_row, anchor_row in zip(candidate, anchor)
                for left, right in zip(candidate_row, anchor_row)
            )
        ):
            continue
        if candidate != anchor:
            selected.add(raw_candidate)
    return frozenset(selected)


def admitted_self_complete_occurrence_alias_spec_orders(
    candidate_orders: Sequence[int] | frozenset[int],
    *,
    eligible_sites_by_spec: Mapping[int, Sequence[int]],
    emissions: Mapping[tuple[int, int], OccurrenceSiteEmission],
) -> frozenset[int]:
    """Admit independent aliases only after complete candidate-side emission."""

    admitted: set[int] = set()
    for raw_order in candidate_orders:
        if type(raw_order) is not int:
            continue
        sites = tuple(eligible_sites_by_spec.get(raw_order, ()))
        if (
            not sites
            or any(type(site) is not int or site < 0 for site in sites)
            or len(set(sites)) != len(sites)
        ):
            continue
        complete = True
        for site in sites:
            emission = emissions.get((raw_order, site))
            if (
                not isinstance(emission, OccurrenceSiteEmission)
                or type(emission.site_index) is not int
                or emission.site_index != site
                or type(emission.row_id) is not int
                or emission.row_id <= 0
                or type(emission.mode_count) is not int
                or emission.mode_count <= 0
                or emission.all_modes_nonzero is not True
                or not isinstance(emission.source_signature, tuple)
                or not emission.source_signature
            ):
                complete = False
                break
        if complete:
            admitted.add(raw_order)
    return frozenset(admitted)
