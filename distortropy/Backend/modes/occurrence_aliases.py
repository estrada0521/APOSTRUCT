"""Admission helpers for replayed dynamic occurrences."""

from __future__ import annotations

from collections.abc import Hashable, Mapping, Sequence
from dataclasses import dataclass


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
    """Admit aliases whose Source emission matches their representative row."""

    admitted: set[int] = set()
    for candidate_order, anchor_order in candidate_to_anchor.items():
        if candidate_order == anchor_order:
            continue
        candidate_sites = tuple(eligible_sites_by_spec.get(candidate_order, ()))
        anchor_sites = tuple(eligible_sites_by_spec.get(anchor_order, ()))
        if not candidate_sites or candidate_sites != anchor_sites:
            continue
        emitted = 0
        for site_index in candidate_sites:
            candidate = emissions.get((candidate_order, site_index))
            anchor = emissions.get((anchor_order, site_index))
            if (
                candidate is None
                or anchor is None
                or candidate.row_id != anchor.row_id
                or candidate.source_signature != anchor.source_signature
                or candidate.atom_operation_records != anchor.atom_operation_records
                or candidate.atom_fractionals != anchor.atom_fractionals
                or candidate.mode_count != anchor.mode_count
                or not candidate.all_modes_nonzero
            ):
                break
            emitted += candidate.mode_count
        else:
            if emitted > 0:
                admitted.add(candidate_order)
    return frozenset(admitted)


def direction_distinct_occurrence_alias_spec_orders(
    render_specs: Sequence[tuple[Mapping[str, object], object, object, object]],
    candidate_to_anchor: Mapping[int, int],
) -> frozenset[int]:
    """Select exact Source aliases whose direction is not the representative's."""

    selected: set[int] = set()
    for candidate_order, anchor_order in candidate_to_anchor.items():
        candidate = render_specs[candidate_order][0].get("direction_matrix")
        anchor = render_specs[anchor_order][0].get("direction_matrix")
        if candidate != anchor:
            selected.add(candidate_order)
    return frozenset(selected)


def admitted_self_complete_occurrence_alias_spec_orders(
    candidate_orders: Sequence[int] | frozenset[int],
    *,
    eligible_sites_by_spec: Mapping[int, Sequence[int]],
    emissions: Mapping[tuple[int, int], OccurrenceSiteEmission],
) -> frozenset[int]:
    """Admit independent aliases only after complete candidate-side emission."""

    admitted: set[int] = set()
    for candidate_order in candidate_orders:
        sites = tuple(eligible_sites_by_spec.get(candidate_order, ()))
        if not sites:
            continue
        candidate_emissions = [
            emissions.get((candidate_order, site)) for site in sites
        ]
        if all(
            emission is not None
            and emission.mode_count > 0
            and emission.all_modes_nonzero
            and emission.source_signature
            for emission in candidate_emissions
        ):
            admitted.add(candidate_order)
    return frozenset(admitted)
