"""Mode-subduction specification helpers."""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import replace
from fractions import Fraction
from numbers import Integral
from typing import Any
import numpy as np
from APOSTRUCT.Backend.source.magnetic import data as magnetic_data
from APOSTRUCT.Backend.modes.engine.decoder import ModeDataDecoder
from APOSTRUCT.Backend.modes.engine.dynamic_subduction import (
    DynamicSubductionRow,
    dynamic_magnetic_subduction_occurrence_alias_rows,
    dynamic_magnetic_subduction_rows,
    dynamic_subduction_occurrence_alias_rows,
    dynamic_subduction_rows,
    kvec_standard_provenance,
    stored_occurrence_alias_candidates,
    _strict_integral_values,
)
from APOSTRUCT.Backend.modes.engine.project.mode_counts import (
    little_records_for_k,
)
from APOSTRUCT.Backend.isotropy.engine.source_data import SourceData as OpdSourceData

from APOSTRUCT.Backend.modes.common import (
    _integer_basis_tuple,
    _isotropy_from_opd_row,
    _isotropy_row_id_from_opd_row,
    _k_label_from_irrep_label,
    _k_params,
    _same_source_kparam,
)
from APOSTRUCT.Backend.modes.request_context import (
    _pml_vector_to_case_k_params,
    _selected_dynamic_gid,
    _source_kparam_record,
)
from APOSTRUCT.Backend.modes.structure_runtime import (
    _selected_subgroup_number,
    _source_split_basis_from_opd_row,
    _source_split_origin_from_opd_row,
    _subgroup_parent_operation_records,
)


# Matches isotropy.engine.get_isotropy._fmt_matrix_value's fixed-six-decimal cell.
_DYNAMIC_SOURCE_MATRIX_ROUNDING_CELL = float(np.nextafter(5e-7, np.inf))


def _opd_source_data(decoder: ModeDataDecoder) -> OpdSourceData:
    """Share immutable Source tables without retaining per-view memo state.

    Both decoders read the same ``data_space``, ``data_little``, and ``const``
    tables from the exact Source directory.  A fresh view keeps case-local
    calculation memoization from leaking across requests.
    """

    return OpdSourceData(decoder.data_dir, tables=decoder.iso)


def _source_kparam_identity(
    values: tuple[int, int, int, int] | list[int],
) -> tuple[Fraction, Fraction, Fraction]:
    return tuple(
        Fraction(int(values[index]), int(values[3]))
        for index in range(3)
    )  # type: ignore[return-value]


def _strict_source_kparam(values: object) -> tuple[int, int, int, int] | None:
    if (
        isinstance(values, (str, bytes, bytearray))
        or not isinstance(values, Sequence)
        or len(values) != 4
        or any(
            isinstance(value, bool) or not isinstance(value, Integral)
            for value in values
        )
    ):
        return None
    normalized = tuple(int(value) for value in values)
    if normalized[3] == 0:
        return None
    return normalized[0], normalized[1], normalized[2], normalized[3]


def _strict_positive_gid(value: object) -> int | None:
    if isinstance(value, bool) or not isinstance(value, Integral) or value <= 0:
        return None
    return int(value)


def _strict_occurrence_alias_anchor(
    value: object,
) -> tuple[int, tuple[int, int, int, int]] | None:
    if not isinstance(value, tuple) or len(value) != 2:
        return None
    gid = _strict_positive_gid(value[0])
    source_kparam = _strict_source_kparam(value[1])
    if gid is None or source_kparam is None:
        return None
    return gid, source_kparam


def _strict_occurrence_alias_heading(
    value: object,
) -> tuple[int, str, tuple[int, int, int, int]] | None:
    if not isinstance(value, tuple) or len(value) != 3:
        return None
    gid = _strict_positive_gid(value[0])
    k_label = value[1]
    source_kparam = _strict_source_kparam(value[2])
    if gid is None or not isinstance(k_label, str) or source_kparam is None:
        return None
    return gid, k_label, source_kparam


def _strict_reciprocal_vector_pml(
    value: object,
) -> tuple[Fraction, Fraction, Fraction] | None:
    if (
        isinstance(value, (str, bytes, bytearray))
        or not isinstance(value, Sequence)
        or len(value) != 3
        or any(not isinstance(component, Fraction) for component in value)
    ):
        return None
    return value[0], value[1], value[2]


def _occurrence_alias_observation_markers_valid(spec: dict[str, Any]) -> bool:
    name = "_occurrence_alias_observation_only"
    return name not in spec or spec.get(name) is True


def _strict_source_kparam_identity(
    values: object,
) -> tuple[Fraction, Fraction, Fraction] | None:
    normalized = _strict_source_kparam(values)
    if normalized is None:
        return None
    return tuple(
        Fraction(normalized[index], normalized[3])
        for index in range(3)
    )  # type: ignore[return-value]


def _selected_magnetic_case_identity(
    slot: object,
) -> tuple[int, tuple[int, int, int, int]] | None:
    """Parse a selected magnetic Source identity without lossy coercion."""

    if not isinstance(slot, dict):
        return None
    irrep = slot.get("irrep")
    if not isinstance(irrep, dict) or not bool(irrep.get("magnetic")):
        return None
    gid = _strict_positive_gid(irrep.get("gid"))
    raw_source_kparam = slot.get("source_kparam")
    if raw_source_kparam is None:
        raw_source_kparam = (0, 0, 0, 1)
    source_kparam = _strict_source_kparam(raw_source_kparam)
    if gid is None or source_kparam is None:
        return None
    return gid, source_kparam


def _strict_case_k_params(values: object) -> tuple[Fraction, ...] | None:
    if (
        isinstance(values, (str, bytes, bytearray))
        or not isinstance(values, Sequence)
        or any(not isinstance(value, Fraction) for value in values)
    ):
        return None
    return tuple(values)


def _row_projects_selected_case_k_params(
    decoder: ModeDataDecoder,
    row: object,
    *,
    selected_gid: int,
    selected_case_k_params: tuple[Fraction, ...],
) -> bool:
    """Keep a canonical Source row already at the selected PML point."""

    if (
        not isinstance(row, DynamicSubductionRow)
        or _strict_positive_gid(getattr(row, "gid", None)) != selected_gid
        or _strict_source_kparam(getattr(row, "source_kparam", None)) is None
    ):
        return False
    try:
        displayed = _pml_vector_to_case_k_params(
            decoder,
            gid=selected_gid,
            reciprocal_vector_pml=row.reciprocal_vector_pml,
        )
    except (
        AttributeError,
        IndexError,
        KeyError,
        TypeError,
        ValueError,
        ZeroDivisionError,
    ):
        return False
    return displayed is not None and tuple(displayed) == selected_case_k_params


def _selected_source_occurrence_matches(
    *,
    gid: int,
    source_kparam: tuple[int, int, int, int] | list[int],
    selected_gid: int | None,
    selected_source_kparam: tuple[int, int, int, int] | list[int] | None,
) -> bool:
    """Match the selected raw Source occurrence, not only its irrep family."""

    return bool(
        selected_gid is not None
        and selected_source_kparam is not None
        and int(gid) == int(selected_gid)
        and _source_kparam_identity(source_kparam)
        == _source_kparam_identity(selected_source_kparam)
    )


def _occurrence_heading_identity(
    row: Any,
) -> tuple[int, str, tuple[int, int, int, int]] | None:
    try:
        gid_values = _strict_integral_values((row.gid,), size=1)
        if gid_values is None or gid_values[0] <= 0:
            return None
        values = _strict_integral_values(row.source_kparam)
        if values is None or values[3] == 0:
            return None
        return gid_values[0], str(row.k_label), values
    except (AttributeError, TypeError, ValueError):
        return None


def _stored_occurrence_alias_specs(
    decoder: ModeDataDecoder,
    raw_rows: Any,
    represented_rows: Any,
    *,
    aliases_for_class: Any,
    spec_for_row: Any,
) -> list[dict[str, Any]]:
    """Build stored-occurrence observations and candidates."""

    raw = tuple(raw_rows)
    represented = tuple(represented_rows)
    represented_by_family: dict[
        tuple[int, str, str, tuple[Any, ...]],
        set[
            tuple[
                tuple[int, int, int, int],
                tuple[Fraction, Fraction, Fraction],
            ]
        ],
    ] = {}
    represented_by_header: dict[
        tuple[int, str, str],
        list[
            tuple[
                tuple[Any, ...],
                tuple[int, int, int, int],
                tuple[Fraction, Fraction, Fraction],
            ]
        ],
    ] = {}
    for row in represented:
        source = _strict_source_kparam(getattr(row, "source_kparam", None))
        pml = _strict_reciprocal_vector_pml(
            getattr(row, "reciprocal_vector_pml", None)
        )
        if source is None or pml is None:
            continue
        header = (
            int(row.gid),
            str(row.irrep_label),
            str(row.k_label),
        )
        occurrences = tuple(row.source_occurrences)
        family = (
            *header,
            occurrences,
        )
        represented_by_family.setdefault(family, set()).add((source, pml))
        represented_by_header.setdefault(header, []).append(
            (occurrences, source, pml)
        )

    def represented_occurrence(row: Any) -> tuple[Any, Any, bool] | None:
        header = (
            int(row.gid),
            str(row.irrep_label),
            str(row.k_label),
        )
        raw_occurrences = tuple(row.source_occurrences)
        family = (*header, raw_occurrences)
        represented_occurrences = represented_by_family.get(family, set())
        if not represented_occurrences and raw_occurrences:
            raw_occurrence_set = set(raw_occurrences)
            represented_occurrences = {
                (source, pml)
                for occurrences, source, pml in represented_by_header.get(header, ())
                if occurrences
                and set(occurrences).issubset(raw_occurrence_set)
            }
        exact = (
            _strict_source_kparam(getattr(row, "source_kparam", None)),
            _strict_reciprocal_vector_pml(
                getattr(row, "reciprocal_vector_pml", None)
            ),
        )
        if None not in exact and exact in represented_occurrences:
            return exact[0], exact[1], False
        if len(represented_occurrences) == 1:
            source, pml = next(iter(represented_occurrences))
            return source, pml, True
        if not represented_occurrences:
            raise ValueError("Source occurrence family has no represented row")
        raise ValueError("Source occurrence family has multiple represented rows")

    candidates = stored_occurrence_alias_candidates(
        decoder,
        raw,
        aliases_for_class=aliases_for_class,
        heading_identity_for_row=_occurrence_heading_identity,
        represented_occurrence_for_row=represented_occurrence,
        represented_heading_identities=(
            *(
                identity
                for row in represented
                if (identity := _occurrence_heading_identity(row)) is not None
            ),
        ),
    )
    represented_exact = {
        (gid_values[0], values)
        for row in represented
        if (gid_values := _strict_integral_values((row.gid,), size=1)) is not None
        and gid_values[0] > 0
        and (values := _strict_integral_values(row.source_kparam)) is not None
    }
    required_anchors = {
        (
            int(candidate.representative_gid),
            _strict_integral_values(candidate.representative_source_kparam),
        )
        for candidate in candidates
    }
    specs: list[dict[str, Any]] = []
    for row in raw:
        values = _strict_integral_values(row.source_kparam)
        if values is None:
            continue
        gid_values = _strict_integral_values((row.gid,), size=1)
        if gid_values is None or gid_values[0] <= 0:
            continue
        key = (gid_values[0], values)
        if key not in required_anchors - represented_exact:
            continue
        spec = spec_for_row(row)
        spec["primary"] = False
        spec["_occurrence_alias_observation_only"] = True
        specs.append(spec)
    for candidate in candidates:
        spec = spec_for_row(candidate.candidate_row)
        spec["primary"] = False
        spec["_occurrence_alias_anchor"] = (
            int(candidate.representative_gid),
            tuple(int(value) for value in candidate.representative_source_kparam),
        )
        spec["_occurrence_alias_candidate"] = _strict_integral_values(
            candidate.candidate_source_kparam
        )
        spec["_occurrence_alias_heading_identity"] = candidate.heading_identity
        specs.append(spec)
    return specs


def _merge_occurrence_alias_specs(
    represented_specs: Sequence[dict[str, Any]],
    alias_specs: Sequence[dict[str, Any]],
) -> list[dict[str, Any]]:
    """Place exact aliases beside their anchor in presentation-PML order."""

    anchors: dict[
        tuple[int, tuple[int, int, int, int]], list[int]
    ] = {}
    for index, spec in enumerate(represented_specs):
        gid = _strict_positive_gid(spec.get("gid"))
        source = _strict_source_kparam(spec.get("source_kparam"))
        if gid is not None and source is not None:
            anchors.setdefault((gid, source), []).append(index)
    grouped: dict[int, list[dict[str, Any]]] = {}
    unplaced: list[dict[str, Any]] = []
    for spec in alias_specs:
        anchor = _strict_occurrence_alias_anchor(
            spec.get("_occurrence_alias_anchor")
        )
        matches = () if anchor is None else anchors.get(anchor, ())
        if len(matches) != 1:
            unplaced.append(spec)
            continue
        grouped.setdefault(matches[0], []).append(spec)

    merged: list[dict[str, Any]] = []

    def pml_key(item: dict[str, Any]) -> tuple[bool, tuple[Fraction, ...]]:
        pml = _strict_reciprocal_vector_pml(item.get("reciprocal_vector_pml"))
        return pml is None, pml or ()

    for index, spec in enumerate(represented_specs):
        cluster = [spec, *grouped.get(index, ())]
        cluster.sort(key=pml_key)
        merged.extend(cluster)
    merged.extend(unplaced)
    return merged


def _canonicalize_exact_pml_alias_specs(
    specs: Sequence[dict[str, Any]],
) -> list[dict[str, Any]]:
    """Collapse one stale heading onto its exact stored Source occurrence."""

    grouped: dict[tuple[int, tuple[Fraction, ...]], list[int]] = {}
    for index, spec in enumerate(specs):
        pml = _strict_reciprocal_vector_pml(spec.get("reciprocal_vector_pml"))
        gid = _strict_positive_gid(spec.get("gid"))
        if gid is None or pml is None:
            continue
        grouped.setdefault((gid, pml), []).append(index)

    replacements: dict[int, tuple[int, dict[str, Any]]] = {}
    for (group_gid, _group_pml), indices in grouped.items():
        candidates = [
            index
            for index in indices
            if "_occurrence_alias_candidate" in specs[index]
        ]
        ordinary = [
            index
            for index in indices
            if index not in candidates
            and "_occurrence_alias_observation_only" not in specs[index]
        ]
        if len(candidates) != 1 or len(ordinary) != 1:
            continue
        candidate_index, ordinary_index = candidates[0], ordinary[0]
        candidate, row = specs[candidate_index], specs[ordinary_index]
        source_kparam = _strict_source_kparam(
            candidate.get("_occurrence_alias_candidate")
        )
        row_source_kparam = _strict_source_kparam(row.get("source_kparam"))
        occurrences = tuple(row.get("source_occurrences") or ())
        candidate_pml = _strict_reciprocal_vector_pml(
            candidate.get("reciprocal_vector_pml")
        )
        k_label = candidate.get("k_label")
        label = candidate.get("label")
        anchor_identity = _strict_occurrence_alias_anchor(
            candidate.get("_occurrence_alias_anchor")
        )
        expected_heading = (group_gid, k_label, source_kparam)
        generic_markers_match = (
            anchor_identity is not None
            and anchor_identity[0] == group_gid
            and _strict_occurrence_alias_heading(
                candidate.get("_occurrence_alias_heading_identity")
            )
            == expected_heading
        )
        if (
            source_kparam is None
            or row_source_kparam is None
            or candidate_pml is None
            or not isinstance(label, str)
            or not isinstance(k_label, str)
            or not generic_markers_match
            or "_occurrence_alias_observation_only" in candidate
            or candidate.get("primary") is not False
            or row.get("primary") is not False
            or label != row.get("label")
            or k_label != row.get("k_label")
            or candidate.get("direction_matrix") != row.get("direction_matrix")
            or _strict_source_kparam(candidate.get("source_kparam")) != source_kparam
            or row_source_kparam == source_kparam
            or not occurrences
            or any(
                _strict_source_kparam(getattr(item, "source_kparam", None))
                != source_kparam
                or _strict_reciprocal_vector_pml(
                    getattr(item, "reciprocal_vector_pml", None)
                )
                != candidate_pml
                for item in occurrences
            )
        ):
            continue
        anchor_rows = [
            spec
            for spec in specs
            if "_occurrence_alias_candidate" not in spec
            and (
                _strict_positive_gid(spec.get("gid")),
                _strict_source_kparam(spec.get("source_kparam")),
            )
            == anchor_identity
            and spec.get("primary") is False
            and spec.get("label") == label
            and spec.get("k_label") == k_label
            and spec.get("direction_matrix") == candidate.get("direction_matrix")
            and _occurrence_alias_observation_markers_valid(spec)
        ]
        if len(anchor_rows) != 1:
            continue
        anchor_occurrences = tuple(anchor_rows[0].get("source_occurrences") or ())
        if not any(
            _strict_source_kparam(getattr(item, "source_kparam", None))
            == source_kparam
            and _strict_reciprocal_vector_pml(
                getattr(item, "reciprocal_vector_pml", None)
            )
            == candidate_pml
            for item in anchor_occurrences
        ):
            continue
        canonical = {
            key: value
            for key, value in candidate.items()
            if not key.startswith("_occurrence_alias_")
        }
        canonical["source_occurrences"] = row["source_occurrences"]
        replacements[ordinary_index] = (candidate_index, canonical)

    canonical_specs: list[dict[str, Any]] = []
    removed = {candidate for candidate, _spec in replacements.values()}
    active_anchors = {
        anchor
        for index, spec in enumerate(specs)
        if index not in removed
        if (
            anchor := _strict_occurrence_alias_anchor(
                spec.get("_occurrence_alias_anchor")
            )
        )
        is not None
    }
    for index, spec in enumerate(specs):
        if index in replacements:
            canonical_specs.append(replacements[index][1])
        elif index not in removed:
            if bool(spec.get("_occurrence_alias_observation_only")):
                key = (
                    _strict_positive_gid(spec.get("gid")),
                    _strict_source_kparam(spec.get("source_kparam")),
                )
                if None not in key and key not in active_anchors:
                    continue
            canonical_specs.append(spec)
    return canonical_specs


def _primary_dynamic_case_is_listed(
    rows: Any,
    *,
    primary_gid: int | None,
    primary_source_kparam: tuple[int, int, int, int] | None,
) -> bool:
    """Test exact Source occurrence identity, not reciprocal-lattice equivalence."""

    if primary_gid is None or primary_source_kparam is None:
        return False
    target = _source_kparam_identity(primary_source_kparam)
    return any(
        int(row.gid) == int(primary_gid)
        and _source_kparam_identity(row.source_kparam) == target
        for row in rows
    )


def _promote_selected_dynamic_occurrences(
    decoder: ModeDataDecoder,
    *,
    sg: int,
    basis: tuple[int, ...],
    operations: tuple[tuple[int, int, int, int, int], ...],
    rows: Any,
    selected_cases: list[tuple[int, tuple[int, int, int, int]]],
    irrep_source: OpdSourceData,
    preserve_discovered_carrier: bool,
) -> list[Any]:
    """Replace deduped formal rows with their selected primary occurrence."""

    promoted_rows = list(rows)
    claimed_indices: set[int] = set()
    promotions: list[tuple[int, Any, Any]] = []
    for selected_gid, selected_source_kparam in selected_cases:
        if _primary_dynamic_case_is_listed(
            promoted_rows,
            primary_gid=selected_gid,
            primary_source_kparam=selected_source_kparam,
        ):
            continue
        for row_index, row in enumerate(tuple(promoted_rows)):
            if row_index in claimed_indices or int(row.gid) != int(selected_gid):
                continue
            selected_alias = next(
                (
                    alias
                    for alias in dynamic_subduction_occurrence_alias_rows(
                        decoder,
                        sg=int(sg),
                        basis=basis,
                        operations=operations,
                        row=row,
                        irrep_source=irrep_source,
                    )
                    if _source_kparam_identity(alias.source_kparam)
                    == _source_kparam_identity(selected_source_kparam)
                ),
                None,
            )
            if selected_alias is None:
                try:
                    standard_kparam = kvec_standard_provenance(
                        decoder,
                        int(sg),
                        int(selected_gid),
                        tuple(int(value) for value in selected_source_kparam),
                    ).standard_kparam
                except (KeyError, TypeError, ValueError):
                    standard_kparam = None
                if (
                    standard_kparam is not None
                    and _source_kparam_identity(row.source_kparam)
                    == _source_kparam_identity(standard_kparam)
                ):
                    selected_alias = replace(
                        row,
                        source_kparam=tuple(int(value) for value in selected_source_kparam),
                    )
            if selected_alias is None:
                continue
            promotions.append((row_index, row, selected_alias))
            promoted_rows[row_index] = selected_alias
            claimed_indices.add(row_index)
            break

    shared_occurrence_families: dict[tuple[Any, ...], set[int]] = {}
    if preserve_discovered_carrier:
        selected_identities = {
            (int(gid), _source_kparam_identity(source_kparam))
            for gid, source_kparam in selected_cases
        }
        for row_index, row in enumerate(rows):
            identity = (int(row.gid), _source_kparam_identity(row.source_kparam))
            family = tuple(row.source_occurrences)
            if identity in selected_identities and family:
                shared_occurrence_families.setdefault(family, set()).add(row_index)
        for row_index, row, _ in promotions:
            family = tuple(row.source_occurrences)
            if family:
                shared_occurrence_families.setdefault(family, set()).add(row_index)
    for row_index, row, selected_alias in promotions:
        family = tuple(row.source_occurrences)
        if len(shared_occurrence_families.get(family, ())) < 2:
            continue
        promoted_rows[row_index] = replace(
            selected_alias,
            direction_matrix=row.direction_matrix,
            source_occurrences=row.source_occurrences,
            carrier_source_kparam=row.source_kparam,
        )
    return promoted_rows


def _promote_selected_magnetic_occurrences(
    decoder: ModeDataDecoder,
    magnetic_source: OpdSourceData,
    *,
    sg: int,
    basis: tuple[int, ...],
    operations: tuple[tuple[int, int, int, int, int], ...],
    rows: Any,
    selected_cases: list[
        tuple[int, tuple[int, int, int, int], tuple[Fraction, ...] | None]
    ],
) -> list[Any]:
    """Replace a deduped row's heading; ``None`` omits the coupled PML guard."""

    promoted_rows = list(rows)
    claimed_indices: set[int] = set()
    for (
        raw_selected_gid,
        raw_selected_source_kparam,
        raw_selected_case_k_params,
    ) in selected_cases:
        selected_gid = _strict_positive_gid(raw_selected_gid)
        selected_source_kparam = _strict_source_kparam(raw_selected_source_kparam)
        selected_case_k_params = (
            None
            if raw_selected_case_k_params is None
            else _strict_case_k_params(raw_selected_case_k_params)
        )
        if (
            selected_gid is None
            or selected_source_kparam is None
            or (
                raw_selected_case_k_params is not None
                and selected_case_k_params is None
            )
        ):
            continue
        target = _strict_source_kparam_identity(selected_source_kparam)
        if any(
            _strict_positive_gid(getattr(row, "gid", None)) == selected_gid
            and _strict_source_kparam_identity(getattr(row, "source_kparam", None))
            == target
            for row in promoted_rows
        ):
            continue
        if selected_case_k_params is not None and any(
            _row_projects_selected_case_k_params(
                decoder,
                row,
                selected_gid=selected_gid,
                selected_case_k_params=selected_case_k_params,
            )
            for row in promoted_rows
        ):
            continue
        for row_index, row in enumerate(tuple(promoted_rows)):
            if (
                row_index in claimed_indices
                or not isinstance(row, DynamicSubductionRow)
                or _strict_positive_gid(getattr(row, "gid", None)) != selected_gid
                or _strict_source_kparam(getattr(row, "source_kparam", None)) is None
            ):
                continue
            try:
                aliases = tuple(
                    dynamic_magnetic_subduction_occurrence_alias_rows(
                        decoder,
                        magnetic_source,
                        sg=int(sg),
                        basis=basis,
                        operations=operations,
                        row=row,
                        requested_source_kparam=selected_source_kparam,
                    )
                )
            except (
                AttributeError,
                IndexError,
                KeyError,
                TypeError,
                ValueError,
                ZeroDivisionError,
            ):
                aliases = ()
            selected_alias = next(
                (
                    alias
                    for alias in aliases
                    if isinstance(alias, DynamicSubductionRow)
                    and _strict_positive_gid(getattr(alias, "gid", None))
                    == selected_gid
                    and _strict_source_kparam_identity(
                        getattr(alias, "source_kparam", None)
                    )
                    == target
                ),
                None,
            )
            if selected_alias is None:
                continue
            promoted_rows[row_index] = replace(
                selected_alias,
                direction_matrix=row.direction_matrix,
                source_occurrences=row.source_occurrences,
                carrier_source_kparam=row.source_kparam,
            )
            claimed_indices.add(row_index)
            break
    return promoted_rows


def _replay_unrepresented_type3_occurrence_aliases(
    decoder: ModeDataDecoder,
    *,
    sg: int,
    basis: tuple[int, ...],
    operations: tuple[tuple[int, int, int, int, int], ...],
    rows: Any,
    selected_cases: list[tuple[int, tuple[int, int, int, int]]],
    irrep_source: OpdSourceData,
) -> list[Any]:
    """Restore distinct type-3 occurrences hidden by family orbit dedupe.

    A selected raw occurrence is promoted first.  Replay is only needed when
    that exact occurrence still has no Source row; type-1 real carriers already
    combine these aliases and must not be duplicated.
    """

    replayed_rows = list(rows)
    for selected_gid, selected_source_kparam in selected_cases:
        if _primary_dynamic_case_is_listed(
            replayed_rows,
            primary_gid=selected_gid,
            primary_source_kparam=selected_source_kparam,
        ):
            continue
        try:
            if int(decoder.little_record_by_gid(int(selected_gid)).irrep_type) != 3:
                continue
        except (IndexError, KeyError, TypeError, ValueError):
            continue
        selected_physical = tuple(
            value % 1
            for value in _source_kparam_identity(selected_source_kparam)
        )
        seen = {
            tuple(value % 1 for value in _source_kparam_identity(row.source_kparam))
            for row in replayed_rows
            if int(row.gid) == int(selected_gid)
        }
        seen.add(selected_physical)
        for row in tuple(replayed_rows):
            if int(row.gid) != int(selected_gid):
                continue
            for alias in dynamic_subduction_occurrence_alias_rows(
                decoder,
                sg=int(sg),
                basis=basis,
                operations=operations,
                row=row,
                irrep_source=irrep_source,
            ):
                physical = tuple(
                    value % 1 for value in _source_kparam_identity(alias.source_kparam)
                )
                if physical in seen:
                    continue
                seen.add(physical)
                replayed_rows.append(alias)
    return replayed_rows


def _order_promoted_dynamic_families_for_presentation(
    decoder: ModeDataDecoder,
    rows: list[Any],
    *,
    promoted_gids: set[int],
) -> list[Any]:
    """Order remaining occurrences by their visible family parameters.

    Coupled primary rows are emitted first in request order.  When the selected
    raw occurrence replaced a deduped formal row, the remaining rows of that
    k-family follow the presentation parameter order rather than reciprocal
    mesh order.  Other families and unresolved parameter rows retain their
    Source order.
    """

    if not promoted_gids:
        return list(rows)
    promoted_labels = {
        str(row.k_label)
        for row in rows
        if int(row.gid) in promoted_gids
    }
    ordered = list(rows)
    for label in promoted_labels:
        indices = [
            index
            for index, row in enumerate(ordered)
            if str(row.k_label) == label
        ]
        if len(indices) <= 1:
            continue
        candidates: list[tuple[tuple[Fraction, ...], int, Any]] = []
        for source_order, index in enumerate(indices):
            row = ordered[index]
            occurrence_vectors = [
                tuple(occurrence.reciprocal_vector_pml)
                for occurrence in tuple(getattr(row, "source_occurrences", ()))
            ] or [tuple(row.reciprocal_vector_pml)]
            occurrence_params = [
                _pml_vector_to_case_k_params(
                    decoder,
                    gid=int(row.gid),
                    reciprocal_vector_pml=vector,
                )
                for vector in occurrence_vectors
            ]
            if any(params is None for params in occurrence_params):
                candidates = []
                break
            canonical_params = min(
                tuple(params)
                for params in occurrence_params
                if params is not None
            )
            candidates.append((canonical_params, source_order, row))
        if not candidates:
            continue
        candidates.sort(key=lambda item: (item[0], item[1]))
        for index, (_params, _source_order, row) in zip(indices, candidates, strict=True):
            ordered[index] = row
    return ordered


def _ordinary_dynamic_subduction_spec(
    decoder: ModeDataDecoder,
    row: DynamicSubductionRow,
    *,
    old_id: int,
    primary: bool,
) -> dict[str, Any]:
    """Expose one computed ordinary subduction through the common mode schema."""

    reciprocal_vector_pml = tuple(row.reciprocal_vector_pml)
    direction_matrix = [list(values) for values in row.direction_matrix]
    return {
        "old_id": int(old_id),
        "gid": int(row.gid),
        "label": str(row.irrep_label),
        "magnetic": False,
        "k_label": str(row.k_label),
        "row_id": None,
        "source_kparam": tuple(int(value) for value in row.source_kparam),
        "carrier_source_kparam": row.carrier_source_kparam,
        "reciprocal_vector_pml": reciprocal_vector_pml,
        "request_k_params": _pml_vector_to_case_k_params(
            decoder,
            gid=int(row.gid),
            reciprocal_vector_pml=reciprocal_vector_pml,
        ),
        "source_occurrences": tuple(row.source_occurrences),
        "direction_matrix": direction_matrix,
        "source_numeric_rows": [
            [
                float(direction_matrix[coordinate][free])
                for coordinate in range(len(direction_matrix))
            ]
            for free in range(len(direction_matrix[0]) if direction_matrix else 0)
        ],
        "frequency": 1,
        "primary": bool(primary),
    }



def _subduced_mode_specs(
    decoder: ModeDataDecoder,
    sg: int,
    selected_irrep: dict[str, Any],
    selected_opd: dict[str, Any] | None,
) -> list[dict[str, Any]]:
    """Return primary + secondary irrep rows for complete mode details."""

    selected_old_id = int(selected_irrep.get("old_id") or 0)
    selected_label = str(
        selected_irrep.get("ordinary_symbol")
        or selected_irrep.get("symbol")
        or selected_irrep.get("label")
        or ""
    )
    selected_row_id = _isotropy_row_id_from_opd_row(selected_opd)
    selected_iso = _isotropy_from_opd_row(selected_opd) or {}
    fallback_source_rows = (
        selected_iso.get("source_numeric_rows")
        or selected_iso.get("source_matrix")
        or []
    )
    fallback_matrix = np.asarray(fallback_source_rows, dtype=float)
    fallback_direction_matrix = (
        fallback_matrix.T.tolist()
        if fallback_matrix.ndim == 2 and fallback_matrix.size
        else None
    )
    fallback = {
        "old_id": selected_old_id,
        "gid": _selected_dynamic_gid(selected_irrep, selected_opd),
        "label": selected_label,
        "k_label": _k_label_from_irrep_label(selected_label),
        "row_id": selected_row_id,
        "frequency": 1,
        "primary": True,
        "source_kparam": _source_kparam_record(selected_opd),
        "source_numeric_rows": fallback_source_rows,
        "source_free_count": int(selected_iso.get("free") or 0),
        "direction_matrix": fallback_direction_matrix,
    }
    if selected_irrep.get("magnetic"):
        return _ordinary_subduced_mode_specs_from_magnetic_embedding(
            decoder,
            sg,
            selected_irrep,
            selected_opd,
        )
    if selected_old_id <= 0:
        child_sg = _selected_subgroup_number(selected_opd)
        basis = _integer_basis_tuple(_source_split_basis_from_opd_row(selected_opd))
        origin = _source_split_origin_from_opd_row(selected_opd)
        if child_sg is None or basis is None or origin is None:
            return [fallback]
        operations = _subgroup_parent_operation_records(
            decoder,
            int(sg),
            int(child_sg),
            basis,
            origin,
        )
        selected_gid = _selected_dynamic_gid(selected_irrep, selected_opd)
        selected_source_kparam = _source_kparam_record(selected_opd)
        irrep_source = _opd_source_data(decoder)
        raw_rows = tuple(
            dynamic_subduction_rows(
                decoder,
                sg=int(sg),
                basis=basis,
                operations=operations,
                irrep_source=irrep_source,
            )
        )
        rows = list(raw_rows)
        if selected_gid is not None and selected_source_kparam is not None:
            rows = _promote_selected_dynamic_occurrences(
                decoder,
                sg=int(sg),
                basis=basis,
                operations=operations,
                rows=rows,
                selected_cases=[(int(selected_gid), selected_source_kparam)],
                irrep_source=irrep_source,
                preserve_discovered_carrier=False,
            )
            rows = _replay_unrepresented_type3_occurrence_aliases(
                decoder,
                sg=int(sg),
                basis=basis,
                operations=operations,
                rows=rows,
                selected_cases=[(int(selected_gid), selected_source_kparam)],
                irrep_source=irrep_source,
            )

        def spec_for_row(row: Any) -> dict[str, Any]:
            return _ordinary_dynamic_subduction_spec(
                decoder,
                row,
                old_id=0,
                primary=int(row.gid) == int(selected_gid or 0)
                and (
                    selected_source_kparam is None
                    or _same_source_kparam(row.source_kparam, selected_source_kparam)
                ),
            )

        specs = [spec_for_row(row) for row in rows]
        specs.extend(
            _stored_occurrence_alias_specs(
                decoder,
                raw_rows,
                rows,
                aliases_for_class=lambda row, _requested: (
                    alias
                    for alias in dynamic_subduction_occurrence_alias_rows(
                        decoder,
                        sg=int(sg),
                        basis=basis,
                        operations=operations,
                        row=row,
                        irrep_source=irrep_source,
                        requested_source_kparam=_requested,
                    )
                    if tuple(int(value) for value in alias.source_kparam)
                    == tuple(int(value) for value in _requested)
                ),
                spec_for_row=spec_for_row,
            )
        )
        specs = _canonicalize_exact_pml_alias_specs(specs)
        for spec in specs:
            if (
                spec["primary"]
                and selected_source_kparam is not None
                and tuple(spec.get("source_kparam") or ()) == tuple(selected_source_kparam)
            ):
                spec["source_free_count"] = int(selected_iso.get("free") or 0)
        if not specs:
            return [fallback]
        if not any(spec["primary"] for spec in specs):
            specs.insert(0, fallback)
        specs.sort(key=lambda item: 0 if item["primary"] else 1)
        return specs
    specs: list[dict[str, Any]] = []
    seen: set[tuple[int, int | None]] = set()
    for entry in decoder.isotropy_subductions_for_row(selected_row_id):
        if int(entry.sg) != int(sg):
            continue
        key = (int(entry.irrep_old_id), int(entry.subgroup_row_id))
        if key in seen:
            continue
        seen.add(key)
        entry_label = str(entry.irrep_label)
        canonical_little = next(
            (
                little
                for little in little_records_for_k(
                    decoder,
                    int(sg),
                    _k_label_from_irrep_label(entry_label),
                )
                if int(little.old_id) == int(entry.irrep_old_id)
            ),
            None,
        )
        specs.append(
            {
                "old_id": int(entry.irrep_old_id),
                "gid": None if canonical_little is None else int(canonical_little.gid),
                "label": entry_label if canonical_little is None else str(canonical_little.label),
                "k_label": _k_label_from_irrep_label(entry_label),
                "row_id": int(entry.subgroup_row_id),
                "opd": str(
                    decoder.iso.isotropy["isotropy_orderparam_label"][
                        int(entry.subgroup_row_id) - 1
                    ]
                ).strip(),
                "domain": int(entry.domain),
                "domain_old": int(entry.domain_old),
                "frequency": int(entry.frequency),
                "primary": (
                    int(entry.irrep_old_id) == selected_old_id
                    and int(entry.subgroup_row_id) == int(selected_row_id)
                ),
            }
        )
    if selected_iso.get("embedding_selected") and fallback_source_rows:
        for spec in specs:
            if not spec.get("primary"):
                continue
            spec["source_numeric_rows"] = fallback_source_rows
            spec["direction_matrix"] = fallback_direction_matrix
    for spec in specs:
        spec["source_free_count"] = int(
            decoder.isotropy_orderparam_freeparam(int(spec["row_id"]))
        )
    child_sg = _selected_subgroup_number(selected_opd)
    basis = _integer_basis_tuple(_source_split_basis_from_opd_row(selected_opd))
    origin = _source_split_origin_from_opd_row(selected_opd)
    if child_sg is not None and basis is not None and origin is not None:
        operations = _subgroup_parent_operation_records(
            decoder,
            int(sg),
            int(child_sg),
            basis,
            origin,
        )
        selected_rows = dynamic_subduction_rows(
            decoder,
            sg=int(sg),
            basis=basis,
            operations=operations,
            irrep_source=_opd_source_data(decoder),
        )
        selected_by_old_id = {
            int(decoder.little_record_by_gid(int(row.gid)).old_id): row
            for row in selected_rows
            if int(decoder.little_record_by_gid(int(row.gid)).old_id) > 0
        }
        for spec in specs:
            if spec.get("primary"):
                continue
            selected_row = selected_by_old_id.get(int(spec.get("old_id") or 0))
            if selected_row is None:
                continue
            spec.update(
                {
                    "gid": int(selected_row.gid),
                    "source_kparam": tuple(int(value) for value in selected_row.source_kparam),
                    "reciprocal_vector_pml": tuple(selected_row.reciprocal_vector_pml),
                    "direction_matrix": [list(values) for values in selected_row.direction_matrix],
                    "source_numeric_rows": [
                        [
                            float(selected_row.direction_matrix[coordinate][free])
                            for coordinate in range(len(selected_row.direction_matrix))
                        ]
                        for free in range(
                            len(selected_row.direction_matrix[0]) if selected_row.direction_matrix else 0
                        )
                    ],
                }
            )
        specs.extend(
            _ordinary_dynamic_subduction_spec(decoder, row, old_id=0, primary=False)
            for row in selected_rows
            if int(decoder.little_record_by_gid(int(row.gid)).old_id) == 0
        )
    specs.sort(key=lambda item: 0 if item["primary"] else 1)
    return specs



def _ordinary_subduced_mode_specs_from_magnetic_embedding(
    decoder: ModeDataDecoder,
    sg: int,
    selected_irrep: dict[str, Any],
    selected_opd: dict[str, Any] | None,
) -> list[dict[str, Any]]:
    """Return time-even ordinary modes fixed by one magnetic embedding."""

    iso = _isotropy_from_opd_row(selected_opd) or {}
    selected_gid = _selected_dynamic_gid(selected_irrep, selected_opd)
    source_kparam = _source_kparam_record(selected_opd)
    source_rows = iso.get("source_numeric_rows") or iso.get("source_matrix") or []
    if selected_gid is None:
        return []
    if source_kparam is None and int(decoder.little_record_by_gid(selected_gid).old_id) > 0:
        source_kparam = (0, 0, 0, 1)
    if source_kparam is None or not source_rows:
        return []
    full_dim = int(decoder.little_record_by_gid(selected_gid).full_dim)
    orderparam: list[float] = []
    for row in source_rows:
        values = [float(value) for value in row]
        orderparam.extend(values[:full_dim] + [0.0] * (48 - full_dim))
    source = _opd_source_data(decoder)
    basis, magnetic_operations = source.orderparam_to_subgroup_magnetic(
        int(selected_gid),
        source_kparam,
        orderparam,
        len(source_rows),
    )
    table = magnetic_data().table
    ordinary_operations = tuple(
        (
            record[0],
            record[1],
            record[2],
            record[3],
            int(table["mag_point_op_mag2nonmag"][record[4] - 1]),
        )
        for record in magnetic_operations
    )
    rows = dynamic_subduction_rows(
        decoder,
        sg=int(sg),
        basis=basis,
        operations=ordinary_operations,
        irrep_source=source,
    )

    def row_spec(row: DynamicSubductionRow) -> dict[str, Any]:
        reciprocal_vector_pml = tuple(row.reciprocal_vector_pml)
        return {
            "old_id": int(decoder.little_record_by_gid(row.gid).old_id),
            "gid": int(row.gid),
            "label": str(row.irrep_label),
            "k_label": str(row.k_label),
            "row_id": None,
            "source_kparam": tuple(int(value) for value in row.source_kparam),
            "reciprocal_vector_pml": reciprocal_vector_pml,
            "request_k_params": _pml_vector_to_case_k_params(
                decoder,
                gid=int(row.gid),
                reciprocal_vector_pml=reciprocal_vector_pml,
            ),
            "source_occurrences": tuple(row.source_occurrences),
            "direction_matrix": [list(values) for values in row.direction_matrix],
            "source_numeric_rows": [
                [float(row.direction_matrix[coordinate][free]) for coordinate in range(len(row.direction_matrix))]
                for free in range(len(row.direction_matrix[0]) if row.direction_matrix else 0)
            ],
            "frequency": 1,
            "primary": False,
        }

    specs = [row_spec(row) for row in rows]
    alias_specs = _stored_occurrence_alias_specs(
        decoder,
        rows,
        rows,
        aliases_for_class=lambda row, requested: (
            alias
            for alias in dynamic_subduction_occurrence_alias_rows(
                decoder,
                sg=int(sg),
                basis=basis,
                operations=ordinary_operations,
                row=row,
                irrep_source=source,
                requested_source_kparam=requested,
            )
            if tuple(int(value) for value in alias.source_kparam)
            == tuple(int(value) for value in requested)
        ),
        spec_for_row=row_spec,
    )
    specs = _merge_occurrence_alias_specs(specs, alias_specs)
    return _canonicalize_exact_pml_alias_specs(specs)


def _selected_primary_stokes_direction(
    *,
    primary: bool,
    magnetic: bool,
    coupled: bool,
    little_type: int,
    selected_source_rows: Any,
    recomputed_direction: list[list[float]],
    carrier_source_kparam: Sequence[int] | None,
    selected_source_rounded: bool = False,
) -> tuple[list[list[float]], list[list[float]] | None]:
    """Retain a selected Source Stokes frame when it spans the same carrier."""

    if not (
        primary
        and magnetic
        and not coupled
        and int(little_type) in {1, 3}
        and isinstance(selected_source_rows, list)
        and len(selected_source_rows) > 1
    ):
        return recomputed_direction, None
    try:
        selected_rows = np.asarray(selected_source_rows, dtype=float)
        recomputed = np.asarray(recomputed_direction, dtype=float)
    except (TypeError, ValueError) as exc:
        raise ValueError("selected Stokes direction is not numeric") from exc
    if (
        selected_rows.ndim != 2
        or recomputed.ndim != 2
        or selected_rows.size == 0
        or recomputed.size == 0
        or not np.all(np.isfinite(selected_rows))
        or not np.all(np.isfinite(recomputed))
    ):
        raise ValueError("selected Stokes direction is not a finite matrix")
    selected = selected_rows.T
    if selected.shape != recomputed.shape:
        raise ValueError("selected and recomputed Stokes directions differ in shape")
    uses_dynamic_rounding_cell = (
        selected_source_rounded and carrier_source_kparam is None
    )
    free_count = int(selected.shape[1])
    if (
        free_count <= 1
        or int(np.linalg.matrix_rank(selected, tol=1e-10)) != free_count
        or int(np.linalg.matrix_rank(recomputed, tol=1e-10)) != free_count
    ):
        raise ValueError("selected or recomputed Stokes direction lost rank")
    if np.array_equal(selected, recomputed):
        return recomputed_direction, None
    if not uses_dynamic_rounding_cell and (
        int(
            np.linalg.matrix_rank(
                np.column_stack((recomputed, selected)), tol=1e-10
            )
        )
        != free_count
    ):
        # Occurrence promotion can attach the selected heading to a conjugate
        # carrier.  Its raw Source columns need an explicit occurrence
        # intertwiner before they can replace that carrier's computed frame.
        if carrier_source_kparam is not None:
            return recomputed_direction, None
        raise ValueError("selected and recomputed Stokes directions differ in span")
    transform, _residuals, transform_rank, _singular = np.linalg.lstsq(
        recomputed, selected, rcond=None
    )
    scale = max(1.0, float(np.max(np.abs(selected))))
    residual_limit = (
        _DYNAMIC_SOURCE_MATRIX_ROUNDING_CELL
        if uses_dynamic_rounding_cell
        else 1e-10 * scale
    )
    transform_is_invalid = (
        int(transform_rank) != free_count
        or int(np.linalg.matrix_rank(transform, tol=1e-10)) != free_count
        or float(np.max(np.abs(recomputed @ transform - selected)))
        > residual_limit
    )
    if transform_is_invalid:
        raise ValueError("selected Stokes direction has no exact carrier transform")
    return selected.tolist(), recomputed.tolist()


def _selected_stokes_rows_use_dynamic_display_precision(
    selected_opd: dict[str, Any] | None,
    selected_iso: dict[str, Any],
) -> bool:
    direction = (
        selected_opd.get("direction") if isinstance(selected_opd, dict) else None
    )
    return bool(
        isinstance(direction, dict)
        and direction.get("dynamic") is True
        and not selected_iso.get("embedding_selected")
    )


def _magnetic_subduced_mode_specs(
    decoder: ModeDataDecoder,
    sg: int,
    selected_irrep: dict[str, Any],
    selected_opd: dict[str, Any] | None,
) -> list[dict[str, Any]]:
    """Return magnetic primary/secondary rows for one selected embedding."""

    iso = _isotropy_from_opd_row(selected_opd) or {}
    selected_gid = _selected_dynamic_gid(selected_irrep, selected_opd)
    source_kparam = _source_kparam_record(selected_opd)
    source_rows = iso.get("source_numeric_rows") or iso.get("source_matrix") or []
    if selected_gid is not None and source_kparam is None and int(decoder.little_record_by_gid(selected_gid).old_id) > 0:
        source_kparam = (0, 0, 0, 1)
    if selected_gid is None or source_kparam is None or not source_rows:
        return []
    selected_k_label = str(
        selected_irrep.get("k_label")
        or _k_label_from_irrep_label(str(selected_irrep.get("ordinary_symbol") or ""))
    )
    full_dim = int(decoder.little_record_by_gid(selected_gid).full_dim)
    orderparam: list[float] = []
    for row in source_rows:
        values = [float(value) for value in row]
        orderparam.extend(values[:full_dim] + [0.0] * (48 - full_dim))
    magnetic_source = _opd_source_data(decoder)
    basis, operations = magnetic_source.orderparam_to_subgroup_magnetic(
        int(selected_gid),
        source_kparam,
        orderparam,
        len(source_rows),
    )
    raw_rows = dynamic_magnetic_subduction_rows(
        decoder,
        magnetic_source,
        sg=int(sg),
        basis=basis,
        operations=operations,
        preferred_kparams={selected_k_label: source_kparam},
    )
    rows = _promote_selected_magnetic_occurrences(
        decoder,
        magnetic_source,
        sg=int(sg),
        basis=basis,
        operations=operations,
        rows=raw_rows,
        selected_cases=[(selected_gid, source_kparam, None)],
    )
    selected_source_rounded = _selected_stokes_rows_use_dynamic_display_precision(
        selected_opd,
        iso,
    )

    def row_spec(row: Any) -> dict[str, Any]:
        reciprocal_vector_pml = tuple(row.reciprocal_vector_pml)
        recomputed_direction = [list(values) for values in row.direction_matrix]
        little = decoder.little_record_by_gid(row.gid)
        primary = _selected_source_occurrence_matches(
            gid=int(row.gid),
            source_kparam=tuple(int(value) for value in row.source_kparam),
            selected_gid=selected_gid,
            selected_source_kparam=source_kparam,
        )
        direction_matrix, retained_recomputed = _selected_primary_stokes_direction(
            primary=primary,
            magnetic=True,
            coupled=False,
            little_type=int(little.irrep_type),
            selected_source_rows=source_rows,
            recomputed_direction=recomputed_direction,
            carrier_source_kparam=row.carrier_source_kparam,
            selected_source_rounded=selected_source_rounded,
        )
        source_free_count = (
            len(direction_matrix[0])
            if direction_matrix
            and all(
                len(values) == len(direction_matrix[0])
                for values in direction_matrix
            )
            else 0
        )
        spec = {
            "old_id": int(little.old_id),
            "gid": int(row.gid),
            "label": str(row.irrep_label),
            "display_label": "m" + str(row.irrep_label),
            "magnetic": True,
            "k_label": str(row.k_label),
            "row_id": None,
            "source_kparam": tuple(int(value) for value in row.source_kparam),
            "carrier_source_kparam": row.carrier_source_kparam,
            "reciprocal_vector_pml": reciprocal_vector_pml,
            "request_k_params": _pml_vector_to_case_k_params(
                decoder,
                gid=int(row.gid),
                reciprocal_vector_pml=reciprocal_vector_pml,
            ),
            "source_occurrences": tuple(row.source_occurrences),
            "direction_matrix": direction_matrix,
            "source_numeric_rows": [
                [
                    float(direction_matrix[coordinate][free])
                    for coordinate in range(len(direction_matrix))
                ]
                for free in range(source_free_count)
            ],
            "source_free_count": source_free_count,
            "frequency": 1,
            "primary": primary,
        }
        if retained_recomputed is not None:
            spec["_selected_stokes_columns"] = True
            spec["_recomputed_direction_matrix"] = retained_recomputed
        return spec

    specs = [row_spec(row) for row in rows]
    alias_specs = _stored_occurrence_alias_specs(
        decoder,
        raw_rows,
        rows,
        aliases_for_class=lambda row, requested: (
            dynamic_magnetic_subduction_occurrence_alias_rows(
                decoder,
                magnetic_source,
                sg=int(sg),
                basis=basis,
                operations=operations,
                row=row,
                requested_source_kparam=requested,
            )
        ),
        spec_for_row=row_spec,
    )
    specs = _merge_occurrence_alias_specs(specs, alias_specs)
    specs = _canonicalize_exact_pml_alias_specs(specs)
    specs.sort(key=lambda item: 0 if item["primary"] else 1)
    return specs



def _coupled_render_specs(
    decoder: ModeDataDecoder,
    sg: int,
    selected_slots: list[dict[str, Any]],
    selected_opd: dict[str, Any] | None,
    *,
    include_displacive: bool,
    include_magnetic: bool,
) -> list[tuple[dict[str, Any], str, int, str]]:
    """Enumerate complete-mode rows from the exact coupled embedding."""

    iso = _isotropy_from_opd_row(selected_opd) or {}
    basis = _integer_basis_tuple(_source_split_basis_from_opd_row(selected_opd))
    raw_operations = iso.get("source_operation_records") or []
    operations = tuple(
        tuple(int(value) for value in record)
        for record in raw_operations
        if isinstance(record, (list, tuple)) and len(record) == 5
    )
    if basis is None or not operations:
        return []
    subgroup = iso.get("subgroup") if isinstance(iso.get("subgroup"), dict) else {}
    magnetic_embedding = bool(subgroup.get("ordinary_number"))
    table = magnetic_data().table if magnetic_embedding else None
    ordinary_operations = (
        tuple(
            (
                record[0],
                record[1],
                record[2],
                record[3],
                int(table["mag_point_op_mag2nonmag"][record[4] - 1]),
            )
            for record in operations
        )
        if table is not None
        else operations
    )
    primary_slots = []
    selected_magnetic_cases: list[
        tuple[int, tuple[int, int, int, int], tuple[Fraction, ...]]
    ] = []
    for index, slot in enumerate(selected_slots):
        magnetic_identity = _selected_magnetic_case_identity(slot)
        mode_kind = "mag" if bool((slot.get("irrep") or {}).get("magnetic")) else "dsp"
        if mode_kind == "mag":
            if magnetic_identity is None:
                continue
            gid, source_kparam = magnetic_identity
            selected_magnetic_cases.append(
                (gid, source_kparam, _k_params(slot.get("k_params") or {}))
            )
        else:
            gid = int((slot.get("irrep") or {}).get("gid") or 0)
            source_kparam = tuple(
                int(value)
                for value in slot.get("source_kparam") or (0, 0, 0, 1)
            )
        primary_slots.append(
            {
                "order": index,
                "mode_kind": mode_kind,
                "gid": gid,
                "k_label": str((slot.get("kpoint") or {}).get("label") or ""),
                "source_kparam": source_kparam,
                "case_k_params": _k_params(slot.get("k_params") or {}),
                "display_case_k_params": _k_params(
                    slot.get("display_k_params")
                    or slot.get("k_params")
                    or {}
                ),
            }
        )
    specs: list[tuple[int, dict[str, Any], str, int, str]] = []
    source_order = 0
    listed_magnetic_primary_orders: set[int] = set()

    def primary_slot_for_row(
        row: DynamicSubductionRow, mode_kind: str
    ) -> dict[str, Any] | None:
        candidates = tuple(
            slot
            for slot in primary_slots
            if slot["mode_kind"] == mode_kind and int(slot["gid"]) == int(row.gid)
        )
        for slot in candidates:
            if _same_source_kparam(slot["source_kparam"], row.source_kparam):
                return slot
            # Promoted dynamic rows store kvec_standard(selected). Match the
            # selected input slot through that standard so coupled headings keep
            # display_case_k_params → Source star[0] (74c C2: Web [1/3,5/6,0]
            # vs Local [1/6,-2/3,0] when standard (5,-3,0,12) was unmatched).
            try:
                standardized = kvec_standard_provenance(
                    decoder,
                    int(sg),
                    int(row.gid),
                    tuple(int(value) for value in slot["source_kparam"]),
                ).standard_kparam
            except (IndexError, KeyError, TypeError, ValueError, ZeroDivisionError):
                continue
            if _same_source_kparam(standardized, row.source_kparam):
                return slot
        if mode_kind == "mag":
            for slot in candidates:
                if (
                    int(slot["order"]) not in listed_magnetic_primary_orders
                    and _row_projects_selected_case_k_params(
                        decoder,
                        row,
                        selected_gid=int(slot["gid"]),
                        selected_case_k_params=tuple(slot["case_k_params"]),
                    )
                ):
                    return slot
        return None

    def spec_for_row(row: Any, mode_kind: str) -> dict[str, Any]:
            gid = int(row.gid)
            primary_slot = primary_slot_for_row(row, mode_kind)
            little_record = decoder.little_record_by_gid(gid)
            old_id = int(little_record.old_id)
            little_type = int(little_record.irrep_type)
            direction_matrix = [list(values) for values in row.direction_matrix]
            source_free_count = (
                len(direction_matrix[0])
                if direction_matrix
                and all(len(values) == len(direction_matrix[0]) for values in direction_matrix)
                else 0
            )
            return {
                "old_id": old_id,
                "little_type": little_type,
                "gid": gid,
                "label": str(row.irrep_label),
                "display_label": ("m" if mode_kind == "mag" else "") + str(row.irrep_label),
                "magnetic": mode_kind == "mag",
                "k_label": str(row.k_label),
                "row_id": None,
                "source_kparam": tuple(int(value) for value in row.source_kparam),
                "carrier_source_kparam": row.carrier_source_kparam,
                "reciprocal_vector_pml": tuple(row.reciprocal_vector_pml),
                "direction_matrix": direction_matrix,
                "source_numeric_rows": [
                    [float(direction_matrix[coordinate][free]) for coordinate in range(len(direction_matrix))]
                    for free in range(source_free_count)
                ],
                "source_free_count": source_free_count,
                "frequency": 1,
                "primary": primary_slot is not None,
                "coupled": True,
                "case_k_params": None if primary_slot is None else primary_slot["case_k_params"],
                "display_case_k_params": (
                    None
                    if primary_slot is None
                    else primary_slot["display_case_k_params"]
                ),
                "_primary_slot_order": (
                    None if primary_slot is None else primary_slot["order"]
                ),
                "source_occurrences": tuple(row.source_occurrences),
            }

    def append_specs(row_specs: Any, mode_kind: str) -> None:
        nonlocal source_order
        for spec in row_specs:
            primary_slot = (
                next(
                    (
                        slot
                        for slot in primary_slots
                        if slot["mode_kind"] == mode_kind
                        and int(slot["gid"]) == int(spec["gid"])
                        and tuple(slot["case_k_params"])
                        == tuple(spec.get("case_k_params") or ())
                    ),
                    None,
                )
                if spec.get("primary") is True
                else None
            )
            spec["primary"] = bool(spec.get("primary") and primary_slot is not None)
            spec["_primary_slot_order"] = (
                None if primary_slot is None else int(primary_slot["order"])
            )
            rank = int(primary_slot["order"]) if primary_slot is not None else len(primary_slots) + source_order
            specs.append(
                (
                    (0 if spec["primary"] else 1) * 1_000_000 + rank,
                    spec,
                    mode_kind,
                    2 if mode_kind == "mag" else 1,
                    str(spec["display_label"]),
                )
            )
            source_order += 1

    if include_displacive:
        ordinary_source = _opd_source_data(decoder)
        preferred: dict[str, list[tuple[int, int, int, int]]] = {}
        for slot in selected_slots:
            if bool((slot.get("irrep") or {}).get("magnetic")):
                continue
            preferred.setdefault(str((slot.get("kpoint") or {}).get("label") or ""), []).append(
                tuple(int(value) for value in slot.get("source_kparam") or (0, 0, 0, 1))
            )
        ordinary_raw_rows = dynamic_subduction_rows(
            decoder,
            sg=int(sg),
            basis=basis,
            operations=ordinary_operations,
            preferred_kparams=preferred,
            irrep_source=ordinary_source,
        )
        selected_cases = [
            (int(slot["gid"]), tuple(slot["source_kparam"]))
            for slot in primary_slots
            if slot["mode_kind"] == "dsp" and int(slot["gid"]) > 0
        ]
        promoted_gids = {
            int(gid)
            for gid, source_kparam in selected_cases
            if not _primary_dynamic_case_is_listed(
                ordinary_raw_rows,
                primary_gid=int(gid),
                primary_source_kparam=tuple(source_kparam),
            )
        }
        ordinary_rows = _promote_selected_dynamic_occurrences(
            decoder,
            sg=int(sg),
            basis=basis,
            operations=ordinary_operations,
            rows=ordinary_raw_rows,
            selected_cases=selected_cases,
            irrep_source=ordinary_source,
            preserve_discovered_carrier=True,
        )
        ordinary_rows = _order_promoted_dynamic_families_for_presentation(
            decoder,
            ordinary_rows,
            promoted_gids=promoted_gids,
        )
        ordinary_specs = [spec_for_row(row, "dsp") for row in ordinary_rows]
        primary_families: dict[tuple[Any, ...], list[dict[str, Any]]] = {}
        for spec in ordinary_specs:
            if spec.get("primary") is not True:
                continue
            occurrences = tuple(spec.get("source_occurrences") or ())
            if occurrences:
                primary_families.setdefault(occurrences, []).append(spec)
        for family_specs in primary_families.values():
            first = min(
                family_specs,
                key=lambda spec: int(spec["_primary_slot_order"]),
            )
            displayed = first["display_case_k_params"]
            for spec in family_specs:
                spec["display_case_k_params"] = displayed
        ordinary_specs.extend(
            _stored_occurrence_alias_specs(
                decoder,
                ordinary_raw_rows,
                ordinary_rows,
                aliases_for_class=lambda row, requested: tuple(
                    alias
                    for alias in dynamic_subduction_occurrence_alias_rows(
                        decoder,
                        sg=int(sg),
                        basis=basis,
                        operations=ordinary_operations,
                        row=row,
                        irrep_source=ordinary_source,
                        requested_source_kparam=requested,
                    )
                    if tuple(int(value) for value in alias.source_kparam)
                    == tuple(int(value) for value in requested)
                ),
                spec_for_row=lambda row: spec_for_row(row, "dsp"),
            )
        )
        append_specs(
            _canonicalize_exact_pml_alias_specs(ordinary_specs),
            "dsp",
        )
    if include_magnetic and magnetic_embedding:
        magnetic_source = _opd_source_data(decoder)
        preferred = {}
        for slot in primary_slots:
            if slot["mode_kind"] != "mag":
                continue
            preferred.setdefault(str(slot["k_label"]), []).append(
                tuple(slot["source_kparam"])
            )
        magnetic_raw_rows = dynamic_magnetic_subduction_rows(
            decoder,
            magnetic_source,
            sg=int(sg),
            basis=basis,
            operations=operations,
            preferred_kparams=preferred,
        )
        magnetic_rows = _promote_selected_magnetic_occurrences(
            decoder,
            magnetic_source,
            sg=int(sg),
            basis=basis,
            operations=operations,
            rows=magnetic_raw_rows,
            selected_cases=selected_magnetic_cases,
        )
        listed_magnetic_primary_orders.update(
            int(slot["order"])
            for slot in primary_slots
            if slot["mode_kind"] == "mag"
            and any(
                int(row.gid) == int(slot["gid"])
                and _same_source_kparam(row.source_kparam, slot["source_kparam"])
                for row in magnetic_raw_rows
            )
        )
        magnetic_specs = [spec_for_row(row, "mag") for row in magnetic_rows]
        magnetic_specs.extend(
            _stored_occurrence_alias_specs(
                decoder,
                magnetic_raw_rows,
                magnetic_rows,
                aliases_for_class=lambda row, requested: (
                    dynamic_magnetic_subduction_occurrence_alias_rows(
                        decoder,
                        magnetic_source,
                        sg=int(sg),
                        basis=basis,
                        operations=operations,
                        row=row,
                        requested_source_kparam=requested,
                    )
                ),
                spec_for_row=lambda row: spec_for_row(row, "mag"),
            )
        )
        append_specs(
            _canonicalize_exact_pml_alias_specs(magnetic_specs),
            "mag",
        )
    specs.sort(key=lambda item: item[0])
    return [(spec, kind, setting, label) for _rank, spec, kind, setting, label in specs]
