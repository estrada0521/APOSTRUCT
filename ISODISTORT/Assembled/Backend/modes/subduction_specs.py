"""Mode-subduction specification helpers.

Extracted mechanically from the former monolithic runtime.
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import replace
from fractions import Fraction
from numbers import Integral
from typing import Any
import numpy as np
from ISODISTORT.Assembled.Backend.source.magnetic import data as magnetic_data
from ISODISTORT.Assembled.Backend.modes.engine.decoder import ModeDataDecoder
from ISODISTORT.Assembled.Backend.modes.engine.dynamic_subduction import (
    DynamicSubductionRow,
    dynamic_magnetic_subduction_occurrence_alias_rows,
    dynamic_magnetic_subduction_rows,
    dynamic_subduction_occurrence_alias_rows,
    dynamic_subduction_rows,
    kvec_standard_provenance,
    stored_occurrence_alias_candidates,
    _strict_integral_values,
)
from ISODISTORT.Assembled.Backend.modes.engine.project.mode_counts import (
    little_records_for_k,
)
from ISODISTORT.Assembled.Backend.isotropy.engine.source_data import SourceData as OpdSourceData

from ISODISTORT.Assembled.Backend.modes.common import (
    _integer_basis_tuple,
    _isotropy_from_opd_row,
    _isotropy_row_id_from_opd_row,
    _k_label_from_irrep_label,
    _k_params,
    _same_source_kparam,
)
from ISODISTORT.Assembled.Backend.modes.request_context import (
    _pml_vector_to_case_k_params,
    _selected_dynamic_gid,
    _source_kparam_record,
)
from ISODISTORT.Assembled.Backend.modes.structure_runtime import (
    _selected_subgroup_number,
    _source_split_basis_from_opd_row,
    _source_split_origin_from_opd_row,
    _subgroup_parent_operation_records,
)


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
    generic_name = "_occurrence_alias_observation_only"
    magnetic_name = "_magnetic_occurrence_alias_observation_only"
    generic_present = generic_name in spec
    magnetic_present = magnetic_name in spec
    if not generic_present and not magnetic_present:
        return True
    return bool(
        generic_present
        and spec.get(generic_name) is True
        and (not magnetic_present or spec.get(magnetic_name) is True)
    )


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
    represented_heading_identities: Sequence[object] = (),
) -> list[dict[str, Any]]:
    """Build fail-closed stored-occurrence observations and candidates."""

    represented = tuple(represented_rows)
    additional_headings = tuple(
        heading
        for value in represented_heading_identities
        if (heading := _strict_occurrence_alias_heading(value)) is not None
    )
    candidates = stored_occurrence_alias_candidates(
        decoder,
        tuple(raw_rows),
        aliases_for_class=aliases_for_class,
        heading_identity_for_row=_occurrence_heading_identity,
        represented_heading_identities=(
            *additional_headings,
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
    for row in raw_rows:
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
            and "_magnetic_occurrence_alias_candidate" not in specs[index]
            and "_occurrence_alias_observation_only" not in specs[index]
            and "_magnetic_occurrence_alias_observation_only" not in specs[index]
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
        magnetic_marker_names = (
            "_magnetic_occurrence_alias_anchor",
            "_magnetic_occurrence_alias_candidate",
            "_magnetic_occurrence_alias_heading_identity",
        )
        has_magnetic_markers = any(
            name in candidate for name in magnetic_marker_names
        )
        magnetic_markers_match = not has_magnetic_markers or (
            _strict_occurrence_alias_anchor(candidate.get(magnetic_marker_names[0]))
            == anchor_identity
            and _strict_source_kparam(candidate.get(magnetic_marker_names[1]))
            == source_kparam
            and _strict_occurrence_alias_heading(
                candidate.get(magnetic_marker_names[2])
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
            or not magnetic_markers_match
            or "_occurrence_alias_observation_only" in candidate
            or "_magnetic_occurrence_alias_observation_only" in candidate
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
            and "_magnetic_occurrence_alias_candidate" not in spec
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
            and not key.startswith("_magnetic_occurrence_alias_")
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
) -> list[Any]:
    """Replace deduped formal rows with their selected primary occurrence."""

    promoted_rows = list(rows)
    claimed_indices: set[int] = set()
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
            promoted_rows[row_index] = replace(
                selected_alias,
                source_occurrences=row.source_occurrences,
            )
            claimed_indices.add(row_index)
            break
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
                source_occurrences=row.source_occurrences,
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


def _replay_distinct_type3_occurrence_aliases(
    decoder: ModeDataDecoder,
    rows: Any,
    *,
    allowed_k_labels: set[str],
    aliases_for_row: Any,
) -> list[Any]:
    """Append invariant type-3 aliases from distinct reciprocal classes.

    Coupled complete-mode rows can consume more than the selected occurrence
    from one k family.  The Source alias builders are the authority for whether
    another raw occurrence remains invariant under the selected subgroup.  We
    only restore physical classes that the deduped family rows do not already
    represent; type-1 carriers and fixed table irreps retain their existing
    combined-row behavior.
    """

    replayed_rows = list(rows)
    seen_by_gid: dict[int, set[tuple[Fraction, Fraction, Fraction]]] = {}
    for row in replayed_rows:
        gid = int(row.gid)
        seen_by_gid.setdefault(gid, set()).add(
            tuple(value % 1 for value in _source_kparam_identity(row.source_kparam))
        )
    for row in tuple(replayed_rows):
        gid = int(row.gid)
        if str(row.k_label) not in allowed_k_labels:
            continue
        occurrence_classes = {
            tuple(
                value % 1
                for value in _source_kparam_identity(occurrence.source_kparam)
            )
            for occurrence in tuple(getattr(row, "source_occurrences", ()))
        }
        if len(occurrence_classes) <= 1:
            continue
        try:
            little_record = decoder.little_record_by_gid(gid)
            if int(little_record.old_id) != 0 or int(little_record.irrep_type) != 3:
                continue
        except (IndexError, KeyError, TypeError, ValueError):
            continue
        seen = seen_by_gid[gid]
        for alias in aliases_for_row(row):
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



def _subduced_mode_specs(
    decoder: ModeDataDecoder,
    sg: int,
    selected_irrep: dict[str, Any],
    selected_opd: dict[str, Any] | None,
) -> list[dict[str, Any]]:
    """Return primary + secondary irrep rows for complete mode details."""

    selected_old_id = int(selected_irrep.get("old_id") or 0)
    raw_selected_label = (
        selected_irrep.get("ordinary_symbol")
        or selected_irrep.get("symbol")
        or selected_irrep.get("label")
        or ""
    )
    selected_label = str(raw_selected_label)
    selected_heading_k_label = (
        _k_label_from_irrep_label(raw_selected_label)
        if isinstance(raw_selected_label, str)
        else None
    )
    selected_row_id = _isotropy_row_id_from_opd_row(selected_opd)
    selected_iso = _isotropy_from_opd_row(selected_opd) or {}
    fallback_source_rows = (
        selected_iso.get("source_numeric_rows")
        or selected_iso.get("source_matrix")
        or []
    )
    try:
        fallback_matrix = np.asarray(fallback_source_rows, dtype=float)
    except (TypeError, ValueError):
        fallback_matrix = np.empty((0, 0), dtype=float)
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
            reciprocal_vector_pml = tuple(row.reciprocal_vector_pml)
            return {
                "old_id": 0,
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
                "primary": int(row.gid) == int(selected_gid or 0)
                and (
                    selected_source_kparam is None
                    or _same_source_kparam(row.source_kparam, selected_source_kparam)
                ),
            }

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
                represented_heading_identities=(
                    ((int(selected_gid), selected_heading_k_label, selected_source_kparam),)
                    if selected_gid is not None
                    and selected_source_kparam is not None
                    and selected_heading_k_label is not None
                    else ()
                ),
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
    if selected_row_id is None:
        raise ValueError("static selected irrep has no exact Source isotropy row")

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
                "frequency": int(entry.frequency),
                "primary": (
                    int(entry.irrep_old_id) == selected_old_id
                    and int(entry.subgroup_row_id) == int(selected_row_id)
                ),
            }
        )
    if not specs:
        raise ValueError("static Source isotropy row has no subduction records")
    if not any(item["primary"] for item in specs):
        raise ValueError("static Source subduction records omit the selected primary")
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
    return [
        {
            "old_id": int(decoder.little_record_by_gid(row.gid).old_id),
            "gid": int(row.gid),
            "label": str(row.irrep_label),
            "k_label": str(row.k_label),
            "row_id": None,
            "source_kparam": tuple(int(value) for value in row.source_kparam),
            "reciprocal_vector_pml": tuple(row.reciprocal_vector_pml),
            "request_k_params": _pml_vector_to_case_k_params(
                decoder,
                gid=int(row.gid),
                reciprocal_vector_pml=tuple(row.reciprocal_vector_pml),
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
        for row in rows
    ]



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
        preferred_kparams={str(selected_irrep.get("k_label") or _k_label_from_irrep_label(str(selected_irrep.get("ordinary_symbol") or ""))): source_kparam},
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
    alias_candidates = stored_occurrence_alias_candidates(
        decoder,
        raw_rows,
        aliases_for_class=lambda row, requested: dynamic_magnetic_subduction_occurrence_alias_rows(
            decoder,
            magnetic_source,
            sg=int(sg),
            basis=basis,
            operations=operations,
            row=row,
            requested_source_kparam=requested,
        ),
        heading_identity_for_row=_occurrence_heading_identity,
        represented_heading_identities=tuple(
            identity
            for row in rows
            if (identity := _occurrence_heading_identity(row)) is not None
        ),
    )

    def row_spec(
        row: Any,
        *,
        alias: Any | None = None,
        observation_only: bool = False,
    ) -> dict[str, Any]:
        reciprocal_vector_pml = tuple(row.reciprocal_vector_pml)
        direction_matrix = [list(values) for values in row.direction_matrix]
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
            "old_id": int(decoder.little_record_by_gid(row.gid).old_id),
            "gid": int(row.gid),
            "label": str(row.irrep_label),
            "display_label": "m" + str(row.irrep_label),
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
            "primary": _selected_source_occurrence_matches(
                gid=int(row.gid),
                source_kparam=tuple(int(value) for value in row.source_kparam),
                selected_gid=selected_gid,
                selected_source_kparam=source_kparam,
            ),
        }
        if alias is not None:
            spec["primary"] = False
            spec["_occurrence_alias_anchor"] = (
                int(alias.representative_gid),
                tuple(int(value) for value in alias.representative_source_kparam),
            )
            spec["_magnetic_occurrence_alias_anchor"] = spec[
                "_occurrence_alias_anchor"
            ]
            spec["_occurrence_alias_candidate"] = tuple(
                int(value) for value in alias.candidate_source_kparam
            )
            spec["_occurrence_alias_heading_identity"] = alias.heading_identity
            spec["_magnetic_occurrence_alias_candidate"] = spec[
                "_occurrence_alias_candidate"
            ]
            spec["_magnetic_occurrence_alias_heading_identity"] = alias.heading_identity
        if observation_only:
            spec["primary"] = False
            spec["_occurrence_alias_observation_only"] = True
            spec["_magnetic_occurrence_alias_observation_only"] = True
        return spec

    specs = [row_spec(row) for row in rows]
    represented_exact = {
        (int(row.gid), tuple(int(value) for value in row.source_kparam))
        for row in rows
    }
    required_anchor_keys = {
        (
            int(candidate.representative_gid),
            tuple(int(value) for value in candidate.representative_source_kparam),
        )
        for candidate in alias_candidates
    }
    specs.extend(
        row_spec(row, observation_only=True)
        for row in raw_rows
        if (
            int(row.gid),
            tuple(int(value) for value in row.source_kparam),
        )
        in required_anchor_keys - represented_exact
    )
    specs.extend(
        row_spec(candidate.candidate_row, alias=candidate)
        for candidate in alias_candidates
    )
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
    slot_source_rows = iso.get("slot_source_numeric_rows") or []
    slot_domain_indices = iso.get("direction_domain_indices") or []
    if len(slot_source_rows) != len(selected_slots) or len(
        slot_domain_indices
    ) != len(selected_slots):
        raise ValueError("coupled OPD does not retain one Source domain per slot")
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
        selected_source_rows = None
        if int(slot_domain_indices[index]) == 0:
            selected_source_rows = [
                [float(value) for value in row] for row in slot_source_rows[index]
            ]
            full_dim = int(decoder.little_record_by_gid(int(gid)).full_dim)
            if not selected_source_rows or any(
                len(row) != full_dim for row in selected_source_rows
            ):
                raise ValueError(
                    f"selected coupled slot {index + 1} Source OPD width differs from its irrep"
                )
        primary_slots.append(
            {
                "order": index,
                "mode_kind": mode_kind,
                "gid": gid,
                "k_label": str((slot.get("kpoint") or {}).get("label") or ""),
                "source_kparam": source_kparam,
                "case_k_params": _k_params(slot.get("k_params") or {}),
                "source_numeric_rows": selected_source_rows,
            }
        )
    specs: list[tuple[int, dict[str, Any], str, int, str]] = []
    source_order = 0

    def spec_for_row(row: Any, mode_kind: str) -> dict[str, Any]:
            gid = int(row.gid)
            primary_slot = next(
                (
                    slot
                    for slot in primary_slots
                    if slot["mode_kind"] == mode_kind
                    and int(slot["gid"]) == gid
                    and _same_source_kparam(slot["source_kparam"], row.source_kparam)
                ),
                None,
            )
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
            if primary_slot is not None and primary_slot["source_numeric_rows"] is not None:
                selected_rows = primary_slot["source_numeric_rows"]
                source_free_count = len(selected_rows)
                direction_matrix = [
                    [selected_rows[free][coordinate] for free in range(source_free_count)]
                    for coordinate in range(len(selected_rows[0]))
                ]
            return {
                "old_id": old_id,
                "little_type": little_type,
                "gid": gid,
                "label": str(row.irrep_label),
                "display_label": ("m" if mode_kind == "mag" else "") + str(row.irrep_label),
                "k_label": str(row.k_label),
                "row_id": None,
                "source_kparam": tuple(int(value) for value in row.source_kparam),
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
                "source_occurrences": tuple(row.source_occurrences),
            }

    def append_specs(row_specs: Any, mode_kind: str) -> None:
        nonlocal source_order
        for spec in row_specs:
            primary_slot = next(
                (
                    slot
                    for slot in primary_slots
                    if slot["mode_kind"] == mode_kind
                    and int(slot["gid"]) == int(spec["gid"])
                    and _same_source_kparam(
                        slot["source_kparam"], spec["source_kparam"]
                    )
                ),
                None,
            )
            spec["primary"] = bool(spec.get("primary") and primary_slot is not None)
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
        )
        ordinary_rows = _order_promoted_dynamic_families_for_presentation(
            decoder,
            ordinary_rows,
            promoted_gids=promoted_gids,
        )
        ordinary_rows = _replay_distinct_type3_occurrence_aliases(
            decoder,
            ordinary_rows,
            allowed_k_labels={str(slot["k_label"]) for slot in primary_slots},
            aliases_for_row=lambda row: dynamic_subduction_occurrence_alias_rows(
                decoder,
                sg=int(sg),
                basis=basis,
                operations=ordinary_operations,
                row=row,
                irrep_source=ordinary_source,
            ),
        )
        ordinary_specs = [spec_for_row(row, "dsp") for row in ordinary_rows]
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
        magnetic_rows = _replay_distinct_type3_occurrence_aliases(
            decoder,
            magnetic_rows,
            allowed_k_labels={str(slot["k_label"]) for slot in primary_slots},
            aliases_for_row=lambda row: dynamic_magnetic_subduction_occurrence_alias_rows(
                decoder,
                magnetic_source,
                sg=int(sg),
                basis=basis,
                operations=operations,
                row=row,
            ),
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
