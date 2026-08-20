"""Group-invariance contradictions in displayed complete-mode vector fields.

Each printed coordinate or vector component denotes its complete rounding
cell, not an exact hidden value.  A definition is refuted only when no
one-to-one row assignment can satisfy one of the declared child-group
operations within those cells.  Passing this necessary condition does not
prove that the unpublished vectors are invariant.
"""

from __future__ import annotations

from collections import Counter, deque
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from fractions import Fraction
from functools import lru_cache
import gzip
import hashlib
import json
from pathlib import Path
import re
import struct
from typing import Any

import numpy as np

from Verification.mathematics.mode_basis import (
    POSITION_SCALE,
    VECTOR_SCALE,
    QuantizedModeDefinition,
    QuantizedStrainDefinition,
    QuantizedStructureOrbit,
)


GROUP_INVARIANCE_THEOREM = "complete_mode_field_group_invariance.v1"
_DATA_SCHEMA = "isodistort.validation.magnetic-group-generators.v3"
_DATA_PATH = Path(__file__).with_name("data") / "magnetic_group_generators.json.gz"
_BNS_LABEL = re.compile(r"(?:[1-9]|[1-9]\d|1\d\d|2[0-2]\d|230)\.[1-9]\d*")
_ORDINARY_LABEL = re.compile(r"(?:[1-9]|[1-9]\d|1\d\d|2[0-2]\d|230)")
_FRACTION = re.compile(r"-?(?:0|[1-9]\d*)(?:/[1-9]\d*)?")
_WEB_SUBGROUP_HEADING = re.compile(r"(?m)^Subgroup details\s*$")
_WEB_SUBGROUP_LINE = re.compile(r"^(?P<label>\d+(?:\.\d+)?)\s+\S.*$")

Matrix = tuple[tuple[int, int, int], tuple[int, int, int], tuple[int, int, int]]
Position = tuple[int, int, int]
Vector = tuple[int, int, int]


@dataclass(frozen=True)
class GroupOperation:
    rotation: Matrix
    translation: tuple[Fraction, Fraction, Fraction]
    time_reversal: int


_IDENTITY = GroupOperation(
    ((1, 0, 0), (0, 1, 0), (0, 0, 1)),
    (Fraction(0), Fraction(0), Fraction(0)),
    1,
)


@dataclass(frozen=True)
class _CarriedRow:
    atom: str
    position: Position
    vector: Vector


@dataclass(frozen=True)
class _GroupAuthority:
    requested_label: str
    bns_label: str
    uni_number: int
    hall_symbol: str
    operation_count: int
    operation_sha256: str
    generators: tuple[GroupOperation, ...]
    nonidentity_operations: tuple[GroupOperation, ...]
    setting_coverage_complete: bool
    settings: tuple["_GroupSettingAuthority", ...]


@dataclass(frozen=True)
class _GroupSettingAuthority:
    setting_key: str
    authority: str
    hall_number: int
    source_setting_ids: tuple[int, ...]
    choice: str
    operation_count: int
    operation_sha256: str
    generators: tuple[GroupOperation, ...]
    nonidentity_operations: tuple[GroupOperation, ...]


@dataclass(frozen=True)
class _PositionGraph:
    edges: tuple[tuple[int, ...], ...]
    bijection: tuple[int, ...] | None


def _strict_object(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
    result: dict[str, Any] = {}
    for key, value in pairs:
        if key in result:
            raise ValueError(f"group authority contains duplicate key {key!r}")
        result[key] = value
    return result


@lru_cache(maxsize=1)
def _load_data() -> Mapping[str, Any]:
    try:
        with gzip.open(_DATA_PATH, "rt", encoding="utf-8") as handle:
            value = json.load(
                handle,
                object_pairs_hook=_strict_object,
                parse_constant=lambda token: (_ for _ in ()).throw(
                    ValueError(f"group authority contains nonfinite {token}")
                ),
            )
    except (EOFError, OSError, json.JSONDecodeError) as exc:
        raise ValueError("invalid magnetic group authority artifact") from exc
    if not isinstance(value, Mapping) or value.get("schema") != _DATA_SCHEMA:
        raise ValueError("invalid magnetic group authority schema")
    if (
        value.get("group_count") != 1651
        or value.get("total_operation_count") != 38_307
        or value.get("setting_count") != 4_551
        or value.get("setting_operation_count") != 78_294
    ):
        raise ValueError("magnetic group authority has incomplete coverage")
    ordinary = value.get("ordinary_groups")
    entries = value.get("entries")
    if (
        not isinstance(ordinary, Mapping)
        or set(ordinary) != {str(number) for number in range(1, 231)}
        or not isinstance(entries, Mapping)
        or len(entries) != 1651
    ):
        raise ValueError("magnetic group authority has an invalid index")
    return value


def _exact_int(value: object, *, field: str) -> int:
    if isinstance(value, bool) or not isinstance(value, int):
        raise ValueError(f"{field} must be an integer")
    return value


def _determinant(matrix: Matrix) -> int:
    first, second, third = matrix
    return (
        first[0] * (second[1] * third[2] - second[2] * third[1])
        - first[1] * (second[0] * third[2] - second[2] * third[0])
        + first[2] * (second[0] * third[1] - second[1] * third[0])
    )


def _operation(value: object, *, label: str, index: int) -> GroupOperation:
    if not isinstance(value, Mapping):
        raise ValueError(f"group {label} generator {index} must be an object")
    if set(value) != {"rotation", "translation", "time_reversal"}:
        raise ValueError(f"group {label} generator {index} has invalid fields")
    raw_rotation = value["rotation"]
    if not isinstance(raw_rotation, list) or len(raw_rotation) != 3:
        raise ValueError(f"group {label} generator {index} needs matrix3 rotation")
    rotation = tuple(
        tuple(_exact_int(component, field="group rotation") for component in row)
        if isinstance(row, list) and len(row) == 3
        else ()
        for row in raw_rotation
    )
    if any(len(row) != 3 for row in rotation):
        raise ValueError(f"group {label} generator {index} needs matrix3 rotation")
    checked_rotation: Matrix = rotation  # type: ignore[assignment]
    if _determinant(checked_rotation) not in {-1, 1}:
        raise ValueError(f"group {label} generator {index} rotation must be unimodular")
    raw_translation = value["translation"]
    if not isinstance(raw_translation, list) or len(raw_translation) != 3:
        raise ValueError(f"group {label} generator {index} needs vector3 translation")
    translation: list[Fraction] = []
    for token in raw_translation:
        if not isinstance(token, str) or _FRACTION.fullmatch(token) is None:
            raise ValueError(f"group {label} generator {index} has invalid translation")
        component = Fraction(token)
        if str(component) != token or not 0 <= component < 1:
            raise ValueError(
                f"group {label} generator {index} translation is not canonical"
            )
        translation.append(component)
    time_reversal = _exact_int(value["time_reversal"], field="group time-reversal sign")
    if time_reversal not in {-1, 1}:
        raise ValueError(f"group {label} generator {index} has invalid time reversal")
    return GroupOperation(
        checked_rotation,
        tuple(translation),  # type: ignore[arg-type]
        time_reversal,
    )


def _matrix_multiply(left: Matrix, right: Matrix) -> Matrix:
    return tuple(
        tuple(
            sum(left[row][axis] * right[axis][column] for axis in range(3))
            for column in range(3)
        )
        for row in range(3)
    )  # type: ignore[return-value]


def _compose(left: GroupOperation, right: GroupOperation) -> GroupOperation:
    translation = tuple(
        (
            left.translation[row]
            + sum(
                Fraction(left.rotation[row][axis]) * right.translation[axis]
                for axis in range(3)
            )
        )
        % 1
        for row in range(3)
    )
    return GroupOperation(
        _matrix_multiply(left.rotation, right.rotation),
        translation,  # type: ignore[arg-type]
        left.time_reversal * right.time_reversal,
    )


def _closure(generators: Sequence[GroupOperation]) -> frozenset[GroupOperation]:
    seen = {_IDENTITY}
    pending = deque([_IDENTITY])
    while pending:
        current = pending.popleft()
        for generator in generators:
            product = _compose(current, generator)
            if product in seen:
                continue
            seen.add(product)
            pending.append(product)
            if len(seen) > 512:
                raise ValueError("group generator closure exceeds the declared bound")
    return frozenset(seen)


def _operation_json(operation: GroupOperation) -> dict[str, Any]:
    return {
        "rotation": [list(row) for row in operation.rotation],
        "translation": [str(component) for component in operation.translation],
        "time_reversal": operation.time_reversal,
    }


def _operation_key(operation: GroupOperation):
    return operation.rotation, operation.translation, operation.time_reversal


def _operation_digest(operations: Sequence[GroupOperation]) -> str:
    rows = [
        _operation_json(operation)
        for operation in sorted(operations, key=_operation_key)
    ]
    encoded = json.dumps(rows, sort_keys=True, separators=(",", ":")).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def _setting_authority(
    value: object,
    *,
    label: str,
    index: int,
) -> _GroupSettingAuthority:
    if not isinstance(value, Mapping) or set(value) != {
        "setting_key",
        "authority",
        "hall_number",
        "source_setting_ids",
        "choice",
        "operation_count",
        "operation_sha256",
        "generators",
    }:
        raise ValueError(f"group authority setting {label}/{index} has invalid fields")
    setting_key = value["setting_key"]
    authority = value["authority"]
    if (
        not isinstance(setting_key, str)
        or re.fullmatch(r"(?:hall|source):[1-9]\d*", setting_key) is None
        or authority not in {"hall", "source"}
        or not setting_key.startswith(authority + ":")
    ):
        raise ValueError(
            f"group authority setting {label}/{index} has invalid identity"
        )
    hall_number = _exact_int(value["hall_number"], field="setting Hall number")
    if not 1 <= hall_number <= 530:
        raise ValueError(
            f"group authority setting {label}/{index} has invalid Hall number"
        )
    raw_source_ids = value["source_setting_ids"]
    if (
        not isinstance(raw_source_ids, list)
        or any(
            isinstance(item, bool) or not isinstance(item, int) or not 1 <= item <= 800
            for item in raw_source_ids
        )
        or raw_source_ids != sorted(set(raw_source_ids))
        or (authority == "hall") != (not raw_source_ids)
    ):
        raise ValueError(
            f"group authority setting {label}/{index} has invalid Source settings"
        )
    choice = value["choice"]
    if not isinstance(choice, str) or choice != choice.strip():
        raise ValueError(f"group authority setting {label}/{index} has invalid choice")
    operation_count = _exact_int(
        value["operation_count"], field="setting operation count"
    )
    if not 1 <= operation_count <= 384:
        raise ValueError(
            f"group authority setting {label}/{index} has invalid operation count"
        )
    operation_sha256 = value["operation_sha256"]
    if (
        not isinstance(operation_sha256, str)
        or re.fullmatch(r"[0-9a-f]{64}", operation_sha256) is None
    ):
        raise ValueError(
            f"group authority setting {label}/{index} has invalid operation digest"
        )
    raw_generators = value["generators"]
    if not isinstance(raw_generators, list) or len(raw_generators) > 12:
        raise ValueError(
            f"group authority setting {label}/{index} has invalid generators"
        )
    generators = tuple(
        _operation(row, label=f"{label}/Hall{hall_number}", index=generator_index)
        for generator_index, row in enumerate(raw_generators)
    )
    operations = _closure(generators)
    if (
        len(operations) != operation_count
        or _operation_digest(operations) != operation_sha256
    ):
        raise ValueError(
            f"group authority setting {label}/{index} fails closure identity"
        )
    return _GroupSettingAuthority(
        setting_key=setting_key,
        authority=authority,
        hall_number=hall_number,
        source_setting_ids=tuple(raw_source_ids),
        choice=choice,
        operation_count=operation_count,
        operation_sha256=operation_sha256,
        generators=generators,
        nonidentity_operations=tuple(
            operation
            for operation in sorted(operations, key=_operation_key)
            if operation != _IDENTITY
        ),
    )


@lru_cache(maxsize=None)
def group_authority(label: str) -> _GroupAuthority:
    """Resolve and revalidate the complete generator closure for one child group."""

    if not isinstance(label, str):
        raise ValueError("subgroup display label must be text")
    data = _load_data()
    ordinary = data["ordinary_groups"]
    entries = data["entries"]
    assert isinstance(ordinary, Mapping)
    assert isinstance(entries, Mapping)
    if _ORDINARY_LABEL.fullmatch(label) is not None:
        bns_label = ordinary.get(label)
    elif _BNS_LABEL.fullmatch(label) is not None:
        bns_label = label
    else:
        raise ValueError(f"invalid subgroup display label {label!r}")
    if not isinstance(bns_label, str) or _BNS_LABEL.fullmatch(bns_label) is None:
        raise ValueError(f"group authority has no exact BNS mapping for {label!r}")
    raw = entries.get(bns_label)
    if not isinstance(raw, Mapping):
        raise ValueError(f"group authority has no entry for {bns_label}")
    if set(raw) != {
        "uni_number",
        "hall_number",
        "hall_symbol",
        "operation_count",
        "operation_sha256",
        "generators",
        "setting_coverage_complete",
        "settings",
    }:
        raise ValueError(f"group authority entry {bns_label} has invalid fields")
    uni_number = _exact_int(raw["uni_number"], field="UNI number")
    if not 1 <= uni_number <= 1651:
        raise ValueError(f"group authority entry {bns_label} has invalid UNI number")
    hall_number = _exact_int(raw["hall_number"], field="Hall number")
    if not 1 <= hall_number <= 530:
        raise ValueError(f"group authority entry {bns_label} has invalid Hall number")
    hall_symbol = raw["hall_symbol"]
    if (
        not isinstance(hall_symbol, str)
        or not hall_symbol
        or hall_symbol != hall_symbol.strip()
    ):
        raise ValueError(f"group authority entry {bns_label} has invalid Hall symbol")
    operation_count = _exact_int(raw["operation_count"], field="group operation count")
    if not 1 <= operation_count <= 384:
        raise ValueError(
            f"group authority entry {bns_label} has invalid operation count"
        )
    operation_sha256 = raw["operation_sha256"]
    if (
        not isinstance(operation_sha256, str)
        or re.fullmatch(r"[0-9a-f]{64}", operation_sha256) is None
    ):
        raise ValueError(
            f"group authority entry {bns_label} has invalid operation digest"
        )
    raw_generators = raw["generators"]
    if not isinstance(raw_generators, list) or len(raw_generators) > 12:
        raise ValueError(f"group authority entry {bns_label} has invalid generators")
    generators = tuple(
        _operation(value, label=bns_label, index=index)
        for index, value in enumerate(raw_generators)
    )
    closure = _closure(generators)
    if (
        len(closure) != operation_count
        or _operation_digest(closure) != operation_sha256
    ):
        raise ValueError(f"group authority entry {bns_label} fails closure identity")
    nonidentity_operations = tuple(
        operation
        for operation in sorted(closure, key=_operation_key)
        if operation != _IDENTITY
    )
    setting_coverage_complete = raw["setting_coverage_complete"]
    if not isinstance(setting_coverage_complete, bool):
        raise ValueError(
            f"group authority entry {bns_label} has invalid setting coverage"
        )
    raw_settings = raw["settings"]
    if not isinstance(raw_settings, list) or not raw_settings:
        raise ValueError(
            f"group authority entry {bns_label} has no coordinate settings"
        )
    settings = tuple(
        _setting_authority(value, label=bns_label, index=index)
        for index, value in enumerate(raw_settings)
    )
    setting_keys = tuple(setting.setting_key for setting in settings)
    hall_numbers = tuple(
        setting.hall_number for setting in settings if setting.authority == "hall"
    )
    if len(set(setting_keys)) != len(setting_keys) or tuple(sorted(hall_numbers)) != hall_numbers:
        raise ValueError(f"group authority entry {bns_label} has invalid setting order")
    standard = next(
        (
            setting
            for setting in settings
            if setting.setting_key == f"hall:{hall_number}"
        ),
        None,
    )
    if standard is None or (
        standard.operation_count != operation_count
        or standard.operation_sha256 != operation_sha256
    ):
        raise ValueError(
            f"group authority entry {bns_label} omits its standard setting"
        )
    incomplete = set(data.get("incomplete_setting_groups") or ())
    if (bns_label in incomplete) == setting_coverage_complete:
        raise ValueError(
            f"group authority entry {bns_label} has inconsistent setting coverage"
        )
    return _GroupAuthority(
        requested_label=label,
        bns_label=bns_label,
        uni_number=uni_number,
        hall_symbol=hall_symbol,
        operation_count=operation_count,
        operation_sha256=operation_sha256,
        generators=generators,
        nonidentity_operations=nonidentity_operations,
        setting_coverage_complete=setting_coverage_complete,
        settings=settings,
    )


def subgroup_label_from_web_text(text: str) -> str:
    if not isinstance(text, str):
        raise ValueError("Web complete-mode payload must be text")
    headings = list(_WEB_SUBGROUP_HEADING.finditer(text))
    if len(headings) != 1:
        raise ValueError(
            "Web complete-mode text must contain one Subgroup details heading"
        )
    tail = text[headings[0].end() :]
    line = next((row.strip() for row in tail.splitlines() if row.strip()), None)
    match = _WEB_SUBGROUP_LINE.fullmatch(line or "")
    if match is None:
        raise ValueError("Web complete-mode text has no exact subgroup display label")
    label = match.group("label")
    if _ORDINARY_LABEL.fullmatch(label) is None and _BNS_LABEL.fullmatch(label) is None:
        raise ValueError(
            f"Web subgroup display label is outside the supported range: {label!r}"
        )
    return label


def _mapping(value: object, *, field: str) -> Mapping[str, Any]:
    if not isinstance(value, Mapping):
        raise ValueError(f"{field} must be an object")
    return value


def subgroup_label_from_local_payload(payload: Mapping[str, Any]) -> str:
    root = _mapping(payload, field="Local payload")
    preview = _mapping(root.get("preview"), field="Local preview")
    selected = _mapping(preview.get("selected"), field="Local selected state")
    details = _mapping(
        selected.get("mode_details"), field="Local complete-mode details"
    )
    subgroup = _mapping(details.get("subgroup_details"), field="Local subgroup details")
    magnetic = _mapping(
        subgroup.get("magnetic_subgroup"), field="Local magnetic subgroup"
    )
    display_label = magnetic.get("display_label")
    number = magnetic.get("number")
    if display_label is not None:
        if not isinstance(display_label, str) or (
            _BNS_LABEL.fullmatch(display_label) is None
            and _ORDINARY_LABEL.fullmatch(display_label) is None
        ):
            raise ValueError("Local subgroup display label is invalid")
        checked_number = _exact_int(number, field="Local subgroup number")
        authority = group_authority(display_label)
        expected_number = (
            authority.uni_number if "." in display_label else int(display_label)
        )
        if checked_number != expected_number:
            raise ValueError("Local subgroup number and display label disagree")
        return display_label
    checked_number = _exact_int(number, field="Local subgroup number")
    if not 1 <= checked_number <= 230:
        raise ValueError("Local subgroup number is outside 1..230")
    return str(checked_number)


def _carried_rows(definition: QuantizedModeDefinition) -> tuple[_CarriedRow, ...]:
    rows: list[_CarriedRow] = []
    current_atom: str | None = None
    keys: set[tuple[str, int, int, int]] = set()
    for row in definition.rows:
        if row.atom_marker is not None:
            current_atom = row.atom_marker
        if current_atom is None:
            raise ValueError(
                f"mode definition starts with an unlabeled atom row: {definition.label!r}"
            )
        key = (current_atom, *row.position)
        if key in keys:
            raise ValueError(
                f"mode definition repeats an atom-position row: {definition.label!r}"
            )
        keys.add(key)
        rows.append(_CarriedRow(current_atom, row.position, row.vector))
    return tuple(rows)


def _ceil(value: Fraction) -> int:
    return -((-value.numerator) // value.denominator)


def _position_candidates(
    rows: Sequence[_CarriedRow], operation: GroupOperation
) -> _PositionGraph:
    by_key = {(row.atom, *row.position): index for index, row in enumerate(rows)}
    result = []
    for row in rows:
        axes = []
        for target_axis in range(3):
            center = (
                sum(
                    operation.rotation[target_axis][source_axis]
                    * row.position[source_axis]
                    for source_axis in range(3)
                )
                + POSITION_SCALE * operation.translation[target_axis]
            )
            bound = Fraction(
                1 + sum(abs(value) for value in operation.rotation[target_axis]), 2
            )
            upper = center + bound
            axes.append(
                tuple(
                    sorted(
                        {
                            integer % POSITION_SCALE
                            for integer in range(
                                _ceil(center - bound),
                                upper.numerator // upper.denominator + 1,
                            )
                        }
                    )
                )
            )
        candidates = tuple(
            index
            for x in axes[0]
            for y in axes[1]
            for z in axes[2]
            if (index := by_key.get((row.atom, x, y, z))) is not None
        )
        result.append(candidates)
    edges = tuple(result)
    targets = tuple(candidates[0] for candidates in edges if len(candidates) == 1)
    bijection = (
        targets
        if len(targets) == len(edges) and len(set(targets)) == len(edges)
        else None
    )
    return _PositionGraph(edges, bijection)


def _vector_action(operation: GroupOperation, kind: str) -> Matrix:
    if kind == "dsp":
        factor = 1
    elif kind == "mag":
        factor = operation.time_reversal * _determinant(operation.rotation)
    else:
        raise ValueError(f"unsupported mode kind {kind!r}")
    return tuple(tuple(factor * value for value in row) for row in operation.rotation)  # type: ignore[return-value]


def _maximum_matching(edges: Sequence[Sequence[int]]) -> tuple[int, ...]:
    source_target = [-1] * len(edges)
    target_source = [-1] * len(edges)
    for root in sorted(
        range(len(edges)), key=lambda source: (len(edges[source]), source)
    ):
        if source_target[root] >= 0:
            continue
        pending = deque([root])
        seen_sources = {root}
        seen_targets: set[int] = set()
        predecessor: dict[int, tuple[int, int]] = {}
        endpoint: tuple[int, int] | None = None
        while pending and endpoint is None:
            source = pending.popleft()
            for target in edges[source]:
                if not 0 <= target < len(edges):
                    raise ValueError(
                        "candidate matching target is outside the row domain"
                    )
                if target in seen_targets:
                    continue
                seen_targets.add(target)
                owner = target_source[target]
                if owner < 0:
                    endpoint = source, target
                    break
                if owner not in seen_sources:
                    seen_sources.add(owner)
                    predecessor[owner] = source, target
                    pending.append(owner)
        if endpoint is None:
            continue
        source, target = endpoint
        while True:
            source_target[source] = target
            target_source[target] = source
            if source == root:
                break
            source, target = predecessor[source]
    return tuple(source_target)


def _graph_digest(edges: Sequence[Sequence[int]]) -> str:
    digest = hashlib.sha256()
    digest.update(b"isodistort.group-invariance.candidate-graph.v1\0")
    digest.update(struct.pack("<I", len(edges)))
    for candidates in edges:
        digest.update(struct.pack("<I", len(candidates)))
        for target in candidates:
            digest.update(struct.pack("<I", target))
    return digest.hexdigest()


def _int64_vectors(
    rows: Sequence[_CarriedRow], operations: Sequence[GroupOperation]
) -> np.ndarray | None:
    max_component = max(abs(value) for row in rows for value in row.vector)
    max_row_norm = max(
        (
            sum(abs(value) for value in transform_row)
            for operation in operations
            for transform_row in _vector_action(operation, "mag")
        ),
        default=1,
    )
    safe_quarter = (2**63 - 1) // 4
    if max_component > safe_quarter // max(1, max_row_norm):
        return None
    return np.asarray([row.vector for row in rows], dtype=np.int64)


def _bijective_vector_edges(
    vectors: np.ndarray,
    targets: Sequence[int],
    transform: Matrix,
) -> tuple[tuple[int, ...], ...]:
    transform_array = np.asarray(transform, dtype=np.int64)
    target_array = np.fromiter(targets, dtype=np.intp, count=len(targets))
    expected = vectors @ transform_array.T
    observed = vectors[target_array]
    bounds_twice = np.asarray(
        [1 + sum(abs(value) for value in row) for row in transform],
        dtype=np.int64,
    )
    admitted = np.all(2 * np.abs(observed - expected) <= bounds_twice, axis=1)
    return tuple(
        (target,) if bool(keep) else ()
        for target, keep in zip(targets, admitted, strict=True)
    )


def _definition_certificate(
    definition: QuantizedModeDefinition,
    operations: Sequence[GroupOperation],
    position_cache: dict[tuple[tuple[tuple[str, Position], ...], int], _PositionGraph],
) -> tuple[dict[str, Any], bytes]:
    rows = _carried_rows(definition)
    pattern = tuple((row.atom, row.position) for row in rows)
    digest = hashlib.sha256()
    digest.update(b"isodistort.group-invariance.definition.v1\0")
    digest.update(struct.pack("<I", definition.ordinal))
    digest.update(definition.kind.encode("ascii") + b"\0")
    digest.update(definition.label.encode("utf-8") + b"\0")
    int64_vectors = _int64_vectors(rows, operations)
    for operation_index, operation in enumerate(operations):
        position_graph = position_cache.get((pattern, operation_index))
        if position_graph is None:
            position_graph = _position_candidates(rows, operation)
            position_cache[(pattern, operation_index)] = position_graph
        transform = _vector_action(operation, definition.kind)
        if int64_vectors is not None and position_graph.bijection is not None:
            vector_edges = _bijective_vector_edges(
                int64_vectors, position_graph.bijection, transform
            )
            matching = tuple(targets[0] if targets else -1 for targets in vector_edges)
        else:
            vector_edges = []
            for source, candidates in enumerate(position_graph.edges):
                expected = tuple(
                    sum(
                        transform[target_axis][source_axis]
                        * rows[source].vector[source_axis]
                        for source_axis in range(3)
                    )
                    for target_axis in range(3)
                )
                bounds = tuple(
                    Fraction(
                        1 + sum(abs(value) for value in transform[target_axis]),
                        2,
                    )
                    for target_axis in range(3)
                )
                vector_edges.append(
                    tuple(
                        target
                        for target in candidates
                        if all(
                            abs(rows[target].vector[axis] - expected[axis])
                            <= bounds[axis]
                            for axis in range(3)
                        )
                    )
                )
            matching = _maximum_matching(vector_edges)
        graph_sha256 = _graph_digest(vector_edges)
        matching_size = sum(target >= 0 for target in matching)
        digest.update(struct.pack("<II", operation_index, matching_size))
        digest.update(bytes.fromhex(graph_sha256))
        for target in matching:
            digest.update(struct.pack("<i", target))
        if matching_size != len(rows):
            unmatched = matching.index(-1)
            return (
                {
                    "definition_index": definition.ordinal,
                    "kind": definition.kind,
                    "label": definition.label,
                    "row_count": len(rows),
                    "operation_index": operation_index,
                    "operation": _operation_json(operation),
                    "candidate_graph_sha256": graph_sha256,
                    "maximum_matching_size": matching_size,
                    "unmatched_source_row": unmatched,
                    "unmatched_source_atom": rows[unmatched].atom,
                    "unmatched_source_position": list(rows[unmatched].position),
                    "unmatched_source_vector": list(rows[unmatched].vector),
                    "unmatched_candidate_targets": list(vector_edges[unmatched]),
                },
                digest.digest(),
            )
    return ({}, digest.digest())


def _position_failure(
    definition: QuantizedModeDefinition,
    operations: Sequence[GroupOperation],
    position_cache: dict[tuple[tuple[tuple[str, Position], ...], int], _PositionGraph],
) -> tuple[dict[str, Any], bytes]:
    rows = _carried_rows(definition)
    pattern = tuple((row.atom, row.position) for row in rows)
    digest = hashlib.sha256()
    digest.update(b"isodistort.group-invariance.position-setting.v1\0")
    digest.update(struct.pack("<I", definition.ordinal))
    for operation_index, operation in enumerate(operations):
        position_graph = position_cache.get((pattern, operation_index))
        if position_graph is None:
            position_graph = _position_candidates(rows, operation)
            position_cache[(pattern, operation_index)] = position_graph
        matching = (
            position_graph.bijection
            if position_graph.bijection is not None
            else _maximum_matching(position_graph.edges)
        )
        graph_sha256 = _graph_digest(position_graph.edges)
        matching_size = sum(target >= 0 for target in matching)
        digest.update(struct.pack("<II", operation_index, matching_size))
        digest.update(bytes.fromhex(graph_sha256))
        for target in matching:
            digest.update(struct.pack("<i", target))
        if matching_size != len(rows):
            unmatched = matching.index(-1)
            return (
                {
                    "stage": "position",
                    "definition_index": definition.ordinal,
                    "kind": definition.kind,
                    "label": definition.label,
                    "row_count": len(rows),
                    "operation_index": operation_index,
                    "operation": _operation_json(operation),
                    "candidate_graph_sha256": graph_sha256,
                    "maximum_matching_size": matching_size,
                    "unmatched_source_row": unmatched,
                    "unmatched_source_atom": rows[unmatched].atom,
                    "unmatched_source_position": list(rows[unmatched].position),
                    "unmatched_source_vector": list(rows[unmatched].vector),
                    "unmatched_candidate_targets": list(
                        position_graph.edges[unmatched]
                    ),
                },
                digest.digest(),
            )
    return {}, digest.digest()


def _periodic_distance(left: Fraction, right: Fraction) -> Fraction:
    difference = (left - right) % POSITION_SCALE
    return min(difference, POSITION_SCALE - difference)


def _representative_orbit_count(
    orbit: QuantizedStructureOrbit,
    operations: Sequence[GroupOperation],
) -> int:
    images: list[
        tuple[
            tuple[Fraction, Fraction, Fraction],
            tuple[tuple[int, int, int], tuple[int, int, int], tuple[int, int, int]],
        ]
    ] = []
    for operation in (_IDENTITY, *operations):
        center = tuple(
            (
                sum(
                    operation.rotation[target][source]
                    * orbit.representative[source]
                    for source in range(3)
                )
                + POSITION_SCALE * operation.translation[target]
            )
            % POSITION_SCALE
            for target in range(3)
        )
        images.append((center, operation.rotation))
    parents = list(range(len(images)))

    def find(index: int) -> int:
        while parents[index] != index:
            parents[index] = parents[parents[index]]
            index = parents[index]
        return index

    def union(left: int, right: int) -> None:
        left_root = find(left)
        right_root = find(right)
        if left_root != right_root:
            parents[right_root] = left_root

    for left in range(len(images)):
        left_center, left_rotation = images[left]
        for right in range(left):
            right_center, right_rotation = images[right]
            if all(
                _periodic_distance(left_center[axis], right_center[axis])
                <= Fraction(
                    sum(
                        abs(
                            left_rotation[axis][source]
                            - right_rotation[axis][source]
                        )
                        for source in range(3)
                    ),
                    2,
                )
                for axis in range(3)
            ):
                union(left, right)
    return len({find(index) for index in range(len(images))})


def _structure_failure(
    structure: Sequence[QuantizedStructureOrbit],
    operations: Sequence[GroupOperation],
) -> dict[str, Any]:
    for orbit_index, orbit in enumerate(structure):
        if orbit.expanded_positions:
            rows = tuple(
                _CarriedRow(orbit.label, position, (0, 0, 0))
                for position in orbit.expanded_positions
            )
            for operation_index, operation in enumerate(operations):
                graph = _position_candidates(rows, operation)
                matching = (
                    graph.bijection
                    if graph.bijection is not None
                    else _maximum_matching(graph.edges)
                )
                if sum(target >= 0 for target in matching) != len(rows):
                    return {
                        "orbit_index": orbit_index,
                        "label": orbit.label,
                        "site": orbit.site,
                        "expected_multiplicity": orbit.multiplicity,
                        "observed_orbit_size": len(rows),
                        "operation_index": operation_index,
                    }
            observed = len(orbit.expanded_positions)
        else:
            observed = _representative_orbit_count(orbit, operations)
        if observed != orbit.multiplicity:
            return {
                "orbit_index": orbit_index,
                "label": orbit.label,
                "site": orbit.site,
                "expected_multiplicity": orbit.multiplicity,
                "observed_orbit_size": observed,
                "operation_index": None,
            }
    return {}


def assess_group_invariance(
    definitions: Sequence[QuantizedModeDefinition | QuantizedStrainDefinition],
    subgroup_label: str,
    structure: Sequence[QuantizedStructureOrbit],
) -> dict[str, Any]:
    """Assess the display-rounding necessary condition for one saved output."""

    if not isinstance(definitions, Sequence):
        raise ValueError("mode definitions must be an ordered sequence")
    authority = group_authority(subgroup_label)
    if not isinstance(structure, Sequence) or not structure:
        raise ValueError("undistorted structure must be a nonempty ordered sequence")
    for orbit in structure:
        if not isinstance(orbit, QuantizedStructureOrbit):
            raise ValueError("undistorted structure has an invalid orbit")
    atomic_definitions: list[QuantizedModeDefinition] = []
    strain_count = 0
    definition_digest = hashlib.sha256()
    definition_digest.update(b"isodistort.group-invariance.output-settings.v2\0")
    kind_counts: Counter[str] = Counter()
    row_count = 0
    for expected_ordinal, definition in enumerate(definitions):
        if isinstance(definition, QuantizedStrainDefinition):
            if definition.ordinal != expected_ordinal:
                raise ValueError("mode definition ordinals must be contiguous and ordered")
            strain_count += 1
            row_count += 1
            continue
        if not isinstance(definition, QuantizedModeDefinition):
            raise ValueError("mode definitions use an unsupported definition type")
        if definition.ordinal != expected_ordinal:
            raise ValueError("mode definition ordinals must be contiguous and ordered")
        if definition.kind not in {"dsp", "mag"}:
            raise ValueError(f"unsupported mode kind {definition.kind!r}")
        kind_counts[definition.kind] += 1
        row_count += len(definition.rows)
        atomic_definitions.append(definition)
    if not atomic_definitions:
        return {
            "theorem": GROUP_INVARIANCE_THEOREM,
            "status": "not_applicable",
            "claim": "atomic_mode_field_group_invariance",
            "reason": "no_atomic_mode_definitions",
            "subgroup_display_label": authority.requested_label,
            "definition_count": len(definitions),
            "definition_row_count": row_count,
            "strain_definition_count": strain_count,
            "atomic_definition_count": 0,
            "undistorted_orbit_count": len(structure),
        }
    compatible_settings: list[str] = []
    satisfying_settings: list[str] = []
    setting_refutations: list[dict[str, Any]] = []
    for setting in authority.settings:
        definition_digest.update(setting.setting_key.encode("ascii") + b"\0")
        definition_digest.update(bytes.fromhex(setting.operation_sha256))
        structure_failure = _structure_failure(
            structure, setting.nonidentity_operations
        )
        if structure_failure:
            setting_refutations.append(
                {
                    "setting_key": setting.setting_key,
                    "authority": setting.authority,
                    "hall_number": setting.hall_number,
                    "source_setting_ids": list(setting.source_setting_ids),
                    "operation_sha256": setting.operation_sha256,
                    "structure_compatible": False,
                    "position_compatible": False,
                    "structure_refutation": structure_failure,
                    "refutations": [],
                }
            )
            continue
        position_cache: dict[
            tuple[tuple[tuple[str, Position], ...], int], _PositionGraph
        ] = {}
        position_failure: dict[str, Any] = {}
        for definition in atomic_definitions:
            position_failure, certificate_digest = _position_failure(
                definition, setting.nonidentity_operations, position_cache
            )
            definition_digest.update(certificate_digest)
            if position_failure:
                break
        if position_failure:
            setting_refutations.append(
                {
                    "setting_key": setting.setting_key,
                    "authority": setting.authority,
                    "hall_number": setting.hall_number,
                    "source_setting_ids": list(setting.source_setting_ids),
                    "operation_sha256": setting.operation_sha256,
                    "structure_compatible": True,
                    "position_compatible": False,
                    "structure_refutation": None,
                    "refutations": [position_failure],
                }
            )
            continue
        compatible_settings.append(setting.setting_key)
        refutations = []
        for definition in atomic_definitions:
            refutation, certificate_digest = _definition_certificate(
                definition, setting.nonidentity_operations, position_cache
            )
            definition_digest.update(certificate_digest)
            if refutation:
                refutations.append({"stage": "vector", **refutation})
                break
        if not refutations:
            satisfying_settings.append(setting.setting_key)
            break
        else:
            setting_refutations.append(
                {
                    "setting_key": setting.setting_key,
                    "authority": setting.authority,
                    "hall_number": setting.hall_number,
                    "source_setting_ids": list(setting.source_setting_ids),
                    "operation_sha256": setting.operation_sha256,
                    "structure_compatible": True,
                    "position_compatible": True,
                    "structure_refutation": None,
                    "refutations": refutations,
                }
            )
    if satisfying_settings:
        status = "satisfied"
        persisted_refutations: list[dict[str, Any]] = []
    elif authority.setting_coverage_complete:
        status = "refuted"
        persisted_refutations = setting_refutations
    else:
        status = "indeterminate"
        persisted_refutations = setting_refutations
    refuted_definitions = {
        int(refutation["definition_index"])
        for setting in persisted_refutations
        for refutation in setting["refutations"]
        if setting["position_compatible"]
    }
    return {
        "theorem": GROUP_INVARIANCE_THEOREM,
        "status": status,
        "claim": "atomic_mode_field_group_invariance_over_declared_presentations",
        "subgroup_display_label": authority.requested_label,
        "subgroup_bns_label": authority.bns_label,
        "subgroup_uni_number": authority.uni_number,
        "hall_symbol": authority.hall_symbol,
        "operation_count": authority.operation_count,
        "operation_sha256": authority.operation_sha256,
        "generator_count": len(authority.generators),
        "nonidentity_operation_count": len(authority.nonidentity_operations),
        "setting_coverage_complete": authority.setting_coverage_complete,
        "setting_count": len(authority.settings),
        "hall_setting_count": sum(
            setting.authority == "hall" for setting in authority.settings
        ),
        "source_presentation_count": sum(
            setting.authority == "source" for setting in authority.settings
        ),
        "position_compatible_setting_count": len(compatible_settings),
        "satisfying_setting_count": len(satisfying_settings),
        "satisfying_setting_keys": satisfying_settings,
        "position_display_scale": POSITION_SCALE,
        "vector_display_scale": VECTOR_SCALE,
        "definition_count": len(definitions),
        "atomic_definition_count": len(atomic_definitions),
        "strain_definition_count": strain_count,
        "definition_row_count": row_count,
        "displacive_definition_count": kind_counts["dsp"],
        "magnetic_definition_count": kind_counts["mag"],
        "site_pattern_count": len(
            {
                tuple((row.atom, row.position) for row in _carried_rows(definition))
                for definition in atomic_definitions
            }
        ),
        "undistorted_orbit_count": len(structure),
        "refuted_definition_count": len(refuted_definitions),
        "definition_certificates_sha256": definition_digest.hexdigest(),
        "setting_refutations": persisted_refutations,
    }


def _nonnegative_int(value: object, *, field: str) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or value < 0:
        raise ValueError(f"{field} must be a nonnegative integer")
    return value


def _sha256(value: object, *, field: str) -> str:
    if not isinstance(value, str) or re.fullmatch(r"[0-9a-f]{64}", value) is None:
        raise ValueError(f"{field} must be one lowercase SHA256 digest")
    return value


def _validate_legacy_group_invariance_certificate(row: Mapping[str, Any]) -> None:
    """Validate a verifier-v2 single-setting certificate."""

    if row.get("theorem") != GROUP_INVARIANCE_THEOREM:
        raise ValueError("invalid group-invariance theorem")
    status = row.get("status")
    if status not in {"satisfied", "refuted"}:
        raise ValueError("invalid group-invariance status")
    if row.get("claim") != "display_rounding_necessary_condition":
        raise ValueError("invalid group-invariance claim")
    display_label = row.get("subgroup_display_label")
    if not isinstance(display_label, str):
        raise ValueError("invalid subgroup display label")
    authority = group_authority(display_label)
    if (
        row.get("subgroup_bns_label") != authority.bns_label
        or row.get("subgroup_uni_number") != authority.uni_number
        or row.get("hall_symbol") != authority.hall_symbol
        or row.get("operation_count") != authority.operation_count
        or row.get("operation_sha256") != authority.operation_sha256
        or row.get("generator_count") != len(authority.generators)
        or row.get("nonidentity_operation_count")
        != len(authority.nonidentity_operations)
    ):
        raise ValueError("invalid subgroup authority identity")
    if (
        row.get("position_display_scale") != POSITION_SCALE
        or row.get("vector_display_scale") != VECTOR_SCALE
    ):
        raise ValueError("invalid group-invariance display scale")
    definition_count = _nonnegative_int(
        row.get("definition_count"), field="group-invariance definition count"
    )
    _nonnegative_int(
        row.get("definition_row_count"), field="group-invariance definition-row count"
    )
    displacive_count = _nonnegative_int(
        row.get("displacive_definition_count"),
        field="group-invariance displacive definition count",
    )
    magnetic_count = _nonnegative_int(
        row.get("magnetic_definition_count"),
        field="group-invariance magnetic definition count",
    )
    site_patterns = _nonnegative_int(
        row.get("site_pattern_count"), field="group-invariance site-pattern count"
    )
    refuted_count = _nonnegative_int(
        row.get("refuted_definition_count"),
        field="group-invariance refuted definition count",
    )
    if (
        displacive_count + magnetic_count != definition_count
        or site_patterns > definition_count
        or refuted_count > definition_count
    ):
        raise ValueError("invalid group-invariance totals")
    _sha256(
        row.get("definition_certificates_sha256"),
        field="group-invariance definition certificate SHA256",
    )
    refutations = row.get("refutations")
    if not isinstance(refutations, list) or len(refutations) != refuted_count:
        raise ValueError("invalid group-invariance refutations")
    prior_index = -1
    for refutation in refutations:
        if not isinstance(refutation, Mapping):
            raise ValueError("invalid group-invariance refutation")
        definition_index = _nonnegative_int(
            refutation.get("definition_index"),
            field="group-invariance refuted definition index",
        )
        row_count = _nonnegative_int(
            refutation.get("row_count"), field="group-invariance row count"
        )
        operation_index = _nonnegative_int(
            refutation.get("operation_index"),
            field="group-invariance operation index",
        )
        matching_size = _nonnegative_int(
            refutation.get("maximum_matching_size"),
            field="group-invariance matching size",
        )
        unmatched = _nonnegative_int(
            refutation.get("unmatched_source_row"),
            field="group-invariance unmatched source row",
        )
        if (
            definition_index <= prior_index
            or definition_index >= definition_count
            or row_count <= 0
            or operation_index >= len(authority.nonidentity_operations)
            or matching_size >= row_count
            or unmatched >= row_count
        ):
            raise ValueError("invalid group-invariance witness")
        prior_index = definition_index
        if refutation.get("kind") not in {"dsp", "mag"}:
            raise ValueError("invalid group-invariance kind")
        for key in ("label", "unmatched_source_atom"):
            value = refutation.get(key)
            if not isinstance(value, str) or not value:
                raise ValueError("invalid group-invariance text")
        for key in ("unmatched_source_position", "unmatched_source_vector"):
            vector = refutation.get(key)
            if (
                not isinstance(vector, list)
                or len(vector) != 3
                or any(
                    isinstance(value, bool) or not isinstance(value, int)
                    for value in vector
                )
            ):
                raise ValueError("invalid group-invariance vector")
        candidates = refutation.get("unmatched_candidate_targets")
        if (
            not isinstance(candidates, list)
            or any(
                isinstance(value, bool)
                or not isinstance(value, int)
                or not 0 <= value < row_count
                for value in candidates
            )
            or candidates != sorted(set(candidates))
        ):
            raise ValueError("invalid group-invariance candidates")
        _sha256(
            refutation.get("candidate_graph_sha256"),
            field="group-invariance candidate-graph SHA256",
        )
        expected_operation = authority.nonidentity_operations[operation_index]
        expected = _operation_json(expected_operation)
        if refutation.get("operation") != expected:
            raise ValueError("invalid group-invariance operation")
    if (status == "refuted") != bool(refutations):
        raise ValueError("invalid group-invariance status")


def _validate_setting_refutation(
    refutation: object,
    *,
    definition_count: int,
    setting: _GroupSettingAuthority,
) -> int:
    if not isinstance(refutation, Mapping):
        raise ValueError("invalid setting-aware group-invariance refutation")
    if refutation.get("stage") not in {"position", "vector"}:
        raise ValueError("invalid group-invariance refutation stage")
    definition_index = _nonnegative_int(
        refutation.get("definition_index"),
        field="group-invariance refuted definition index",
    )
    row_count = _nonnegative_int(
        refutation.get("row_count"), field="group-invariance row count"
    )
    operation_index = _nonnegative_int(
        refutation.get("operation_index"),
        field="group-invariance operation index",
    )
    matching_size = _nonnegative_int(
        refutation.get("maximum_matching_size"),
        field="group-invariance matching size",
    )
    unmatched = _nonnegative_int(
        refutation.get("unmatched_source_row"),
        field="group-invariance unmatched source row",
    )
    if (
        definition_index >= definition_count
        or row_count <= 0
        or operation_index >= len(setting.nonidentity_operations)
        or matching_size >= row_count
        or unmatched >= row_count
    ):
        raise ValueError("invalid setting-aware group-invariance witness")
    if refutation.get("kind") not in {"dsp", "mag"}:
        raise ValueError("invalid group-invariance kind")
    for key in ("label", "unmatched_source_atom"):
        value = refutation.get(key)
        if not isinstance(value, str) or not value:
            raise ValueError("invalid group-invariance text")
    for key in ("unmatched_source_position", "unmatched_source_vector"):
        vector = refutation.get(key)
        if (
            not isinstance(vector, list)
            or len(vector) != 3
            or any(
                isinstance(value, bool) or not isinstance(value, int)
                for value in vector
            )
        ):
            raise ValueError("invalid group-invariance vector")
    candidates = refutation.get("unmatched_candidate_targets")
    if (
        not isinstance(candidates, list)
        or any(
            isinstance(value, bool)
            or not isinstance(value, int)
            or not 0 <= value < row_count
            for value in candidates
        )
        or candidates != sorted(set(candidates))
    ):
        raise ValueError("invalid group-invariance candidates")
    _sha256(
        refutation.get("candidate_graph_sha256"),
        field="group-invariance candidate-graph SHA256",
    )
    expected = _operation_json(setting.nonidentity_operations[operation_index])
    if refutation.get("operation") != expected:
        raise ValueError("invalid setting-aware group-invariance operation")
    return definition_index


def _validate_presentation_group_invariance_certificate(
    row: Mapping[str, Any]
) -> None:
    status = row.get("status")
    display_label = row.get("subgroup_display_label")
    if not isinstance(display_label, str):
        raise ValueError("invalid subgroup display label")
    authority = group_authority(display_label)
    if status == "not_applicable":
        definition_rows = row.get("definition_row_count")
        if (
            row.get("reason") != "no_atomic_mode_definitions"
            or row.get("atomic_definition_count") != 0
            or row.get("definition_count") != row.get("strain_definition_count")
            or definition_rows not in {None, row.get("definition_count")}
            or not isinstance(row.get("undistorted_orbit_count"), int)
            or int(row["undistorted_orbit_count"]) <= 0
        ):
            raise ValueError("invalid not-applicable group-invariance result")
        return
    if status not in {"satisfied", "refuted", "indeterminate"}:
        raise ValueError("invalid presentation-aware group-invariance status")
    if (
        row.get("subgroup_bns_label") != authority.bns_label
        or row.get("subgroup_uni_number") != authority.uni_number
        or row.get("hall_symbol") != authority.hall_symbol
        or row.get("operation_count") != authority.operation_count
        or row.get("operation_sha256") != authority.operation_sha256
        or row.get("generator_count") != len(authority.generators)
        or row.get("nonidentity_operation_count")
        != len(authority.nonidentity_operations)
        or row.get("setting_coverage_complete")
        is not authority.setting_coverage_complete
        or row.get("setting_count") != len(authority.settings)
        or row.get("hall_setting_count")
        != sum(setting.authority == "hall" for setting in authority.settings)
        or row.get("source_presentation_count")
        != sum(setting.authority == "source" for setting in authority.settings)
    ):
        raise ValueError("invalid presentation-aware subgroup authority identity")
    if (
        row.get("position_display_scale") != POSITION_SCALE
        or row.get("vector_display_scale") != VECTOR_SCALE
    ):
        raise ValueError("invalid group-invariance display scale")
    definition_count = _nonnegative_int(
        row.get("definition_count"), field="group-invariance definition count"
    )
    atomic_count = _nonnegative_int(
        row.get("atomic_definition_count"),
        field="group-invariance atomic definition count",
    )
    strain_count = _nonnegative_int(
        row.get("strain_definition_count"),
        field="group-invariance strain definition count",
    )
    displacive_count = _nonnegative_int(
        row.get("displacive_definition_count"),
        field="group-invariance displacive definition count",
    )
    magnetic_count = _nonnegative_int(
        row.get("magnetic_definition_count"),
        field="group-invariance magnetic definition count",
    )
    compatible_count = _nonnegative_int(
        row.get("position_compatible_setting_count"),
        field="group-invariance compatible setting count",
    )
    satisfying_count = _nonnegative_int(
        row.get("satisfying_setting_count"),
        field="group-invariance satisfying setting count",
    )
    site_patterns = _nonnegative_int(
        row.get("site_pattern_count"), field="group-invariance site-pattern count"
    )
    refuted_count = _nonnegative_int(
        row.get("refuted_definition_count"),
        field="group-invariance refuted definition count",
    )
    orbit_count = _nonnegative_int(
        row.get("undistorted_orbit_count"),
        field="group-invariance undistorted orbit count",
    )
    _nonnegative_int(
        row.get("definition_row_count"),
        field="group-invariance definition-row count",
    )
    if (
        atomic_count <= 0
        or atomic_count + strain_count != definition_count
        or displacive_count + magnetic_count != atomic_count
        or site_patterns > atomic_count
        or refuted_count > atomic_count
        or satisfying_count > compatible_count
        or compatible_count > len(authority.settings)
        or orbit_count <= 0
    ):
        raise ValueError("invalid presentation-aware group-invariance totals")
    _sha256(
        row.get("definition_certificates_sha256"),
        field="group-invariance definition certificate SHA256",
    )
    available = {setting.setting_key: setting for setting in authority.settings}
    satisfying = row.get("satisfying_setting_keys")
    if (
        not isinstance(satisfying, list)
        or any(not isinstance(key, str) or key not in available for key in satisfying)
        or satisfying != sorted(set(satisfying))
        or len(satisfying) != satisfying_count
    ):
        raise ValueError("invalid satisfying presentation settings")
    setting_rows = row.get("setting_refutations")
    if not isinstance(setting_rows, list):
        raise ValueError("invalid presentation-aware refutations")
    if status == "satisfied":
        if satisfying_count != 1 or setting_rows or refuted_count:
            raise ValueError("invalid satisfied presentation-aware proof")
        return
    if satisfying_count or len(setting_rows) != len(authority.settings):
        raise ValueError("incomplete presentation-aware refutation coverage")
    if (status == "refuted") != authority.setting_coverage_complete:
        raise ValueError("invalid presentation-aware coverage verdict")
    observed: set[str] = set()
    observed_compatible = 0
    refuted_definitions: set[int] = set()
    for setting_row in setting_rows:
        if not isinstance(setting_row, Mapping):
            raise ValueError("invalid presentation-aware setting row")
        key = setting_row.get("setting_key")
        setting = available.get(key) if isinstance(key, str) else None
        if (
            setting is None
            or key in observed
            or setting_row.get("authority") != setting.authority
            or setting_row.get("hall_number") != setting.hall_number
            or setting_row.get("source_setting_ids")
            != list(setting.source_setting_ids)
            or setting_row.get("operation_sha256") != setting.operation_sha256
            or not isinstance(setting_row.get("structure_compatible"), bool)
            or not isinstance(setting_row.get("position_compatible"), bool)
        ):
            raise ValueError("invalid presentation-aware setting identity")
        observed.add(key)
        structure_compatible = bool(setting_row["structure_compatible"])
        position_compatible = bool(setting_row["position_compatible"])
        structure_refutation = setting_row.get("structure_refutation")
        refutations = setting_row.get("refutations")
        if not isinstance(refutations, list):
            raise ValueError("invalid presentation-aware setting refutations")
        if not structure_compatible:
            if position_compatible or refutations or not isinstance(
                structure_refutation, Mapping
            ):
                raise ValueError("invalid undistorted-structure refutation")
            for field in (
                "orbit_index",
                "expected_multiplicity",
                "observed_orbit_size",
            ):
                _nonnegative_int(
                    structure_refutation.get(field),
                    field=f"undistorted structure {field}",
                )
            for field in ("label", "site"):
                if not isinstance(structure_refutation.get(field), str) or not str(
                    structure_refutation[field]
                ):
                    raise ValueError("invalid undistorted-structure text")
            operation_index = structure_refutation.get("operation_index")
            if operation_index is not None:
                checked = _nonnegative_int(
                    operation_index, field="undistorted structure operation index"
                )
                if checked >= len(setting.nonidentity_operations):
                    raise ValueError("invalid undistorted-structure operation")
            continue
        if structure_refutation is not None or not refutations:
            raise ValueError("compatible presentation needs an atomic refutation")
        observed_compatible += position_compatible
        for refutation in refutations:
            definition_index = _validate_setting_refutation(
                refutation,
                definition_count=atomic_count,
                setting=setting,
            )
            if position_compatible:
                refuted_definitions.add(definition_index)
    if (
        observed != set(available)
        or observed_compatible != compatible_count
        or len(refuted_definitions) != refuted_count
    ):
        raise ValueError("invalid presentation-aware refutation summary")


def validate_group_invariance_certificate(row: Mapping[str, Any]) -> None:
    """Validate a persisted setting-aware or legacy proof certificate."""

    if row.get("claim") == "display_rounding_necessary_condition":
        _validate_legacy_group_invariance_certificate(row)
        return
    if row.get("claim") in {
        "atomic_mode_field_group_invariance",
        "atomic_mode_field_group_invariance_over_declared_presentations",
    }:
        _validate_presentation_group_invariance_certificate(row)
        return
    if row.get("theorem") != GROUP_INVARIANCE_THEOREM or row.get("claim") != (
        "display_rounding_necessary_condition_over_declared_hall_settings"
    ):
        raise ValueError("invalid setting-aware group-invariance claim")
    status = row.get("status")
    if status not in {"satisfied", "refuted", "indeterminate"}:
        raise ValueError("invalid setting-aware group-invariance status")
    display_label = row.get("subgroup_display_label")
    if not isinstance(display_label, str):
        raise ValueError("invalid subgroup display label")
    authority = group_authority(display_label)
    hall_settings = tuple(
        setting for setting in authority.settings if setting.authority == "hall"
    )
    if (
        row.get("subgroup_bns_label") != authority.bns_label
        or row.get("subgroup_uni_number") != authority.uni_number
        or row.get("hall_symbol") != authority.hall_symbol
        or row.get("operation_count") != authority.operation_count
        or row.get("operation_sha256") != authority.operation_sha256
        or row.get("generator_count") != len(authority.generators)
        or row.get("nonidentity_operation_count")
        != len(authority.nonidentity_operations)
        or row.get("setting_coverage_complete")
        is not authority.setting_coverage_complete
        or row.get("setting_count") != len(hall_settings)
    ):
        raise ValueError("invalid setting-aware subgroup authority identity")
    if (
        row.get("position_display_scale") != POSITION_SCALE
        or row.get("vector_display_scale") != VECTOR_SCALE
    ):
        raise ValueError("invalid group-invariance display scale")
    definition_count = _nonnegative_int(
        row.get("definition_count"), field="group-invariance definition count"
    )
    _nonnegative_int(
        row.get("definition_row_count"), field="group-invariance definition-row count"
    )
    displacive_count = _nonnegative_int(
        row.get("displacive_definition_count"),
        field="group-invariance displacive definition count",
    )
    magnetic_count = _nonnegative_int(
        row.get("magnetic_definition_count"),
        field="group-invariance magnetic definition count",
    )
    site_patterns = _nonnegative_int(
        row.get("site_pattern_count"), field="group-invariance site-pattern count"
    )
    refuted_count = _nonnegative_int(
        row.get("refuted_definition_count"),
        field="group-invariance refuted definition count",
    )
    compatible_count = _nonnegative_int(
        row.get("position_compatible_setting_count"),
        field="group-invariance compatible setting count",
    )
    satisfying_count = _nonnegative_int(
        row.get("satisfying_setting_count"),
        field="group-invariance satisfying setting count",
    )
    if (
        displacive_count + magnetic_count != definition_count
        or site_patterns > definition_count
        or refuted_count > definition_count
        or satisfying_count > compatible_count
        or compatible_count > len(authority.settings)
    ):
        raise ValueError("invalid setting-aware group-invariance totals")
    _sha256(
        row.get("definition_certificates_sha256"),
        field="group-invariance definition certificate SHA256",
    )
    satisfying_halls = row.get("satisfying_hall_numbers")
    available = {setting.hall_number: setting for setting in hall_settings}
    if (
        not isinstance(satisfying_halls, list)
        or satisfying_halls != sorted(set(satisfying_halls))
        or len(satisfying_halls) != satisfying_count
        or any(hall not in available for hall in satisfying_halls)
    ):
        raise ValueError("invalid satisfying Hall settings")
    setting_refutations = row.get("setting_refutations")
    if not isinstance(setting_refutations, list):
        raise ValueError("invalid setting-aware refutations")
    if status == "satisfied":
        if satisfying_count <= 0 or setting_refutations or refuted_count:
            raise ValueError("invalid satisfied setting-aware proof")
        return
    if satisfying_count or len(setting_refutations) != len(hall_settings):
        raise ValueError("incomplete setting-aware refutation coverage")
    if (status == "refuted") != authority.setting_coverage_complete:
        raise ValueError("invalid setting-aware coverage verdict")
    observed_halls: set[int] = set()
    observed_compatible = 0
    refuted_definitions: set[int] = set()
    for setting_row in setting_refutations:
        if not isinstance(setting_row, Mapping):
            raise ValueError("invalid setting-aware setting row")
        hall_number = _nonnegative_int(
            setting_row.get("hall_number"), field="setting-aware Hall number"
        )
        setting = available.get(hall_number)
        if (
            setting is None
            or hall_number in observed_halls
            or setting_row.get("operation_sha256") != setting.operation_sha256
            or not isinstance(setting_row.get("position_compatible"), bool)
        ):
            raise ValueError("invalid setting-aware setting identity")
        observed_halls.add(hall_number)
        position_compatible = bool(setting_row["position_compatible"])
        observed_compatible += position_compatible
        refutations = setting_row.get("refutations")
        if not isinstance(refutations, list) or not refutations:
            raise ValueError("setting-aware setting needs a refutation")
        prior_index = -1
        for refutation in refutations:
            definition_index = _validate_setting_refutation(
                refutation,
                definition_count=definition_count,
                setting=setting,
            )
            if definition_index <= prior_index:
                raise ValueError("setting-aware refutations are not ordered")
            prior_index = definition_index
            if position_compatible:
                if refutation["stage"] != "vector":
                    raise ValueError("compatible setting has a position refutation")
                refuted_definitions.add(definition_index)
            elif refutation["stage"] != "position" or len(refutations) != 1:
                raise ValueError("incompatible setting needs one position refutation")
    if (
        observed_halls != set(available)
        or observed_compatible != compatible_count
        or len(refuted_definitions) != refuted_count
    ):
        raise ValueError("invalid setting-aware refutation summary")


__all__ = [
    "GROUP_INVARIANCE_THEOREM",
    "GroupOperation",
    "assess_group_invariance",
    "group_authority",
    "subgroup_label_from_local_payload",
    "subgroup_label_from_web_text",
    "validate_group_invariance_certificate",
]
