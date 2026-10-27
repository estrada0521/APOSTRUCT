"""Fixed-space dimension of the complete mode list.

For each atom block (a kind and the set of child atoms its definitions occupy), the number of definitions equals the
dimension of the subspace of that block's displacement (magnetic, occupancy) field left invariant by the child group:
the average over the group's operations of the trace of their action, a permutation of the atoms times the vector
action on each fixed atom. With linear independence and group invariance, this makes the definitions of a block a
basis of its whole invariant subspace. The irrep each definition is attributed to is not part of the claim."""

from __future__ import annotations

from collections import defaultdict
from collections.abc import Sequence
from fractions import Fraction
from typing import Any

from .group_invariance import (
    _DECLARED_CLAIM,
    _DERIVED_CLAIM,
    _carried_rows,
    _maximum_matching,
    _position_candidates,
    _structure_failure,
    _vector_action,
    presentation_authority,
)
from .mode_basis import QuantizedModeDefinition, QuantizedStrainDefinition, QuantizedStructureOrbit

FIXED_SPACE_THEOREM = "complete_mode_fixed_space_dimension.v1"


def _block_dimension(kind: str, rows, operations) -> Fraction | None:
    width = len(rows) if kind == "occ" else 3 * len(rows)
    total = Fraction(width)
    for operation in operations:
        graph = _position_candidates(rows, operation)
        matching = graph.bijection if graph.bijection is not None else _maximum_matching(graph.edges)
        if any(target < 0 for target in matching):
            return None
        action = _vector_action(operation, kind)
        trace = 1 if kind == "occ" else sum(action[axis][axis] for axis in range(3))
        total += trace * sum(1 for source, target in enumerate(matching) if source == target)
    return total / (len(operations) + 1)


def assess_fixed_space_dimension(
    definitions: Sequence[QuantizedModeDefinition | QuantizedStrainDefinition],
    subgroup_label: str,
    structure: Sequence[QuantizedStructureOrbit],
    satisfying_setting_keys: Sequence[str] = (),
) -> dict[str, Any]:
    """Count against dimension per atom block, in the presentation where the definitions were found invariant, else
    in the first presentation whose operations permute every block's atoms."""

    atomic = [definition for definition in definitions if isinstance(definition, QuantizedModeDefinition)]
    base = {
        "theorem": FIXED_SPACE_THEOREM,
        "definition_count": len(definitions),
        "definition_row_count": sum(len(definition.rows) for definition in atomic) + (len(definitions) - len(atomic)),
    }
    if not atomic:
        return {**base, "status": "not_applicable", "reason": "no_atomic_mode_definitions", "group_count": 0,
                "refuted_group_count": 0, "groups": []}
    blocks: dict[tuple[str, tuple], list[QuantizedModeDefinition]] = defaultdict(list)
    for definition in atomic:
        rows = _carried_rows(definition)
        blocks[(definition.kind, tuple((row.atom, row.position) for row in rows))].append(definition)
    wanted = list(satisfying_setting_keys)
    settings = []
    for claim in (_DECLARED_CLAIM, _DERIVED_CLAIM):
        for setting in presentation_authority(subgroup_label, claim).settings:
            if setting.setting_key not in {s.setting_key for s in settings}:
                settings.append(setting)
    settings.sort(key=lambda setting: (setting.setting_key not in wanted, wanted.index(setting.setting_key) if setting.setting_key in wanted else 0))
    for setting in settings:
        operations = tuple(setting.nonidentity_operations)
        if _structure_failure(structure, operations):
            continue
        groups = []
        for (kind, pattern), members in blocks.items():
            dimension = _block_dimension(kind, _carried_rows(members[0]), operations)
            if dimension is None:
                break
            groups.append({"kind": kind, "atom": pattern[0][0], "atoms": len(pattern), "definitions": len(members),
                           "dimension": str(dimension)})
        else:
            refuted = [group for group in groups if Fraction(group["dimension"]) != group["definitions"]]
            return {**base, "status": "refuted" if refuted else "satisfied", "setting_key": setting.setting_key,
                    "group_count": len(groups), "refuted_group_count": len(refuted), "groups": groups}
    return {**base, "status": "indeterminate", "reason": "no_presentation_permutes_the_atoms",
            "group_count": len(blocks), "refuted_group_count": 0, "groups": []}
