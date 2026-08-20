"""Compare the requested and locally selected subgroup frame."""

from __future__ import annotations

from fractions import Fraction
import re
from typing import Any, Mapping, Sequence

from Verification.comparison.basis_lattice import (
    unimodular_basis_change as _unimodular_basis_change,
)


def _request_slot_identity(
    selections: Sequence[Mapping[str, Any]],
) -> tuple[Any, ...] | None:
    if not selections:
        return None
    return (
        "slots",
        tuple(
            (
                str(selection.get("k_label") or ""),
                tuple(
                    sorted(
                        (str(key), str(Fraction(str(value))))
                        for key, value in (selection.get("k_params") or {}).items()
                    )
                ),
                str(selection.get("irrep") or ""),
            )
            for selection in selections
        ),
    )


def _local_slot_identity(preview: dict[str, Any]) -> tuple[Any, ...] | None:
    slots = ((preview.get("selected") or {}).get("selection_slots") or [])
    if not slots:
        return None
    return (
        "slots",
        tuple(
            (
                str((slot.get("kpoint") or {}).get("label") or ""),
                tuple(
                    sorted(
                        (str(key), str(Fraction(str(value))))
                        for key, value in (slot.get("k_params") or {}).items()
                    )
                ),
                str((slot.get("irrep") or {}).get("symbol") or ""),
            )
            for slot in slots
            if isinstance(slot, dict)
        ),
    )


def _rational_tokens(value: object) -> tuple[str, ...]:
    pattern = r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?(?:/\d+)?"
    values = []
    for token in re.findall(pattern, str(value)):
        fraction = Fraction(token)
        if "." in token or "e" in token.lower():
            fraction = fraction.limit_denominator(1_000_000)
        values.append(str(fraction))
    return tuple(values)


def _request_selected_state(
    selections: Sequence[Mapping[str, Any]], subgroup: Mapping[str, Any]
) -> tuple[Any, ...] | None:
    number = subgroup.get("number")
    basis = subgroup.get("basis")
    origin = subgroup.get("origin")
    if (
        isinstance(number, bool)
        or not isinstance(number, int)
        or not isinstance(basis, str)
        or not basis
        or not isinstance(origin, str)
        or not origin
    ):
        return None
    return (
        _request_slot_identity(selections),
        number,
        _rational_tokens(basis),
        _rational_tokens(origin),
    )


def _local_selected_state(
    preview: dict[str, Any], selections: Sequence[Mapping[str, Any]]
) -> tuple[Any, ...] | None:
    opd = ((preview.get("selected") or {}).get("orderparam") or {})
    isotropy = opd.get("isotropy") if isinstance(opd, dict) else None
    if not isinstance(isotropy, dict):
        return None
    subgroup = isotropy.get("subgroup")
    if isinstance(subgroup, dict):
        subgroup = subgroup.get("ordinary_number") or subgroup.get("number")
    basis = isotropy.get("basis")
    origin = isotropy.get("origin")
    if subgroup is None or basis is None or origin is None:
        return None
    return (
        _local_slot_identity(preview)
        or _request_slot_identity(selections),
        int(subgroup),
        _rational_tokens(basis),
        _rational_tokens(origin),
    )


def compare_selected_state(
    preview: dict[str, Any],
    *,
    selections: Sequence[Mapping[str, Any]],
    web_subgroup: Mapping[str, Any],
    operation_set_equivalent: bool = False,
) -> dict[str, Any]:
    expected = _request_selected_state(selections, web_subgroup)
    actual = _local_selected_state(preview, selections)
    basis_change = None
    if (
        expected is not None
        and actual is not None
        and expected[0] == actual[0]
        and expected[1] == actual[1]
        and (expected[3] == actual[3] or operation_set_equivalent)
    ):
        basis_change = _unimodular_basis_change(expected[2], actual[2])
    status = (
        "missing"
        if expected is None or actual is None
        else "exact"
        if expected == actual
        else "equivalent_basis"
        if basis_change is not None and expected[3] == actual[3]
        else "equivalent_embedding"
        if basis_change is not None and operation_set_equivalent
        else "diff"
    )
    return {
        "status": status,
        "request": expected,
        "local": actual,
        **(
            {"basis_change": basis_change}
            if status in {"equivalent_basis", "equivalent_embedding"}
            else {}
        ),
    }


__all__ = ["compare_selected_state"]
