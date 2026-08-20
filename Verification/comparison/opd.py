"""Compare a selected Web OPD row with the corresponding local OPD row."""

from __future__ import annotations

from fractions import Fraction
import re
from typing import Any, Mapping, Sequence

from Verification.comparison.basis_lattice import unimodular_basis_change
from APOSTRUCT.Backend.modes.engine.decoder import ModeDataDecoder
from APOSTRUCT.Backend.modes.engine.subgroup_structure.presentation_transport import (
    embedded_child_operation_records,
    embedded_magnetic_child_operation_records,
)
from APOSTRUCT.Backend.source.tables import source_tables


_OPD_ROW_RE = re.compile(
    r"^(?P<token>\S+)\s+"
    r"(?P<direction>\([^)]*\))\s+"
    r"(?P<child>\d+(?:\.\d+)?)\s+"
    r"(?P<symbol>[^,]+),\s*"
    r"basis=\{(?P<basis>[^}]*)\},\s*"
    r"origin=(?P<origin>\([^)]*\)),\s*"
    r"s=(?P<s>-?\d+),\s*i=(?P<i>-?\d+),\s*"
    r"k-active=\s*(?P<k_active>.*)$"
)

_K_ACTIVE_RE = re.compile(r"\s*\([^()]*\)(?:\s*[,;]\s*\([^()]*\))*\s*\Z")


def row_label(row: object) -> str:
    if not isinstance(row, Mapping):
        return ""
    isotropy = row.get("isotropy")
    direction = row.get("direction")
    if not isinstance(isotropy, Mapping):
        isotropy = {}
    if not isinstance(direction, Mapping):
        direction = {}
    return str(isotropy.get("opd_label") or direction.get("label") or "").strip()


def parse_web_opd(raw: str, label: str) -> dict[str, Any]:
    visible = raw.split("*", 1)[1].strip() if "*" in raw else raw.strip()
    match = _OPD_ROW_RE.fullmatch(visible)
    if match is None:
        raise ValueError(f"cannot parse selected Web OPD row: {raw!r}")
    return {
        "label": label,
        "display_token": match.group("token"),
        "direction": match.group("direction"),
        "child_sg": match.group("child"),
        "symbol": match.group("symbol").strip(),
        "basis": match.group("basis").strip(),
        "origin": match.group("origin").strip(),
        "s": int(match.group("s")),
        "i": int(match.group("i")),
        "k_active": match.group("k_active").strip(),
    }


def _basis_text(value: object) -> str:
    if isinstance(value, str):
        return value.strip().removeprefix("{").removesuffix("}")
    if not isinstance(value, Sequence) or isinstance(value, (str, bytes)):
        return ""
    rows: list[str] = []
    for row in value:
        if not isinstance(row, Sequence) or isinstance(row, (str, bytes)):
            return ""
        rows.append("(" + ",".join(str(item) for item in row) + ")")
    return ",".join(rows)


def local_opd(row: Mapping[str, Any], *, magnetic: bool) -> dict[str, Any]:
    isotropy = row.get("isotropy")
    direction = row.get("direction")
    if not isinstance(isotropy, Mapping):
        raise ValueError("selected local OPD has no isotropy payload")
    if not isinstance(direction, Mapping):
        direction = {}
    subgroup = isotropy.get("subgroup")
    if not isinstance(subgroup, Mapping):
        subgroup = {}
    display_label = subgroup.get("display_label")
    if magnetic and display_label:
        child_sg = str(display_label)
    else:
        child_sg = str(subgroup.get("number") or subgroup.get("ordinary_number") or "")
    basis = isotropy.get("basis_text") or _basis_text(isotropy.get("basis"))
    return {
        "label": row_label(row),
        "direction": str(
            isotropy.get("display_opd")
            or isotropy.get("source_opd")
            or direction.get("opd")
            or isotropy.get("direction_opd")
            or ""
        ).strip(),
        "child_sg": child_sg,
        "symbol": str(subgroup.get("symbol") or ""),
        "basis": str(basis or "").strip(),
        "origin": str(isotropy.get("origin") or "").strip(),
        "s": isotropy.get("s", isotropy.get("arms")),
        "i": isotropy.get("i"),
        "k_active": str(
            isotropy.get("k_active") or isotropy.get("k_active_text") or ""
        ).strip(),
    }


def _parenthesized_rows(value: object) -> tuple[str, ...]:
    text = str(value or "").strip()
    rows: list[str] = []
    depth = 0
    start: int | None = None
    for index, character in enumerate(text):
        if character == "(":
            if depth == 0:
                start = index + 1
            depth += 1
        elif character == ")":
            depth -= 1
            if depth < 0:
                return ()
            if depth == 0 and start is not None:
                rows.append(text[start:index])
                start = None
    return tuple(rows) if depth == 0 else ()


def _fraction_tuple(value: str, *, count: int) -> tuple[Fraction, ...] | None:
    parts = tuple(item.strip().replace("−", "-") for item in value.split(","))
    if len(parts) != count:
        return None
    try:
        return tuple(Fraction(item) for item in parts)
    except (ValueError, ZeroDivisionError):
        return None


def _basis(value: object) -> tuple[tuple[Fraction, Fraction, Fraction], ...] | None:
    rows = _parenthesized_rows(value)
    parsed = tuple(_fraction_tuple(row, count=3) for row in rows)
    if len(parsed) != 3 or any(row is None for row in parsed):
        return None
    return parsed  # type: ignore[return-value]


def _origin(value: object) -> tuple[Fraction, Fraction, Fraction] | None:
    rows = _parenthesized_rows(value)
    if len(rows) != 1:
        return None
    return _fraction_tuple(rows[0], count=3)  # type: ignore[return-value]


def _determinant(matrix: Sequence[Sequence[Fraction]]) -> Fraction:
    return (
        matrix[0][0] * (matrix[1][1] * matrix[2][2] - matrix[1][2] * matrix[2][1])
        - matrix[0][1]
        * (matrix[1][0] * matrix[2][2] - matrix[1][2] * matrix[2][0])
        + matrix[0][2]
        * (matrix[1][0] * matrix[2][1] - matrix[1][1] * matrix[2][0])
    )


def _inverse(
    matrix: tuple[tuple[Fraction, Fraction, Fraction], ...],
) -> tuple[tuple[Fraction, Fraction, Fraction], ...] | None:
    det = _determinant(matrix)
    if det == 0:
        return None
    return (
        (
            (matrix[1][1] * matrix[2][2] - matrix[1][2] * matrix[2][1]) / det,
            (matrix[0][2] * matrix[2][1] - matrix[0][1] * matrix[2][2]) / det,
            (matrix[0][1] * matrix[1][2] - matrix[0][2] * matrix[1][1]) / det,
        ),
        (
            (matrix[1][2] * matrix[2][0] - matrix[1][0] * matrix[2][2]) / det,
            (matrix[0][0] * matrix[2][2] - matrix[0][2] * matrix[2][0]) / det,
            (matrix[0][2] * matrix[1][0] - matrix[0][0] * matrix[1][2]) / det,
        ),
        (
            (matrix[1][0] * matrix[2][1] - matrix[1][1] * matrix[2][0]) / det,
            (matrix[0][1] * matrix[2][0] - matrix[0][0] * matrix[2][1]) / det,
            (matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0]) / det,
        ),
    )


def _origin_coset_equal(web: Mapping[str, Any], local: Mapping[str, Any]) -> bool:
    web_basis = _basis(web.get("basis"))
    web_origin = _origin(web.get("origin"))
    local_origin = _origin(local.get("origin"))
    if web_basis is None or web_origin is None or local_origin is None:
        return False
    inverse = _inverse(web_basis)
    if inverse is None:
        return False
    delta = tuple(local_origin[index] - web_origin[index] for index in range(3))
    coefficients = tuple(
        sum(delta[axis] * inverse[axis][column] for axis in range(3))
        for column in range(3)
    )
    return all(value.denominator == 1 for value in coefficients)


def _operation_record_set(
    records: object,
) -> frozenset[tuple[Fraction, Fraction, Fraction, int]] | None:
    if not isinstance(records, Sequence) or isinstance(records, (str, bytes)):
        return None
    normalized: set[tuple[Fraction, Fraction, Fraction, int]] = set()
    for record in records:
        if (
            not isinstance(record, Sequence)
            or isinstance(record, (str, bytes))
            or len(record) != 5
        ):
            return None
        try:
            denominator = int(record[3])
            if denominator == 0:
                return None
            normalized.add(
                (
                    Fraction(int(record[0]), denominator) % 1,
                    Fraction(int(record[1]), denominator) % 1,
                    Fraction(int(record[2]), denominator) % 1,
                    int(record[4]),
                )
            )
        except (TypeError, ValueError, ZeroDivisionError):
            return None
    return frozenset(normalized) if normalized else None


def _embedded_operation_sets_equal(
    web: Mapping[str, Any],
    *,
    parent_sg: int | None,
    parent_setting_id: int | None,
    subgroup_setting_id: int | None,
    child_ordinary_sg: int | None,
    child_magnetic_group: int | None,
    local_operation_records: object,
) -> bool:
    web_basis = _basis(web.get("basis"))
    web_origin = _origin(web.get("origin"))
    local_records = _operation_record_set(local_operation_records)
    if (
        web_basis is None
        or web_origin is None
        or local_records is None
        or parent_sg is None
        or child_ordinary_sg is None
    ):
        return False
    try:
        tables = source_tables()
        source_basis, source_origin = tables.subgroup_embedding_from_cinter(
            int(parent_sg),
            int(child_ordinary_sg),
            web_basis,
            web_origin,
            parent_setting_id=(
                int(parent_setting_id) if parent_setting_id is not None else None
            ),
            subgroup_setting_id=(
                int(subgroup_setting_id) if subgroup_setting_id is not None else None
            ),
        )
        decoder = ModeDataDecoder(tables=tables)
        if child_magnetic_group is None:
            web_records = embedded_child_operation_records(
                decoder,
                parent_sg=int(parent_sg),
                child_sg=int(child_ordinary_sg),
                subgroup_basis=source_basis,
                subgroup_origin=source_origin,
            )
        else:
            web_records = embedded_magnetic_child_operation_records(
                decoder,
                parent_sg=int(parent_sg),
                child_magnetic_group=int(child_magnetic_group),
                subgroup_basis=source_basis,
                subgroup_origin=source_origin,
            )
    except (IndexError, KeyError, TypeError, ValueError, ZeroDivisionError):
        return False
    return _operation_record_set(web_records) == local_records


def _basis_relation(web: object, local: object) -> tuple[bool, object | None]:
    web_basis = _basis(web)
    local_basis = _basis(local)
    if web_basis is None or local_basis is None:
        return False, None
    web_flat = tuple(value for row in web_basis for value in row)
    local_flat = tuple(value for row in local_basis for value in row)
    change = unimodular_basis_change(web_flat, local_flat)
    return change is not None, change


def _compact(value: object) -> str:
    return re.sub(r"\s+", "", str(value or "")).replace("−", "-")


def _direction(value: object) -> str:
    return _compact(value).replace(";", ",")


def _k_active(value: object) -> frozenset[tuple[Fraction, Fraction, Fraction]] | None:
    if not isinstance(value, str) or _K_ACTIVE_RE.fullmatch(value) is None:
        return None
    rows = _parenthesized_rows(value)
    parsed = tuple(_fraction_tuple(row, count=3) for row in rows)
    if not parsed or any(row is None for row in parsed):
        return None
    return frozenset(parsed)  # type: ignore[arg-type]


def _k_active_relation(web: object, local: object) -> tuple[bool, bool]:
    web_active = _k_active(web)
    local_active = _k_active(local)
    if web_active is None or local_active is None:
        return False, False
    exact = web_active == local_active
    def physical_representative(
        vector: tuple[Fraction, Fraction, Fraction],
    ) -> tuple[Fraction, Fraction, Fraction]:
        quotient = tuple(component % 1 for component in vector)
        conjugate = tuple((-component) % 1 for component in vector)
        return min(quotient, conjugate)  # type: ignore[return-value]

    web_quotient = frozenset(physical_representative(vector) for vector in web_active)
    local_quotient = frozenset(
        physical_representative(vector) for vector in local_active
    )
    return web_quotient == local_quotient, exact


def compare_opd(
    web: Mapping[str, Any],
    local: Mapping[str, Any],
    *,
    source: str,
    direction_matched_by_label: bool = False,
    parent_sg: int | None = None,
    parent_setting_id: int | None = None,
    subgroup_setting_id: int | None = None,
    child_ordinary_sg: int | None = None,
    child_magnetic_group: int | None = None,
    local_operation_records: object = None,
) -> dict[str, Any]:
    same_lattice, basis_change = _basis_relation(web.get("basis"), local.get("basis"))
    physical_k_active, strict_k_active = _k_active_relation(
        web.get("k_active"), local.get("k_active")
    )
    lattice_origin_equivalent = same_lattice and _origin_coset_equal(web, local)
    operation_set_equivalent = same_lattice and (
        not lattice_origin_equivalent
        and _embedded_operation_sets_equal(
            web,
            parent_sg=parent_sg,
            parent_setting_id=parent_setting_id,
            subgroup_setting_id=subgroup_setting_id,
            child_ordinary_sg=child_ordinary_sg,
            child_magnetic_group=child_magnetic_group,
            local_operation_records=local_operation_records,
        )
    )
    origin_equivalent = lattice_origin_equivalent or operation_set_equivalent
    physical = {
        "child_sg": str(web.get("child_sg")) == str(local.get("child_sg")),
        "k_active": physical_k_active,
        "dir": (
            True
            if direction_matched_by_label
            else _direction(web.get("direction")) == _direction(local.get("direction"))
        ),
        "sublattice": same_lattice,
        "origin_coset": origin_equivalent,
        "s": web.get("s") == local.get("s"),
        "i": web.get("i") == local.get("i"),
    }
    strict = {
        "basis": _compact(web.get("basis")) == _compact(local.get("basis")),
        "origin": _compact(web.get("origin")) == _compact(local.get("origin")),
        "k_active": strict_k_active,
    }
    physical_pass = all(physical.values())
    strict_pass = physical_pass and all(strict.values())
    return {
        "exist": True,
        "level": (
            "strict"
            if strict_pass
            else "physical"
            if physical_pass
            else "present_nonphysical"
        ),
        "source": source,
        "label": web.get("label"),
        "physical": physical,
        "strict": strict,
        "origin_equivalence": (
            "lattice_coset"
            if lattice_origin_equivalent
            else "operation_set"
            if operation_set_equivalent
            else None
        ),
        "web": dict(web),
        "local": dict(local),
        **({"basis_change": basis_change} if basis_change is not None else {}),
    }


__all__ = ["compare_opd", "local_opd", "parse_web_opd", "row_label"]
