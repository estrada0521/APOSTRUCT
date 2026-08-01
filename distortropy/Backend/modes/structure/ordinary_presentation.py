"""Source-only matching of Formula15 rows to presentation-grid orbits."""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from numbers import Integral
import re
from typing import Any, Mapping, Sequence

import gemmi

from distortropy.Backend.exactmath import (
    fraction_matrix_inverse3,
    fraction_row_multiply3 as _row_multiply,
)


ExactPoint = tuple[Fraction, Fraction, Fraction]
ExactMatrix3 = tuple[ExactPoint, ExactPoint, ExactPoint]


def _exact_integer(value: Any) -> int:
    if isinstance(value, bool) or not isinstance(value, Integral):
        raise TypeError("Source occurrence identity must contain exact integers")
    return int(value)


def _exact_integer_record(value: Any, length: int) -> tuple[int, ...]:
    if not isinstance(value, (list, tuple)) or len(value) != length:
        raise ValueError(f"Source occurrence identity requires {length} integers")
    return tuple(_exact_integer(item) for item in value)


@dataclass(frozen=True)
class FormulaChildAtom:
    """One Source occurrence expressed in the selected child setting."""

    source_flat_index: int
    kernel_fraction_ordinal: int
    kernel_fraction: tuple[int, int, int, int]
    parent_branch_ordinal: int
    parent_coset_record: tuple[int, int, int, int, int]
    centering_ordinal: int
    parent_xyz: ExactPoint
    xyz: ExactPoint


@dataclass(frozen=True)
class FormulaChildSite:
    """One Formula15 child site with its complete atom membership."""

    formula_site: str
    atoms: tuple[FormulaChildAtom, ...]
    representative_atom_index: int


def _mod1(value: Fraction) -> Fraction:
    return value % 1


def _input_fraction(value: Any) -> Fraction:
    """Recover the exact decimal carried through presentation float math."""

    if isinstance(value, float):
        decimal = Fraction(str(round(value, 9)))
        crystallographic = decimal.limit_denominator(384)
        if abs(crystallographic - decimal) <= Fraction(1, 200_000_000):
            return crystallographic
        return decimal
    return Fraction(str(value))


def _unfolded_point(values: Any) -> ExactPoint | None:
    if isinstance(values, str):
        text = values.strip()
        if text.startswith("(") and text.endswith(")"):
            text = text[1:-1]
        values = tuple(part.strip() for part in text.split(","))
    if not isinstance(values, (list, tuple)):
        return None
    if len(values) == 4:
        try:
            record = _exact_integer_record(values, 4)
            denominator = record[3]
            if denominator <= 0:
                return None
            return tuple(Fraction(record[axis], denominator) for axis in range(3))  # type: ignore[return-value]
        except (TypeError, ValueError, ZeroDivisionError):
            return None
    if len(values) != 3:
        return None
    try:
        return tuple(_input_fraction(values[axis]) for axis in range(3))  # type: ignore[return-value]
    except (TypeError, ValueError, ZeroDivisionError):
        return None


def _matrix3(values: Any) -> ExactMatrix3 | None:
    if not isinstance(values, (list, tuple)) or len(values) != 3:
        return None
    rows = tuple(_unfolded_point(row) for row in values)
    if any(row is None for row in rows):
        return None
    return tuple(row for row in rows if row is not None)  # type: ignore[return-value]


def _inverse3(matrix: ExactMatrix3) -> ExactMatrix3 | None:
    try:
        return fraction_matrix_inverse3(matrix)
    except ValueError:
        return None


def _site_multiplicity(site: str) -> int | None:
    match = re.fullmatch(r"\s*(\d+)\s*[^\s]+\s*", str(site or ""))
    return int(match.group(1)) if match else None


def _centering_translations(child_sg: int) -> tuple[ExactPoint, ...] | None:
    if (
        isinstance(child_sg, bool)
        or not isinstance(child_sg, Integral)
        or not 1 <= int(child_sg) <= 230
    ):
        return None
    try:
        space_group = gemmi.find_spacegroup_by_number(int(child_sg))
        if space_group is None:
            return None
        operations = tuple(space_group.operations())
    except (RuntimeError, TypeError, ValueError):
        return None
    denominator = int(gemmi.Op.DEN)
    identity = (
        (denominator, 0, 0),
        (0, denominator, 0),
        (0, 0, denominator),
    )
    translations: list[ExactPoint] = []
    for operation in operations:
        rotation = tuple(tuple(int(value) for value in row) for row in operation.rot)
        if rotation != identity:
            continue
        translation = tuple(
            _mod1(Fraction(int(operation.tran[axis]), denominator)) for axis in range(3)
        )
        if translation not in translations:
            translations.append(translation)  # type: ignore[arg-type]
    zero = (Fraction(0), Fraction(0), Fraction(0))
    if zero not in translations:
        return None
    return (zero, *(translation for translation in translations if translation != zero))


def formula_child_sites_in_presentation(
    *,
    child_sg: int | None,
    formula_rows: Sequence[Mapping[str, Any]],
    subgroup_basis: Sequence[Sequence[Any]],
    subgroup_origin: Sequence[Any],
    child_pml_to_cinter_matrix: Sequence[Sequence[Any]],
    child_pml_to_cinter_origin: Sequence[Any],
) -> tuple[FormulaChildSite, ...] | None:
    """Expand Formula15 sites directly into selected-cell atom members."""

    try:
        exact_child_sg = _exact_integer(child_sg)
    except (TypeError, ValueError):
        return None
    if not 1 <= exact_child_sg <= 230 or not formula_rows:
        return None
    result_basis = _matrix3(subgroup_basis)
    result_origin = _unfolded_point(subgroup_origin)
    result_inverse = None if result_basis is None else _inverse3(result_basis)
    child_pml_to_cinter = _matrix3(child_pml_to_cinter_matrix)
    child_cinter_origin = _unfolded_point(child_pml_to_cinter_origin)
    child_cinter_to_pml = (
        None if child_pml_to_cinter is None else _inverse3(child_pml_to_cinter)
    )
    if (
        result_basis is None
        or result_origin is None
        or result_inverse is None
        or child_pml_to_cinter is None
        or child_cinter_origin is None
        or child_cinter_to_pml is None
    ):
        return None

    sites: list[FormulaChildSite] = []
    used_identities: set[tuple[int, int]] = set()
    used_points: set[ExactPoint] = set()
    for row in formula_rows:
        raw_site = row.get("site")
        if not isinstance(raw_site, str) or not raw_site:
            return None
        site = raw_site
        multiplicity = _site_multiplicity(site)
        occurrences = row.get("_wyckoff_source_occurrences")
        if (
            multiplicity is None
            or not isinstance(occurrences, (list, tuple))
            or not occurrences
        ):
            return None
        occurrence_count = len(occurrences)
        if multiplicity <= 0 or multiplicity % occurrence_count:
            return None
        expansion = multiplicity // occurrence_count
        if expansion == 1:
            translations = ((Fraction(0), Fraction(0), Fraction(0)),)
        else:
            all_translations = _centering_translations(exact_child_sg)
            if all_translations is None or len(all_translations) != expansion:
                return None
            translations = all_translations

        atoms: list[FormulaChildAtom] = []
        representative_atom_index: int | None = None
        try:
            representative_flat_index = _exact_integer(
                row["_wyckoff_source_representative_flat_index"]
            )
        except (KeyError, TypeError, ValueError):
            return None
        if representative_flat_index < 0:
            return None
        for occurrence in occurrences:
            if not isinstance(occurrence, Mapping):
                return None
            parent_pml = _unfolded_point(occurrence.get("parent_point_fraction"))
            try:
                flat_index = _exact_integer(occurrence["flat_index"])
                kernel_fraction_ordinal = _exact_integer(
                    occurrence["kernel_fraction_ordinal"]
                )
                kernel_fraction = _exact_integer_record(
                    occurrence["kernel_fraction"], 4
                )
                parent_branch_ordinal = _exact_integer(
                    occurrence["parent_branch_ordinal"]
                )
                parent_coset_record = _exact_integer_record(
                    occurrence["parent_coset_record"], 5
                )
            except (KeyError, TypeError, ValueError):
                return None
            if (
                parent_pml is None
                or flat_index < 0
                or kernel_fraction_ordinal < 0
                or kernel_fraction[3] <= 0
                or parent_branch_ordinal < 0
                or parent_coset_record[3] <= 0
                or parent_coset_record[4] <= 0
            ):
                return None
            child_pml = _row_multiply(
                tuple(parent_pml[axis] - result_origin[axis] for axis in range(3)),
                result_inverse,
            )
            child_cinter_linear = _row_multiply(child_pml, child_pml_to_cinter)
            child_cinter = tuple(
                child_cinter_linear[axis] + child_cinter_origin[axis]
                for axis in range(3)
            )
            for centering_ordinal, translation in enumerate(translations):
                identity = (flat_index, centering_ordinal)
                xyz = tuple(
                    _mod1(child_cinter[axis] + translation[axis]) for axis in range(3)
                )
                centered_cinter = tuple(
                    child_cinter[axis] + translation[axis] for axis in range(3)
                )
                centered_child_pml = _row_multiply(
                    tuple(
                        centered_cinter[axis] - child_cinter_origin[axis]
                        for axis in range(3)
                    ),
                    child_cinter_to_pml,
                )
                parent_xyz_linear = _row_multiply(centered_child_pml, result_basis)
                parent_xyz = tuple(
                    parent_xyz_linear[axis] + result_origin[axis] for axis in range(3)
                )
                if identity in used_identities or xyz in used_points:
                    return None
                used_identities.add(identity)
                used_points.add(xyz)
                if identity == (representative_flat_index, 0):
                    representative_atom_index = len(atoms)
                atoms.append(
                    FormulaChildAtom(
                        source_flat_index=flat_index,
                        kernel_fraction_ordinal=kernel_fraction_ordinal,
                        kernel_fraction=kernel_fraction,
                        parent_branch_ordinal=parent_branch_ordinal,
                        parent_coset_record=parent_coset_record,
                        centering_ordinal=centering_ordinal,
                        parent_xyz=parent_xyz,
                        xyz=xyz,
                    )
                )
        if len(atoms) != multiplicity or representative_atom_index is None:
            return None
        sites.append(
            FormulaChildSite(
                formula_site=site,
                atoms=tuple(atoms),
                representative_atom_index=representative_atom_index,
            )
        )
    return tuple(sites)
