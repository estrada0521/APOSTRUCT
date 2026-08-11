"""Canonical child-site and atom identity for complete mode output."""

from __future__ import annotations

from dataclasses import dataclass, replace
from fractions import Fraction
from numbers import Integral
import re
from typing import Mapping, Sequence

from APOSTRUCT.Backend.exactmath import fraction_matrix_inverse3

from APOSTRUCT.Backend.modes.structure.ordinary_presentation import (
    FormulaChildAtom,
    FormulaChildSite,
    _input_fraction,
)


ExactPoint = tuple[Fraction, Fraction, Fraction]
OperationRecord = tuple[int, int, int, int, int]


def _exact_nonempty_string(value: object, field: str) -> str:
    if not isinstance(value, str):
        raise TypeError(f"{field} must be an exact string")
    if not value:
        raise ValueError(f"{field} must not be empty")
    return value


def _exact_nonnegative_integer(value: object, field: str) -> int:
    if isinstance(value, bool) or not isinstance(value, Integral):
        raise TypeError(f"{field} must be an exact integer")
    result = int(value)
    if result < 0:
        raise ValueError(f"{field} must not be negative")
    return result


def _exact_integer_record(value: object, length: int, field: str) -> tuple[int, ...]:
    if not isinstance(value, (list, tuple)) or len(value) != length:
        raise ValueError(f"{field} requires {length} integers")
    if any(isinstance(item, bool) or not isinstance(item, Integral) for item in value):
        raise TypeError(f"{field} must contain exact integers")
    return tuple(int(item) for item in value)


@dataclass(frozen=True)
class ChildAtom:
    """One physical atom member of a displayed child Wyckoff site."""

    atom_id: str
    child_site_id: str
    child_site_order: int
    member_order: int
    source_raw_index: int
    source_flat_index: int
    source_kernel_fraction: tuple[int, int, int, int]
    source_parent_coset_record: tuple[int, int, int, int, int]
    centering_ordinal: int
    parent_xyz: ExactPoint
    structure_xyz: ExactPoint


@dataclass(frozen=True)
class ChildSite:
    """One displayed child Wyckoff site and its complete atom membership."""

    child_site_id: str
    child_site_order: int
    wyckoff_site: str
    representative_atom_id: str
    representative_xyz: ExactPoint
    atom_ids: tuple[str, ...]


@dataclass(frozen=True)
class ChildAtomPresentationRow:
    """One Source display row bound once to a canonical child atom."""

    atom_id: str
    source_raw_index: int
    centering_ordinal: int
    operation_record: OperationRecord


@dataclass(frozen=True)
class ChildAtomLayout:
    """The single atom/site identity consumed by structure, modes, and viewer."""

    sites: tuple[ChildSite, ...]
    atoms: tuple[ChildAtom, ...]
    presentation_rows: tuple[ChildAtomPresentationRow, ...] = ()
    point_operation_count: int = 0


def exact_operation_record(
    value: Sequence[int],
    *,
    point_operation_count: int,
) -> OperationRecord:
    if isinstance(point_operation_count, bool) or not isinstance(
        point_operation_count, Integral
    ):
        raise TypeError("point-operation count must be an exact integer")
    if point_operation_count <= 0:
        raise ValueError("point-operation count must be positive")
    operation_count = int(point_operation_count)
    if len(value) != 5:
        raise ValueError(f"invalid atom operation record: {value!r}")
    if any(isinstance(item, bool) or not isinstance(item, Integral) for item in value):
        raise TypeError(f"atom operation record must contain exact integers: {value!r}")
    record = tuple(int(item) for item in value)
    if record[3] <= 0 or not 1 <= record[4] <= operation_count:
        raise ValueError(f"invalid atom operation record: {value!r}")
    return record  # type: ignore[return-value]


def operation_record_row_assignments(
    canonical_records: Sequence[Sequence[int]],
    row_records: Sequence[Sequence[int]],
    *,
    presentation_basis_pml: Sequence[int],
    point_operation_count: int,
) -> tuple[tuple[int, ...], ...]:
    """Assign every supplied mode row to one canonical physical atom."""

    if isinstance(point_operation_count, bool) or not isinstance(
        point_operation_count, Integral
    ):
        raise TypeError("point-operation count must be an exact integer")
    if point_operation_count <= 0:
        raise ValueError("point-operation count must be positive")
    operation_count = int(point_operation_count)
    canonical = tuple(
        exact_operation_record(record, point_operation_count=operation_count)
        for record in canonical_records
    )
    rows = tuple(
        exact_operation_record(record, point_operation_count=operation_count)
        for record in row_records
    )
    if len(presentation_basis_pml) != 9:
        raise ValueError(
            f"presentation basis requires 9 entries, got {len(presentation_basis_pml)}"
        )
    if any(
        isinstance(value, bool) or not isinstance(value, Integral)
        for value in presentation_basis_pml
    ):
        raise TypeError("presentation basis must contain exact integers")
    basis_values = tuple(int(value) for value in presentation_basis_pml)
    basis = tuple(
        tuple(Fraction(basis_values[3 * row + column]) for column in range(3))
        for row in range(3)
    )
    basis_inverse = fraction_matrix_inverse3(basis)

    def occurrence_key(
        record: OperationRecord,
    ) -> tuple[int, tuple[Fraction, Fraction, Fraction]]:
        translation = tuple(Fraction(record[axis], record[3]) for axis in range(3))
        coordinates = tuple(
            sum(translation[row] * basis_inverse[row][column] for row in range(3))
            for column in range(3)
        )
        return record[4], tuple(value % 1 for value in coordinates)  # type: ignore[return-value]

    canonical_by_key: dict[tuple[int, tuple[Fraction, Fraction, Fraction]], int] = {}
    for canonical_index, record in enumerate(canonical):
        key = occurrence_key(record)
        previous = canonical_by_key.setdefault(key, canonical_index)
        if previous != canonical_index:
            raise ValueError(
                "canonical child atoms share one affine occurrence identity: "
                f"{previous}, {canonical_index}"
            )

    assignments: list[list[int]] = [[] for _record in canonical]
    for row_index, record in enumerate(rows):
        canonical_index = canonical_by_key.get(occurrence_key(record))
        if canonical_index is None:
            raise ValueError(
                "mode row did not identify one canonical child atom: "
                f"row={row_index}, record={record}"
            )
        assignments[canonical_index].append(row_index)
    return tuple(tuple(indices) for indices in assignments)


def child_atom_layout_from_formula_sites(
    formula_sites: Sequence[FormulaChildSite],
    *,
    label_prefix: str,
) -> ChildAtomLayout:
    """Build one orbit-major layout from a complete Formula15 partition."""

    label_prefix = _exact_nonempty_string(label_prefix, "child atom label prefix")
    if not formula_sites:
        raise ValueError("Formula15 produced no child sites")
    for formula_site in formula_sites:
        if not isinstance(formula_site, FormulaChildSite):
            raise TypeError("Formula15 child site has an invalid type")
        _exact_nonempty_string(formula_site.formula_site, "Formula15 Wyckoff site")
        representative_index = _exact_nonnegative_integer(
            formula_site.representative_atom_index,
            "Formula15 representative atom index",
        )
        if not formula_site.atoms or representative_index >= len(formula_site.atoms):
            raise ValueError("Formula15 representative atom index is out of range")
        for source_atom in formula_site.atoms:
            if not isinstance(source_atom, FormulaChildAtom):
                raise TypeError("Formula15 child atom has an invalid type")
            _exact_nonnegative_integer(
                source_atom.source_flat_index, "Formula15 flat index"
            )
            _exact_nonnegative_integer(
                source_atom.kernel_fraction_ordinal,
                "Formula15 kernel fraction ordinal",
            )
            kernel_fraction = _exact_integer_record(
                source_atom.kernel_fraction, 4, "Formula15 kernel fraction"
            )
            _exact_nonnegative_integer(
                source_atom.parent_branch_ordinal,
                "Formula15 parent branch ordinal",
            )
            parent_coset = _exact_integer_record(
                source_atom.parent_coset_record, 5, "Formula15 parent coset"
            )
            _exact_nonnegative_integer(
                source_atom.centering_ordinal, "Formula15 centering ordinal"
            )
            if kernel_fraction[3] <= 0:
                raise ValueError(
                    "Formula15 kernel fraction denominator must be positive"
                )
            if parent_coset[3] <= 0 or parent_coset[4] <= 0:
                raise ValueError("Formula15 parent coset record is invalid")
            for field, point in (
                ("Formula15 parent point", source_atom.parent_xyz),
                ("Formula15 child point", source_atom.xyz),
            ):
                if (
                    not isinstance(point, tuple)
                    or len(point) != 3
                    or any(not isinstance(value, Fraction) for value in point)
                ):
                    raise TypeError(f"{field} must contain exact fractions")
    base_atoms = [
        atom
        for site in formula_sites
        for atom in site.atoms
        if atom.centering_ordinal == 0
    ]
    kernel_ordinals = {atom.kernel_fraction_ordinal for atom in base_atoms}
    parent_ordinals = {atom.parent_branch_ordinal for atom in base_atoms}
    if not kernel_ordinals or not parent_ordinals:
        raise ValueError("Formula15 child sites contain no base occurrences")
    kernel_count = max(kernel_ordinals) + 1
    parent_count = max(parent_ordinals) + 1
    if kernel_ordinals != set(range(kernel_count)) or parent_ordinals != set(
        range(parent_count)
    ):
        raise ValueError("Formula15 occurrence ordinals are not contiguous")
    occurrence_pairs = {
        (atom.kernel_fraction_ordinal, atom.parent_branch_ordinal)
        for atom in base_atoms
    }
    if occurrence_pairs != {
        (kernel_ordinal, parent_ordinal)
        for kernel_ordinal in range(kernel_count)
        for parent_ordinal in range(parent_count)
    }:
        raise ValueError("Formula15 does not cover the occurrence product")
    if any(
        atom.source_flat_index
        != atom.kernel_fraction_ordinal * parent_count + atom.parent_branch_ordinal
        for atom in base_atoms
    ):
        raise ValueError("Formula15 flat occurrence order is inconsistent")
    cosets_by_parent: dict[int, tuple[int, int, int, int, int]] = {}
    for atom in base_atoms:
        previous = cosets_by_parent.setdefault(
            atom.parent_branch_ordinal,
            atom.parent_coset_record,
        )
        if previous != atom.parent_coset_record:
            raise ValueError("Formula15 parent coset identity is inconsistent")

    atoms: list[ChildAtom] = []
    sites: list[ChildSite] = []
    seen_atom_ids: set[str] = set()
    for child_site_order, formula_site in enumerate(formula_sites):
        multiplicity = re.fullmatch(r"\s*(\d+)\s*\S+\s*", formula_site.formula_site)
        if multiplicity is None or int(multiplicity.group(1)) != len(
            formula_site.atoms
        ):
            raise ValueError(
                "Formula15 child-site membership differs from its multiplicity"
            )
        child_site_id = f"{label_prefix}_{child_site_order + 1}"
        site_atoms: list[ChildAtom] = []
        representative_atom_id: str | None = None
        for member_order, source_atom in enumerate(formula_site.atoms):
            kernel_ordinal = int(source_atom.kernel_fraction_ordinal)
            parent_ordinal = int(source_atom.parent_branch_ordinal)
            centering_ordinal = int(source_atom.centering_ordinal)
            source_raw_index = parent_ordinal * kernel_count + kernel_ordinal
            atom_id = f"{label_prefix}:{source_raw_index}:{centering_ordinal}"
            if atom_id in seen_atom_ids:
                raise ValueError("Formula15 produced duplicate child atom identity")
            seen_atom_ids.add(atom_id)
            atom = ChildAtom(
                atom_id=atom_id,
                child_site_id=child_site_id,
                child_site_order=child_site_order,
                member_order=member_order,
                source_raw_index=source_raw_index,
                source_flat_index=int(source_atom.source_flat_index),
                source_kernel_fraction=tuple(
                    int(value) for value in source_atom.kernel_fraction
                ),
                source_parent_coset_record=tuple(
                    int(value) for value in source_atom.parent_coset_record
                ),
                centering_ordinal=centering_ordinal,
                parent_xyz=source_atom.parent_xyz,
                structure_xyz=source_atom.xyz,
            )
            if member_order == int(formula_site.representative_atom_index):
                representative_atom_id = atom_id
            site_atoms.append(atom)
        if representative_atom_id is None or not site_atoms:
            raise ValueError("Formula15 child site has no representative atom")
        atoms.extend(site_atoms)
        sites.append(
            ChildSite(
                child_site_id=child_site_id,
                child_site_order=child_site_order,
                wyckoff_site=formula_site.formula_site,
                representative_atom_id=representative_atom_id,
                representative_xyz=site_atoms[
                    int(formula_site.representative_atom_index)
                ].structure_xyz,
                atom_ids=tuple(atom.atom_id for atom in site_atoms),
            )
        )
    if len(atoms) != sum(len(site.atom_ids) for site in sites):
        raise ValueError("Formula15 child sites do not partition child atoms")
    return ChildAtomLayout(sites=tuple(sites), atoms=tuple(atoms))


def child_atom_layout_in_presentation_order(
    layout: ChildAtomLayout,
    *,
    label_prefix: str,
    atom_ids: Sequence[str],
    atom_positions: Mapping[str, Sequence[Fraction | int | float]],
) -> ChildAtomLayout:
    """Order and name complete child sites from the Source display atom table."""

    label_prefix = _exact_nonempty_string(label_prefix, "child atom label prefix")
    atom_by_id = {atom.atom_id: atom for atom in layout.atoms}
    ordered_ids = tuple(
        _exact_nonempty_string(atom_id, "Source display atom identity")
        for atom_id in atom_ids
    )
    if (
        len(atom_by_id) != len(layout.atoms)
        or len(ordered_ids) != len(layout.atoms)
        or len(set(ordered_ids)) != len(ordered_ids)
        or set(ordered_ids) != set(atom_by_id)
    ):
        raise ValueError(
            "Source display atom order is not a permutation of the child layout"
        )
    order_by_id = {atom_id: order for order, atom_id in enumerate(ordered_ids)}
    positions = {}
    for atom_id, xyz in atom_positions.items():
        exact_atom_id = _exact_nonempty_string(
            atom_id, "Source display atom-position identity"
        )
        if not isinstance(xyz, (list, tuple)):
            raise TypeError("Source display atom position must be a sequence")
        positions[exact_atom_id] = tuple(_input_fraction(value) % 1 for value in xyz)
    if set(positions) != set(atom_by_id) or any(
        len(xyz) != 3 for xyz in positions.values()
    ):
        raise ValueError("Source display atom positions do not cover the child layout")
    ordered_sites = sorted(
        layout.sites,
        key=lambda site: min(order_by_id[atom_id] for atom_id in site.atom_ids),
    )

    sites: list[ChildSite] = []
    atoms: list[ChildAtom] = []
    for child_site_order, source_site in enumerate(ordered_sites):
        member_ids = tuple(sorted(source_site.atom_ids, key=order_by_id.__getitem__))
        child_site_id = f"{label_prefix}_{child_site_order + 1}"
        site_atoms = tuple(
            replace(
                atom_by_id[atom_id],
                child_site_id=child_site_id,
                child_site_order=child_site_order,
                member_order=member_order,
                structure_xyz=positions[atom_id],
            )
            for member_order, atom_id in enumerate(member_ids)
        )
        atoms.extend(site_atoms)
        sites.append(
            replace(
                source_site,
                child_site_id=child_site_id,
                child_site_order=child_site_order,
                representative_xyz=positions[source_site.representative_atom_id],
                atom_ids=member_ids,
            )
        )
    return ChildAtomLayout(
        sites=tuple(sites),
        atoms=tuple(atoms),
        presentation_rows=layout.presentation_rows,
        point_operation_count=layout.point_operation_count,
    )


def regroup_child_atom_layout(
    layout: ChildAtomLayout,
    *,
    label_prefix: str,
    groups: Sequence[tuple[Sequence[int], str, Sequence[Fraction | int]]],
) -> ChildAtomLayout:
    """Regroup complete ordinary sites while preserving every atom identity."""

    label_prefix = _exact_nonempty_string(label_prefix, "child atom label prefix")
    if not layout.sites or not layout.atoms or not groups:
        raise ValueError("child atom regrouping requires a complete layout")
    atom_by_id = {atom.atom_id: atom for atom in layout.atoms}
    if len(atom_by_id) != len(layout.atoms):
        raise ValueError("child atom layout contains duplicate atom identities")
    used_sites: set[int] = set()
    used_atoms: set[str] = set()
    sites: list[ChildSite] = []
    atoms: list[ChildAtom] = []
    for child_site_order, (raw_indices, wyckoff_site, raw_xyz) in enumerate(groups):
        if not isinstance(raw_indices, (list, tuple)):
            raise TypeError("ordinary-site regrouping indices must be a sequence")
        indices = tuple(
            _exact_nonnegative_integer(index, "ordinary-site regrouping index")
            for index in raw_indices
        )
        exact_wyckoff_site = _exact_nonempty_string(
            wyckoff_site, "magnetic Wyckoff site"
        )
        if (
            not indices
            or len(set(indices)) != len(indices)
            or any(not 0 <= index < len(layout.sites) for index in indices)
            or any(index in used_sites for index in indices)
        ):
            raise ValueError(f"invalid ordinary-site regrouping: {raw_indices!r}")
        if not isinstance(raw_xyz, (list, tuple)) or len(raw_xyz) != 3:
            raise ValueError(f"invalid child-site representative: {raw_xyz!r}")
        if any(
            isinstance(value, bool) or not isinstance(value, (Fraction, Integral))
            for value in raw_xyz
        ):
            raise TypeError("child-site representative must contain exact numbers")
        representative_xyz = tuple(Fraction(value) % 1 for value in raw_xyz)
        source_sites = tuple(layout.sites[index] for index in indices)
        member_ids = tuple(
            atom_id for site in source_sites for atom_id in site.atom_ids
        )
        if (
            not member_ids
            or len(set(member_ids)) != len(member_ids)
            or any(atom_id not in atom_by_id for atom_id in member_ids)
            or any(atom_id in used_atoms for atom_id in member_ids)
        ):
            raise ValueError("magnetic child site does not partition canonical atoms")
        representative_ids = tuple(
            atom_id
            for atom_id in member_ids
            if tuple(value % 1 for value in atom_by_id[atom_id].structure_xyz)
            == representative_xyz
        )
        if len(representative_ids) != 1:
            raise ValueError(
                "magnetic child-site representative is not one canonical atom"
            )
        child_site_id = f"{label_prefix}_{child_site_order + 1}"
        regrouped_atoms = tuple(
            replace(
                atom_by_id[atom_id],
                child_site_id=child_site_id,
                child_site_order=child_site_order,
                member_order=member_order,
            )
            for member_order, atom_id in enumerate(member_ids)
        )
        used_sites.update(indices)
        used_atoms.update(member_ids)
        atoms.extend(regrouped_atoms)
        sites.append(
            ChildSite(
                child_site_id=child_site_id,
                child_site_order=child_site_order,
                wyckoff_site=exact_wyckoff_site,
                representative_atom_id=representative_ids[0],
                representative_xyz=representative_xyz,  # type: ignore[arg-type]
                atom_ids=member_ids,
            )
        )
    if used_sites != set(range(len(layout.sites))) or used_atoms != set(atom_by_id):
        raise ValueError("child-site regrouping is not a complete partition")
    return ChildAtomLayout(
        sites=tuple(sites),
        atoms=tuple(atoms),
        presentation_rows=layout.presentation_rows,
        point_operation_count=layout.point_operation_count,
    )
