"""Supercell atom expansion and final atom ordering helpers."""

from __future__ import annotations

from fractions import Fraction
from typing import Iterable

from ISODISTORT.Assembled.Backend.exactmath import (
    fraction_identity3,
    fraction_matrix_inverse3,
)
from ISODISTORT.Assembled.Backend.modes.engine.input import Case
from ISODISTORT.Assembled.Backend.modes.engine.records import WyckoffRow

class IsotropyAtomExpansionMixin:
    def wyckoff_multiplicity(self, sg: int, row: WyckoffRow) -> int:
        """Return the generated conventional Wyckoff orbit multiplicity."""

        site_records = self.wyc_pg_elements_records(sg, row)
        return len(self.wyc_pg_cosets_records(sg, site_records))

    def supercell_atom_count(self, case: Case, row: WyckoffRow) -> int:
        """Return the conventional-orbit multiplicity times the k-cell index."""

        return self.wyckoff_multiplicity(case.sg, row) * self.k_supercell_index(case)

    def display_distortion_atom_fractionals(
        self,
        sg: int,
        row: WyckoffRow,
        site_params: Iterable[object] | None,
        basis_pml: Iterable[int],
    ) -> tuple[tuple[Fraction, Fraction, Fraction], ...]:
        """Construct DISPLAY DISTORTION's printed Point column.

        The operation/base part is an affine coordinate and the free Wyckoff
        parameter part is a direct vector.  They therefore cross the setting
        boundary by different Source transforms and are added only at the
        display edge.
        """

        matrix = tuple(int(value) for value in basis_pml)
        translations = self.get_new_fractionals_from_basis_columns(matrix)
        base_vectors = self.wyckoff_fraction_vectors(row)
        site_records = self.wyc_pg_elements_records(int(sg), row)
        cosets = self.wyc_pg_cosets_records(int(sg), site_records)
        params = tuple(Fraction(str(value)) for value in (site_params or ()))
        pml_to_cml = self.pml_to_cml_matrix(int(sg))
        cml_to_cinter = self.cml_to_cinter_matrix(int(sg))
        out: list[tuple[Fraction, Fraction, Fraction]] = []
        for operation in cosets:
            branch = tuple(
                self.vrot_fraction(int(sg), int(operation[4]), vector)
                for vector in base_vectors
            )
            internal_parameter = tuple(
                sum(
                    (params[index] if index < len(params) else Fraction(0))
                    * branch[index + 1][axis]
                    for index in range(3)
                )
                for axis in range(3)
            )
            cml_parameter = tuple(
                sum(internal_parameter[row_index] * pml_to_cml[row_index][axis] for row_index in range(3))
                for axis in range(3)
            )
            display_parameter = tuple(
                sum(cml_parameter[row_index] * cml_to_cinter[row_index][axis] for row_index in range(3))
                for axis in range(3)
            )
            branch_base = self.iso.fraction_operation_record(
                branch[0], int(operation[4])
            )
            for translation in translations:
                combined = self.add_translation_to_operation_record(translation, operation)
                position_record = self.add_translation_to_operation_record(branch_base[:4], combined)
                internal_base = tuple(
                    Fraction(int(position_record[axis]), int(position_record[3]))
                    for axis in range(3)
                )
                display_base = self.xyz_change_setting_point(int(sg), "pml", "cinter", internal_base)
                out.append(
                    tuple(display_base[axis] + display_parameter[axis] for axis in range(3))
                )
        return tuple(out)

    def supercell_atom_operation_records(
        self,
        sg: int,
        row: WyckoffRow,
        basis: Iterable[int],
        site_params: Iterable[object] | None = None,
    ) -> tuple[tuple[int, int, int, int, int], ...]:
        """Return operation records in the same order as atom-table positions."""

        matrix = tuple(int(x) for x in basis)
        site_records = self.wyc_pg_elements_records(sg, row)
        cosets = self.wyc_pg_cosets_records(sg, site_records)
        translations = self.get_new_fractionals_from_basis_columns(matrix)
        out: list[tuple[int, int, int, int, int]] = []
        for coset in cosets:
            for tx, ty, tz, tden in translations:
                translation_record = (int(tx), int(ty), int(tz), int(tden))
                out.append(self.add_translation_to_operation_record(translation_record, coset))
        return tuple(out)

    def _atom_positions_duplicate_after_fold(
        self,
        representative: tuple[Fraction, Fraction, Fraction],
        cosets: Iterable[tuple[int, int, int, int, int]],
        matrix: tuple[int, ...],
        inverse: tuple[tuple[Fraction, Fraction, Fraction], ...],
        translations: Iterable[tuple[int, int, int, int]],
    ) -> bool:
        for coset in cosets:
            seen: set[tuple[Fraction, Fraction, Fraction]] = set()
            denominator = int(coset[3])
            tau = tuple(Fraction(int(coset[axis]), denominator) for axis in range(3))
            point = self.iso.point_ops[int(coset[4]) - 1]
            for tx, ty, tz, tden in translations:
                translation = (Fraction(tx, tden), Fraction(ty, tden), Fraction(tz, tden))
                position = tuple(
                    sum(Fraction(int(point[axis, col])) * representative[col] for col in range(3))
                    + tau[axis]
                    + translation[axis]
                    for axis in range(3)
                )
                raw_super_coords = tuple(
                    sum(position[col] * inverse[col][axis] for col in range(3))
                    for axis in range(3)
                )
                if all(Fraction(0) <= value < Fraction(1) for value in raw_super_coords):
                    displayed = position
                else:
                    super_coords = tuple(value % 1 for value in raw_super_coords)
                    displayed = tuple(
                        sum(super_coords[col] * Fraction(matrix[3 * col + axis]) for col in range(3))
                        for axis in range(3)
                    )
                if displayed in seen:
                    return True
                seen.add(displayed)
        return False

    def atom_expansion_uses_column_fallback(
        self,
        sg: int,
        row: WyckoffRow,
        site_params: Iterable[object] | None,
        basis: Iterable[int],
    ) -> bool:
        matrix = tuple(int(x) for x in basis)
        if self.pml_to_cinter_matrix(sg) != fraction_identity3():
            return False
        inverse = fraction_matrix_inverse3(
            tuple(
                tuple(Fraction(matrix[3 * row_index + col]) for col in range(3))
                for row_index in range(3)
            )
        )
        representative = self.wyckoff_representative(row, site_params)
        site_records = self.wyc_pg_elements_records(sg, row)
        cosets = self.wyc_pg_cosets_records(sg, site_records)
        translations = self.get_new_fractionals_from_basis(matrix)
        return self._atom_positions_duplicate_after_fold(representative, cosets, matrix, inverse, translations)

    @staticmethod
    def final_atom_order_indices(
        operation_records: Iterable[tuple[int, int, int, int, int]],
    ) -> tuple[int, ...]:
        """Return the display-order permutation for final atom/mode tables.

        Keep this separate from ``wyc_pg_cosets_`` and
        ``supercell_atom_operation_records`` because those are routine-boundary
        data paths. This implementation preserves the input operation order;
        SG126/e/DT and SG223/l/SM are controls for that behavior, not a general
        Source ordering theorem.
        """

        records = tuple(operation_records)
        return tuple(range(len(records)))

    @staticmethod
    def final_atom_position_order_indices(
        positions: Iterable[tuple[Fraction, Fraction, Fraction]],
        order: Iterable[int],
    ) -> tuple[int, ...]:
        """Apply final display-order corrections that depend on folded positions."""

        return tuple(int(item) for item in order)
