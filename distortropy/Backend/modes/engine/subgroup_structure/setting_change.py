"""Setting-change and vector-record conversion helpers."""

from __future__ import annotations

from fractions import Fraction
from functools import lru_cache
import math
from typing import Iterable, Sequence

from distortropy.Backend.exactmath import (
    fraction_identity3,
    fraction_matrix_inverse3,
    integer_determinant3,
)
from distortropy.Backend.lattice_quotient import (
    integer_inverse_denominator,
    integral_row_images_source_order,
)
from distortropy.Backend.source.iso_data import (
    pml_to_cml_matrix_from_table,
    setting_change_matrix_from_table,
    setting_to_cinter_affine_from_table,
)


@lru_cache(maxsize=256)
def _new_fractionals_from_basis_columns_cached(
    matrix: tuple[int, ...],
) -> tuple[tuple[int, int, int, int], ...]:
    det = abs(IsotropySettingChangeMixin._det3_row_major(matrix))
    if det <= 0:
        raise ValueError(f"singular basis matrix: {matrix}")
    return integral_row_images_source_order(
        matrix,
        integer_inverse_denominator(matrix),
    )


class IsotropySettingChangeMixin:
    @staticmethod
    def _det3_row_major(matrix: Iterable[int]) -> int:
        m = tuple(int(x) for x in matrix)
        if len(m) != 9:
            raise ValueError(f"expected 9 matrix entries, got {len(m)}")
        return integer_determinant3(m)

    @staticmethod
    def get_new_fractionals_from_basis(matrix: Iterable[int]) -> tuple[tuple[int, int, int, int], ...]:
        """Enumerate `get_new_fractionals_` representatives for a basis matrix.

        `data_isotropy:isotropy_basis` is stored in the same row-major
        convention as the printed superlattice basis. The enumeration scans
        triples over the reduced ``matinv_`` denominator, applies ``P``, reduces the
        rational vector, and keeps the integral cases. The shared exact
        quotient kernel produces the same ordered translations directly.
        """

        m = tuple(int(x) for x in matrix)
        det = abs(IsotropySettingChangeMixin._det3_row_major(m))
        if det <= 0:
            raise ValueError(f"singular basis matrix: {m}")
        transposed = tuple(
            m[3 * col + row]
            for row in range(3)
            for col in range(3)
        )
        return integral_row_images_source_order(
            transposed,
            integer_inverse_denominator(m),
        )

    @staticmethod
    def get_new_fractionals_from_basis_columns(matrix: Iterable[int]) -> tuple[tuple[int, int, int, int], ...]:
        """Return internal translations using the transpose multiplication convention.

        Public superlattice bases are row-major, but final atom/mode expansion
        for identity-setting rows follows the column multiplication
        convention. The row helper above remains available to its row-vector
        callers; this helper is only used at the final expansion boundary.
        """

        return _new_fractionals_from_basis_columns_cached(
            tuple(int(value) for value in matrix)
        )

    def cml_to_cinter_matrix(self, sg: int) -> tuple[tuple[Fraction, Fraction, Fraction], ...]:
        """Return the row-vector transform from cml to cinter for an SG."""

        return setting_to_cinter_affine_from_table(self.iso.space, sg, "cml")[0]

    def pml_to_cml_matrix(self, sg: int) -> tuple[tuple[Fraction, Fraction, Fraction], ...]:
        """Return the row-vector transform from pml to cml for an SG lattice."""

        return pml_to_cml_matrix_from_table(self.iso.space, sg)

    def pml_to_cinter_matrix(self, sg: int) -> tuple[tuple[Fraction, Fraction, Fraction], ...]:
        """Return the row-vector transform used by `vector_change_setting_(pml,cinter)`."""

        return setting_to_cinter_affine_from_table(self.iso.space, sg, "pml")[0]

    def cinter_to_pml_matrix(self, sg: int) -> tuple[tuple[Fraction, Fraction, Fraction], ...]:
        """Return the inverse direct-space transform from cinter to PML."""

        return fraction_matrix_inverse3(self.pml_to_cinter_matrix(int(sg)))

    def cinter_translation_to_pml_record(
        self,
        sg: int,
        translation: Sequence[Fraction | int | float],
    ) -> tuple[int, int, int, int, int]:
        """Represent a cinter lattice translation as an identity PML record."""

        values = tuple(
            Fraction(str(value)) if isinstance(value, float) else Fraction(value)
            for value in translation
        )
        if len(values) != 3:
            raise ValueError(f"expected cinter translation triplet, got {translation!r}")
        matrix = self.cinter_to_pml_matrix(int(sg))
        transformed = tuple(
            sum(values[row] * matrix[row][col] for row in range(3))
            for col in range(3)
        )
        denominator = 1
        for value in transformed:
            denominator = math.lcm(denominator, value.denominator)
        numerators = tuple(int(value * denominator) for value in transformed)
        identity_point_op = int(self.generate_space_group_records(int(sg))[0][4])
        return numerators[0], numerators[1], numerators[2], denominator, identity_point_op

    def setting_change_matrix(
        self,
        sg: int,
        from_setting: str,
        to_setting: str,
    ) -> tuple[tuple[Fraction, Fraction, Fraction], ...]:
        """Return the row-vector matrix used by ``vector_change_setting_``.

        The Source setting catalog supports Miller-Love primitive (``pml``),
        Miller-Love conventional (``cml``), and internal conventional
        (``cinter``) coordinates.  All matrices satisfy ``v_to = v_from * M``.
        """

        return setting_change_matrix_from_table(
            self.iso.space,
            sg,
            from_setting,
            to_setting,
        )

    def reciprocal_setting_change_vector(
        self,
        sg: int,
        from_setting: str,
        to_setting: str,
        vector: Iterable[Fraction],
    ) -> tuple[Fraction, Fraction, Fraction]:
        """Transform reciprocal-coordinate row vectors between settings.

        ``vector_change_setting_`` acts on direct-space row coordinates as
        ``v_to = v_from * M``.  mode-kernel prints k vectors in the reciprocal basis
        of the conventional setting, so the corresponding coordinate transform
        is ``k_to = k_from * inverse(M).T``.
        """

        values = tuple(Fraction(value) for value in vector)
        if len(values) != 3:
            raise ValueError(f"expected 3 reciprocal coordinates, got {values}")
        direct = self.setting_change_matrix(sg, from_setting, to_setting)
        inverse = fraction_matrix_inverse3(direct)
        inverse_transpose = tuple(
            tuple(inverse[col][row] for col in range(3))
            for row in range(3)
        )
        return tuple(
            sum(values[row] * inverse_transpose[row][col] for row in range(3))
            for col in range(3)
        )  # type: ignore[return-value]

    def _cml_to_cinter_origin(self, sg: int) -> tuple[Fraction, Fraction, Fraction]:
        return setting_to_cinter_affine_from_table(self.iso.space, sg, "cml")[1]

    def _point_to_cinter_affine(
        self,
        sg: int,
        setting: str,
    ) -> tuple[tuple[tuple[Fraction, Fraction, Fraction], ...], tuple[Fraction, Fraction, Fraction]]:
        setting = setting.strip().lower()
        identity = fraction_identity3()
        zero = (Fraction(0), Fraction(0), Fraction(0))
        if setting == "cinter":
            return identity, zero
        if setting == "cml":
            return self.cml_to_cinter_matrix(sg), self._cml_to_cinter_origin(sg)
        if setting == "pml":
            return self.pml_to_cinter_matrix(sg), self._cml_to_cinter_origin(sg)
        raise KeyError(f"unsupported point setting {setting!r}")

    def xyz_change_setting_point(
        self,
        sg: int,
        from_setting: str,
        to_setting: str,
        point: Iterable[Fraction],
    ) -> tuple[Fraction, Fraction, Fraction]:
        """Apply the affine point-coordinate part of ``xyz_change_setting_``."""

        values = tuple(point)
        if len(values) != 3:
            raise ValueError(f"expected 3 point coordinates, got {values}")
        from_matrix, from_origin = self._point_to_cinter_affine(sg, from_setting)
        to_matrix, to_origin = self._point_to_cinter_affine(sg, to_setting)
        cinter = tuple(
            sum(values[row] * from_matrix[row][col] for row in range(3)) + from_origin[col]
            for col in range(3)
        )
        shifted = tuple(cinter[axis] - to_origin[axis] for axis in range(3))
        inverse_to = fraction_matrix_inverse3(to_matrix)
        return tuple(
            sum(shifted[row] * inverse_to[row][col] for row in range(3))
            for col in range(3)
        )  # type: ignore[return-value]

    @staticmethod
    def _fraction_vector_to_record(
        vector: Iterable[Fraction],
    ) -> tuple[int, int, int, int]:
        values = tuple(vector)
        if len(values) != 3:
            raise ValueError(f"expected 3-vector, got {values}")
        den = 1
        for value in values:
            den = math.lcm(den, value.denominator)
        nums = [int(value * den) for value in values]
        gcd = den
        for num in nums:
            gcd = math.gcd(gcd, abs(num))
        if gcd > 1:
            nums = [num // gcd for num in nums]
            den //= gcd
        if den < 0:
            nums = [-num for num in nums]
            den = -den
        return (nums[0], nums[1], nums[2], den)

    def vector_change_setting_record(
        self,
        sg: int,
        from_setting: str,
        to_setting: str,
        record: Iterable[int],
    ) -> tuple[int, int, int, int]:
        """Apply ``vector_change_setting_`` to one 4-int rational vector."""

        raw = tuple(int(x) for x in record)
        if len(raw) != 4:
            raise ValueError(f"expected 4-int vector record, got {raw}")
        if raw[3] == 0:
            raise ValueError(f"zero vector denominator in setting transform: {raw}")
        vector = tuple(Fraction(raw[axis], raw[3]) for axis in range(3))
        matrix = self.setting_change_matrix(sg, from_setting, to_setting)
        transformed = tuple(
            sum(vector[row] * matrix[row][col] for row in range(3))
            for col in range(3)
        )
        return self._fraction_vector_to_record(transformed)
