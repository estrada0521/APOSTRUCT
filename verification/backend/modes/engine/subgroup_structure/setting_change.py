"""Setting-change and vector-record conversion helpers."""

from __future__ import annotations

from fractions import Fraction

from verification.backend.source.iso_data import (
    pml_to_cml_matrix_from_table,
    setting_change_matrix_from_table,
    setting_to_cinter_affine_from_table,
)


class IsotropySettingChangeMixin:

    def cml_to_cinter_matrix(self, sg: int) -> tuple[tuple[Fraction, Fraction, Fraction], ...]:
        """Return the row-vector transform from cml to cinter for an SG."""

        memo = self.__dict__.setdefault("_cml_to_cinter_matrix", {})
        if sg not in memo:
            memo[sg] = setting_to_cinter_affine_from_table(self.iso.space, sg, "cml")[0]
        return memo[sg]

    def pml_to_cml_matrix(self, sg: int) -> tuple[tuple[Fraction, Fraction, Fraction], ...]:
        """Return the row-vector transform from pml to cml for an SG lattice."""

        return pml_to_cml_matrix_from_table(self.iso.space, sg)

    def pml_to_cinter_matrix(self, sg: int) -> tuple[tuple[Fraction, Fraction, Fraction], ...]:
        """Return the row-vector transform used by `vector_change_setting_(pml,cinter)`."""

        memo = self.__dict__.setdefault("_pml_to_cinter_matrix", {})
        if sg not in memo:
            memo[sg] = setting_to_cinter_affine_from_table(self.iso.space, sg, "pml")[0]
        return memo[sg]

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

        cache = self.__dict__.setdefault("_setting_change_matrix_cache", {})
        key = (int(sg), str(from_setting), str(to_setting))
        if key not in cache:
            cache[key] = setting_change_matrix_from_table(self.iso.space, *key)
        return cache[key]
