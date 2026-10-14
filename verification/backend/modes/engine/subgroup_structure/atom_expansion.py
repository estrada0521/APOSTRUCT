"""Supercell atom expansion and final atom ordering helpers."""

from __future__ import annotations

from fractions import Fraction
from functools import lru_cache
from typing import Iterable

from verification.backend import native
from verification.backend.modes.engine.records import WyckoffRow


@lru_cache(maxsize=1024)
def _site_atom_table(
    sg: int,
    row_id: int,
    site_pg: int,
    params: tuple[Fraction, ...],
    basis: tuple[int, ...],
    representative: tuple[Fraction, ...] | None,
):
    return native.site_atom_table(sg, row_id, site_pg, params, basis, representative)


class IsotropyAtomExpansionMixin:
    def wyckoff_multiplicity(self, sg: int, row: WyckoffRow) -> int:
        """Return the generated conventional Wyckoff orbit multiplicity."""

        site_records = self.wyc_pg_elements_records(sg, row)
        return len(self.wyc_pg_cosets_records(sg, site_records))

    def site_atom_table(
        self,
        sg: int,
        row: WyckoffRow,
        site_params: Iterable[object] | None,
        basis_pml: Iterable[int],
        representative_pml: tuple[Fraction, Fraction, Fraction] | None = None,
    ) -> tuple[
        tuple[tuple[Fraction, Fraction, Fraction], ...],
        tuple[tuple[int, int, int, int, int], ...],
        tuple[int, ...] | None,
    ]:
        """Points, operation records and (given the parent representative)
        the rebased display order of one parent site's supercell atoms."""

        return _site_atom_table(
            int(sg),
            int(row.row_id),
            int(row.site_pg),
            tuple(Fraction(str(value)) for value in (site_params or ())),
            tuple(int(value) for value in basis_pml),
            None if representative_pml is None else tuple(Fraction(value) for value in representative_pml),
        )
