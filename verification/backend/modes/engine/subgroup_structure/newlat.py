"""Integer superlattice basis construction for ``newlat_`` records."""

from __future__ import annotations

from functools import lru_cache
from typing import Iterable

from verification.backend import native

class NewlatMixin:
    @staticmethod
    def newlat_basis(
        mode_count: int,
        vector_records: Iterable[Iterable[int]],
    ) -> tuple[int, ...]:
        """Construct the integer lattice basis defined by ``newlat_`` records.

        `newlat_` receives a small list of 4-int rational vectors and constructs
        a 3x3 integer basis whose rows have integral dot product with the first
        `mode_count` vectors.  It then performs a short max-norm reduction by
        adding or subtracting other basis rows.  This routine is used by
        `find_isotropy_subgroup_` for parametric-k superlattices.
        """

        return _newlat_basis(
            int(mode_count),
            tuple(tuple(int(x) for x in record) for record in vector_records),
        )

@lru_cache(maxsize=4096)
def _newlat_basis(
    mode_count: int,
    records: tuple[tuple[int, ...], ...],
) -> tuple[int, ...]:
    return native.newlat_basis(mode_count, records)
