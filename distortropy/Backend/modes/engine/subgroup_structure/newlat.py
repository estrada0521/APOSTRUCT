"""Integer superlattice basis construction for ``newlat_`` records."""

from __future__ import annotations

from fractions import Fraction
from typing import Iterable

from distortropy.Backend.exactmath import fraction_matrix_multiply3

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

        records = tuple(tuple(int(x) for x in record) for record in vector_records)
        if mode_count <= 0:
            raise ValueError(f"newlat mode_count must be positive, got {mode_count}")
        if len(records) < mode_count:
            raise ValueError(f"newlat needs {mode_count} vector records, got {len(records)}")
        max_den = records[0][3]
        if max_den <= 0:
            raise ValueError(f"newlat first denominator must be positive, got {max_den}")

        rows: list[list[int]] = []
        for pivot_index in range(3):
            found: list[int] | None = None
            for pivot in range(1, max_den + 1):
                combo_count = (max_den + 1) ** pivot_index
                for combo_serial in range(combo_count):
                    combo: list[int] = []
                    value = combo_serial
                    for _axis in range(pivot_index):
                        combo.append(value % (max_den + 1))
                        value //= max_den + 1
                    candidate = [0, 0, 0]
                    candidate[pivot_index] = pivot
                    for axis, component in enumerate(combo):
                        candidate[axis] = component
                    ok = True
                    for record in records[:mode_count]:
                        den = record[3]
                        if den == 0:
                            raise ValueError(f"newlat vector has zero denominator: {record}")
                        dot = sum(candidate[axis] * record[axis] for axis in range(pivot_index + 1))
                        if dot % den != 0:
                            ok = False
                            break
                    if ok:
                        found = candidate
                        break
                if found is not None:
                    break
            if found is None:
                raise ValueError(f"newlat failed to find basis row {pivot_index + 1}")
            rows.append(found)

        for _pass_index in range(3):
            for target_index in range(3):
                for other_index in range(3):
                    if other_index == target_index:
                        continue
                    for sign in (-1, 1):
                        candidate = [
                            rows[target_index][axis] + sign * rows[other_index][axis]
                            for axis in range(3)
                        ]
                        if max(abs(value) for value in candidate) < max(abs(value) for value in rows[target_index]):
                            rows[target_index] = candidate

        return tuple(value for row in rows for value in row)

    @staticmethod
    def transform_basis_rows(
        basis: Iterable[int],
        matrix: tuple[tuple[Fraction, Fraction, Fraction], ...],
    ) -> tuple[tuple[Fraction, Fraction, Fraction], ...]:
        rows = tuple(tuple(Fraction(int(x)) for x in tuple(basis)[3 * row:3 * row + 3]) for row in range(3))
        return fraction_matrix_multiply3(rows, matrix)

    @staticmethod
    def format_fraction_matrix(
        matrix: Iterable[Iterable[Fraction]],
    ) -> list[list[str]]:
        return [[str(value) for value in row] for row in matrix]
