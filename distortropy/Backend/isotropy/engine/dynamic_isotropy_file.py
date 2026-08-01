"""Dynamic isotropy rows and their Source emission order."""

from __future__ import annotations

from dataclasses import dataclass, field
from fractions import Fraction
import re
from typing import Iterable

LABEL_RE = re.compile(r"^(?P<prefix>P|C|S|(?P<dim>\d+)D)(?P<index>\d+)$")


@dataclass(frozen=True)
class DynamicIsotropyRow:
    subgroup_number: int
    free: int
    basis_values: tuple[int, ...]
    origin_values: tuple[int, ...]
    direction: str
    matrix: tuple[tuple[str, ...], ...]
    magnetic_subgroup_selection: object | None = field(
        default=None,
        compare=False,
        repr=False,
    )
    magnetic_operation_records: tuple[tuple[int, int, int, int, int], ...] = field(
        default=(),
        compare=False,
        repr=False,
    )

    @property
    def basis(self) -> str:
        return _basis(list(self.basis_values))

    @property
    def origin(self) -> str:
        return _origin(list(self.origin_values))

    @property
    def opd(self) -> str:
        return _opd_from_matrix([list(row) for row in self.matrix])

def dynamic_row_write_order_key(row: DynamicIsotropyRow) -> tuple[int, int]:
    """Return the final ``*.iso`` ordering key used by ``get_isotropy_``.

    Candidate discovery order is normalized to the file emission order:
    rows are grouped by OPD dimension (`P`, `C`, `S`, `nD`) and occurrence
    index.
    """

    label = row.direction.strip()
    match = LABEL_RE.match(label)
    if match:
        return (row.free, int(match.group("index")))
    return (row.free, 0)


def sort_dynamic_rows_for_file(rows: Iterable[DynamicIsotropyRow]) -> tuple[DynamicIsotropyRow, ...]:
    return tuple(sorted(rows, key=dynamic_row_write_order_key))


def _opd_from_matrix(matrix: list[list[str]]) -> str:
    params = "abcdefghijklmnopqrstuvwxyz"
    if not matrix:
        return "()"
    width = max(len(row) for row in matrix)
    comps: list[str] = []
    for col in range(width):
        expr = ""
        for row_index, row in enumerate(matrix):
            coef = row[col] if col < len(row) else "0"
            if coef == "0":
                continue
            negative = coef.startswith("-")
            mag = coef[1:] if negative else coef
            body = params[row_index] if mag == "1" else f"{mag}{params[row_index]}"
            if negative:
                term = f"-{body}"
            else:
                term = body if not expr else f"+{body}"
            expr += term
        comps.append(expr or "0")
    return "(" + ",".join(comps) + ")"


def _origin(values: list[int]) -> str:
    den = values[3] or 1
    out = []
    for value in values[:3]:
        frac = Fraction(value, den)
        out.append(str(frac.numerator) if frac.denominator == 1 else f"{frac.numerator}/{frac.denominator}")
    return "(" + ",".join(out) + ")"


def _basis(values: list[int]) -> str:
    rows = [values[index:index + 3] for index in range(0, 9, 3)]
    return ",".join("(" + ",".join(str(item) for item in row) + ")" for row in rows)
