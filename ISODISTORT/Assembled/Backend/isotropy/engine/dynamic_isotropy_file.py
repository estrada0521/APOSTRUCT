#!/usr/bin/env python3
"""Parse upstream dynamic isotropy subgroup ``*.iso`` database files.

These files are materialized by ``Source/iso`` when ``DISPLAY ISOTROPY`` needs
OPD rows that are not already present in static tables.  They are the output of
the ``get_isotropy_`` generation layer before ``isotropy_display_isotropy_``
applies display-setting conversion and public formatting.
"""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass, field
from fractions import Fraction
import json
import re
from pathlib import Path
from typing import Iterable

HEADER_RE = re.compile(r"Space Group\s+(\d+),\s+Irrep\s+([^,]+?)(?:,\s+kparam\s+(.+))?$")
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

    def to_jsonable(self) -> dict[str, object]:
        data = asdict(self)
        data.pop("magnetic_subgroup_selection", None)
        data.pop("magnetic_operation_records", None)
        data["basis"] = self.basis
        data["origin"] = self.origin
        data["opd"] = self.opd
        return data


@dataclass(frozen=True)
class DynamicIsotropyFile:
    path: str
    space_group: int
    irrep: str
    kparam: str
    raw_kparam: tuple[int, ...]
    rows: tuple[DynamicIsotropyRow, ...]

    def to_jsonable(self) -> dict[str, object]:
        return {
            "schema": "isotropy.disassembled.dynamic_isotropy_file.v1",
            "path": self.path,
            "space_group": self.space_group,
            "irrep": self.irrep,
            "kparam": self.kparam,
            "raw_kparam": list(self.raw_kparam),
            "rows": [row.to_jsonable() for row in self.rows],
        }


def dynamic_row_write_order_key(row: DynamicIsotropyRow) -> tuple[int, int]:
    """Return the final ``*.iso`` ordering key used by ``get_isotropy_``.

    GDB row-write probes observe internal search order.  The binary writes
    those candidates to a temporary file, then rewinds it and emits final rows
    grouped by OPD dimension (`P`, `C`, `S`, `nD`) and occurrence index.
    """

    label = row.direction.strip()
    match = LABEL_RE.match(label)
    if match:
        return (row.free, int(match.group("index")))
    return (row.free, 0)


def sort_dynamic_rows_for_file(rows: Iterable[DynamicIsotropyRow]) -> tuple[DynamicIsotropyRow, ...]:
    return tuple(sorted(rows, key=dynamic_row_write_order_key))


def _fmt_number(text: str) -> str:
    value = float(text)
    if abs(value) < 5e-7:
        value = 0.0
    if abs(value - round(value)) < 5e-7:
        return str(int(round(value)))
    return f"{value:.6f}".rstrip("0").rstrip(".")


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


def _header_values(line: str) -> list[int] | None:
    parts = line.split()
    if len(parts) != 15:
        return None
    try:
        return [int(round(float(value))) for value in parts]
    except ValueError:
        return None


def parse_dynamic_isotropy_file(path: Path) -> DynamicIsotropyFile:
    lines = [line.rstrip() for line in Path(path).read_text(errors="replace").splitlines()]
    if not lines:
        raise ValueError(f"empty isotropy file: {path}")
    header = HEADER_RE.match(lines[0].strip())
    if not header:
        raise ValueError(f"unexpected isotropy header in {path}: {lines[0]!r}")
    params: list[int] = []
    index = 1
    if index < len(lines):
        maybe_params = lines[index].split()
        try:
            maybe_ints = [int(value) for value in maybe_params]
        except ValueError:
            maybe_ints = []
        # Parametric-k files carry one metadata line before row headers.  Its
        # length depends on k-parameter dimensionality, while row headers are
        # always the 15-int records parsed by _header_values.
        if maybe_ints and _header_values(lines[index].strip()) is None:
            params = maybe_ints
            index += 1
    rows: list[DynamicIsotropyRow] = []
    while index < len(lines):
        line = lines[index].strip()
        if not line or line == "done":
            index += 1
            continue
        values = _header_values(line)
        if values is None:
            break
        subgroup = values[0]
        free = values[1]
        basis_values = values[2:11]
        origin_values = values[11:15]
        index += 1
        direction = lines[index].strip()
        index += 1
        matrix_values: list[str] = []
        next_header: list[int] | None = None
        while index < len(lines):
            next_line = lines[index].strip()
            if not next_line:
                index += 1
                continue
            if next_line == "done":
                break
            next_header = _header_values(next_line)
            if next_header is not None:
                break
            matrix_values.extend(_fmt_number(value) for value in next_line.split())
            index += 1
        width = len(matrix_values) // free if free else 0
        matrix: list[tuple[str, ...]] = []
        for row_index in range(free):
            start = row_index * width
            matrix.append(tuple(matrix_values[start:start + width]))
        rows.append(DynamicIsotropyRow(
            subgroup_number=subgroup,
            free=free,
            basis_values=tuple(basis_values),
            origin_values=tuple(origin_values),
            direction=direction,
            matrix=tuple(matrix),
        ))
    return DynamicIsotropyFile(
        path=str(path),
        space_group=int(header.group(1)),
        irrep=header.group(2).strip(),
        kparam=(header.group(3) or "").strip(),
        raw_kparam=tuple(params),
        rows=tuple(rows),
    )


def _main(argv: Iterable[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("path", type=Path)
    parser.add_argument("--indent", type=int, default=2)
    args = parser.parse_args(list(argv) if argv is not None else None)
    parsed = parse_dynamic_isotropy_file(args.path)
    print(json.dumps(parsed.to_jsonable(), ensure_ascii=False, indent=args.indent, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(_main())
