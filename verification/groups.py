"""Space and magnetic group operations from the Source tables, computed here from their definitions.

The verifier reads the same convention data as apo (the Source blob, whose sections the gate checks against the BYU
tables) but builds group elements, products and types with its own code, so a mistake in apo's group code cannot
cancel out in a comparison.
"""

from __future__ import annotations

from fractions import Fraction
from functools import lru_cache
from pathlib import Path

from verification.backend.source.compiled import SourceBlob

_BLOB = Path(__file__).resolve().parents[1] / "bin/source.blob"


@lru_cache(maxsize=None)
def _section(name: str) -> tuple:
    return _blob()._section(name)


@lru_cache(maxsize=1)
def _blob() -> SourceBlob:
    return SourceBlob(_BLOB)


def point_matrix(op: int) -> tuple[tuple[int, int, int], ...]:
    """Point operation op (1-72) in the stored order; it acts on a column vector as its transpose (_column_apply)."""

    m = _section("group.point_ops.matrix")[9 * (op - 1):9 * op]
    return tuple(tuple(int(m[3 * i + j]) for j in range(3)) for i in range(3))


def _mul(a, b):
    return tuple(tuple(sum(a[i][k] * b[k][j] for k in range(3)) for j in range(3)) for i in range(3))


@lru_cache(maxsize=None)
def _op_of_matrix(hexagonal: bool) -> dict:
    return {point_matrix(op): op for op in (range(49, 73) if hexagonal else range(1, 49))}


def point_product(left: int, right: int) -> int:
    """The point operation whose matrix is left's times right's (0 across the cubic and hexagonal families)."""

    if (left > 48) != (right > 48):
        return 0
    return _op_of_matrix(left > 48).get(_mul(point_matrix(left), point_matrix(right)), 0)


def _column_apply(m, v):
    return tuple(sum(m[k][axis] * v[k] for k in range(3)) for axis in range(3))


def _ml(lattice: int) -> tuple[tuple[tuple[int, ...], ...], tuple[tuple[int, ...], ...], int]:
    record = _section("group.lattices.ml")[36 * (lattice - 1):36 * lattice]
    den = int(_section("group.lattices.ml_denom")[2 * (lattice - 1)])
    to_ml = tuple(tuple(int(record[i + 3 * j]) for j in range(3)) for i in range(3))
    back = tuple(tuple(int(record[9 + i + 3 * j]) for j in range(3)) for i in range(3))
    return to_ml, back, den


def _generators(sg: int) -> list[tuple[int, int, int, int, int]]:
    point_group = _section("group.space_groups.point_group")
    counts = _section("group.point_groups.generator_count")
    start = sum(int(counts[int(point_group[s - 1]) - 1]) for s in range(1, sg))
    n = int(counts[int(point_group[sg - 1]) - 1])
    flat = _section("group.space_groups.generators")
    return [tuple(int(v) for v in flat[5 * (start + g):5 * (start + g + 1)]) for g in range(n)]


@lru_cache(maxsize=None)
def space_group_records(sg: int) -> tuple[tuple[int, int, int, int, int], ...]:
    """The coset representatives (x, y, z, den, point op) of a space group: the closure of its generators, the
    translations in the ml basis reduced into [0, 1), by ascending point operation."""

    gens = [(g[4], tuple(Fraction(g[a], g[3]) % 1 for a in range(3))) for g in _generators(sg)]
    identity = 49 if any(op > 48 for op, _ in gens) else 1
    elements = {identity: (Fraction(0),) * 3}
    frontier = [identity]
    while frontier:
        op = frontier.pop()
        t = elements[op]
        rotation = point_matrix(op)
        for gop, gt in gens:
            product = point_product(gop, op)
            if product in elements:
                continue
            moved = _column_apply(rotation, gt)
            elements[product] = tuple((moved[a] + t[a]) % 1 for a in range(3))
            frontier.append(product)
    to_ml, _back, _den = _ml(int(_section("group.space_groups.lattice")[sg - 1]))
    out = []
    for op in sorted(elements):
        v = [x % 1 for x in _column_apply(tuple(zip(*to_ml)), elements[op])]
        den = 1
        for x in v:
            den = den * x.denominator // _gcd(den, x.denominator)
        out.append((*(int(x * den) for x in v), den, op))
    return tuple(out)


def _gcd(a: int, b: int) -> int:
    while b:
        a, b = b, a % b
    return a


@lru_cache(maxsize=None)
def _lattice_ops(lattice: int) -> frozenset[int]:
    lattices = _section("group.space_groups.lattice")
    return frozenset(r[4] for sg in range(1, 231) if int(lattices[sg - 1]) == lattice for r in space_group_records(sg))


def pml_point_matrix(sg: int, op: int) -> tuple[int, ...] | None:
    """Point operation op written in the ml basis of sg's lattice (9 entries, column-major), None when no space group
    of that lattice uses it."""

    lattice = int(_section("group.space_groups.lattice")[sg - 1])
    if op not in _lattice_ops(lattice):
        return None
    to_ml, back, den = _ml(lattice)
    m = tuple(zip(*point_matrix(op)))
    product = _mul(_mul(to_ml, m), back)
    out = [0] * 9
    for i in range(3):
        for j in range(3):
            if product[i][j] % den:
                raise ValueError(f"point operation {op} is not integral in the ml basis of SG{sg}")
            out[i + 3 * j] = product[i][j] // den
    return tuple(out)


def _magnetic_parts(op: int) -> tuple[int, bool]:
    return int(_section("group.magnetic_ops.ordinary")[op - 1]), bool(_section("group.magnetic_ops.reversed")[op - 1])


@lru_cache(maxsize=None)
def _magnetic_op_of() -> dict:
    return {_magnetic_parts(op): op for op in range(1, 145)}


def magnetic_point_op_product(left: int, right: int) -> int:
    """The magnetic point operation of right's ordinary part times left's, time reversed when exactly one is."""

    (lo, lr), (ro, rr) = _magnetic_parts(left), _magnetic_parts(right)
    ordinary = point_product(ro, lo)
    return _magnetic_op_of().get((ordinary, lr != rr), 0) if ordinary else 0


def magnetic_ordinary_group(group: int) -> int:
    """The space group a magnetic group is listed under (the families are stored by space group)."""

    start = _section("group.magnetic_groups.by_sg_start")
    for sg in range(230, 0, -1):
        if int(start[sg - 1]) <= group:
            return sg
    raise ValueError(f"magnetic group {group} precedes every family")


def _magnetic_generators(group: int) -> list[tuple[int, int, int, int, int]]:
    count = int(_section("group.magnetic_groups.generator_count")[group - 1])
    flat = _section("group.magnetic_groups.generators")
    return [tuple(int(v) for v in flat[20 * (group - 1) + 5 * i:20 * (group - 1) + 5 * i + 5]) for i in range(count)]


def magnetic_group_type(group: int) -> int:
    """1 without time reversal, 2 when time reversal itself is an element, 4 when it comes only with a translation that
    is not a lattice vector, 3 otherwise; from the closure of the generators (translations in the ml basis)."""

    gens = [(g[4], tuple(Fraction(g[a], g[3]) for a in range(3))) for g in _magnetic_generators(group)]
    if not any(_magnetic_parts(op)[1] for op, _ in gens):
        return 1
    identity = _magnetic_op_of()[(49 if any(_magnetic_parts(op)[0] > 48 for op, _ in gens) else 1, False)]
    elements = {(identity, (Fraction(0),) * 3)}
    frontier = list(elements)
    sg = magnetic_ordinary_group(group)
    while frontier:
        op, t = frontier.pop()
        m = pml_point_matrix(sg, _magnetic_parts(op)[0])
        for gop, gt in gens:
            moved = tuple(sum(m[axis + 3 * k] * gt[k] for k in range(3)) for axis in range(3))
            element = (magnetic_point_op_product(op, gop), tuple((moved[a] + t[a]) % 1 for a in range(3)))
            if element not in elements:
                elements.add(element)
                frontier.append(element)
    reversed_identity = [t for op, t in elements if _magnetic_parts(op) in ((1, True), (49, True))]
    if any(all(x == 0 for x in t) for t in reversed_identity):
        return 2
    return 4 if reversed_identity else 3


__all__ = [
    "magnetic_group_type",
    "magnetic_ordinary_group",
    "magnetic_point_op_product",
    "pml_point_matrix",
    "point_matrix",
    "point_product",
    "space_group_records",
]
