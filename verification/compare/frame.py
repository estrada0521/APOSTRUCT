"""Child frames related by a parent symmetry: a local subgroup in another domain than the Web's.

A candidate maps local child coordinates to Web child coordinates, y_web = y_local U + s, with
U = B_local Q B_web^-1 integral for a rotation Q of the parent lattice metric and s fixed by the undistorted atoms
(every Web site representative is the image of a local orbit point of the same parent site). The mode definitions
then decide among the candidates (Definitions.transport)."""

from __future__ import annotations

from fractions import Fraction
import itertools
import math
import re
from typing import Any, Mapping, Sequence

Matrix = tuple[tuple[Fraction, Fraction, Fraction], ...]


def _det(m: Sequence[Sequence[Fraction]]) -> Fraction:
    return (m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1]) - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
            + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0]))


def _mul(a: Sequence[Sequence[Any]], b: Sequence[Sequence[Any]]) -> tuple:
    return tuple(tuple(sum(a[i][k] * b[k][j] for k in range(3)) for j in range(3)) for i in range(3))


def _inverse(m: Matrix) -> Matrix | None:
    d = _det(m)
    if d == 0:
        return None
    c = lambda i, j: m[(j + 1) % 3][(i + 1) % 3] * m[(j + 2) % 3][(i + 2) % 3] - m[(j + 1) % 3][(i + 2) % 3] * m[(j + 2) % 3][(i + 1) % 3]
    return tuple(tuple(Fraction(c(i, j)) / d for j in range(3)) for i in range(3))


_UNIMODULAR = tuple(
    m for m in (tuple(tuple(v[3 * i:3 * i + 3]) for i in range(3)) for v in itertools.product((-1, 0, 1), repeat=9))
    if abs(_det(m)) == 1
)


def _metric(cell: Sequence[float]) -> tuple:
    a, b, c = cell[:3]
    ca, cb, cg = (math.cos(math.radians(x)) for x in cell[3:])
    return ((a * a, a * b * cg, a * c * cb), (a * b * cg, b * b, b * c * ca), (a * c * cb, b * c * ca, c * c))


def parent_rotations(cell: Sequence[float]) -> list[tuple]:
    """Rotations x -> x Q (row vectors, entries -1..1) keeping the parent metric, the identity first."""

    M = _metric(cell)
    out = []
    for Q in _UNIMODULAR:
        QMQ = _mul(_mul(Q, M), tuple(zip(*Q)))
        if all(abs(QMQ[i][j] - M[i][j]) <= 1e-3 * math.sqrt(M[i][i] * M[j][j]) for i in range(3) for j in range(3)):
            out.append(Q)
    out.sort(key=lambda q: q != ((1, 0, 0), (0, 1, 0), (0, 0, 1)))
    return out


def _rows(text: str) -> Matrix:
    return tuple(tuple(Fraction(x) for x in r.split(",")) for r in re.findall(r"\(([^()]*)\)", text))


def _kind(label: str) -> str:
    return re.sub(r"_\d+$", "", str(label))


def _wrap(x: float) -> float:
    return x - math.floor(x)


def _close(a: Sequence[float], b: Sequence[float], tol: float) -> bool:
    return all(min(abs(_wrap(x - y)), 1 - abs(_wrap(x - y))) <= tol for x, y in zip(a, b))


def domain_frames(
    *,
    parent_cell: Sequence[float],
    web_basis: str,
    web_origin: str,
    local_basis: Sequence[Sequence[float]],
    local_origin: str,
    local_atoms: Sequence[Mapping[str, Any]],
    web_atoms: Sequence[Any],
    tol: float = 2e-3,
    limit: int = 512,
) -> list[dict[str, Any]]:
    """Candidate frames (U, s, Q), ordered by rotation (the identity first) and by local orbit point."""

    web_b = _rows(web_basis)
    local_b = tuple(tuple(Fraction(v).limit_denominator(1000) for v in row) for row in local_basis)
    web_inv = _inverse(web_b) if len(web_b) == 3 else None
    if web_inv is None or len(local_b) != 3:
        return []
    local_points: dict[str, list[tuple[float, float, float]]] = {}
    for atom in local_atoms:
        points = atom.get("_presentation_orbit_points") or [atom.get("xyz")]
        local_points.setdefault(_kind(atom.get("label")), []).extend(tuple(float(v) for v in p) for p in points if p)
    web_reps: dict[str, list[tuple[float, float, float]]] = {}
    web_count: dict[str, int] = {}
    for atom in web_atoms:
        web_reps.setdefault(_kind(atom.label), []).append(tuple(atom.xyz))
        match = re.match(r"\d+", str(atom.site))
        web_count[_kind(atom.label)] = web_count.get(_kind(atom.label), 0) + (int(match.group()) if match else 0)
    if set(local_points) != set(web_reps) or any(len(local_points[k]) != web_count[k] for k in local_points):
        return []
    k0 = min(local_points, key=lambda k: len(local_points[k]))
    w0 = web_reps[k0][0]
    frames: list[dict[str, Any]] = []
    for Q in parent_rotations(parent_cell):
        U = _mul(_mul(local_b, Q), web_inv)
        if any(v.denominator != 1 for row in U for v in row):
            continue
        Uf = tuple(tuple(float(v) for v in row) for row in U)
        moved = {k: [tuple(sum(p[i] * Uf[i][j] for i in range(3)) for j in range(3)) for p in pts] for k, pts in local_points.items()}
        for a in moved[k0]:
            s = tuple(w0[j] - a[j] for j in range(3))
            if all(any(_close(w, tuple(p[j] + s[j] for j in range(3)), tol) for p in moved[k]) for k, reps in web_reps.items() for w in reps):
                frames.append({"basis_change": tuple(tuple(int(v) for v in row) for row in U),
                               "shift": tuple(_wrap(x) for x in s), "rotation": Q})
                if len(frames) >= limit:
                    return frames
    return frames


def transported_k(vectors: Sequence[Sequence[Fraction]], Q: Sequence[Sequence[int]]) -> frozenset:
    """Arms k Q^-T, each up to sign and modulo 1."""

    Qinv = _inverse(tuple(tuple(Fraction(v) for v in row) for row in Q))
    out = set()
    for k in vectors:
        kk = tuple(sum(k[i] * Qinv[j][i] for i in range(3)) for j in range(3))
        a, b = tuple(v % 1 for v in kk), tuple((-v) % 1 for v in kk)
        out.add(min(a, b))
    return frozenset(out)
