"""FINDSYM-derived lattice presentation helpers for local ISODISTORT.

These helpers are runtime-local ports of small FINDSYM routines.  They read the
shared Source/data_* tables through ``SourceData`` and never call FINDSYM
binaries or import FINDSYM/Disassembled modules.
"""

from __future__ import annotations

from fractions import Fraction
from functools import lru_cache
import math
from typing import Iterable, Sequence

from ISODISTORT.Assembled.Backend.exactmath import (
    integer_adjugate3,
    integer_determinant3,
)
from ISODISTORT.Assembled.Backend.isotropy.engine.source_data import SourceData


ATOL = 1e-5


def _matinv_unimodular(matrix: Iterable[int]) -> tuple[int, ...]:
    m = tuple(int(value) for value in matrix)
    if len(m) < 9:
        raise IndexError("tuple index out of range")
    det = integer_determinant3(m[:9])
    if abs(det) != 1:
        raise ValueError(f"expected unimodular matrix, det={det}")
    return tuple(value // det for value in integer_adjugate3(m[:9]))


def _matmlt(left: Iterable[int], right: Iterable[int]) -> tuple[int, ...]:
    """Port upstream ``matmlt_`` storage semantics for flat 3x3 buffers."""

    l = tuple(int(value) for value in left)
    r = tuple(int(value) for value in right)
    out = [0] * 9
    for outer in range(3):
        right_offset = 0
        for group in range(3):
            total = 0
            for k in range(3):
                total += l[outer + 3 * k] * r[right_offset + k]
            out[outer + 3 * group] = total
            right_offset += 3
    return tuple(out)


def _vmlt(matrix: Iterable[int], vector: tuple[int, int, int, int]) -> tuple[int, int, int, int]:
    m = tuple(int(value) for value in matrix)
    x, y, z, den = (int(item) for item in vector)
    return (
        x * m[0] + y * m[3] + z * m[6],
        x * m[1] + y * m[4] + z * m[7],
        x * m[2] + y * m[5] + z * m[8],
        den,
    )


def _rows(matrix: Sequence[int]) -> list[list[int]]:
    values = tuple(int(value) for value in matrix)
    return [list(values[offset:offset + 3]) for offset in range(0, 9, 3)]


def _flat(matrix: Sequence[Sequence[int]] | Sequence[int]) -> tuple[int, ...]:
    if len(matrix) == 9 and not isinstance(matrix[0], (list, tuple)):  # type: ignore[index]
        return tuple(int(value) for value in matrix)  # type: ignore[arg-type]
    return tuple(int(value) for row in matrix for value in row)  # type: ignore[union-attr]


def _get_gmat(lattice: Sequence[float]) -> list[list[float]]:
    a, b, c, alpha, beta, gamma = (float(value) for value in lattice)
    return [
        [a * a, a * b * math.cos(math.radians(gamma)), a * c * math.cos(math.radians(beta))],
        [a * b * math.cos(math.radians(gamma)), b * b, b * c * math.cos(math.radians(alpha))],
        [a * c * math.cos(math.radians(beta)), b * c * math.cos(math.radians(alpha)), c * c],
    ]


def _transform_gram(transform: Sequence[Sequence[int]], gram: Sequence[Sequence[float]]) -> list[list[float]]:
    return [
        [
            sum(
                int(transform[row][p]) * float(gram[p][q]) * int(transform[col][q])
                for p in range(3)
                for q in range(3)
            )
            for col in range(3)
        ]
        for row in range(3)
    ]


def _fraction_integral_delta(
    left: tuple[int, int, int, int],
    right: tuple[int, int, int, int],
) -> bool:
    denominator = math.lcm(int(left[3]), int(right[3]))
    for axis in range(3):
        delta = int(left[axis]) * (denominator // int(left[3])) - int(right[axis]) * (
            denominator // int(right[3])
        )
        if delta % denominator:
            return False
    return True


def _monoclinic_pml_candidate_matrix(
    sg: int,
    cinter_transform: Sequence[Sequence[int]] | Sequence[int],
    *,
    data: SourceData,
) -> tuple[int, ...] | None:
    inverse = _matinv_unimodular(_flat(cinter_transform))
    rows: list[int] = []
    for vector in ((1, 0, 0, 1), (0, 1, 0, 1), (0, 0, 1, 1)):
        changed = data.vector_change_setting_record(int(sg), "pml", "cinter", vector)
        changed = _vmlt(inverse, changed)
        changed = data.vector_change_setting_record(int(sg), "cinter", "pml", changed)
        if changed[3] != 1:
            return None
        rows.extend(changed[:3])
    return tuple(rows)


def _monoclinic_candidate_preserves_translations(
    sg: int,
    cinter_transform: Sequence[Sequence[int]] | Sequence[int],
    *,
    data: SourceData,
) -> bool:
    pml_matrix = _monoclinic_pml_candidate_matrix(int(sg), cinter_transform, data=data)
    if pml_matrix is None:
        return False
    for record in data.generate_space_group_records(int(sg))[1:]:
        transformed = _vmlt(pml_matrix, record[:4])
        if not _fraction_integral_delta(record[:4], transformed):
            return False
    return True


def _axis_tuple_for_monoclinic_sg(sg: int, data: SourceData) -> tuple[int, int, int]:
    vector = data.vector_change_setting_record(int(sg), "cml", "cinter", (0, 0, 1, 1))
    if vector[0] != 0:
        return 0, 1, 2
    if vector[1] != 0:
        return 1, 0, 2
    if vector[2] != 0:
        return 2, 0, 1
    raise ValueError(f"failed to determine monoclinic unique axis for SG{sg}: {vector}")


def monoclinic_lattparam_transform(
    sg: int,
    lattice: Sequence[float],
    *,
    data: SourceData | None = None,
) -> tuple[tuple[int, int, int], ...]:
    """Return FINDSYM ``monoclinic_lattparam_`` basis transform.

    The returned matrix follows the OPD presentation convention:

    ``web_basis = transform * local_basis``.
    """

    source = data or SourceData()
    gram = _get_gmat(lattice)
    transform = [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
    identity = [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
    unique_axis, axis_a, axis_c = _axis_tuple_for_monoclinic_sg(int(sg), source)

    def accept(candidate: list[list[int]]) -> bool:
        return _monoclinic_candidate_preserves_translations(int(sg), candidate, data=source)

    def accumulate(candidate: list[list[int]]) -> None:
        nonlocal transform
        transform = _rows(_matmlt(_flat(transform), _flat(candidate)))

    axis_pair = [axis_a, axis_c]
    for _ in range(2):
        other, target = axis_pair
        step = -1 if gram[target][other] > 0 else 1
        offset = 0
        while True:
            offset += step
            candidate = [row[:] for row in identity]
            candidate[target][other] = offset
            if not accept(candidate):
                continue
            next_gram = _transform_gram(candidate, gram)
            if next_gram[target][target] < gram[target][target] - ATOL:
                gram = next_gram
                accumulate(candidate)
                offset = 0
                continue
            break
        axis_pair = [axis_pair[1], axis_pair[0]]

    low_axis, high_axis = sorted((axis_a, axis_c))
    if gram[high_axis][high_axis] < gram[low_axis][low_axis] - ATOL:
        candidate = [[0, 0, 0], [0, 0, 0], [0, 0, 0]]
        candidate[unique_axis][unique_axis] = 1
        candidate[low_axis][high_axis] = -1
        candidate[high_axis][low_axis] = 1
        if accept(candidate):
            gram = _transform_gram(candidate, gram)
            accumulate(candidate)

    if gram[high_axis][low_axis] > ATOL:
        candidate = [row[:] for row in identity]
        candidate[unique_axis][unique_axis] = -1
        candidate[high_axis][high_axis] = -1
        if accept(candidate):
            accumulate(candidate)

    return tuple(tuple(int(value) for value in row) for row in transform)


def _cell_cartesian_basis(
    params: tuple[float, float, float, float, float, float],
) -> list[list[float]]:
    a, b, c, alpha, beta, gamma = params
    alpha_r = math.radians(alpha)
    beta_r = math.radians(beta)
    gamma_r = math.radians(gamma)
    va = [a, 0.0, 0.0]
    vb = [b * math.cos(gamma_r), b * math.sin(gamma_r), 0.0]
    cx = c * math.cos(beta_r)
    cy = c * (math.cos(alpha_r) - math.cos(beta_r) * math.cos(gamma_r)) / math.sin(gamma_r)
    cz = math.sqrt(max(0.0, c * c - cx * cx - cy * cy))
    return [va, vb, [cx, cy, cz]]


def _float_matmul(left: list[list[float]], right: list[list[float]]) -> list[list[float]]:
    return [
        [sum(left[row][k] * right[k][col] for k in range(3)) for col in range(3)]
        for row in range(3)
    ]


def _cartesian_vector(row: list[Fraction], parent: list[list[float]]) -> list[float]:
    return [sum(float(row[idx]) * parent[idx][col] for idx in range(3)) for col in range(3)]


@lru_cache(maxsize=1)
def _default_source_data() -> SourceData:
    return SourceData()


def _child_lattice_params_from_basis(
    basis_rows: list[list[Fraction]],
    parent_cell: tuple[float, float, float, float, float, float],
) -> tuple[float, float, float, float, float, float]:
    parent = _cell_cartesian_basis(parent_cell)
    vectors = [_cartesian_vector(row, parent) for row in basis_rows]

    def norm(vector: list[float]) -> float:
        return math.sqrt(sum(value * value for value in vector))

    def angle(vec_a: list[float], vec_b: list[float]) -> float:
        denom = norm(vec_a) * norm(vec_b)
        if denom <= 0:
            return 90.0
        value = max(-1.0, min(1.0, sum(vec_a[i] * vec_b[i] for i in range(3)) / denom))
        return math.degrees(math.acos(value))

    return (
        norm(vectors[0]),
        norm(vectors[1]),
        norm(vectors[2]),
        angle(vectors[1], vectors[2]),
        angle(vectors[0], vectors[2]),
        angle(vectors[0], vectors[1]),
    )


def _apply_int_basis_transform(
    transform: tuple[tuple[int, int, int], ...],
    basis_rows: list[list[Fraction]],
) -> list[list[Fraction]]:
    return [
        [
            sum(Fraction(transform[row][k]) * basis_rows[k][col] for k in range(3))
            for col in range(3)
        ]
        for row in range(3)
    ]


def _findsym_monoclinic_lattparam_basis_candidate(
    subgroup: int,
    basis_rows: list[list[Fraction]],
    parent_cell: tuple[float, float, float, float, float, float] | None,
    *,
    data: SourceData | None = None,
) -> list[list[Fraction]]:
    if not (3 <= int(subgroup) <= 15) or parent_cell is None:
        return basis_rows
    try:
        transform = monoclinic_lattparam_transform(
            int(subgroup),
            _child_lattice_params_from_basis(basis_rows, parent_cell),
            data=data or _default_source_data(),
        )
    except Exception:
        return basis_rows
    return _apply_int_basis_transform(transform, basis_rows)


def _metric_for_fraction_basis(
    basis_rows: list[list[Fraction]],
    parent_cell: tuple[float, float, float, float, float, float],
) -> list[list[float]]:
    parent = _cell_cartesian_basis(parent_cell)
    basis_float = [[float(value) for value in row] for row in basis_rows]
    child = _float_matmul(basis_float, parent)
    return [
        [sum(child[row][axis] * child[col][axis] for axis in range(3)) for col in range(3)]
        for row in range(3)
    ]


_HEX_BASAL_TRANSFORMS: tuple[tuple[tuple[int, int, int], ...], ...] = (
    ((1, 0, 0), (0, 1, 0), (0, 0, 1)),
    ((0, -1, 0), (1, 1, 0), (0, 0, 1)),
    ((-1, -1, 0), (1, 0, 0), (0, 0, 1)),
    ((-1, 0, 0), (0, -1, 0), (0, 0, 1)),
    ((0, 1, 0), (-1, -1, 0), (0, 0, 1)),
    ((1, 1, 0), (-1, 0, 0), (0, 0, 1)),
)


_TRICLINIC_NCMP = 1e-6
_TRICLINIC_OFFDIAG_SNAP = 1.7e-5


def _triclinic_ncmp(value: float) -> bool:
    """Port upstream ``ncmp_``: significant when ``|value| > 1e-6``."""

    return abs(value) > _TRICLINIC_NCMP


def _triclinic_nint(value: float) -> int:
    """Fortran ``NINT``: round half away from zero."""

    return int(math.floor(value + 0.5)) if value >= 0.0 else int(math.ceil(value - 0.5))


def _float_congruence(transform: list[list[float]], gram: list[list[float]]) -> list[list[float]]:
    """``G' = M G M^T`` for row-vector basis transforms ``rows' = M rows``."""

    return [
        [
            sum(transform[row][p] * gram[p][q] * transform[col][q] for p in range(3) for q in range(3))
            for col in range(3)
        ]
        for row in range(3)
    ]


def _gram_lattparam(gram: list[list[float]]) -> tuple[float, float, float, float, float, float]:
    a = math.sqrt(gram[0][0])
    b = math.sqrt(gram[1][1])
    c = math.sqrt(gram[2][2])
    alpha = math.degrees(math.acos(max(-1.0, min(1.0, gram[1][2] / b / c))))
    beta = math.degrees(math.acos(max(-1.0, min(1.0, gram[0][2] / a / c))))
    gamma = math.degrees(math.acos(max(-1.0, min(1.0, gram[0][1] / a / b))))
    return (a, b, c, alpha, beta, gamma)


def _triclinic_lattparam_sort(gram: list[list[float]]) -> tuple[list[list[float]], list[list[float]] | None]:
    """Ascending stable diagonal sort; ``det < 0`` negates the last-placed axis."""

    used = [False, False, False]
    maxdiag = max(gram[0][0], gram[1][1], gram[2][2])
    perm = [[0.0] * 3 for _ in range(3)]
    changed = False
    picked = 0
    for position in range(3):
        best = maxdiag + 1.0
        for index in range(3):
            if used[index]:
                continue
            value = gram[index][index]
            if value < best and _triclinic_ncmp(value - best):
                picked = index
                best = value
        used[picked] = True
        perm[position][picked] = 1.0
        if picked != position:
            changed = True
    det = (
        perm[0][0] * (perm[1][1] * perm[2][2] - perm[1][2] * perm[2][1])
        - perm[0][1] * (perm[1][0] * perm[2][2] - perm[1][2] * perm[2][0])
        + perm[0][2] * (perm[1][0] * perm[2][1] - perm[1][1] * perm[2][0])
    )
    if det < 0:
        perm[2] = [-value for value in perm[2]]
    if not changed:
        return gram, None
    return _float_congruence(perm, gram), perm


def _triclinic_lattparam_sign(
    gram: list[list[float]],
    zero_default: int,
) -> tuple[list[list[float]], list[list[float]] | None]:
    """Mixed-sign sector fix toward all-negative off-diagonals.

    Zero off-diagonals (below the upstream 1.7e-5 snap) count as
    ``zero_default``: -1 in the first pass, +1 after the sum transform.
    """

    for row, col in ((0, 1), (0, 2), (1, 2)):
        if abs(gram[row][col]) < _TRICLINIC_OFFDIAG_SNAP:
            gram[row][col] = 0.0
            gram[col][row] = 0.0

    def sign_of(value: float) -> int:
        if not _triclinic_ncmp(value):
            return zero_default
        return 1 if value > 0 else -1

    sign_zeta = sign_of(gram[0][1])
    sign_eta = sign_of(gram[0][2])
    sign_xi = sign_of(gram[1][2])
    total = sign_zeta + sign_eta + sign_xi
    if abs(total) == 3:
        return gram, None
    flips = [[1.0 if i == j else 0.0 for j in range(3)] for i in range(3)]
    if total == sign_xi:
        flips[0][0] = -1.0
    if total == sign_eta:
        flips[1][1] = -1.0
    if total == sign_zeta:
        flips[2][2] = -1.0
    return _float_congruence(flips, gram), flips


def _triclinic_lattparam_reduce(gram: list[list[float]], max_outer: int = 10) -> list[list[float]]:
    """Faithful float port of FINDSYM ``triclinic_lattparam_`` on a Gram matrix.

    Returns the accumulated row transform.  Tie policy is the binary's: NINT
    pair shears accepted only on a >1e-6 diagonal change, tolerance-stable
    sort, all-negative sign normalization with asymmetric zero defaults, and
    the ``c' = -sign(eta) a - sign(xi) b + c`` sum transform, iterated to a
    lattice-parameter fixpoint.
    """

    total = [[1.0 if i == j else 0.0 for j in range(3)] for i in range(3)]

    def apply(step: list[list[float]], current: list[list[float]]) -> list[list[float]]:
        nonlocal total
        total = _float_matmul(step, total)
        return _float_congruence(step, current)

    previous = _gram_lattparam(gram)
    for _ in range(max_outer):
        order = [0, 1, 2]
        failures = 0
        while True:
            pivot, second, third = order
            shear_second = _triclinic_nint(-(gram[second][pivot] / gram[pivot][pivot]))
            shear_third = _triclinic_nint(-(gram[third][pivot] / gram[pivot][pivot]))
            if shear_second == 0 and shear_third == 0:
                failures += 1
            else:
                step = [[1.0 if i == j else 0.0 for j in range(3)] for i in range(3)]
                step[second][pivot] = float(shear_second)
                step[third][pivot] = float(shear_third)
                candidate = _float_congruence(step, gram)
                if any(_triclinic_ncmp(candidate[d][d] - gram[d][d]) for d in range(3)):
                    gram = apply(step, gram)
                else:
                    failures += 1
            if failures > 2:
                break
            order = [(value + 1) % 3 for value in order]
        gram, sort_step = _triclinic_lattparam_sort(gram)
        if sort_step is not None:
            total = _float_matmul(sort_step, total)
        gram, sign_step = _triclinic_lattparam_sign(gram, -1)
        if sign_step is not None:
            total = _float_matmul(sign_step, total)
        product = gram[0][1] * gram[0][2] * gram[1][2]
        slack = (
            gram[0][0]
            + gram[1][1]
            - abs(2.0 * gram[0][1])
            - abs(2.0 * gram[0][2])
            - abs(2.0 * gram[1][2])
        )
        if product < 0.0 and slack < 0.0 and _triclinic_ncmp(product) and _triclinic_ncmp(slack):
            step = [[1.0 if i == j else 0.0 for j in range(3)] for i in range(3)]
            step[2][0] = -1.0 if gram[0][2] > 0.0 else 1.0
            step[2][1] = -1.0 if gram[1][2] > 0.0 else 1.0
            gram = apply(step, gram)
            gram, sort_step = _triclinic_lattparam_sort(gram)
            if sort_step is not None:
                total = _float_matmul(sort_step, total)
            gram, sign_step = _triclinic_lattparam_sign(gram, 1)
            if sign_step is not None:
                total = _float_matmul(sign_step, total)
        params = _gram_lattparam(gram)
        if all(not _triclinic_ncmp(params[index] - previous[index]) for index in range(6)):
            break
        previous = params
    return total


def _triclinic_lattparam_basis_candidate(
    subgroup: int,
    basis_rows: list[list[Fraction]],
    parent_cell: tuple[float, float, float, float, float, float] | None,
) -> list[list[Fraction]]:
    """Present a triclinic child basis with the FINDSYM ``triclinic_lattparam_`` policy.

    The upstream ``param_6`` contract is kept: when the reduced lattice
    parameters equal the input parameters within 1e-6, the caller discards the
    transform and the input basis is returned unchanged.
    """

    if int(subgroup) > 2 or parent_cell is None:
        return basis_rows
    try:
        rows = [[Fraction(value) for value in row] for row in basis_rows]
        gram = _metric_for_fraction_basis(rows, parent_cell)
    except Exception:
        return basis_rows
    cell = list(_gram_lattparam(gram))
    for index in range(3, 6):
        if abs(cell[index] - 90.0) < 0.001:
            cell[index] = 90.0
    gram = _get_gmat(cell)
    transform = _triclinic_lattparam_reduce(gram)
    reduced = _float_congruence(transform, _get_gmat(cell))
    final_params = _gram_lattparam(reduced)
    if all(not _triclinic_ncmp(final_params[index] - cell[index]) for index in range(6)):
        return basis_rows
    step = [[Fraction(_triclinic_nint(value)) for value in row] for row in transform]
    return [
        [sum(step[i][k] * rows[k][j] for k in range(3)) for j in range(3)]
        for i in range(3)
    ]


def _hex_trig_basal_basis_candidate(
    subgroup: int,
    basis_rows: list[list[Fraction]],
    *,
    parametric: bool,
    coupled: bool = False,
) -> list[list[Fraction]]:
    """Apply the anchored fixed-coupled trigonal basal convention."""

    power: int | None = None
    if coupled and not parametric:
        power = 4 if int(subgroup) == 160 else None
    if power is None or power == 0:
        return basis_rows
    transform = _HEX_BASAL_TRANSFORMS[power]
    return [
        [sum(Fraction(transform[row][k]) * basis_rows[k][col] for k in range(3)) for col in range(3)]
        for row in range(3)
    ]


def present_opd_basis_rows(
    subgroup: int,
    basis_rows: list[list[Fraction]],
    parent_cell: tuple[float, float, float, float, float, float] | None,
    *,
    data: SourceData | None = None,
    parametric: bool = False,
    coupled: bool = False,
) -> list[list[Fraction]]:
    """Apply OPD-table lattice presentation transforms for a subgroup basis."""

    rows = _triclinic_lattparam_basis_candidate(int(subgroup), basis_rows, parent_cell)
    rows = _hex_trig_basal_basis_candidate(
        int(subgroup),
        rows,
        parametric=parametric,
        coupled=coupled,
    )
    return _findsym_monoclinic_lattparam_basis_candidate(int(subgroup), rows, parent_cell, data=data)
