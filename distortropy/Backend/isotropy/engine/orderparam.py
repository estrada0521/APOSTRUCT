"""Order-parameter labeling, reduction, and invariance kernels."""

from __future__ import annotations

from collections.abc import Sequence

from distortropy.Backend.isotropy.engine.numerics import ncmp, xrowop2


def orderparam_label(free: int, index: int) -> str:
    """Return the public label shape produced by ``orderparam_label_``."""

    if free == 1:
        prefix = "P"
    elif free == 2:
        prefix = "C"
    elif free == 3:
        prefix = "S"
    elif 4 <= free <= 99:
        prefix = f"{free}D"
    else:
        raise ValueError(f"invalid order-parameter dimension: {free}")
    return f"{prefix}{index}"


def _buffer_value(buffer: Sequence[float], row: int, col: int, *, stride: int) -> float:
    index = row + col * stride
    if index >= len(buffer):
        return 0.0
    return float(buffer[index])


def _nice_score(rows: int, cols: int, buffer: Sequence[float], *, stride: int = 48) -> tuple[int, int, int, int, int]:
    """Score one OPD matrix exactly like ``nice_orderparam_``.

    The upstream buffer is a fixed-stride Fortran work array.  ``rows`` and
    ``cols`` are the two integer arguments passed to ``nice_orderparam_``.
    """

    score = [0, 0, 0, 0, 0]
    for row in range(rows):
        seen_nonzero_in_row = False
        for col in range(cols):
            value = _buffer_value(buffer, row, col, stride=stride)
            if not ncmp(value):
                continue
            if ncmp(abs(value) - 1.0):
                score[1] += 1
            if seen_nonzero_in_row:
                score[0] += 1
            else:
                score[2] += 1
            if value < 0.0:
                score[3] += 1
            score[4] += row + 1
            seen_nonzero_in_row = True
    return tuple(score)  # type: ignore[return-value]


def nice_orderparam(rows: int, cols: int, candidate_a: Sequence[float], candidate_b: Sequence[float], *, stride: int = 48) -> int:
    """Compare two candidates with the ``nice_orderparam_`` ordering.

    Return values identify the preferred candidate:

    - ``1``: candidate A is preferred
    - ``2``: candidate B is preferred
    - ``0``: candidates tie
    """

    first_a = float(candidate_a[0]) if candidate_a else 0.0
    first_b = float(candidate_b[0]) if candidate_b else 0.0
    if ncmp(first_a) and not ncmp(first_b):
        return 1
    if not ncmp(first_a) and ncmp(first_b):
        return 2

    score_a = _nice_score(rows, cols, candidate_a, stride=stride)
    score_b = _nice_score(rows, cols, candidate_b, stride=stride)
    for item_a, item_b in zip(score_a, score_b):
        if item_a < item_b:
            return 1
        if item_b < item_a:
            return 2
    return 0


def orthogonal_orderparam(rows: int, cols: int, buffer: Sequence[float], *, stride: int = 48) -> tuple[int, tuple[float, ...]]:
    """Apply the ``orthogonal_orderparam_`` row-orthogonalization contract.

    Return ``(flag, buffer)`` where ``flag`` reports whether rows changed:

    - ``1``: rows were already orthogonal
    - ``0``: a reverse Gram-Schmidt pass was applied
    """

    out = [float(value) for value in buffer]
    required = max(0, (rows - 1) * stride + cols)
    if len(out) < required:
        out.extend([0.0] * (required - len(out)))

    def value(row: int, col: int) -> float:
        return out[row * stride + col]

    already_orthogonal = True
    for row in range(1, rows):
        for prior in range(row):
            dot = sum(value(row, col) * value(prior, col) for col in range(cols))
            if ncmp(dot):
                already_orthogonal = False
                break
        if not already_orthogonal:
            break
    if already_orthogonal:
        return 1, tuple(out)

    for row in range(rows - 2, -1, -1):
        for other in range(row + 1, rows):
            dot = sum(value(other, col) * value(row, col) for col in range(cols))
            norm = sum(value(other, col) * value(other, col) for col in range(cols))
            if ncmp(norm):
                for col in range(cols):
                    out[row * stride + col] -= (dot / norm) * value(other, col)
    return 0, tuple(out)


def eqs_to_orderparam(cols: int, eq_count: int, equations: Sequence[float], *, eq_stride: int = 50, op_stride: int = 48) -> tuple[int, tuple[float, ...]]:
    """Convert equations to the ``eqs_to_orderparam_`` OPD layout.

    ``equations`` is the column-stride work matrix where equation
    ``r`` and variable ``c`` live at ``r + c * eq_stride``.  The returned OPD
    buffer uses the usual stride-48 row layout.
    """

    def eq(row: int, col: int) -> float:
        index = row + col * eq_stride
        if index >= len(equations):
            return 0.0
        return float(equations[index])

    pivots: list[int] = []
    for row in range(eq_count):
        for col in range(cols):
            if ncmp(eq(row, col)):
                pivots.append(col)
                break
        else:
            raise ValueError("eqs_to_orderparam_: equation row has no pivot")

    free_count = 0
    out = [0.0] * (max(cols - eq_count, 0) * op_stride)
    pivot_set = set(pivots)
    for source_col in range(cols - 1, -1, -1):
        if source_col in pivot_set:
            continue
        free_count += 1
        row_index = free_count - 1
        required = row_index * op_stride + cols
        if len(out) < required:
            out.extend([0.0] * (required - len(out)))
        # Variables are scanned from the last source column back to the
        # first, but writes them into the OPD buffer in reversed display order.
        out[row_index * op_stride + (cols - source_col - 1)] = 1.0
        for eq_row, pivot_col in enumerate(pivots):
            pivot_value = eq(eq_row, pivot_col)
            if not ncmp(pivot_value):
                raise ValueError("eqs_to_orderparam_: zero pivot")
            out[row_index * op_stride + (cols - pivot_col - 1)] = -(eq(eq_row, source_col) / pivot_value)

    if cols - eq_count != free_count:
        raise ValueError("eqs_to_orderparam_: inconsistent free dimension")
    return free_count, tuple(out)


def orderparam_to_eqs(cols: int, free: int, orderparam: Sequence[float], *, op_stride: int = 48, eq_stride: int = 50) -> tuple[int, tuple[float, ...]]:
    """Convert OPD rows to the ``orderparam_to_eqs_`` equation layout.

    The input OPD rows use the stride-48 layout. The returned
    equation matrix uses stride 50 and the same reversed display-column order
    as ``eqs_to_orderparam_``.
    """

    # Work on a dense copy and reduce the OPD rows to RREF-like form.
    reduced = [
        [
            float(orderparam[index]) if index < len(orderparam) else 0.0
            for col in range(cols)
            for index in (row * op_stride + col,)
        ]
        for row in range(max(free, 0))
    ]
    for row in range(free):
        pivot = None
        for col in range(cols):
            if ncmp(reduced[row][col]):
                pivot = col
                break
        if pivot is None:
            raise ValueError("orderparam_to_eqs_: OPD row has no pivot")
        pivot_value = reduced[row][pivot]
        for col in range(pivot, cols):
            reduced[row][col] /= pivot_value
        for other in range(free):
            if other == row:
                continue
            factor = reduced[other][pivot]
            if ncmp(factor):
                for col in range(cols):
                    reduced[other][col] -= factor * reduced[row][col]

    column_state = [0] * cols
    eq_count = 0
    assignments: list[tuple[int, int, float]] = []

    for row in range(free):
        pivot = None
        for col in range(cols):
            if ncmp(reduced[row][col]):
                pivot = col
                break
        if pivot is None:
            raise ValueError("orderparam_to_eqs_: reduced row has no pivot")
        if column_state[pivot] != 0:
            raise ValueError("orderparam_to_eqs_: duplicate pivot")
        column_state[pivot] = -1
        for col in range(pivot + 1, cols):
            value = reduced[row][col]
            if not ncmp(value):
                continue
            existing = column_state[col]
            if existing < 1:
                if existing != 0:
                    raise ValueError("orderparam_to_eqs_: invalid column state")
                eq_count += 1
                column_state[col] = eq_count
                assignments.append((eq_count - 1, col, 1.0))
            assignments.append((column_state[col] - 1, pivot, -value))

    # Columns untouched by the OPD span become simple unit constraints.  The
    # The loop walks source columns left-to-right while its display index
    # counts down, which is equivalent to set_eq's reversed mapping.
    for col in range(cols):
        if column_state[col] == 0:
            eq_count += 1
            assignments.append((eq_count - 1, col, 1.0))

    if cols - free != eq_count:
        raise ValueError("orderparam_to_eqs_: inconsistent equation count")
    final_length = max(
        max(cols - free, 0) * eq_stride,
        (cols - 1) * eq_stride + eq_count if eq_count else 0,
    )
    for eq_index, source_col, _value in assignments:
        display_col = cols - source_col - 1
        final_length = max(final_length, eq_index + display_col * eq_stride + 1)
    out = [0.0] * final_length
    for eq_index, source_col, value in assignments:
        display_col = cols - source_col - 1
        out[eq_index + display_col * eq_stride] = value
    return eq_count, tuple(out)


def orderparam_add_eqs(
    cols: int,
    eq_count: int,
    equations: Sequence[float],
    irrep_matrix: Sequence[float],
    *,
    eq_stride: int = 50,
    matrix_stride: int = 48,
) -> tuple[int, tuple[float, ...]]:
    """Append one group element's ``orderparam_add_eqs_`` constraints.

    ``equations`` is the constraint matrix with column stride 50.
    ``irrep_matrix`` is the stride-48 matrix buffer passed in from
    ``get_irreps_``.  The routine appends the equations ``D - I`` for one
    group element, periodically row-reduces them with ``xrowop2_``, and trims
    trailing zero rows.
    """

    out = [float(value) for value in equations]
    required = max(0, (cols - 1) * eq_stride + max(eq_count, 1))
    if len(out) < required:
        out.extend([0.0] * (required - len(out)))

    def set_eq(row: int, display_col: int, value: float) -> None:
        index = row + display_col * eq_stride
        if index >= len(out):
            out.extend([0.0] * (index + 1 - len(out)))
        out[index] = value

    def matrix_value(source_row: int, source_col: int) -> float:
        index = source_row * matrix_stride + source_col
        if index >= len(irrep_matrix):
            return 0.0
        return float(irrep_matrix[index])

    def row_is_nonzero(row: int) -> bool:
        for display_col in range(cols):
            index = row + display_col * eq_stride
            if index < len(out) and ncmp(out[index]):
                return True
        return False

    active_rows = int(eq_count)
    for source_col in range(cols):
        row = active_rows
        for source_row in range(cols - 1, -1, -1):
            display_col = cols - source_row - 1
            value = matrix_value(source_row, source_col)
            if source_row == source_col:
                value -= 1.0
            set_eq(row, display_col, value)
        active_rows += 1

        if source_col == cols - 1 or active_rows == eq_stride:
            reduced = list(xrowop2(out, active_rows, cols, stride=eq_stride))
            if len(reduced) > len(out):
                out.extend([0.0] * (len(reduced) - len(out)))
            out[: len(reduced)] = reduced

            while active_rows > 0 and not row_is_nonzero(active_rows - 1):
                active_rows -= 1
            if active_rows == eq_stride:
                raise ValueError("orderparam_add_eqs_: equation buffer full")

    return active_rows, tuple(out)


def orderparam_check(
    dim: int,
    row_count: int,
    irrep_matrix: Sequence[float],
    orderparam: Sequence[float],
    *,
    matrix_stride: int = 48,
    op_stride: int = 48,
) -> int:
    """Evaluate the ``orderparam_check_`` invariance condition.

    Return ``1`` when every active OPD row is invariant under the supplied
    irrep matrix, otherwise ``0``. Buffers use the fixed stride-48 layout
    consumed by ``orderparam_to_subgroup_``.
    """

    def matrix_value(row: int, col: int) -> float:
        index = row * matrix_stride + col
        if index >= len(irrep_matrix):
            return 0.0
        return float(irrep_matrix[index])

    def op_value(row: int, col: int) -> float:
        index = row * op_stride + col
        if index >= len(orderparam):
            return 0.0
        return float(orderparam[index])

    for row in range(row_count):
        for col in range(dim):
            transformed = 0.0
            for source_row in range(dim):
                transformed += matrix_value(source_row, col) * op_value(row, source_row)
            if ncmp(transformed - op_value(row, col)):
                return 0
    return 1
