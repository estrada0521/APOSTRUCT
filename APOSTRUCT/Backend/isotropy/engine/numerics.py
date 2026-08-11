"""Fixed-stride numeric kernels used by the isotropy engine."""

from __future__ import annotations


def ncmp(value: float, *, eps: float = 1e-6) -> bool:
    """Return whether a scalar is nonzero at the ``ncmp_`` tolerance."""

    return eps < abs(float(value))


def xrowop2(
    buffer: list[float] | tuple[float, ...],
    rows: int,
    cols: int,
    *,
    stride: int = 50,
    in_place: bool = False,
) -> tuple[float, ...] | list[float]:
    """Reduce a real matrix in the fixed-stride ``xrowop2_`` layout.

    The work matrix is column-major: row ``r`` and column
    ``c`` live at ``r + c * stride``.  The routine performs Gauss-Jordan
    reduction over the active ``rows`` x ``cols`` window, preserving the fixed
    stride buffer layout.
    """

    out = buffer if in_place else [float(value) for value in buffer]
    required = max(0, (cols - 1) * stride + rows)
    if len(out) < required:
        out.extend([0.0] * (required - len(out)))

    active = [
        [out[row + col * stride] for col in range(cols)]
        for row in range(rows)
    ]

    pivot_col = 0
    for pivot_row in range(rows):
        while pivot_col < cols:
            swap_row = None
            for row in range(pivot_row, rows):
                if abs(active[row][pivot_col]) > 1e-6:
                    swap_row = row
                    break
            if swap_row is None:
                pivot_col += 1
                continue

            if swap_row != pivot_row:
                active[pivot_row], active[swap_row] = active[swap_row], active[pivot_row]

            pivot_values = active[pivot_row]
            pivot = pivot_values[pivot_col]
            for col in range(pivot_col, cols):
                pivot_values[col] /= pivot

            for row in range(rows):
                if row == pivot_row:
                    continue
                row_values = active[row]
                factor = row_values[pivot_col]
                if abs(factor) <= 1e-6:
                    continue
                for col in range(cols):
                    if col == pivot_col:
                        continue
                    row_values[col] -= factor * pivot_values[col]
                row_values[pivot_col] = 0.0

            pivot_col += 1
            break
        else:
            break
    for row, values in enumerate(active):
        for col, value in enumerate(values):
            out[row + col * stride] = value
    return out if in_place else tuple(out)
