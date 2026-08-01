"""Shared numerical kernels for Source little-irrep representations."""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from fractions import Fraction
import math
from typing import Any

import numpy as np


def decode_little_sparse_matrix(
    little: Mapping[str, Sequence[Any]],
    constants: Sequence[float],
    *,
    gid: int,
    dim: int,
    position: int,
) -> np.ndarray:
    """Decode one already-selected Source sparse-matrix position."""

    base = int(little["little_irr_full_matrices_irr_pointer"][gid - 1])
    matrix_index = base + position - 1
    start = int(little["little_irr_full_matrices_pointer"][matrix_index - 1])
    end = int(little["little_irr_full_matrices_pointer"][matrix_index])
    raw = little["little_irr_full_matrices"][(start - 1) * 3 : (end - 1) * 3]
    out = np.zeros((dim, dim), dtype=float)
    for offset in range(0, len(raw), 3):
        row, col, code = (int(value) for value in raw[offset : offset + 3])
        out[row - 1, col - 1] = constants[code]
    return out


def real_phase_operator(
    dim: int,
    phases: tuple[Fraction, ...],
    *,
    gid: int | None = None,
) -> np.ndarray:
    """Return the real block representation of Source star-arm phases."""

    context = "" if gid is None else f" for gid={gid}"
    if not phases:
        return np.eye(dim)
    if dim % len(phases):
        raise ValueError(f"phase count {len(phases)} does not divide dim {dim}{context}")
    block_size = dim // len(phases)
    out = np.zeros((dim, dim), dtype=float)
    for arm, phase in enumerate(phases):
        start = arm * block_size
        angle = 2.0 * math.pi * float(phase)
        c = math.cos(angle)
        s = math.sin(angle)
        if abs(c) < 1e-12:
            c = 0.0
        if abs(s) < 1e-12:
            s = 0.0
        if block_size == 1:
            if abs(s) > 1e-10:
                nphase = "" if gid is None else f" nphase={len(phases)}"
                raise ValueError(
                    f"non-real scalar phase{context}: phase={phase} dim={dim}{nphase}"
                )
            out[start, start] = c
            continue
        if block_size % 2:
            if abs(s) > 1e-10:
                block = "block" if gid is None else "block_size"
                raise ValueError(
                    f"odd real block cannot hold complex phase{context}: "
                    f"phase={phase} {block}={block_size}"
                )
            out[start : start + block_size, start : start + block_size] = c * np.eye(block_size)
            continue
        half = block_size // 2
        out[start : start + half, start : start + half] = c * np.eye(half)
        out[start : start + half, start + half : start + block_size] = s * np.eye(half)
        out[start + half : start + block_size, start : start + half] = -s * np.eye(half)
        out[start + half : start + block_size, start + half : start + block_size] = c * np.eye(half)
    return out
