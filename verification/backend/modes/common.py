"""Shared mode-construction helpers."""

from __future__ import annotations

from fractions import Fraction
from functools import lru_cache
import math
from typing import Any
from verification.backend.exactmath import (
    fraction_row_multiply3,
)
from verification.backend.source.tables import (
    SOURCE as _SOURCE_ROOT,
    SourceTables,
    source_tables,
)
from verification.backend.modes.engine.decoder import ModeDataDecoder
from verification.backend.modes.engine.input import parse_fraction_text

_fraction_row_multiply = fraction_row_multiply3

def _k_params(k_params: dict[str, str] | None) -> tuple[Fraction, ...]:
    if not k_params:
        return ()
    values: list[Fraction] = []
    for key in ("a", "b", "g", "alpha", "beta", "gamma"):
        value = k_params.get(key)
        if value is not None and str(value).strip():
            values.append(parse_fraction_text(str(value).strip()))
    return tuple(values)



def _site_label(site: dict[str, Any]) -> str:
    wyckoff = str(site.get("wyckoff") or "")
    multiplicity = str(site.get("wyckoff_multiplicity") or site.get("multiplicity") or "")
    if wyckoff and wyckoff[0].isdigit():
        return wyckoff
    return f"{multiplicity}{wyckoff}" if multiplicity or wyckoff else ""



@lru_cache(maxsize=1)
def _source_tables() -> SourceTables:
    """Reuse the immutable Source tables within the current process."""

    return source_tables()



@lru_cache(maxsize=1)
def _mode_decoder() -> ModeDataDecoder:
    """Reuse the immutable Source mode tables within one process."""

    return ModeDataDecoder(str(_SOURCE_ROOT), tables=_source_tables())

def _row_multiply(row: list[float], matrix: list[list[float]]) -> list[float]:
    return [
        sum(row[col] * matrix[col][axis] for col in range(3))
        for axis in range(3)
    ]



def _metric_from_cell(params: tuple[float, float, float, float, float, float]) -> list[list[float]]:
    a, b, c, alpha, beta, gamma = params
    ca = math.cos(math.radians(alpha))
    cb = math.cos(math.radians(beta))
    cg = math.cos(math.radians(gamma))
    return [
        [a * a, a * b * cg, a * c * cb],
        [a * b * cg, b * b, b * c * ca],
        [a * c * cb, b * c * ca, c * c],
    ]


