from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction

from verification.backend.source.tables import LittleIrrep, WyckoffRow as WyckoffRow

LittleIrrepRecord = LittleIrrep


@dataclass(frozen=True)
class LittleKVectorRecord:
    gid: int
    sg: int
    kslot: int
    nmod: int
    pointer: int
    records: tuple[tuple[int, int, int, int], ...]
    vectors: tuple[tuple[Fraction, Fraction, Fraction], ...]
