from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction

from verification.backend.source.tables import LittleIrrep, WyckoffRow as WyckoffRow

LittleIrrepRecord = LittleIrrep


@dataclass(frozen=True)
class ImageRecord:
    image_id: int
    label: str
    dimension: int
    order: int
    landau: int
    subgroup_count: int
    subgroup_pointer: int
    label2: str | None
    image_type: int | None
    generator_count: int | None
    diagonalize_pointer: int | None


@dataclass(frozen=True)
class LittleKVectorRecord:
    gid: int
    sg: int
    kslot: int
    nmod: int
    pointer: int
    records: tuple[tuple[int, int, int, int], ...]
    vectors: tuple[tuple[Fraction, Fraction, Fraction], ...]
