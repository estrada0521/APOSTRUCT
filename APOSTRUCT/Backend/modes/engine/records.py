from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction

import numpy as np

from APOSTRUCT.Backend.source.tables import (
    IsotropySubduction as IsotropySubduction,
    LittleIrrep,
    WyckoffRow as WyckoffRow,
)

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
class WyckoffSubduction:
    wyckoff: WyckoffRow
    subduction_entry: int
    pointer: int
    count: int
    pairs: tuple[tuple[int, int], ...]  # (site point-group irrep, frequency)


@dataclass(frozen=True)
class IsotropySubgroupRow:
    row_id: int
    subgroup: int
    arms: int
    direction: int
    basis: tuple[int, ...]
    origin: tuple[int, int, int, int]


@dataclass(frozen=True)
class ProjectCandidateRows:
    parent_ops: tuple[int, ...]
    site_ops: tuple[int, ...]
    pg_irrep: int
    rows: np.ndarray


@dataclass(frozen=True)
class ProjectSelection:
    parent_ops: tuple[int, ...]
    site_ops: tuple[int, ...]
    pg_irrep: int
    selected_indices: tuple[int, ...]
    real_blocks: np.ndarray
    reduced_rows: np.ndarray


@dataclass(frozen=True)
class LittleKVectorRecord:
    gid: int
    sg: int
    kslot: int
    nmod: int
    pointer: int
    records: tuple[tuple[int, int, int, int], ...]
    vectors: tuple[tuple[Fraction, Fraction, Fraction], ...]
