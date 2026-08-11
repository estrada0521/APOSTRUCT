"""Runtime input records shared by mode construction."""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
import re


@dataclass(frozen=True)
class Case:
    sg: int
    wyckoff: str
    k_label: str
    params: tuple[float, float, float, float, float, float]
    title: str = ""
    atom_label: str = "Fe"
    k_direction: str | None = None
    site_params: tuple[float, ...] | None = None
    k_params: tuple[Fraction, ...] = ()


def k_label_from_irrep(label: str) -> str:
    match = re.match(r"([A-Z]+)", label.strip())
    if not match:
        raise ValueError(f"cannot infer k label from irrep label {label!r}")
    return match.group(1)


def parse_fraction_text(value: str | int | float) -> Fraction:
    return Fraction(value)
