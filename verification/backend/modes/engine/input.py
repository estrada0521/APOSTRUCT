"""Runtime input records shared by mode construction."""

from __future__ import annotations

from fractions import Fraction
import re


def k_label_from_irrep(label: str) -> str:
    match = re.match(r"([A-Z]+)", label.strip())
    if not match:
        raise ValueError(f"cannot infer k label from irrep label {label!r}")
    return match.group(1)


def parse_fraction_text(value: str | int | float) -> Fraction:
    return Fraction(value)
