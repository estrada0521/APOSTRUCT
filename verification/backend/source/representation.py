"""Shared numerical kernels for Source little-irrep representations."""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from fractions import Fraction
import math
from typing import Any

import numpy as np


def _translation_mod_one(x: int, y: int, z: int, den: int) -> tuple[tuple[int, int], ...]:
    """Reduced (numerator, denominator) of each coordinate of (x, y, z)/den mod 1."""

    out = []
    for value in (x, y, z):
        numerator, denominator = Fraction(value, den).as_integer_ratio()
        numerator %= denominator
        divisor = math.gcd(numerator, denominator)
        out.append((numerator // divisor, denominator // divisor))
    return tuple(out)
