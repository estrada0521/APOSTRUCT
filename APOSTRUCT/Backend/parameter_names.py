"""Stable order-parameter variable names shared across backend projections."""

from __future__ import annotations


def parameter_name(index: int) -> str:
    """Return Web-style a..z,aa.. parameter names in stable order."""

    if type(index) is not int or index < 0:
        raise ValueError("parameter index must be a nonnegative integer")
    value = index + 1
    parts = []
    while value:
        value, remainder = divmod(value - 1, 26)
        parts.append(chr(ord("a") + remainder))
    return "".join(reversed(parts))


__all__ = ["parameter_name"]
