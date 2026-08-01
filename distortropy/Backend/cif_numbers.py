"""Numeric conversion for ordinary CIF values."""

from __future__ import annotations

import re
from typing import Any


def float_cif_number(value: Any) -> float | None:
    """Return a CIF number without its estimated-standard-deviation suffix."""

    text = str(value or "").strip()
    if not text or text in {".", "?"}:
        return None
    text = re.sub(r"\([^)]*\)", "", text)
    try:
        return float(text)
    except ValueError:
        return None


__all__ = ["float_cif_number"]
