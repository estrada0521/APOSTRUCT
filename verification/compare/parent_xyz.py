"""The displayed parent representative of a mapped CIF site."""

from __future__ import annotations

import re
from typing import Any

from verification.symmetry.fraction_expression import (
    evaluate_fraction_expression,
    split_coordinate_expression3,
)


def _safe_eval_fraction_expr(expr: str, params: dict[str, Any]) -> float:
    expr = re.sub(
        r"(?<![A-Za-z0-9_])([+-]?(?:\d+(?:/\d+)?|\d*\.\d+))([xyz])\b", r"\1*\2", expr
    )
    return float(
        evaluate_fraction_expression(
            expr,
            params,
            "unknown Wyckoff parameter {name}",
            f"unsupported Wyckoff expression {expr!r}",
        )
    )


def parent_display_xyz(site: dict[str, Any]) -> list[float]:
    """Return the Source-formula representative for a mapped CIF site.

    The CIF may contain any symmetry-equivalent representative, while
    the local display uses the representative implied by the matched Source
    Wyckoff formula. Use the stored formula and fitted parameters when
    available, then fall back to the raw CIF fraction.
    """

    parts = split_coordinate_expression3(site.get("wyckoff_formula"))
    params = site.get("wyckoff_params") or {}
    if parts and isinstance(params, dict):
        try:
            return [_safe_eval_fraction_expr(part, params) % 1.0 for part in parts]
        except Exception:
            pass
    xyz = site.get("fract") or (0, 0, 0)
    return [float(value) % 1.0 for value in xyz]
