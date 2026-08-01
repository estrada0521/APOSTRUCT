"""Presentation transforms shared by complete-mode definition builders."""

from __future__ import annotations

import math
from typing import Any, Iterable, Mapping, Sequence

from distortropy.Backend.modes.periodic import periodic_float_close3


Vector3 = tuple[float, float, float]
ModeRow = dict[str, Vector3]
BucketKey = tuple[int, int, int]


def _inverse3(matrix: Sequence[Sequence[float]]) -> tuple[Vector3, Vector3, Vector3]:
    a = tuple(tuple(float(value) for value in row) for row in matrix)
    det = (
        a[0][0] * (a[1][1] * a[2][2] - a[1][2] * a[2][1])
        - a[0][1] * (a[1][0] * a[2][2] - a[1][2] * a[2][0])
        + a[0][2] * (a[1][0] * a[2][1] - a[1][1] * a[2][0])
    )
    if abs(det) <= 1e-15:
        raise ValueError("mode presentation basis is singular")
    return (
        (
            (a[1][1] * a[2][2] - a[1][2] * a[2][1]) / det,
            (a[0][2] * a[2][1] - a[0][1] * a[2][2]) / det,
            (a[0][1] * a[1][2] - a[0][2] * a[1][1]) / det,
        ),
        (
            (a[1][2] * a[2][0] - a[1][0] * a[2][2]) / det,
            (a[0][0] * a[2][2] - a[0][2] * a[2][0]) / det,
            (a[0][2] * a[1][0] - a[0][0] * a[1][2]) / det,
        ),
        (
            (a[1][0] * a[2][1] - a[1][1] * a[2][0]) / det,
            (a[0][1] * a[2][0] - a[0][0] * a[2][1]) / det,
            (a[0][0] * a[1][1] - a[0][1] * a[1][0]) / det,
        ),
    )


def _row_multiply(row: Sequence[float], matrix: Sequence[Sequence[float]]) -> Vector3:
    return tuple(
        sum(float(row[col]) * float(matrix[col][axis]) for col in range(3))
        for axis in range(3)
    )  # type: ignore[return-value]


def centering_translations(symbol: str) -> tuple[Vector3, ...]:
    """Return conventional centering translations for a Hermann-Mauguin symbol."""
    letter = str(symbol).strip()[:1].upper()
    extra = {
        "A": ((0.0, 0.5, 0.5),),
        "B": ((0.5, 0.0, 0.5),),
        "C": ((0.5, 0.5, 0.0),),
        "I": ((0.5, 0.5, 0.5),),
        "F": ((0.0, 0.5, 0.5), (0.5, 0.0, 0.5), (0.5, 0.5, 0.0)),
        "R": ((2 / 3, 1 / 3, 1 / 3), (1 / 3, 2 / 3, 2 / 3)),
    }.get(letter, ())
    return ((0.0, 0.0, 0.0),) + extra


def _same_vector(left: Sequence[float], right: Sequence[float], tol: float) -> bool:
    return all(abs(float(left[axis]) - float(right[axis])) <= tol for axis in range(3))


def _bucket_width(tol: float) -> float | None:
    if not math.isfinite(tol) or tol <= 0.0:
        return None
    width = max(2.0 * tol, math.sqrt(tol))
    if width <= 0.0 or not math.isfinite(1.0 / width):
        return None
    return width


def _periodic_bucket_key(values: Sequence[float], count: int) -> BucketKey | None:
    folded = tuple(float(value) % 1.0 for value in values)
    if len(folded) != 3 or not all(math.isfinite(value) for value in folded):
        return None
    return tuple(min(count - 1, int(math.floor(value * count))) for value in folded)  # type: ignore[return-value]


def _periodic_bucket_candidates(
    buckets: Mapping[BucketKey, Sequence[int]],
    values: Sequence[float],
    center: BucketKey,
    *,
    periodic_count: int,
    tol: float,
) -> list[int]:
    folded = tuple(float(value) % 1.0 for value in values)
    axis_keys: list[tuple[int, ...]] = []
    for axis, key in enumerate(center):
        keys = [key]
        lower = key / periodic_count
        upper = (key + 1) / periodic_count
        if folded[axis] - lower <= tol:
            keys.append((key - 1) % periodic_count)
        if upper - folded[axis] <= tol:
            keys.append((key + 1) % periodic_count)
        axis_keys.append(tuple(dict.fromkeys(keys)))
    candidates = {
        index
        for first in axis_keys[0]
        for second in axis_keys[1]
        for third in axis_keys[2]
        for index in buckets.get((first, second, third), ())
    }
    return sorted(candidates)


def present_mode_rows(
    parent_rows: Iterable[Mapping[str, Sequence[float]]],
    *,
    basis: Sequence[Sequence[float]],
    origin: Sequence[float],
    centering_symbol: str,
    tol: float = 1e-8,
    include_centering_ordinal: bool = False,
) -> dict[str, Any]:
    """Transform parent rows to the selected child setting and Web normalization."""
    inverse = _inverse3(basis)
    shift = tuple(float(value) for value in origin)
    primitive: list[dict[str, Any]] = []
    for row in parent_rows:
        parent_xyz = tuple(float(value) for value in row["xyz"])
        parent_dxyz = tuple(float(value) for value in row["dxyz"])
        xyz = _row_multiply(tuple(parent_xyz[i] - shift[i] for i in range(3)), inverse)
        metadata = {
            key: value for key, value in row.items() if key not in {"xyz", "dxyz"}
        }
        primitive.append(
            {
                **metadata,
                "xyz": tuple(value % 1.0 for value in xyz),
                "dxyz": _row_multiply(parent_dxyz, inverse),
            }
        )

    scale = max((abs(value) for row in primitive for value in row["dxyz"]), default=0.0)
    if scale > 0:
        primitive = [
            {**row, "dxyz": tuple(value / scale for value in row["dxyz"])}
            for row in primitive
        ]

    translations = centering_translations(centering_symbol)
    quotient: list[ModeRow] = []
    width = _bucket_width(tol)
    periodic_count = max(1, int(math.floor(1.0 / width))) if width is not None else None
    orbit_position_buckets: dict[BucketKey, list[int]] = {}
    position_index_complete = True
    for row in primitive:
        position_key = (
            _periodic_bucket_key(row["xyz"], periodic_count)
            if width is not None and periodic_count is not None
            else None
        )
        position_candidates = (
            _periodic_bucket_candidates(
                orbit_position_buckets,
                row["xyz"],
                position_key,
                periodic_count=periodic_count,
                tol=tol,
            )
            if position_key is not None and position_index_complete
            else list(range(len(quotient)))
        )
        # Position proximity is a complete prefilter; vector equality remains final.
        candidate_indices = position_candidates
        if any(
            _same_vector(row["dxyz"], quotient[index]["dxyz"], tol)
            and any(
                periodic_float_close3(
                    row["xyz"],
                    tuple(quotient[index]["xyz"][i] + translation[i] for i in range(3)),
                    tol,
                )
                for translation in translations
            )
            for index in candidate_indices
        ):
            continue
        quotient.append(row)
        if width is not None and periodic_count is not None:
            translated_keys = tuple(
                _periodic_bucket_key(
                    tuple(
                        (row["xyz"][axis] + translation[axis]) % 1.0
                        for axis in range(3)
                    ),
                    periodic_count,
                )
                for translation in translations
            )
            if all(key is not None for key in translated_keys):
                for translated_key in translated_keys:
                    orbit_position_buckets.setdefault(translated_key, []).append(
                        len(quotient) - 1
                    )  # type: ignore[arg-type]
            else:
                position_index_complete = False

    conventional: list[ModeRow] = []
    position_buckets: dict[BucketKey, list[int]] = {}
    for row in quotient:
        for centering_ordinal, translation in enumerate(translations):
            candidate = {
                **row,
                "xyz": tuple((row["xyz"][i] + translation[i]) % 1.0 for i in range(3)),
                "dxyz": row["dxyz"],
                **(
                    {"_presentation_centering_ordinal": centering_ordinal}
                    if include_centering_ordinal
                    else {}
                ),
            }
            position_key = (
                _periodic_bucket_key(candidate["xyz"], periodic_count)
                if width is not None and periodic_count is not None
                else None
            )
            candidate_indices = (
                _periodic_bucket_candidates(
                    position_buckets,
                    candidate["xyz"],
                    position_key,
                    periodic_count=periodic_count,
                    tol=tol,
                )
                if position_key is not None
                else range(len(conventional))
            )
            if any(
                periodic_float_close3(candidate["xyz"], conventional[index]["xyz"], tol)
                and _same_vector(candidate["dxyz"], conventional[index]["dxyz"], tol)
                for index in candidate_indices
            ):
                continue
            conventional.append(candidate)
            if position_key is not None:
                position_buckets.setdefault(position_key, []).append(
                    len(conventional) - 1
                )
    return {"primitive_rows": quotient, "rows": conventional, "scale": scale}


def child_lattice_cartesian(
    parent_lattice: Sequence[float], basis: Sequence[Sequence[float]]
) -> tuple[Vector3, Vector3, Vector3]:
    """Build child Cartesian row vectors from six parent lattice parameters."""
    a, b, c, alpha, beta, gamma = (float(value) for value in parent_lattice)
    ar, br, gr = (math.radians(value) for value in (alpha, beta, gamma))
    parent = (
        (a, 0.0, 0.0),
        (b * math.cos(gr), b * math.sin(gr), 0.0),
        (
            c * math.cos(br),
            c * (math.cos(ar) - math.cos(br) * math.cos(gr)) / math.sin(gr),
            c
            * math.sqrt(
                max(
                    0.0,
                    1.0
                    - math.cos(br) ** 2
                    - ((math.cos(ar) - math.cos(br) * math.cos(gr)) / math.sin(gr))
                    ** 2,
                )
            ),
        ),
    )
    return tuple(_row_multiply(row, parent) for row in basis)  # type: ignore[return-value]


def mode_normfactor(
    primitive_rows: Iterable[Mapping[str, Sequence[float]]],
    child_cartesian: Sequence[Sequence[float]],
) -> float | None:
    """Return the complete-mode normfactor from primitive displacement rows."""
    total = 0.0
    for row in primitive_rows:
        cartesian = _row_multiply(row["dxyz"], child_cartesian)
        total += sum(value * value for value in cartesian)
    return None if total <= 0 else 1.0 / math.sqrt(total)
