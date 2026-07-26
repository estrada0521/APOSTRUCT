"""Shared mode-construction helpers.

Extracted mechanically from the former monolithic runtime.
"""

from __future__ import annotations

from fractions import Fraction
from functools import lru_cache
import math
from pathlib import Path
import re
from typing import Any
from ISODISTORT.Assembled.Backend.exactmath import (
    fraction_matrix_inverse3,
    fraction_matrix_multiply3,
    fraction_row_multiply3,
)
from ISODISTORT.Assembled.Backend.source.tables import (
    SOURCE as _SOURCE_ROOT,
    SourceTables,
    source_tables,
)
from ISODISTORT.Assembled.Backend.modes.engine.decoder import ModeDataDecoder
from ISODISTORT.Assembled.Backend.modes.engine.input import parse_fraction_text

_fraction_matmul = fraction_matrix_multiply3
_fraction_row_multiply = fraction_row_multiply3

def _cell_params(cif_info: dict[str, Any]) -> tuple[float, float, float, float, float, float]:
    lattice = cif_info.get("lattice") or {}
    keys = ("a", "b", "c", "alpha", "beta", "gamma")
    try:
        return tuple(float(lattice[key]) for key in keys)  # type: ignore[return-value]
    except (KeyError, TypeError, ValueError) as exc:
        raise ValueError(f"cannot build lattice parameters from CIF summary: {lattice!r}") from exc



def _site_params(site: dict[str, Any]) -> tuple[float, ...] | None:
    params = site.get("wyckoff_params") or {}
    if not isinstance(params, dict) or not params:
        return None
    values: list[float] = []
    ordered_keys = ("x", "y", "z")
    for index, key in enumerate(ordered_keys):
        if key in params:
            values.append(float(params[key]))
        elif any(later in params for later in ordered_keys[index + 1:]):
            # The Source mode engine consumes positional x/y/z slots. A
            # Wyckoff form such as (0,y,z) therefore still needs the x=0
            # placeholder so y and z land on the correct parameter vectors.
            values.append(0.0)
    while values and abs(values[-1]) < 1e-15:
        values.pop()
    return tuple(values) or None



def _k_params(k_params: dict[str, str] | None) -> tuple[Fraction, ...]:
    if not k_params:
        return ()
    values: list[Fraction] = []
    for key in ("a", "b", "g", "alpha", "beta", "gamma"):
        value = k_params.get(key)
        if value is not None and str(value).strip():
            values.append(parse_fraction_text(str(value).strip()))
    return tuple(values)



def _mapped_sites(cif_info: dict[str, Any]) -> list[dict[str, Any]]:
    return [
        site
        for site in (cif_info.get("atom_sites") or [])
        if isinstance(site, dict) and site.get("wyckoff")
    ]



def _site_label(site: dict[str, Any]) -> str:
    wyckoff = str(site.get("wyckoff") or "")
    multiplicity = str(site.get("wyckoff_multiplicity") or site.get("multiplicity") or "")
    if wyckoff and wyckoff[0].isdigit():
        return wyckoff
    return f"{multiplicity}{wyckoff}" if multiplicity or wyckoff else ""



@lru_cache(maxsize=1)
def _assembled_data() -> SourceTables:
    """Reuse the immutable Source tables within the current process."""

    return source_tables()



@lru_cache(maxsize=4)
def _mode_decoder(source_dir: str | None) -> ModeDataDecoder:
    """Reuse the immutable Source mode tables within one process."""

    data_dir = source_dir or "Source"
    if Path(data_dir).resolve() == Path(_SOURCE_ROOT).resolve():
        assembled = _assembled_data()
        return ModeDataDecoder(data_dir, tables=assembled)
    return ModeDataDecoder(data_dir)



def _fold_fractional_xyz(xyz: Any) -> tuple[float, float, float] | None:
    try:
        values = tuple(float(Fraction(str(value))) % 1.0 for value in xyz)
    except (TypeError, ValueError, ZeroDivisionError):
        return None
    if len(values) != 3:
        return None
    return tuple(0.0 if abs(value) < 1e-10 or abs(value - 1.0) < 1e-10 else value for value in values)  # type: ignore[return-value]



def _mat4_multiply_fraction(left: list[list[Fraction]], right: list[list[Fraction]]) -> list[list[Fraction]]:
    return [
        [sum(left[row][idx] * right[idx][col] for idx in range(4)) for col in range(4)]
        for row in range(4)
    ]



def _normalize_setting_matrix(raw: list[list[int]]) -> list[list[Fraction]]:
    den = Fraction(int(raw[3][3] or 1))
    out = [[Fraction(int(raw[row][col]), 1) / den for col in range(4)] for row in range(4)]
    out[0][3] = out[1][3] = out[2][3] = Fraction(0)
    out[3][3] = Fraction(1)
    return out



def _fraction_floor(value: Fraction) -> int:
    return value.numerator // value.denominator



def _fraction_mod01(value: Fraction) -> Fraction:
    return value - _fraction_floor(value)



def _fraction_matrix_inverse_3(
    matrix: list[list[float]],
) -> tuple[tuple[Fraction, Fraction, Fraction], ...] | None:
    rows = tuple(
        tuple(Fraction(str(value)).limit_denominator(1000000) for value in row)
        for row in matrix
    )
    try:
        return fraction_matrix_inverse3(rows)
    except (TypeError, ValueError):
        return None



def _basis_matrix(row: dict[str, Any] | None) -> list[list[float]] | None:
    if not isinstance(row, dict):
        return None
    basis = row.get("basis_pml_to_cinter") or row.get("basis")
    if isinstance(basis, list) and len(basis) == 9:
        basis = [basis[0:3], basis[3:6], basis[6:9]]
    if not (
        isinstance(basis, list)
        and len(basis) == 3
        and all(isinstance(item, list) and len(item) == 3 for item in basis)
    ):
        return None
    try:
        return [[float(value) for value in item] for item in basis]
    except (TypeError, ValueError):
        return None



def _transform_dxyz(vector: list[Any], matrix: list[list[float]] | None) -> list[float]:
    values = [float(value) for value in vector]
    if matrix is None:
        return values
    return [sum(matrix[row][col] * values[col] for col in range(3)) for row in range(3)]



def _origin_record_from_any(origin: Any) -> tuple[int, int, int, int] | None:
    if isinstance(origin, tuple) and len(origin) == 4:
        try:
            return tuple(int(value) for value in origin)  # type: ignore[return-value]
        except (TypeError, ValueError):
            return None
    values = _origin_vector(origin)
    den = 1
    fracs = [Fraction(str(value)).limit_denominator(1000000) for value in values]
    for value in fracs:
        den = math.lcm(den, value.denominator)
    return tuple(int(value * den) for value in fracs) + (den,)  # type: ignore[return-value]



def _integer_basis_tuple(basis: list[list[float]] | None) -> tuple[int, ...] | None:
    if basis is None:
        return None
    out: list[int] = []
    for row in basis:
        for value in row:
            frac = Fraction(str(value)).limit_denominator(1000000)
            if frac.denominator != 1:
                return None
            out.append(int(frac))
    return tuple(out) if len(out) == 9 else None



def _matrix_from_basis_tuple(basis: tuple[int, ...]) -> tuple[tuple[Fraction, Fraction, Fraction], ...]:
    return tuple(
        tuple(Fraction(int(basis[row * 3 + col]), 1) for col in range(3))
        for row in range(3)
    )  # type: ignore[return-value]



def _fraction_vecadd(
    left: tuple[Fraction, Fraction, Fraction],
    right: tuple[Fraction, Fraction, Fraction],
) -> tuple[Fraction, Fraction, Fraction]:
    return tuple(left[index] + right[index] for index in range(3))  # type: ignore[return-value]



def _fraction_vecsub(
    left: tuple[Fraction, Fraction, Fraction],
    right: tuple[Fraction, Fraction, Fraction],
) -> tuple[Fraction, Fraction, Fraction]:
    return tuple(left[index] - right[index] for index in range(3))  # type: ignore[return-value]



def _origin_record_vector(origin: tuple[int, int, int, int]) -> tuple[Fraction, Fraction, Fraction]:
    den = int(origin[3])
    return (Fraction(int(origin[0]), den), Fraction(int(origin[1]), den), Fraction(int(origin[2]), den))



def _det3_float(matrix: list[list[float]]) -> float:
    return (
        matrix[0][0] * (matrix[1][1] * matrix[2][2] - matrix[1][2] * matrix[2][1])
        - matrix[0][1] * (matrix[1][0] * matrix[2][2] - matrix[1][2] * matrix[2][0])
        + matrix[0][2] * (matrix[1][0] * matrix[2][1] - matrix[1][1] * matrix[2][0])
    )



def _float_matrix_inverse_3(matrix: list[list[float]]) -> list[list[float]] | None:
    det = _det3_float(matrix)
    if abs(det) < 1e-12:
        return None
    return [
        [
            (matrix[1][1] * matrix[2][2] - matrix[1][2] * matrix[2][1]) / det,
            (matrix[0][2] * matrix[2][1] - matrix[0][1] * matrix[2][2]) / det,
            (matrix[0][1] * matrix[1][2] - matrix[0][2] * matrix[1][1]) / det,
        ],
        [
            (matrix[1][2] * matrix[2][0] - matrix[1][0] * matrix[2][2]) / det,
            (matrix[0][0] * matrix[2][2] - matrix[0][2] * matrix[2][0]) / det,
            (matrix[0][2] * matrix[1][0] - matrix[0][0] * matrix[1][2]) / det,
        ],
        [
            (matrix[1][0] * matrix[2][1] - matrix[1][1] * matrix[2][0]) / det,
            (matrix[0][1] * matrix[2][0] - matrix[0][0] * matrix[2][1]) / det,
            (matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0]) / det,
        ],
    ]



def _row_multiply(row: list[float], matrix: list[list[float]]) -> list[float]:
    return [
        sum(row[col] * matrix[col][axis] for col in range(3))
        for axis in range(3)
    ]



def _origin_vector(origin: Any) -> list[float]:
    if isinstance(origin, (list, tuple)) and len(origin) >= 4:
        try:
            denominator = float(origin[3])
            if denominator:
                return [float(origin[index]) / denominator for index in range(3)]
        except (TypeError, ValueError):
            pass
    if isinstance(origin, (list, tuple)) and len(origin) == 3:
        try:
            return [float(Fraction(str(origin[index]))) for index in range(3)]
        except (TypeError, ValueError, ZeroDivisionError):
            pass
    if isinstance(origin, str):
        text = origin.strip()
        if text.startswith("(") and text.endswith(")"):
            parts = [part.strip() for part in text[1:-1].split(",")]
            if len(parts) >= 3:
                try:
                    return [float(Fraction(part)) for part in parts[:3]]
                except (ValueError, ZeroDivisionError):
                    pass
    return [0.0, 0.0, 0.0]



def _fold01(value: float) -> float:
    folded = value % 1.0
    if abs(folded - 1.0) < 1e-10 or abs(folded) < 1e-10:
        return 0.0
    return folded
def _isotropy_row_id_from_opd_row(selected_opd: dict[str, Any] | None) -> int | None:
    if not isinstance(selected_opd, dict):
        return None
    iso = selected_opd.get("isotropy") or selected_opd
    if not isinstance(iso, dict):
        return None
    value = iso.get("row_id")
    try:
        return int(value)
    except (TypeError, ValueError):
        return None



def _origin_from_opd_row(selected_opd: dict[str, Any] | None) -> tuple[int, int, int, int] | str | None:
    if not isinstance(selected_opd, dict):
        return None
    iso = selected_opd.get("isotropy") or selected_opd
    if not isinstance(iso, dict):
        return None
    origin = iso.get("origin")
    if isinstance(origin, str) and origin.strip():
        return origin
    if not (isinstance(origin, list) and len(origin) >= 4):
        return None
    try:
        return tuple(int(origin[index]) for index in range(4))  # type: ignore[return-value]
    except (TypeError, ValueError):
        return None



def _isotropy_from_opd_row(selected_opd: dict[str, Any] | None) -> dict[str, Any] | None:
    if not isinstance(selected_opd, dict):
        return None
    iso = selected_opd.get("isotropy") or selected_opd
    return iso if isinstance(iso, dict) else None



def _is_full_parameter_opd(selected_opd: dict[str, Any] | None, freeparam: int | None) -> bool:
    if not isinstance(selected_opd, dict) or freeparam is None:
        return False
    direction = selected_opd.get("direction")
    if not isinstance(direction, dict):
        return False
    try:
        dimension = int(direction.get("dimension"))
        subduction = int(direction.get("subduction"))
    except (TypeError, ValueError):
        return False
    return int(freeparam) == dimension and subduction == dimension



def _freeparam_from_opd_row(
    decoder: ModeDataDecoder,
    selected_opd: dict[str, Any] | None,
) -> int | None:
    iso = _isotropy_from_opd_row(selected_opd)
    if isinstance(iso, dict) and iso.get("free") is not None:
        try:
            return int(iso["free"])
        except (TypeError, ValueError):
            return None
    row_id = _isotropy_row_id_from_opd_row(selected_opd)
    if row_id is None:
        return None
    try:
        return decoder.isotropy_orderparam_freeparam(int(row_id))
    except (IndexError, TypeError, ValueError):
        return None



def _k_label_from_irrep_label(label: str) -> str:
    text = str(label).strip()
    if len(text) > 1 and text.startswith("m") and text[1].isupper():
        text = text[1:]
    match = re.match(r"^[A-Za-z]+", text)
    return match.group(0) if match else text



def _same_source_kparam(
    left: tuple[int, int, int, int] | list[int] | None,
    right: tuple[int, int, int, int] | list[int] | None,
) -> bool:
    if left is None or right is None or len(left) < 4 or len(right) < 4:
        return False
    return all(
        Fraction(int(left[index]), int(left[3])) % 1 == Fraction(int(right[index]), int(right[3])) % 1
        for index in range(3)
    )



def _metric_from_cell(params: tuple[float, float, float, float, float, float]) -> list[list[float]]:
    a, b, c, alpha, beta, gamma = params
    ca = math.cos(math.radians(alpha))
    cb = math.cos(math.radians(beta))
    cg = math.cos(math.radians(gamma))
    return [
        [a * a, a * b * cg, a * c * cb],
        [a * b * cg, b * b, b * c * ca],
        [a * c * cb, b * c * ca, c * c],
    ]



def _dot_metric(u: list[float], v: list[float], metric: list[list[float]]) -> float:
    return sum(u[i] * metric[i][j] * v[j] for i in range(3) for j in range(3))



def _cell_from_basis(
    params: tuple[float, float, float, float, float, float],
    basis: list[list[float]] | None,
) -> dict[str, float]:
    keys = ("a", "b", "c", "alpha", "beta", "gamma")
    if basis is None:
        return {key: value for key, value in zip(keys, params)}
    metric = _metric_from_cell(params)
    lengths = [math.sqrt(max(_dot_metric(row, row, metric), 0.0)) for row in basis]

    def angle(i: int, j: int) -> float:
        denom = lengths[i] * lengths[j]
        if denom == 0:
            return 0.0
        value = max(-1.0, min(1.0, _dot_metric(basis[i], basis[j], metric) / denom))
        return math.degrees(math.acos(value))

    # alpha is angle(b, c), beta is angle(a, c), gamma is angle(a, b).
    return {
        "a": lengths[0],
        "b": lengths[1],
        "c": lengths[2],
        "alpha": angle(1, 2),
        "beta": angle(0, 2),
        "gamma": angle(0, 1),
    }
