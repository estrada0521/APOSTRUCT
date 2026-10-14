"""ctypes bindings to the C kernels of APOSTRUCT (bin/_apostruct.so, built by build.sh)."""

from __future__ import annotations

import ctypes
import functools
import hashlib
import json
import os
import inspect
import threading
from pathlib import Path
from typing import Sequence

from fractions import Fraction

import numpy as np

_ERRORS: dict[int, tuple[type[Exception], str]] = {
    1: (MemoryError, "native isotropy traversal ran out of memory"),
    2: (ValueError, "orderparam_to_eqs_: OPD row has no pivot"),
    3: (ValueError, "orderparam_to_eqs_: reduced row has no pivot"),
    4: (ValueError, "orderparam_to_eqs_: duplicate pivot"),
    5: (ValueError, "orderparam_to_eqs_: invalid column state"),
    6: (ValueError, "orderparam_to_eqs_: inconsistent equation count"),
    7: (ValueError, "eqs_to_orderparam_: equation row has no pivot"),
    8: (ValueError, "eqs_to_orderparam_: zero pivot"),
    9: (ValueError, "eqs_to_orderparam_: inconsistent free dimension"),
    10: (ValueError, "fixed-space intersection changed free dimension"),
    11: (KeyError, "kernel fraction coset not found"),
}

PUBLIC = Path(__file__).resolve().parents[3]


def _sources() -> list[Path]:
    paths = [*(PUBLIC / "core").glob("*.c"), *(PUBLIC / "core").glob("*.h"), PUBLIC / "source" / "blob_build.c"]
    return sorted(paths, key=lambda path: path.relative_to(PUBLIC).as_posix())


def _load() -> ctypes.CDLL:
    path = PUBLIC / "bin" / "_apostruct.so"
    stamp = PUBLIC / "bin" / "_apostruct.source"
    source = hashlib.sha256(b"".join(item.read_bytes() for item in _sources())).hexdigest()
    if not path.is_file() or not stamp.is_file() or stamp.read_text().strip() != source:
        raise ImportError(f"native kernels are missing or stale; run sh {PUBLIC / 'build.sh'}")
    library = ctypes.CDLL(str(path))
    i64 = ctypes.POINTER(ctypes.c_int64)
    i32 = ctypes.POINTER(ctypes.c_int)
    f64 = ctypes.POINTER(ctypes.c_double)
    irrep = (ctypes.c_int, ctypes.c_int, i64)
    library.orderparam_to_subgroup_rows.argtypes = irrep + (
        f64, ctypes.c_int64, ctypes.c_int, ctypes.c_int, i64, i32, i32, ctypes.c_int, i64, i64, i32,
    )
    library.orderparam_to_subgroup_rows.restype = ctypes.c_int
    blob = PUBLIC / "bin" / "source.blob"
    library.tables_open.argtypes = (ctypes.c_char_p,)
    library.tables_open.restype = ctypes.c_int
    if library.tables_open(str(blob).encode()):
        raise ImportError(f"the Source blob is missing or unreadable; run sh {PUBLIC / 'build.sh'}")
    library.magnetic_records.argtypes = (ctypes.c_int, ctypes.c_int, ctypes.c_int, i64)
    library.magnetic_records.restype = ctypes.c_int
    library.magnetic_records_iso.argtypes = (ctypes.c_int, ctypes.c_int, i64)
    library.magnetic_records_iso.restype = ctypes.c_int
    library.newlat_basis.argtypes = (ctypes.c_int, ctypes.c_int, i64, i64)
    library.newlat_basis.restype = ctypes.c_int
    library.opd_basis_presentation.argtypes = (ctypes.c_int, i64, f64, ctypes.c_int, ctypes.c_int, ctypes.c_int, i64)
    library.opd_basis_presentation.restype = ctypes.c_int
    library.op_new_fractionals.argtypes = (i64, i64, ctypes.c_int, i32)
    library.op_new_fractionals.restype = ctypes.c_int
    library.isotropy_rows_for.argtypes = irrep + (ctypes.c_int, ctypes.c_int)
    library.isotropy_rows_for.restype = ctypes.c_int
    library.isotropy_rows_ints.argtypes = (ctypes.POINTER(i64),)
    library.isotropy_rows_ints.restype = ctypes.c_int64
    library.isotropy_rows_doubles.argtypes = (ctypes.POINTER(f64),)
    library.isotropy_rows_doubles.restype = ctypes.c_int64
    library.site_atom_table.argtypes = (
        ctypes.c_int, ctypes.c_int, ctypes.c_int, i64, i64, i64, ctypes.c_int64, i64, i64, i64,
    )
    library.site_atom_table.restype = ctypes.c_int64
    library.mapped_subgroup_records.argtypes = (ctypes.c_int, ctypes.c_int, i64, i64, ctypes.c_int64, i64)
    library.mapped_subgroup_records.restype = ctypes.c_int
    library.id_subgroup_identify.argtypes = (ctypes.c_int, i64, ctypes.c_int, i64, ctypes.c_int, i64)
    library.id_subgroup_identify.restype = ctypes.c_int
    library.match_inter_wyckoff.argtypes = (
        ctypes.c_int, ctypes.c_int, f64, f64, ctypes.c_int, ctypes.c_double, i32, f64, i32, f64,
    )
    library.match_inter_wyckoff.restype = ctypes.c_int
    return library


_LIBRARY = _load()


_OPD_BASIS_ERRORS = {
    1: (ArithmeticError, "OPD basis presentation value out of range"),
    2: (ZeroDivisionError, "float division by zero"),
    3: (ValueError, "OPD basis presentation Source tables are inconsistent"),
}


def opd_basis_presentation(
    subgroup: int,
    rows: Sequence[Sequence[Fraction]],
    cell: Sequence[float] | None,
    *,
    parametric: bool,
    coupled: bool,
    magnetic_group: int | None,
) -> list[list[Fraction]]:
    """lattice_conventions.present_opd_basis_rows of one subgroup basis."""

    words = np.array(
        [part for row in rows for value in row for part in (*_i128_words(Fraction(value).numerator), *_i128_words(Fraction(value).denominator))],
        dtype=np.int64,
    )
    cell_values = None if cell is None else np.ascontiguousarray([float(value) for value in cell], dtype=float)
    out = np.zeros(36, dtype=np.int64)
    status = _LIBRARY.opd_basis_presentation(
        int(subgroup), _pointer(words, ctypes.c_int64),
        None if cell_values is None else _pointer(cell_values, ctypes.c_double), int(bool(parametric)),
        int(bool(coupled)), int(magnetic_group or 0), _pointer(out, ctypes.c_int64),
    )
    if status:
        kind, message = _OPD_BASIS_ERRORS.get(status, (RuntimeError, f"OPD basis presentation failed ({status})"))
        raise kind(message)
    values = _words_to_fractions(out.tolist())
    return [list(values[3 * row : 3 * row + 3]) for row in range(3)]


def new_fractionals(basis: Sequence[int]) -> tuple[tuple[int, int, int, int], ...]:
    """get_new_fractionals_: the parent lattice translations x M / m for the x in
    [0, m)^3 (lexicographic) with x M = 0 (mod m), m the inverse denominator."""

    matrix = _basis_array(basis)
    det = abs(_determinant3(basis))
    if not det:
        raise ValueError(f"singular basis matrix: {tuple(int(value) for value in basis[:9])}")
    out = np.zeros((det + 1, 4), dtype=np.int64)
    count = np.zeros(1, dtype=np.intc)
    status = _LIBRARY.op_new_fractionals(
        _pointer(matrix, ctypes.c_int64), _pointer(out, ctypes.c_int64), det + 1, _pointer(count, ctypes.c_int)
    )
    if status:
        _raise_orderparam(status)
    return tuple(map(tuple, out[: int(count[0])].tolist()))


def isotropy_rows(gid: int, kparam: Sequence[int], *, magnetic: bool, occurrence: int) -> tuple[list[int], list[float]]:
    """get_isotropy_ dynamic rows of one irrep as (row records, orthogonal matrices).
    occurrence > 0 stops at the occurrence-th one-parameter row."""

    status = _LIBRARY.isotropy_rows_for(*_irrep_args(gid, kparam), int(bool(magnetic)), int(occurrence))
    if status:
        _raise_isotropy_rows(status)
    ints = ctypes.POINTER(ctypes.c_int64)()
    doubles = ctypes.POINTER(ctypes.c_double)()
    int_count = _LIBRARY.isotropy_rows_ints(ctypes.byref(ints))
    double_count = _LIBRARY.isotropy_rows_doubles(ctypes.byref(doubles))
    return ints[:int_count], doubles[:double_count]


def _raise_isotropy_rows(status: int) -> None:
    if status >= 3000:
        kind, message = ValueError, "magnetic Source tables are inconsistent"
    elif status >= 2000:
        if status - 2000 >= 100:
            _raise_orderparam(status - 2100)
        kind, message = _BRIDGE_WEIGHT_ERRORS.get(status - 2000, (RuntimeError, f"bridge weights failed ({status - 2000})"))
    elif status >= 1000:
        kind, message = _TABLE_ERRORS.get(status - 1000, (RuntimeError, f"traversal tables failed ({status - 1000})"))
    elif status > 200:
        kind, message = _ID_SUBGROUP_ERRORS.get(200 - status, (RuntimeError, f"id_subgroup_ status {status}"))
    elif status > 100:
        kind, message = _ORDERPARAM_ERRORS.get(status, (RuntimeError, f"orderparam_to_subgroup_ failed ({status})"))
    else:
        kind, message = _ERRORS.get(status, (RuntimeError, f"native isotropy traversal failed ({status})"))
    raise kind(message)


def _pointer(array: np.ndarray, kind):
    return array.ctypes.data_as(ctypes.POINTER(kind))


_BRIDGE_WEIGHT_ERRORS = {
    1: (KeyError, "operation record not found in generated space-group records"),
    2: (ValueError, "phase count does not divide the irrep dimension"),
    3: (ValueError, "non-real scalar phase"),
    4: (ValueError, "odd real block cannot hold complex phase"),
    5: (ArithmeticError, "bridge weight phase out of range"),
}


_PROJECTION_ENTRY_ERRORS = {
    1: (ValueError, "projection Source tables are inconsistent"),
    2: (KeyError, "projection Source entry not found"),
    3: (ValueError, "projection Source entry is not integral"),
    4: (ValueError, "projection entry capacity exceeded"),
}


_MODE_BLOCK_ERRORS = {
    1: (MemoryError, "mode block projection ran out of memory"),
    2: (ValueError, "recordwise projection produced no Source modes"),
    3: (TypeError, "int() argument must be a string, a bytes-like object or a real number, not 'NoneType'"),
}





_K_RECORD_ERRORS = {
    1: (ValueError, "k-vector Source tables are inconsistent"),
    2: (ValueError, "k label expects more parameters"),
    3: (ValueError, "missing little_k_star_conv2ml pointer"),
    4: (ArithmeticError, "k-vector value out of range"),
    5: (ValueError, "newlat failed to find a basis"),
}


def newlat_basis(mode_count: int, records: Sequence[Sequence[int]]) -> tuple[int, ...]:
    """newlat_: the integer superlattice basis of the first mode_count records."""

    array = np.ascontiguousarray([[int(value) for value in record[:4]] for record in records] or [[0, 0, 0, 1]], dtype=np.int64)
    out = np.zeros(9, dtype=np.int64)
    status = _LIBRARY.newlat_basis(int(mode_count), len(records), _pointer(array, ctypes.c_int64), _pointer(out, ctypes.c_int64))
    if status:
        kind, message = _K_RECORD_ERRORS.get(status, (RuntimeError, f"newlat failed ({status})"))
        raise kind(message)
    return tuple(out.tolist())


_RECORDWISE_ERRORS = {
    1: (MemoryError, "recordwise projection ran out of memory"),
    2: (ValueError, "bridge Source family has no active basis values"),
    3: (ValueError, "recordwise projection Source tables are missing"),
    4: (IndexError, "tuple index out of range"),
    5: (ValueError, "recordwise projection capacity exceeded"),
    6: (ValueError, "recordwise projection has no irrep matrices"),
    7: (RuntimeError, "native bridge source modes failed"),
}


def _irrep_args(gid: int, kparam: Sequence[int]) -> tuple:
    params = np.array([int(value) for value in kparam] or [0], dtype=np.int64)
    return int(gid), len(kparam), _pointer(params, ctypes.c_int64)


_ORDERPARAM_ERRORS = {
    101: (KeyError, "operation record not found in generated space-group records"),
    102: (ValueError, "orderparam_to_subgroup_ has no Source answer"),
    103: (ArithmeticError, "orderparam_to_subgroup_ value out of range"),
    104: (ZeroDivisionError, "operation record has a zero denominator"),
    105: (ValueError, "Source little irrep tables are inconsistent"),
    106: (ValueError, "the little irrep has no native irrep tables"),
}


def _raise_orderparam(status: int) -> None:
    kind, message = _ORDERPARAM_ERRORS.get(status, (RuntimeError, f"orderparam_to_subgroup_ failed ({status})"))
    raise kind(message)


def orderparam_to_subgroup_rows(
    gid: int,
    kparam: Sequence[int],
    candidates: Sequence[Sequence[int]],
    matrix_ops: Sequence[int],
    signs: Sequence[int],
    values: Sequence[float],
    row_count: int,
    max_iter: int = 256,
) -> tuple[tuple[int, ...], tuple[tuple[int, int, int, int, int], ...]]:
    """orderparam_to_subgroup_ basis and kept records."""

    array = np.array([float(value) for value in values], dtype=float)
    try:
        records = np.array(candidates, dtype=np.int64).reshape(len(candidates), 5)
    except OverflowError:
        raise ArithmeticError("orderparam_to_subgroup_ value out of range") from None
    ops = np.array(matrix_ops, dtype=np.intc)
    sign_values = np.array(signs, dtype=np.intc)
    basis = np.zeros(9, dtype=np.int64)
    kept = np.zeros((max(1, len(candidates)), 5), dtype=np.int64)
    count = np.zeros(1, dtype=np.intc)
    status = _LIBRARY.orderparam_to_subgroup_rows(
        *_irrep_args(gid, kparam), _pointer(array, ctypes.c_double),
        len(array), int(row_count), len(records), _pointer(records, ctypes.c_int64),
        _pointer(ops, ctypes.c_int), _pointer(sign_values, ctypes.c_int), int(max_iter),
        _pointer(basis, ctypes.c_int64), _pointer(kept, ctypes.c_int64), _pointer(count, ctypes.c_int),
    )
    if status:
        _raise_orderparam(status)
    return (
        tuple(int(value) for value in basis.tolist()),
        tuple(tuple(int(value) for value in record) for record in kept[: int(count[0])].tolist()),
    )


_TABLE_ERRORS = {
    1: (ArithmeticError, "traversal table value out of range"),
    2: (KeyError, "kernel fraction coset not found"),
    3: (KeyError, "operation record not represented by ispace/kernel pair"),
    4: (KeyError, "point-operation product not found"),
    5: (ValueError, "kernel fractions are not distinct modulo the kernel lattice"),
}


_WYCKOFF_SPLIT_ERRORS = {
    1: (ArithmeticError, "Wyckoff split value out of range"),
    2: (KeyError, "child operation image not found"),
    3: (ValueError, "child Wyckoff orbit collision"),
    4: (KeyError, "no child Wyckoff candidate"),
    5: (ValueError, "could not solve child formula"),
    6: (MemoryError, "Wyckoff split capacity exceeded"),
    7: (TypeError, "child Wyckoff orbit has no site point group"),
}


_K_FAMILY_ERRORS: dict[int, tuple[type[Exception], str]] = {
    1: (ArithmeticError, "k family arithmetic out of range"),
    2: (ValueError, "k family Source tables are inconsistent"),
    3: (MemoryError, "k family search ran out of memory"),
    4: (ValueError, "k family output capacity exceeded"),
    5: (ValueError, "lattice holohedry has no identity operation"),
    6: (ZeroDivisionError, "zero denominator in k family Source tables"),
    7: (IndexError, "special k family has no star vector"),
}

def _k_family_error(status: int) -> Exception:
    kind, message = _K_FAMILY_ERRORS.get(-status, (RuntimeError, f"k family search failed ({status})"))
    return kind(message)


def _fraction_triplet(values: Sequence[int]) -> tuple[Fraction, Fraction, Fraction]:
    return tuple(Fraction(values[2 * axis], values[2 * axis + 1]) for axis in range(3))  # type: ignore[return-value]


_CHILD_SITE_ERRORS: dict[int, tuple[type[Exception], str]] = {
    1: (ValueError, "Source Wyckoff geometry is inconsistent"),
    2: (ArithmeticError, "child site arithmetic out of range"),
    3: (ValueError, "parent Wyckoff row id does not identify one Source row"),
    4: (ValueError, "Source Wyckoff split produced no Formula15 rows"),
    5: (ValueError, "Formula15 did not partition the presentation atoms"),
    6: (MemoryError, "child site expansion ran out of memory"),
    7: (ValueError, "child site output capacity exceeded"),
    8: (ValueError, "zero denominator"),
}


def _i128_words(value: int) -> tuple[int, int]:
    if not -(1 << 127) <= value < (1 << 127):
        raise ArithmeticError("child site arithmetic out of range")
    low = value & 0xFFFFFFFFFFFFFFFF
    return value >> 64, low - (1 << 64) if low >= (1 << 63) else low


_ROW_ASSIGNMENT_ERRORS: dict[int, tuple[type[Exception], str]] = {
    1: (ValueError, "invalid atom operation record"),
    2: (ValueError, "child centering ordinal out of range"),
    3: (ValueError, "presentation basis is singular"),
    4: (ValueError, "canonical child atoms share one affine occurrence identity"),
    5: (ValueError, "mode row did not identify one canonical child atom"),
    6: (MemoryError, "row assignment ran out of memory"),
    7: (ArithmeticError, "row assignment value out of range"),
}


_MAGNETIC_ERRORS: dict[int, tuple[type[Exception], str]] = {
    1: (ValueError, "magnetic Source tables are inconsistent"),
    2: (IndexError, "magnetic value out of range"),
    3: (MemoryError, "magnetic kernel ran out of memory"),
    4: (KeyError, "magnetic point operation not found"),
    5: (ValueError, "dependent magnetic Wyckoff parameter slots"),
    6: (KeyError, "no magnetic Wyckoff formula"),
    7: (ValueError, "canonical atoms share one child-cinter point"),
    8: (ValueError, "ordinary child sites do not partition magnetic atoms"),
    9: (KeyError, "magnetic Source operation leaves the canonical atom population"),
    10: (ValueError, "magnetic operation is not an atom permutation"),
    11: (ValueError, "magnetic atom component splits an ordinary child site"),
    12: (ValueError, "canonical magnetic orbit does not equal atom membership"),
    13: (ValueError, "magnetic representative does not identify one canonical atom"),
    14: (ValueError, "singular magnetic basis"),
    15: (ValueError, "magnetic kernel capacity exceeded"),
}

_MAGNETIC_SETTINGS = {"closure": 0, "binary": 1, "pml": 1, "cinter": 2}


def _magnetic_error(status: int) -> Exception:
    kind, message = _MAGNETIC_ERRORS.get(status, (RuntimeError, f"magnetic kernel failed ({status})"))
    return kind(message)


def magnetic_records(group: int, setting: str = "binary") -> tuple[tuple[int, int, int, int, int], ...]:
    """The generated Source records of a magnetic group in one coordinate setting."""

    out = np.zeros((192, 5), dtype=np.int64)
    count = _LIBRARY.magnetic_records(int(group), _MAGNETIC_SETTINGS[setting], 192, _pointer(out, ctypes.c_int64))
    if count < 0:
        raise _magnetic_error(-count)
    return tuple(map(tuple, out[:count].tolist()))


def magnetic_records_isotropy(group: int) -> tuple[tuple[int, int, int, int, int], ...]:
    """generate_spacegroup_magnetic_'s closure as the isotropy tables build it."""

    out = np.zeros((192, 5), dtype=np.int64)
    count = _LIBRARY.magnetic_records_iso(int(group), 192, _pointer(out, ctypes.c_int64))
    if count < 0:
        raise _magnetic_error(-count)
    return tuple(map(tuple, out[:count].tolist()))


def _rational_words(points: Sequence[Sequence[Fraction]]) -> np.ndarray:
    return np.array(
        [
            word
            for point in points
            for value in point
            for part in (value.numerator, value.denominator)
            for word in _i128_words(int(part))
        ],
        dtype=np.int64,
    )


def _raise_coupled_domain_orbit(status: int) -> None:
    if status == 20:
        raise KeyError("domain action leaves the retained domain subspaces")
    if status == 21:
        raise ValueError("magnetic point-operation tables are inconsistent")
    if status >= 2000:
        _subgroup_ops_records(-(status - 2000), np.zeros((0, 5), dtype=np.int64))
    if status >= 1100:
        _raise_orderparam(status - 1100)
    if status >= 1000:
        kind, message = _BRIDGE_WEIGHT_ERRORS.get(status - 1000, (RuntimeError, f"bridge weights failed ({status - 1000})"))
        raise kind(message)
    if status > 100:
        _raise_orderparam(status)
    kind, message = _ERRORS.get(status, (RuntimeError, f"coupled domain orbit failed ({status})"))
    raise kind(message)


_SITE_ATOM_ERRORS = {
    1: (ValueError, "site atom table Source tables are inconsistent"),
    2: (ArithmeticError, "site atom table value out of range"),
    3: (ValueError, "singular basis matrix"),
    4: (ValueError, "parent representative has no Source orbit transport"),
    5: (ValueError, "Source site cosets overlap"),
    6: (ValueError, "Source atom operations are not unique in the supercell"),
    7: (ValueError, "Source representative order is not an atom permutation"),
    8: (KeyError, "rebased atom operation has no Source atom"),
}


def site_atom_table(
    sg: int,
    row_id: int,
    site_pg: int,
    params: Sequence[Fraction],
    basis: Sequence[int],
    representative: Sequence[Fraction] | None,
) -> tuple[
    tuple[tuple[Fraction, Fraction, Fraction], ...],
    tuple[tuple[int, int, int, int, int], ...],
    tuple[int, ...] | None,
]:
    """One parent site's supercell atoms: displayed points, operation records, rebased order."""

    padded = (tuple(params) + (Fraction(0),) * 3)[:3]
    param_words = _rational_words([padded])
    basis_array = np.ascontiguousarray([int(value) for value in basis[:9]], dtype=np.int64)
    b = basis_array.tolist()
    determinant = abs(
        b[0] * (b[4] * b[8] - b[5] * b[7]) - b[1] * (b[3] * b[8] - b[5] * b[6]) + b[2] * (b[3] * b[7] - b[4] * b[6])
    )
    capacity = 48 * max(1, determinant)
    points = np.zeros(12 * capacity, dtype=np.int64)
    records = np.zeros(5 * capacity, dtype=np.int64)
    order = np.zeros(capacity, dtype=np.int64)
    representative_words = (
        None if representative is None else _rational_words([representative])
    )
    count = _LIBRARY.site_atom_table(
        int(sg), int(row_id), int(site_pg), _pointer(param_words, ctypes.c_int64),
        _pointer(basis_array, ctypes.c_int64),
        None if representative_words is None else _pointer(representative_words, ctypes.c_int64),
        capacity, _pointer(points, ctypes.c_int64), _pointer(records, ctypes.c_int64), _pointer(order, ctypes.c_int64),
    )
    if count < 0:
        kind, message = _SITE_ATOM_ERRORS.get(-count, (RuntimeError, f"site atom table failed ({-count})"))
        raise kind(message)
    words = points[: 12 * count].tolist()
    values = [(words[i] << 64) | (words[i + 1] & 0xFFFFFFFFFFFFFFFF) for i in range(0, len(words), 2)]
    return (
        tuple(
            (Fraction(values[i], values[i + 1]), Fraction(values[i + 2], values[i + 3]), Fraction(values[i + 4], values[i + 5]))
            for i in range(0, len(values), 6)
        ),
        tuple(map(tuple, records[: 5 * count].reshape(count, 5).tolist())),
        None if representative is None else tuple(order[:count].tolist()),
    )


_CHILD_SETTING_ERRORS = {
    1: (ValueError, "child setting Source tables are inconsistent"),
    2: (ArithmeticError, "child setting value out of range"),
    3: (ValueError, "child presentation transform is singular"),
    4: (ValueError, "child presentation transform does not define a conventional lattice"),
    5: (ValueError, "child space-group action contains duplicate affine operations"),
    6: (ValueError, "child setting transform collapses affine operations"),
    7: (ValueError, "displayed child action matches no Source inter setting"),
    8: (ValueError, "Source Wyckoff parameters are linearly dependent"),
    9: (ValueError, "Formula15 child site contains duplicate points"),
    10: (ValueError, "Source child-setting Wyckoff identity is not unique"),
    11: (ValueError, "equivalent Source child settings disagree on one Wyckoff identity"),
    12: (ValueError, "child setting capacity exceeded"),
}


def _raise_child_setting(status: int) -> None:
    kind, message = _CHILD_SETTING_ERRORS.get(status, (RuntimeError, f"child setting failed ({status})"))
    raise kind(message)


def _words_to_fractions(words: Sequence[int]) -> tuple[Fraction, ...]:
    values = [(words[i] << 64) | (words[i + 1] & 0xFFFFFFFFFFFFFFFF) for i in range(0, len(words), 2)]
    return tuple(Fraction(values[i], values[i + 1]) for i in range(0, len(values), 2))


_CHILD_LAYOUT_ERRORS = {
    1: (MemoryError, "child layout ran out of memory"),
    2: (ArithmeticError, "child layout value out of range"),
    3: (ValueError, "Formula15 child sites are inconsistent"),
    4: (ValueError, "Formula15 occurrence ordinals are not contiguous"),
    5: (ValueError, "Formula15 does not cover the occurrence product"),
    6: (ValueError, "Formula15 flat occurrence order is inconsistent"),
    7: (ValueError, "Formula15 parent coset identity is inconsistent"),
    8: (ValueError, "Formula15 child-site membership differs from its multiplicity"),
    9: (ValueError, "Formula15 produced duplicate child atom identity"),
    10: (ValueError, "Source atom table does not bijectively cover canonical child atoms"),
    11: (ValueError, "Source atom table contains a nonfinite position"),
    12: (ValueError, "invalid atom operation record"),
    13: (ValueError, "child layout capacity exceeded"),
    14: (ValueError, "Source atom positions and operation records differ"),
    15: (ValueError, "magnetic presentation atoms share one position"),
    16: (ValueError, "magnetic child site does not partition canonical atoms"),
    17: (ValueError, "magnetic child-site representative is not one canonical atom"),
}


def _raise_child_layout(status: int) -> None:
    family, code = divmod(status, 1000)
    if family == 1:
        if code >= 100:
            kind, message = _WYCKOFF_SPLIT_ERRORS.get(code - 100, (RuntimeError, f"Wyckoff split failed ({code - 100})"))
        else:
            kind, message = _CHILD_SITE_ERRORS.get(code, (RuntimeError, f"child sites failed ({code})"))
    elif family == 2:
        kind, message = _CHILD_SETTING_ERRORS.get(code, (RuntimeError, f"child setting failed ({code})"))
    elif family == 3:
        kind, message = _SITE_ATOM_ERRORS.get(code, (RuntimeError, f"site atom table failed ({code})"))
    elif family == 4:
        kind, message = (ValueError, "mode presentation basis is singular" if code == 1 else f"mode presentation failed ({code})")
    elif family == 5:
        kind, message = _ROW_ASSIGNMENT_ERRORS.get(code, (RuntimeError, f"row assignment failed ({code})"))
    elif family == 6:
        raise _magnetic_error(code)
    else:
        kind, message = _CHILD_LAYOUT_ERRORS.get(code, (RuntimeError, f"child layout failed ({status})"))
    raise kind(message)


def _determinant3(values: Sequence[int]) -> int:
    b = [int(value) for value in values[:9]]
    return abs(
        b[0] * (b[4] * b[8] - b[5] * b[7]) - b[1] * (b[3] * b[8] - b[5] * b[6]) + b[2] * (b[3] * b[7] - b[4] * b[6])
    )


_SUBGROUP_OPS_ERRORS = {
    1: (KeyError, "operation record not represented in the parent space group"),
    2: (ValueError, "singular subgroup basis"),
    3: (RuntimeError, "generate_domains2 point-coset scan did not close"),
    4: (ArithmeticError, "subgroup operation value out of range"),
    5: (ValueError, "subgroup operation capacity exceeded"),
    6: (ZeroDivisionError, "operation record has a zero denominator"),
    7: (ValueError, "subgroup operation Source tables are inconsistent"),
}

Record = tuple[int, int, int, int, int]


def _subgroup_ops_records(status: int, out: np.ndarray) -> tuple[Record, ...]:
    if status < 0:
        kind, message = _SUBGROUP_OPS_ERRORS.get(-status, (RuntimeError, f"subgroup operation failed ({-status})"))
        raise kind(message)
    return tuple(map(tuple, out[:status].tolist()))


def _basis_array(basis: Sequence[int]) -> np.ndarray:
    return np.array([int(value) for value in basis[:9]], dtype=np.int64)


def _operation_capacity(basis: Sequence[int]) -> int:
    return 192 * max(1, _determinant3(basis))


def mapped_subgroup_records(
    parent_sg: int, child_sg: int, basis: Sequence[int], origin: Sequence[int]
) -> tuple[Record, ...]:
    """map_subgroup2_: the child's operations as parent operation records."""

    out = np.zeros((192, 5), dtype=np.int64)
    status = _LIBRARY.mapped_subgroup_records(
        int(parent_sg), int(child_sg), _pointer(_basis_array(basis), ctypes.c_int64),
        _pointer(np.array([int(value) for value in origin[:4]], dtype=np.int64), ctypes.c_int64), len(out),
        _pointer(out, ctypes.c_int64),
    )
    return _subgroup_ops_records(status, out)


_ID_SUBGROUP_ERRORS = {
    -1: (ValueError, "rowop matrix dump is too narrow"),
    -2: (ArithmeticError, "id_subgroup_ value out of range"),
    -3: (ValueError, "id_subgroup_ matrix did not reduce to integer"),
    -4: (ValueError, "id_subgroup_ Source tables are inconsistent"),
    -5: (KeyError, "id_subgroup_ operation has no Source record"),
    -6: (ValueError, "id_subgroup_ input exceeds capacity"),
    -7: (ZeroDivisionError, "id_subgroup_ operation with zero denominator"),
}


def _id_subgroup_call(function, parent: int, basis: Sequence[int], operations: Sequence[Sequence[int]], flag: int, *extra):
    basis_array = np.ascontiguousarray([int(value) for value in basis[:9]], dtype=np.int64)
    records = [tuple(int(value) for value in record[:5]) for record in operations if len(record) >= 5]
    ops = np.ascontiguousarray(records if records else np.zeros((0, 5)), dtype=np.int64).reshape(-1, 5)
    status = function(
        int(parent), _pointer(basis_array, ctypes.c_int64), len(ops), _pointer(ops, ctypes.c_int64), int(flag), *extra
    )
    if status < 0:
        kind, message = _ID_SUBGROUP_ERRORS.get(status, (RuntimeError, f"native id_subgroup_ status {status}"))
        raise kind(message)
    return status


def id_subgroup_identify(
    parent: int, basis: Sequence[int], operations: Sequence[Sequence[int]], flag: int
) -> tuple[int, tuple[int, ...], tuple[int, int, int, int]] | None:
    """id_subgroup_: (subgroup, basis, origin) of the first accepted candidate, or None."""

    out = np.zeros(14, dtype=np.int64)
    if not _id_subgroup_call(_LIBRARY.id_subgroup_identify, parent, basis, operations, flag, _pointer(out, ctypes.c_int64)):
        return None
    values = out.tolist()
    return values[0], tuple(values[1:10]), tuple(values[10:14])


def match_inter_wyckoff(
    sg: int,
    setting_id: int,
    fract: Sequence[float],
    source_values: Sequence[float],
    wanted_multiplicity: int | None,
    tol: float,
) -> tuple[int, tuple[float, float, float], tuple[bool, bool, bool], tuple[float, float, float]] | None:
    """First Wyckoff row offset of an inter setting holding fract, with its parameters."""

    row = np.zeros(1, dtype=np.intc)
    params = np.zeros(3, dtype=float)
    mask = np.zeros(3, dtype=np.intc)
    representative = np.zeros(3, dtype=float)
    status = _LIBRARY.match_inter_wyckoff(
        int(sg), int(setting_id), _pointer(np.array(fract, dtype=float), ctypes.c_double),
        _pointer(np.array(source_values, dtype=float), ctypes.c_double), int(wanted_multiplicity or 0), float(tol),
        _pointer(row, ctypes.c_int), _pointer(params, ctypes.c_double), _pointer(mask, ctypes.c_int),
        _pointer(representative, ctypes.c_double),
    )
    if status < 0:
        raise ValueError(f"Wyckoff site match failed ({status})")
    if not status:
        return None
    return int(row[0]), tuple(params.tolist()), tuple(bool(value) for value in mask), tuple(representative.tolist())


def k_catalog(sg: int) -> dict:
    """The k-point catalog of a space group (k_catalog.c)."""

    _LIBRARY.k_catalog_json.restype = ctypes.c_char_p
    text = _LIBRARY.k_catalog_json(int(sg))
    if text is None:
        raise ValueError("k-point catalog Source tables are inconsistent")
    return json.loads(text.decode())


def k_coordinate_map(sg: int, kslot: int) -> tuple[tuple[str, ...], tuple[Fraction, ...], tuple[tuple[Fraction, ...], ...]]:
    """The affine map behind a K-family's displayed expression: (names, origin, columns)."""

    out = np.zeros(29, dtype=np.int64)
    status = _LIBRARY.k_coordinate_map_words(int(sg), int(kslot), _pointer(out, ctypes.c_int64))
    if status:
        raise (ArithmeticError if status == 3 else ValueError)("K-family display expression is not affine rational")
    o = out.tolist()
    count = o[1]
    return (
        tuple(chr(o[2 + i]) for i in range(count)),
        tuple(Fraction(o[5 + 2 * axis], o[6 + 2 * axis]) for axis in range(3)),
        tuple(tuple(Fraction(o[11 + 6 * i + 2 * axis], o[12 + 6 * i + 2 * axis]) for axis in range(3)) for i in range(count)),
    )


def k_source_kparam_of_visible(sg: int, kslot: int, visible: Sequence[Fraction]) -> tuple[int, int, int, int] | None:
    """Source's internal K record of a visible K vector, None when the slot has none."""

    words = np.array([part for value in visible[:3] for part in (Fraction(value).numerator, Fraction(value).denominator)], dtype=np.int64)
    out = np.zeros(4, dtype=np.int64)
    found = _LIBRARY.k_source_kparam_of_visible(int(sg), int(kslot), _pointer(words, ctypes.c_int64), _pointer(out, ctypes.c_int64))
    if found < 0:
        raise (ArithmeticError if found == -2 else ValueError)("k-vector Source tables are inconsistent")
    return tuple(out.tolist()) if found else None  # type: ignore[return-value]


def k_numeric_star(
    parent_sg: int, gid: int, visible: Sequence[Fraction] | None, source_visible: Sequence[Fraction] | None
) -> list[str]:
    """The numeric displayed star arms of a parametric gid at a visible k (k_catalog.c)."""

    def words(values: Sequence[Fraction] | None) -> np.ndarray:
        return np.array([part for value in (values or ())[:3] for part in (Fraction(value).numerator, Fraction(value).denominator)] or [0], dtype=np.int64)

    _LIBRARY.k_numeric_star_json.restype = ctypes.c_char_p
    text = _LIBRARY.k_numeric_star_json(
        int(parent_sg), int(gid), int(visible is not None), _pointer(words(visible), ctypes.c_int64),
        int(source_visible is not None), _pointer(words(source_visible), ctypes.c_int64),
    )
    if text is None:
        raise ValueError("k-vector Source tables are inconsistent")
    return json.loads(text.decode())


def k_display_kvector(sg: int, kslot: int) -> str:
    _LIBRARY.k_display_kvector.restype = ctypes.c_char_p
    text = _LIBRARY.k_display_kvector(int(sg), int(kslot))
    if text is None:
        raise ValueError("k-vector Source tables are inconsistent")
    return text.decode()


def k_display_star(sg: int, kslot: int, dim: int) -> list[str]:
    _LIBRARY.k_display_star_json.restype = ctypes.c_char_p
    text = _LIBRARY.k_display_star_json(int(sg), int(kslot), int(dim))
    if text is None:
        raise ValueError("k-vector Source tables are inconsistent")
    return json.loads(text.decode())


_BLOCK_PLAN_ERRORS = {
    1: "mode topology requires Source rows and centerings",
    2: "mode topology does not cover every canonical child atom",
    3: "mode presentation basis is singular",
    4: "mode rows do not cover the canonical presentation layout",
    5: "mode rows do not bijectively cover canonical child atoms",
}


_MODE_TEXT_ERRORS = {
    1: (ValueError, "mode label Source tables are inconsistent"),
    2: (ArithmeticError, "k-vector value out of range"),
    3: (ArithmeticError, "K-family display expression is not affine rational"),
    4: (IndexError, "index out of range"),
    5: (KeyError, "direction_matrix"),
    7: (MemoryError, "mode label text ran out of memory"),
    8: (np.linalg.LinAlgError, "SVD did not converge"),
}


def _raise_mode_text(status: int) -> None:
    kind, message = _MODE_TEXT_ERRORS.get(status, (RuntimeError, f"mode label failed ({status})"))
    raise kind(message)


_MODE_DETAIL_ERRORS: dict[int, type[Exception]] = {
    1: ValueError, 2: KeyError, 3: IndexError, 4: TypeError, 5: ZeroDivisionError, 6: ArithmeticError,
    7: MemoryError, 8: RuntimeError, 9: np.linalg.LinAlgError,
}


def _raise_mode_detail(kind: int, status: int, message: str) -> None:
    if kind == 100:
        _raise_child_layout(status)
    if kind == 101:
        _raise_child_setting(status)
    if kind == 102:
        _raise_orderparam(status)
    if kind == 103:
        raise _magnetic_error(status)
    if kind in (110, 111, 112):
        table = {110: _RECORDWISE_ERRORS, 111: _PROJECTION_ENTRY_ERRORS, 112: _MODE_BLOCK_ERRORS}[kind]
        error, text = table.get(status, (RuntimeError, f"mode block failed ({status})"))
        raise error(text)
    if kind == 113:
        if status >= 100:
            error, text = _ROW_ASSIGNMENT_ERRORS.get(status - 100, (RuntimeError, f"row assignment failed ({status - 100})"))
            raise error(text)
        if status == 6:
            raise MemoryError("block plan ran out of memory")
        raise ValueError(_BLOCK_PLAN_ERRORS.get(status, f"block plan failed ({status})"))
    error = _MODE_DETAIL_ERRORS.get(kind, RuntimeError)
    if error is KeyError and len(message) >= 2 and message[0] == message[-1] == "'":
        message = message[1:-1]
    raise error(message)


class OpdNotFoundError(ValueError):
    """The requested OPD is absent from the computed candidate rows."""


def _raise_selection(kind: int, status: int, message: str) -> None:
    if kind == 120:
        _raise_isotropy_rows(status)
    if kind == 102:
        _raise_orderparam(status)
    if kind == 103:
        raise _magnetic_error(status)
    tables = {111: _PROJECTION_ENTRY_ERRORS, 121: _ID_SUBGROUP_ERRORS, 122: _SUBGROUP_OPS_ERRORS, 124: _BRIDGE_WEIGHT_ERRORS}
    if kind in tables:
        error, text = tables[kind].get(status, (RuntimeError, f"selection kernel failed ({kind}, {status})"))
        raise error(text)
    if kind == 123:
        _raise_coupled_domain_orbit(status)
    if kind == 125:
        raise OverflowError(message)
    if kind == 126:
        raise OpdNotFoundError(message)
    _raise_mode_detail(kind, status, message)


def _selection_call(request: dict) -> dict:
    data = json.dumps(request, default=list).encode()
    _LIBRARY.selection_json.restype = ctypes.c_char_p
    text = _LIBRARY.selection_json(ctypes.c_char_p(data), ctypes.c_int64(len(data)))
    if text is None:
        raise MemoryError("selection stage ran out of memory")
    out = json.loads(text)
    if len(out) == 1 and "error" in out:
        error = out["error"]
        _raise_selection(int(error["kind"]), int(error["status"]), str(error["message"]))
    return out


def k_points_state(parent_state: dict, **options) -> dict:
    """select_k_points of one parent state (sel_kpoints.c)."""

    return _selection_call({"call": "k_points", "state": parent_state, "options": options})


def _parent_call(request: dict) -> dict:
    data = json.dumps(request, default=list).encode()
    _LIBRARY.parent_json.restype = ctypes.c_char_p
    text = _LIBRARY.parent_json(ctypes.c_char_p(data), ctypes.c_int64(len(data)))
    if text is None:
        raise MemoryError("parent stage ran out of memory")
    out = json.loads(text)
    if len(out) == 1 and "error" in out:
        error = out["error"]
        _raise_selection(int(error["kind"]), int(error["status"]), str(error["message"]))
    return out


def _cif_summary(data: bytes, from_file: bool) -> dict:
    _LIBRARY.cif_summary_json.restype = ctypes.c_char_p
    text = _LIBRARY.cif_summary_json(ctypes.c_char_p(data), ctypes.c_int64(len(data)), int(from_file))
    if text is None:
        raise MemoryError("CIF reader ran out of memory")
    out = json.loads(text)
    if len(out) == 1 and "error" in out:
        error = out["error"]
        raise {1: ValueError, 3: IndexError}.get(int(error["kind"]), RuntimeError)(str(error["message"]))
    return out


def read_cif_summary(cif: os.PathLike | str) -> dict:
    """The parent-stage CIF summary of one file (par_cif.c)."""

    return _cif_summary(str(cif).encode(), True)


def read_cif_summary_from_text(cif_text: str) -> dict:
    return _cif_summary(cif_text.encode(), False)


def build_parent_state_from_cif_info(cif_info: dict) -> dict:
    """The parent International setting and every site's Wyckoff position (par_state.c)."""

    return _parent_call({"call": "parent_state", "cif_info": cif_info})


def _uses_library(code) -> bool:
    return "_LIBRARY" in code.co_names or any(
        _uses_library(const) for const in code.co_consts if inspect.iscode(const)
    )


def _serialized(function):
    @functools.wraps(function)
    def call(*args, **kwargs):
        with _LIBRARY_LOCK:
            return function(*args, **kwargs)

    return call


# The C core keeps process-wide caches and result buffers, so callers on
# different threads (the GUI server) take turns.
_LIBRARY_LOCK = threading.RLock()
for _name, _value in list(globals().items()):
    if inspect.isfunction(_value) and _value.__module__ == __name__ and _uses_library(_value.__code__):
        globals()[_name] = _serialized(_value)
