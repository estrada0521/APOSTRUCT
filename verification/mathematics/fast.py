"""The Local theorems of mode_basis, group_invariance and fixed_space through native/checks.c: the same statuses,
from one mode-details payload, without building the Python definition objects. Group presentations still come from
group_invariance."""

from __future__ import annotations

import ctypes
import json
from functools import lru_cache
from pathlib import Path

from .group_invariance import (
    _DECLARED_CLAIM,
    _DERIVED_CLAIM,
    group_authority,
    presentation_authority,
    subgroup_label_from_local_payload,
)

PUBLIC = Path(__file__).resolve().parents[2]


@lru_cache(maxsize=1)
def _library() -> ctypes.CDLL:
    path = PUBLIC / "bin" / "_mathcheck.so"
    if not path.is_file():
        raise ImportError(f"the mathematics checks are not built; run sh {PUBLIC / 'build.sh'}")
    lib = ctypes.CDLL(str(path))
    lib.mc_open.argtypes = (ctypes.c_char_p, ctypes.c_int64)
    lib.mc_open.restype = ctypes.c_void_p
    for name in ("mc_error", "mc_subgroup"):
        getattr(lib, name).argtypes = (ctypes.c_void_p,)
        getattr(lib, name).restype = ctypes.c_char_p
    for name in ("mc_basis", "mc_definitions"):
        getattr(lib, name).argtypes = (ctypes.c_void_p,)
        getattr(lib, name).restype = ctypes.c_int
    words = ctypes.POINTER(ctypes.c_int64)
    lib.mc_invariance.argtypes = (ctypes.c_void_p, words)
    lib.mc_invariance.restype = ctypes.c_int
    lib.mc_fixed.argtypes = (ctypes.c_void_p, words, ctypes.POINTER(ctypes.c_int))
    lib.mc_fixed.restype = ctypes.c_int
    lib.mc_close.argtypes = (ctypes.c_void_p,)
    return lib


def _encode(settings) -> ctypes.Array:
    words = [len(settings)]
    for setting in settings:
        operations = setting.nonidentity_operations
        words.append(len(operations))
        for operation in operations:
            words.extend(value for row in operation.rotation for value in row)
            words.extend(t.numerator for t in operation.translation)
            words.extend(t.denominator for t in operation.translation)
            words.append(operation.time_reversal)
    return (ctypes.c_int64 * len(words))(*words)


def assess_payload(text: bytes) -> dict[str, object]:
    """Statuses of the three theorems for one payload holding `preview.selected.mode_details` or
    `selected.mode_details`; raises ValueError where the Python checks reject the payload."""

    lib = _library()
    handle = lib.mc_open(text, len(text))
    try:
        error = lib.mc_error(handle)
        if error:
            raise ValueError(error.decode())
        subgroup = json.loads(lib.mc_subgroup(handle))
        label = subgroup_label_from_local_payload(
            {"preview": {"selected": {"mode_details": {"subgroup_details": {"magnetic_subgroup": subgroup}}}}}
        )
        basis = lib.mc_basis(handle)
        if basis < 0:
            raise ValueError(lib.mc_error(handle).decode())
        authority = group_authority(label)
        declared = authority.settings
        derived = lambda: presentation_authority(label, _DERIVED_CLAIM).settings[len(declared):]
        invariance, keys = "refuted" if authority.setting_coverage_complete else "indeterminate", []
        for settings in (declared, None):
            settings = settings if settings is not None else derived()
            found = lib.mc_invariance(handle, _encode(settings))
            if found == -3:
                raise ValueError(lib.mc_error(handle).decode())
            if found == -2:
                invariance = "not_applicable"
                break
            if found >= 0:
                invariance, keys = "satisfied", [settings[found]]
                break
        refuted, fixed = ctypes.c_int(0), "indeterminate"
        done = {setting.setting_key for setting in keys}
        for settings in (keys + [s for s in declared if s.setting_key not in done], None):
            if settings is None:
                done |= {setting.setting_key for setting in declared}
                settings = [s for s in derived() if s.setting_key not in done]
            found = lib.mc_fixed(handle, _encode(settings), ctypes.byref(refuted))
            if found == -3:
                raise ValueError(lib.mc_error(handle).decode())
            if found == -2:
                fixed = "not_applicable"
                break
            if found >= 0:
                fixed = "refuted" if refuted.value else "satisfied"
                break
        return {"basis": ("satisfied", "refuted")[basis], "invariance": invariance, "fixed": fixed,
                "definitions": lib.mc_definitions(handle)}
    finally:
        lib.mc_close(handle)


if __name__ == "__main__":
    import gzip
    import sys

    for name in sys.argv[1:]:
        data = Path(name).read_bytes()
        print(name, json.dumps(assess_payload(gzip.decompress(data) if data[:2] == b"\x1f\x8b" else data)))
