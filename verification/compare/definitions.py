"""The mode-definition comparison (native/definitions.c): Web definition tables against one apo output's definitions."""

from __future__ import annotations

import ctypes
import json
from functools import lru_cache
from pathlib import Path
from typing import Any, Sequence

PUBLIC = Path(__file__).resolve().parents[2]
KINDS = ("dsp", "mag", "occ", "strain")


@lru_cache(maxsize=1)
def _library() -> ctypes.CDLL:
    path = PUBLIC / "bin" / "_definitions.so"
    if not path.is_file():
        raise ImportError(f"the definition comparison is not built; run sh {PUBLIC / 'build.sh'}")
    lib = ctypes.CDLL(str(path))
    lib.vd_open.argtypes = (ctypes.c_char_p, ctypes.c_int64, ctypes.c_char_p, ctypes.c_int64, ctypes.c_double)
    lib.vd_open.restype = ctypes.c_void_p
    for name in ("vd_error", "vd_result"):
        getattr(lib, name).argtypes = (ctypes.c_void_p,)
        getattr(lib, name).restype = ctypes.c_char_p
    for name in ("vd_status", "vd_strict_status"):
        getattr(lib, name).argtypes = (ctypes.c_void_p, ctypes.c_int)
        getattr(lib, name).restype = ctypes.c_char_p
    for name in ("vd_web_count", "vd_accepted"):
        getattr(lib, name).argtypes = (ctypes.c_void_p, ctypes.c_int)
        getattr(lib, name).restype = ctypes.c_int
    lib.vd_transport.argtypes = (
        ctypes.c_void_p,
        ctypes.c_int,
        ctypes.POINTER(ctypes.c_int64),
        ctypes.POINTER(ctypes.c_double),
    )
    lib.vd_transport.restype = ctypes.c_int
    lib.vd_close.argtypes = (ctypes.c_void_p,)
    return lib


class Definitions:
    """One comparison; the local text is a mode-details object, or an apo output or Local artifact holding one."""

    def __init__(self, web_text: str, local_json: bytes, tol: float) -> None:
        self._lib = _library()
        web = web_text.encode("utf-8")
        self._handle = self._lib.vd_open(web, len(web), local_json, len(local_json), float(tol))
        error = self._lib.vd_error(self._handle)
        if error is not None:
            self.close()
            raise ValueError(error.decode())
        self._result: dict[str, Any] | None = None

    def close(self) -> None:
        if self._handle:
            self._lib.vd_close(self._handle)
            self._handle = None

    def __del__(self) -> None:
        self.close()

    def web_count(self, kind: str) -> int:
        return self._lib.vd_web_count(self._handle, KINDS.index(kind))

    def accepted(self, kind: str) -> bool:
        return bool(self._lib.vd_accepted(self._handle, KINDS.index(kind)))

    def transport(self, kind: str, basis: Sequence[Sequence[int]], shift: Sequence[Any]) -> bool:
        """Adopt the definitions carried into the request child frame when they close there."""

        words = (ctypes.c_int64 * 9)(*(int(value) for row in basis for value in row))
        offsets = (ctypes.c_double * 3)(*(float(value) for value in shift))
        return bool(self._lib.vd_transport(self._handle, KINDS.index(kind), words, offsets))

    def payload(self, kind: str) -> dict[str, Any]:
        index = KINDS.index(kind)
        return {
            "status": self._lib.vd_status(self._handle, index).decode(),
            "strict_status": self._lib.vd_strict_status(self._handle, index).decode(),
        }

    def mode_definition(self) -> tuple[dict[str, Any], str | None, bool, bool]:
        """The mode_definition layers, first divergence, physical and strict verdicts."""

        if self._result is None:
            self._result = json.loads(self._lib.vd_result(self._handle))
        if "error" in self._result:
            raise ValueError(self._result["error"])
        return (
            self._result["layers"],
            self._result["divergence"],
            self._result["physical"],
            self._result["strict"],
        )


__all__ = ["Definitions", "KINDS"]
