"""Read-only, key-lazy access to the baked Source blob (native/source.blob)."""

from __future__ import annotations

from collections.abc import Iterator, Sequence
import hashlib
import mmap
from pathlib import Path
import struct
from typing import Any

import numpy as np


TABLE_NAMES = ("space", "irreps", "little", "isotropy", "wyckoff", "images", "magnetic")
AUTHORITY_FILES = (
    "data_space.txt",
    "data_irreps.txt",
    "data_little.txt",
    "data_isotropy.txt",
    "data_wyckoff.txt",
    "data_images.txt",
    "data_magnetic.txt",
    "const.dat",
)
_KINDS = {1: "<i1", 2: "<i2", 3: "<i4", 4: "<i8", 5: "<f8", 6: "|b1"}


def file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


class CompiledSequence(Sequence[Any]):
    """One immutable Source section of the blob with cached Python scalar values."""

    def __init__(self, values: tuple[Any, ...]):
        self._values = values

    def _load(self) -> tuple[Any, ...]:
        return self._values

    def __len__(self) -> int:
        return len(self._values)

    def __getitem__(self, index: int | slice) -> Any:
        if isinstance(index, slice):
            return list(self._values[index])
        return self._values[index]

    def __iter__(self) -> Iterator[Any]:
        return iter(self._values)


class CompiledTable(dict[str, Sequence[Any]]):
    """Decode individual table sections only when first requested."""

    def __init__(self, blob: "SourceBlob", table: str, names: Sequence[str]):
        self._blob = blob
        self._table = table
        self._names = tuple(names)
        super().__init__()

    def __missing__(self, key: str) -> Sequence[Any]:
        if key not in self._names:
            raise KeyError(key)
        sequence = CompiledSequence(self._blob.section(f"{self._table}.{key}"))
        dict.__setitem__(self, key, sequence)
        return sequence

    def __len__(self) -> int:
        return len(self._names)

    def __iter__(self) -> Iterator[str]:
        return iter(self._names)

    def __contains__(self, key: object) -> bool:
        return key in self._names

    def get(self, key: str, default: Any = None) -> Any:
        return self[key] if key in self._names else default

    def keys(self):
        return self._names

    def items(self):
        return ((key, self[key]) for key in self._names)

    def values(self):
        return (self[key] for key in self._names)

    @staticmethod
    def _immutable(*_args: Any, **_kwargs: Any) -> None:
        raise TypeError("compiled Source tables are read-only")

    __setitem__ = _immutable
    __delitem__ = _immutable
    clear = _immutable
    pop = _immutable
    popitem = _immutable
    setdefault = _immutable
    update = _immutable
    __ior__ = _immutable

class SourceBlob:
    """The mmapped blob baked from Source by native/generator/source_blob_build."""

    def __init__(self, path: Path):
        self.path = Path(path)
        try:
            with self.path.open("rb") as stream:
                self._data = mmap.mmap(stream.fileno(), 0, access=mmap.ACCESS_READ)
        except OSError as exc:
            raise ValueError(f"cannot open the Source blob {self.path}; run build.sh") from exc
        magic, version, count, index_offset, size = struct.unpack_from("<8sIIQQ", self._data, 0)
        if magic != b"APSTBLB1" or version != 1 or size != len(self._data):
            raise ValueError(f"unsupported Source blob: {self.path}")
        self._index: dict[str, tuple[int, int, int]] = {}
        for position in range(count):
            name, kind, _reserved, length, offset, _bytes = struct.unpack_from(
                "<96sIIQQQ", self._data, index_offset + 128 * position
            )
            self._index[name.split(b"\0")[0].decode()] = (kind, length, offset)

    def names(self, table: str) -> tuple[str, ...]:
        prefix = table + "."
        return tuple(name[len(prefix):] for name in self._index if name.startswith(prefix))

    def section(self, name: str) -> tuple[Any, ...]:
        kind, length, offset = self._index[name]
        if kind == 7:
            offsets = np.frombuffer(self._data, "<i8", length + 1, offset).tolist()
            pool = offset + 8 * (length + 1)
            return tuple(
                self._data[pool + offsets[index]:pool + offsets[index + 1]].decode()
                for index in range(length)
            )
        return tuple(np.frombuffer(self._data, _KINDS[kind], length, offset).tolist())


class CompiledSource:
    """Verify the Source authority against the blob and expose its six tables lazily."""

    def __init__(self, source: Path, blob: Path):
        self.source = Path(source)
        self.blob = SourceBlob(blob)
        recorded = dict(line.split("=", 1) for line in self.blob.section("source.sha256"))
        if set(recorded) != set(AUTHORITY_FILES):
            raise ValueError("Source blob records an unexpected authority set")
        for name, digest in recorded.items():
            if file_sha256(self.source / name) != digest:
                raise ValueError(f"Source authority hash mismatch: {name}; run build.sh")
        self.tables = {name: CompiledTable(self.blob, name, self.blob.names(name)) for name in TABLE_NAMES}

    def table(self, name: str) -> CompiledTable:
        return self.tables[name]
