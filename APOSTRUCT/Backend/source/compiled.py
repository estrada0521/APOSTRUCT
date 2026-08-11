"""Read-only, key-lazy access to compiled Source table arrays."""

from __future__ import annotations

from collections.abc import Iterator, Mapping, Sequence
import hashlib
import json
from pathlib import Path
from typing import Any

import numpy as np


SCHEMA = "APOSTRUCT.source-compiled.v1"
TABLE_NAMES = ("space", "irreps", "little", "isotropy", "wyckoff", "images")


def file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _python_value(value: Any) -> Any:
    return value.item() if isinstance(value, np.generic) else value


class CompiledSequence(Sequence[Any]):
    """A key-lazy immutable Source section with cached Python scalar values."""

    def __init__(self, root: Path, entry: Mapping[str, Any]):
        self._root = root
        self._entry = entry
        self._mapped_values: np.ndarray | None = None
        self._values: tuple[Any, ...] | None = None

    def _load(self) -> tuple[Any, ...]:
        if self._values is not None:
            return self._values
        path = self._root / str(self._entry["file"])
        if file_sha256(path) != self._entry["sha256"]:
            raise ValueError(f"compiled Source section hash mismatch: {path.name}")
        mapped = np.load(path, mmap_mode="r", allow_pickle=False)
        if mapped.dtype.str != self._entry["dtype"] or list(mapped.shape) != self._entry["shape"]:
            raise ValueError(f"compiled Source section metadata mismatch: {path.name}")
        mapped.flags.writeable = False
        self._mapped_values = mapped
        self._values = (
            tuple(mapped.tolist())
            if mapped.ndim == 1
            else tuple(_python_value(value) for value in mapped)
        )
        return self._values

    def __len__(self) -> int:
        return len(self._load())

    def __getitem__(self, index: int | slice) -> Any:
        values = self._load()
        if isinstance(index, slice):
            return list(values[index])
        return values[index]

    def __iter__(self) -> Iterator[Any]:
        return iter(self._load())


class CompiledTable(dict[str, Sequence[Any]]):
    """Load and verify individual table sections only when first requested."""

    def __init__(self, root: Path, entries: Mapping[str, Mapping[str, Any]]):
        self._root = root
        self._entries = dict(entries)
        super().__init__()

    def __missing__(self, key: str) -> Sequence[Any]:
        try:
            entry = self._entries[key]
        except KeyError:
            raise KeyError(key) from None
        sequence = CompiledSequence(self._root, entry)
        sequence._load()
        dict.__setitem__(self, key, sequence)
        return sequence

    def __len__(self) -> int:
        return len(self._entries)

    def __iter__(self) -> Iterator[str]:
        return iter(self._entries)

    def __contains__(self, key: object) -> bool:
        return key in self._entries

    def get(self, key: str, default: Any = None) -> Any:
        return self[key] if key in self._entries else default

    def keys(self):
        return self._entries.keys()

    def items(self):
        return ((key, self[key]) for key in self._entries)

    def values(self):
        return (self[key] for key in self._entries)

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

    @property
    def loaded_keys(self) -> tuple[str, ...]:
        return tuple(dict.__iter__(self))


class CompiledSource:
    """Verify Source authority and expose its six compiled tables lazily."""

    def __init__(self, source: Path, compiled: Path):
        self.source = Path(source)
        self.compiled = Path(compiled)
        manifest_path = self.compiled / "manifest.json"
        try:
            manifest = json.loads(manifest_path.read_text())
        except (OSError, json.JSONDecodeError) as exc:
            raise ValueError(f"cannot read compiled Source manifest: {manifest_path}") from exc
        if manifest.get("schema") != SCHEMA:
            raise ValueError(f"unsupported compiled Source manifest: {manifest_path}")
        source_entries = manifest.get("source", {})
        if set(source_entries) != {
            "data_space.txt",
            "data_irreps.txt",
            "data_little.txt",
            "data_isotropy.txt",
            "data_wyckoff.txt",
            "data_images.txt",
            "const.dat",
        }:
            raise ValueError("compiled Source manifest has an unexpected authority set")
        for name, entry in source_entries.items():
            path = self.source / name
            if file_sha256(path) != entry["sha256"]:
                raise ValueError(f"Source authority hash mismatch: {name}")
        tables = manifest.get("tables", {})
        if set(tables) != set(TABLE_NAMES):
            raise ValueError("compiled Source manifest has an unexpected table set")
        self.tables = {
            name: CompiledTable(self.compiled, tables[name]) for name in TABLE_NAMES
        }

    def table(self, name: str) -> CompiledTable:
        return self.tables[name]
