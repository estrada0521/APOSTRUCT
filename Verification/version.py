"""Content identity for the standalone single-branch comparator."""

from __future__ import annotations

import hashlib
from importlib import metadata
from pathlib import Path
import sys

from Verification.result import COMPARATOR_VERSION as COMPARATOR_SCHEMA


_DIRECT_JUDGE_DEPENDENCIES = (
    "APOSTRUCT/Backend/exactmath.py",
    "APOSTRUCT/Backend/fraction_expression.py",
    "APOSTRUCT/Backend/source/iso_data.py",
    "APOSTRUCT/Backend/source/tables.py",
    "APOSTRUCT/Backend/source/magnetic.py",
    "APOSTRUCT/Backend/reciprocal/catalog.py",
    "APOSTRUCT/Backend/isotropy/catalog.py",
    "APOSTRUCT/Backend/modes/mode_detail_text.py",
    "APOSTRUCT/Backend/modes/engine/decoder.py",
    "APOSTRUCT/Backend/modes/engine/subgroup_structure/presentation_transport.py",
    "APOSTRUCT/Backend/modes/structure/magnetic_wyckoff.py",
)

_JUDGE_DISTRIBUTIONS = ("gemmi", "numpy", "scipy", "sympy")


def _comparator_sources(project_root: Path) -> tuple[Path, ...]:
    root = Path(project_root).resolve()
    verification = root / "Verification"
    candidates = [verification / "compare.py", verification / "result.py"]
    for directory in ("comparison", "parsers"):
        candidates.extend((verification / directory).rglob("*.py"))
    candidates.extend(root / relative for relative in _DIRECT_JUDGE_DEPENDENCIES)
    missing = tuple(path for path in candidates if not path.is_file())
    if missing:
        names = ", ".join(path.relative_to(root).as_posix() for path in missing)
        raise ValueError(f"missing comparator source files: {names}")
    return tuple(
        sorted(
            candidates,
            key=lambda path: path.relative_to(root).as_posix(),
        )
    )


def _environment_identity() -> tuple[str, ...]:
    versions = [
        f"python={sys.version_info.major}.{sys.version_info.minor}.{sys.version_info.micro}"
    ]
    for distribution in _JUDGE_DISTRIBUTIONS:
        try:
            value = metadata.version(distribution)
        except metadata.PackageNotFoundError:
            value = "missing"
        versions.append(f"{distribution}={value}")
    return tuple(versions)


def comparator_version(project_root: Path) -> str:
    """Return a path-independent identity for all executable judge sources."""

    root = Path(project_root).resolve()
    verification = root / "Verification"
    sources = _comparator_sources(project_root)
    if not sources:
        raise ValueError(f"no comparator sources under {verification}")
    digest = hashlib.sha256()
    digest.update(b"isodistort.validation.comparator.v1\0")
    for value in _environment_identity():
        digest.update(value.encode("utf-8"))
        digest.update(b"\0")
    for path in sources:
        relative = path.relative_to(root).as_posix().encode("utf-8")
        digest.update(relative)
        digest.update(b"\0")
        digest.update(path.read_bytes())
        digest.update(b"\0")
    return f"{COMPARATOR_SCHEMA}+{digest.hexdigest()[:12]}"


def require_unchanged_comparator(project_root: Path, expected: str) -> None:
    """Reject a run whose judge sources changed before persistence."""

    observed = comparator_version(project_root)
    if observed != expected:
        raise RuntimeError(
            "Validation comparator sources changed during execution; "
            "no results were persisted"
        )


__all__ = [
    "COMPARATOR_SCHEMA",
    "comparator_version",
    "require_unchanged_comparator",
]
