"""Content identity for the independent mathematics verifier."""

from __future__ import annotations

import hashlib
from importlib import metadata
from pathlib import Path
import sys


MATH_VERIFIER_SCHEMA = "mathematics.v4"
_DATA_DEPENDENCIES = (
    "Verification/mathematics/data/README.md",
    "Verification/mathematics/data/SPGLIB-LICENSE.txt",
    "Verification/mathematics/data/magnetic_group_generators.json.gz",
    "Verification/mathematics/data/upstream/spglib-2.7.0/magnetic_hall_symbols.yaml",
    "Verification/mathematics/data/upstream/spglib-2.7.0/msg_numbers.csv",
)
_DISTRIBUTIONS = ("numpy", "sympy")


def _verifier_sources(project_root: Path) -> tuple[Path, ...]:
    root = Path(project_root).resolve()
    directory = root / "Verification" / "mathematics"
    candidates = []
    candidates.extend(directory.rglob("*.py"))
    candidates.extend(root / relative for relative in _DATA_DEPENDENCIES)
    missing = tuple(path for path in candidates if not path.is_file())
    if missing:
        names = ", ".join(path.relative_to(root).as_posix() for path in missing)
        raise ValueError(f"missing mathematics verifier source files: {names}")
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
    for distribution in _DISTRIBUTIONS:
        try:
            value = metadata.version(distribution)
        except metadata.PackageNotFoundError:
            value = "missing"
        versions.append(f"{distribution}={value}")
    return tuple(versions)


def math_verifier_version(project_root: Path) -> str:
    """Return a path-independent identity for the executable proof boundary."""

    root = Path(project_root).resolve()
    digest = hashlib.sha256()
    digest.update(b"isodistort.validation.mathematics.verifier.v4\0")
    for value in _environment_identity():
        digest.update(value.encode("utf-8"))
        digest.update(b"\0")
    for path in _verifier_sources(root):
        digest.update(path.relative_to(root).as_posix().encode("utf-8"))
        digest.update(b"\0")
        digest.update(path.read_bytes())
        digest.update(b"\0")
    return f"{MATH_VERIFIER_SCHEMA}+{digest.hexdigest()[:12]}"


def require_unchanged_math_verifier(project_root: Path, expected: str) -> None:
    """Stop before persistence when proof code or its environment changed."""

    observed = math_verifier_version(project_root)
    if observed != expected:
        raise RuntimeError(
            "mathematics verifier sources changed during execution; "
            "no further results were persisted"
        )


__all__ = [
    "MATH_VERIFIER_SCHEMA",
    "math_verifier_version",
    "require_unchanged_math_verifier",
]
