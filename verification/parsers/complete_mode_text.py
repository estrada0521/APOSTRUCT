#!/usr/bin/env python3
"""Parse visible complete-mode-detail text for validation.

The complete-modes detail page has two displacive sections with different
roles:

- "Displacive mode definitions" is the mode-basis object: labels,
  normfactors, atom positions, and displacement vectors.
- "Displacive mode amplitudes" is the current amplitude state entered on the
  previous distortion page.  It is useful for UI/state checks, but it is not
  the kernel mode-basis comparison target.
"""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import gzip
from html import unescape
import json
from pathlib import Path
import re
from typing import Iterable


FLOAT_RE = re.compile(r"[-+]?(?:\d+\.\d+|\d+)(?:[eE][-+]?\d+)?")


@dataclass(frozen=True)
class Lattice:
    a: float
    b: float
    c: float
    alpha: float
    beta: float
    gamma: float


@dataclass(frozen=True)
class SpaceGroupHeader:
    number: int | None
    symbol: str


@dataclass(frozen=True)
class SubgroupDetails:
    number: int | None
    symbol: str
    basis: str
    origin: str
    s: int | None
    i: int | None
    k_active: str
    display_label: str | None = None


@dataclass(frozen=True)
class StructureAtom:
    label: str
    site: str
    xyz: tuple[float, float, float]
    occ: float | None
    displ: float | None
    # magnetic ページの追加列(mx my mz mag)。非磁気ページでは None
    mxyz: tuple[float, float, float] | None = None
    mag: float | None = None


@dataclass(frozen=True)
class DistortedAtom:
    label: str
    site: str
    xyz: tuple[float, float, float]
    occ: float
    displ: float
    mxyz: tuple[float, float, float] | None = None
    mag: float | None = None


@dataclass(frozen=True)
class ModeVectorRow:
    atom_label: str | None
    xyz: tuple[float, float, float]
    dxyz: tuple[float, float, float]


@dataclass(frozen=True)
class DisplaciveModeDefinition:
    label: str
    normfactor: float | None
    rows: tuple[ModeVectorRow, ...]


@dataclass(frozen=True)
class StrainModeDefinition:
    label: str
    normfactor: float | None
    components: tuple[float, float, float, float, float, float]


@dataclass(frozen=True)
class DisplaciveAmplitudeRow:
    label: str
    values: tuple[float, ...]


@dataclass(frozen=True)
class CompleteModeDetails:
    parent: SpaceGroupHeader | None
    parent_lattice: Lattice | None
    parent_atoms: tuple[StructureAtom, ...]
    subgroup: SubgroupDetails | None
    undistorted_lattice: Lattice | None
    undistorted_atoms: tuple[StructureAtom, ...]
    lattice: Lattice | None
    distorted_atoms: tuple[DistortedAtom, ...]
    displacive_definitions: tuple[DisplaciveModeDefinition, ...]
    displacive_amplitudes: tuple[DisplaciveAmplitudeRow, ...]
    magnetic_definitions: tuple[DisplaciveModeDefinition, ...] = ()
    strain_definitions: tuple[StrainModeDefinition, ...] = ()
    magnetic_amplitudes: tuple[DisplaciveAmplitudeRow, ...] = ()


def _plain_text(html: str) -> str:
    text = re.sub(r"(?is)<br\s*/?>", "\n", html)
    text = re.sub(r"(?is)</(?:p|div|pre|body|html|h[1-6]|tr)>", "\n", text)
    text = re.sub(r"(?is)<(?:p|div|pre|body|html|h[1-6]|tr|table|tbody)[^>]*>", "\n", text)
    text = re.sub(r"(?is)</?(?:b|i|u|span|td|th|font)[^>]*>", "", text)
    text = re.sub(r"(?is)<script[^>]*>.*?</script>", "\n", text)
    text = re.sub(r"(?is)<style[^>]*>.*?</style>", "\n", text)
    text = re.sub(r"(?is)<[^>]+>", " ", text)
    return unescape(text).replace("\r\n", "\n").replace("\r", "\n")


def _section(text: str, start: str, end: str | None = None) -> str:
    start_index = text.find(start)
    if start_index < 0:
        return ""
    start_index += len(start)
    if end is None:
        return text[start_index:]
    end_index = text.find(end, start_index)
    if end_index < 0:
        return text[start_index:]
    return text[start_index:end_index]


def _numbers(line: str) -> list[float]:
    return [float(item) for item in FLOAT_RE.findall(line)]


def parse_distorted_lattice(text: str) -> Lattice | None:
    match = re.search(
        r"a=\s*(?P<a>[-+\d.]+),\s*b=\s*(?P<b>[-+\d.]+),\s*c=\s*(?P<c>[-+\d.]+),?\s*"
        r"alpha=\s*(?P<alpha>[-+\d.]+),\s*beta=\s*(?P<beta>[-+\d.]+),\s*gamma=\s*(?P<gamma>[-+\d.]+)",
        text,
    )
    if not match:
        return None
    values = {key: float(value) for key, value in match.groupdict().items()}
    return Lattice(**values)


def parse_parent_header(text: str) -> SpaceGroupHeader | None:
    match = re.search(r"Parent structure\s*\((?P<number>\d+)\s+(?P<symbol>[^)]*)\)", text)
    if match is None:
        match = re.search(
            r"Parent structure\s*\n\s*(?P<number>\d+)\s+(?P<symbol>[^\n]+)",
            text,
        )
    if not match:
        return None
    return SpaceGroupHeader(number=int(match.group("number")), symbol=match.group("symbol").strip())


def parse_subgroup_details(text: str) -> SubgroupDetails | None:
    body = _section(text, "Subgroup details", "Undistorted superstructure")
    for raw in body.splitlines():
        line = raw.strip()
        if not line:
            continue
        match = re.match(r"^(?P<number>\d+(?:\.\d+)?)\s+(?P<symbol>.*?),\s*basis=\{(?P<basis>.*?)\},\s*origin=(?P<origin>\([^)]*\)),\s*s=(?P<s>-?\d+),\s*i=(?P<i>-?\d+),\s*k-active=\s*(?P<k_active>.*)$", line)
        if not match:
            continue
        display_label = match.group("number")
        return SubgroupDetails(
            number=int(display_label.split(".", 1)[0]),
            symbol=match.group("symbol").strip(),
            basis=match.group("basis").strip(),
            origin=match.group("origin").strip(),
            s=int(match.group("s")),
            i=int(match.group("i")),
            k_active=match.group("k_active").strip(),
            display_label=display_label,
        )
    return None


def parse_structure_atoms(text: str, start: str, end: str) -> tuple[StructureAtom, ...]:
    body = _cut_at_headings(
        _section(text, start, end),
        "Displacive mode definitions",
        "Magnetic mode definitions",
        "Parent-cell strain mode definitions",
    )
    rows: list[StructureAtom] = []
    for raw in body.splitlines():
        line = raw.strip()
        if not line or line.startswith("a=") or line.startswith("atom site"):
            continue
        parts = line.split()
        if len(parts) < 5:
            continue
        try:
            xyz = tuple(float(value) for value in parts[2:5])
        except ValueError:
            continue
        occ = None
        displ = None
        if len(parts) >= 6:
            try:
                occ = float(parts[5])
            except ValueError:
                occ = None
        if len(parts) >= 7:
            try:
                displ = float(parts[6])
            except ValueError:
                displ = None
        mxyz, mag = _optional_moments(parts)
        rows.append(StructureAtom(label=parts[0], site=parts[1], xyz=xyz, occ=occ,
                                  displ=displ, mxyz=mxyz, mag=mag))
    return tuple(rows)


def _optional_moments(parts: list[str]) -> tuple[tuple[float, float, float] | None, float | None]:
    mxyz = None
    mag = None
    if len(parts) >= 10:
        try:
            mxyz = tuple(float(value) for value in parts[7:10])
        except ValueError:
            mxyz = None
    if mxyz is not None and len(parts) >= 11:
        try:
            mag = float(parts[10])
        except ValueError:
            mag = None
    return mxyz, mag


def _cut_at_headings(body: str, *stops: str) -> str:
    """後続節見出しの最初の出現で打ち切る(見出し欠落ページの over-capture 防止)。"""
    for stop in stops:
        idx = body.find(stop)
        if idx >= 0:
            body = body[:idx]
    return body


_DEFINITION_HEADINGS = (
    "Displacive mode definitions",
    "Magnetic mode definitions",
    "Parent-cell strain mode definitions",
)


def parse_distorted_atoms(text: str) -> tuple[DistortedAtom, ...]:
    # 全磁気ページには "Displacive mode definitions" が無い。end をその見出しに
    # 固定すると magnetic mode 行が偽 distorted 原子として捕捉されるため、
    # どの definitions 見出しでも打ち切る。
    body = _cut_at_headings(_section(text, "Distorted superstructure"), *_DEFINITION_HEADINGS)
    rows: list[DistortedAtom] = []
    for raw in body.splitlines():
        line = raw.strip()
        if not line or line.startswith("a=") or line.startswith("atom site"):
            continue
        parts = line.split()
        if len(parts) < 7:
            continue
        try:
            values = [float(value) for value in parts[2:7]]
        except ValueError:
            continue
        mxyz, mag = _optional_moments(parts)
        rows.append(
            DistortedAtom(
                label=parts[0],
                site=parts[1],
                xyz=(values[0], values[1], values[2]),
                occ=values[3],
                displ=values[4],
                mxyz=mxyz,
                mag=mag,
            )
        )
    return tuple(rows)


def parse_displacive_definitions(text: str) -> tuple[DisplaciveModeDefinition, ...]:
    body = _cut_at_headings(
        _section(text, "Displacive mode definitions"),
        "Displacive mode amplitudes",
        "Magnetic mode definitions",
        "Magnetic mode amplitudes",
        "Parent-cell strain mode definitions",
        "Parent-cell strain mode amplitudes",
    )
    modes: list[DisplaciveModeDefinition] = []
    current_label: str | None = None
    current_norm: float | None = None
    current_rows: list[ModeVectorRow] = []

    def flush() -> None:
        nonlocal current_label, current_norm, current_rows
        if current_label is None:
            return
        modes.append(
            DisplaciveModeDefinition(
                label=current_label,
                normfactor=current_norm,
                rows=tuple(current_rows),
            )
        )
        current_label = None
        current_norm = None
        current_rows = []

    for raw in body.splitlines():
        line = raw.strip()
        if not line or line.startswith("atom "):
            continue
        if "normfactor" in line:
            flush()
            label, tail = line.split("normfactor", 1)
            current_label = label.strip()
            if (
                current_label.count("(") == current_label.count(")") + 1
                and re.search(r"\([^()]*$", current_label) is not None
            ):
                current_label += ")"
            nums = _numbers(tail)
            current_norm = nums[0] if nums else None
            continue
        nums = _numbers(line)
        if current_label is None or len(nums) < 6:
            continue
        parts = line.split()
        atom_label: str | None = None
        # Continuation rows contain only numeric columns.  A leading atom token
        # is present only for the first row of an atom orbit block.
        if parts and not FLOAT_RE.fullmatch(parts[0]):
            atom_label = parts[0]
        values = nums[-6:]
        current_rows.append(
            ModeVectorRow(
                atom_label=atom_label,
                xyz=(values[0], values[1], values[2]),
                dxyz=(values[3], values[4], values[5]),
            )
        )
    flush()
    return tuple(modes)


def parse_magnetic_definitions(text: str) -> tuple[DisplaciveModeDefinition, ...]:
    body = _cut_at_headings(
        _section(text, "Magnetic mode definitions"),
        "Magnetic mode amplitudes",
        "Parent-cell strain mode definitions",
        "Parent-cell strain mode amplitudes",
    )
    # The magnetic table has the same row grammar as the displacive table;
    # only the vector headings and label kind differ.
    return parse_displacive_definitions(
        "Displacive mode definitions\n" + body + "\nDisplacive mode amplitudes"
    )


def parse_strain_definitions(text: str) -> tuple[StrainModeDefinition, ...]:
    body = _section(text, "Parent-cell strain mode definitions", "Parent-cell strain mode amplitudes")
    modes: list[StrainModeDefinition] = []
    label: str | None = None
    normfactor: float | None = None
    for raw in body.splitlines():
        line = raw.strip()
        if not line or line.startswith("e1"):
            continue
        if "normfactor" in line:
            label_text, tail = line.split("normfactor", 1)
            label = label_text.strip()
            values = _numbers(tail)
            normfactor = values[0] if values else None
            continue
        values = _numbers(line)
        if label is None or len(values) != 6:
            continue
        modes.append(
            StrainModeDefinition(
                label=label,
                normfactor=normfactor,
                components=tuple(values),  # type: ignore[arg-type]
            )
        )
        label = None
        normfactor = None
    return tuple(modes)


def parse_displacive_amplitudes(text: str) -> tuple[DisplaciveAmplitudeRow, ...]:
    return _parse_amplitude_section(text, "Displacive mode amplitudes")


def parse_magnetic_amplitudes(text: str) -> tuple[DisplaciveAmplitudeRow, ...]:
    return _parse_amplitude_section(text, "Magnetic mode amplitudes")


def _parse_amplitude_section(text: str, start: str) -> tuple[DisplaciveAmplitudeRow, ...]:
    # 後続節まで読み進めると、6 数値の mode 行から末尾 2-3 個を吸った
    # 偽 amplitude 行が量産される。最初の後続見出しで打ち切る。
    body = _cut_at_headings(_section(text, start), *_DEFINITION_HEADINGS)
    rows: list[DisplaciveAmplitudeRow] = []
    for raw in body.splitlines():
        line = raw.strip()
        if not line or line.startswith("mode ") or line.startswith("</"):
            continue
        match = re.match(
            r"^(?P<label>.*?)(?P<values>(?:\s+[-+]?(?:\d+\.\d+|\d+)(?:[eE][-+]?\d+)?){2,3})\s*$",
            line,
        )
        if not match:
            continue
        values = tuple(float(value) for value in FLOAT_RE.findall(match.group("values")))
        rows.append(DisplaciveAmplitudeRow(label=match.group("label").rstrip(), values=values))
    return tuple(rows)


def parse_complete_mode_details(html: str) -> CompleteModeDetails:
    return parse_complete_mode_details_from_text(_plain_text(html))


def parse_complete_mode_details_from_text(text: str) -> CompleteModeDetails:
    """可視本文テキスト(observation の正本)から直接パースする。"""
    parent_section = _section(text, "Parent structure", "Subgroup details")
    undistorted_section = _section(text, "Undistorted superstructure", "Distorted superstructure")
    distorted_section = _section(text, "Distorted superstructure", "Displacive mode definitions")
    return CompleteModeDetails(
        parent=parse_parent_header(text),
        parent_lattice=parse_distorted_lattice(parent_section),
        parent_atoms=parse_structure_atoms(text, "Parent structure", "Subgroup details"),
        subgroup=parse_subgroup_details(text),
        undistorted_lattice=parse_distorted_lattice(undistorted_section),
        undistorted_atoms=parse_structure_atoms(text, "Undistorted superstructure", "Distorted superstructure"),
        lattice=parse_distorted_lattice(distorted_section) or parse_distorted_lattice(text),
        distorted_atoms=parse_distorted_atoms(text),
        displacive_definitions=parse_displacive_definitions(text),
        displacive_amplitudes=parse_displacive_amplitudes(text),
        magnetic_definitions=parse_magnetic_definitions(text),
        strain_definitions=parse_strain_definitions(text),
        magnetic_amplitudes=parse_magnetic_amplitudes(text),
    )


def parse_complete_mode_structure(html: str) -> CompleteModeDetails:
    """Parse only Parent/Subgroup/Undistorted structure sections.

    This intentionally stops before the Mode definition section so structure
    validation cannot spend time parsing or reporting vector rows.
    """

    text = _plain_text(html)
    parent_section = _section(text, "Parent structure", "Subgroup details")
    undistorted_section = _section(text, "Undistorted superstructure", "Distorted superstructure")
    return CompleteModeDetails(
        parent=parse_parent_header(text),
        parent_lattice=parse_distorted_lattice(parent_section),
        parent_atoms=parse_structure_atoms(text, "Parent structure", "Subgroup details"),
        subgroup=parse_subgroup_details(text),
        undistorted_lattice=parse_distorted_lattice(undistorted_section),
        undistorted_atoms=parse_structure_atoms(text, "Undistorted superstructure", "Distorted superstructure"),
        lattice=None,
        distorted_atoms=(),
        displacive_definitions=(),
        displacive_amplitudes=(),
    )


def definition_signature(details: CompleteModeDetails, ndigits: int = 5) -> list[dict[str, object]]:
    """Return the strict kernel-comparison payload.

    This intentionally excludes displacive amplitude rows.  The rows are part
    of the page state, not the mode-basis definition.
    """

    out: list[dict[str, object]] = []
    for mode in details.displacive_definitions:
        out.append(
            {
                "label": mode.label,
                "normfactor": None if mode.normfactor is None else round(mode.normfactor, ndigits),
                "rows": [
                    {
                        "atom": row.atom_label,
                        "xyz": [round(value, ndigits) for value in row.xyz],
                        "dxyz": [round(value, ndigits) for value in row.dxyz],
                    }
                    for row in mode.rows
                ],
            }
        )
    return out


def target_signature(details: CompleteModeDetails, ndigits: int = 5) -> dict[str, object]:
    """Return the complete-mode-detail surface currently worth comparing.

    Distorted-superstructure rows and displacive amplitudes are deliberately
    excluded.  They are downstream state, while the current oracle boundary is
    Parent / Subgroup / Undistorted superstructure / Displacive definitions.
    """

    def lattice_payload(lattice: Lattice | None) -> dict[str, float] | None:
        if lattice is None:
            return None
        return {
            "a": round(lattice.a, ndigits),
            "b": round(lattice.b, ndigits),
            "c": round(lattice.c, ndigits),
            "alpha": round(lattice.alpha, ndigits),
            "beta": round(lattice.beta, ndigits),
            "gamma": round(lattice.gamma, ndigits),
        }

    def atoms_payload(rows: tuple[StructureAtom, ...]) -> list[dict[str, object]]:
        return [
            {
                "label": row.label,
                "site": row.site,
                "xyz": [round(value, ndigits) for value in row.xyz],
                "occ": None if row.occ is None else round(row.occ, ndigits),
            }
            for row in rows
        ]

    return {
        "parent": None if details.parent is None else asdict(details.parent),
        "parent_lattice": lattice_payload(details.parent_lattice),
        "parent_atoms": atoms_payload(details.parent_atoms),
        "subgroup": None if details.subgroup is None else asdict(details.subgroup),
        "undistorted_lattice": lattice_payload(details.undistorted_lattice),
        "undistorted_atoms": atoms_payload(details.undistorted_atoms),
        "displacive_definitions": definition_signature(details, ndigits=ndigits),
    }


def _json_default(obj: object) -> object:
    if hasattr(obj, "__dataclass_fields__"):
        return asdict(obj)
    if isinstance(obj, tuple):
        return list(obj)
    raise TypeError(f"object of type {type(obj).__name__} is not JSON serializable")


def summarize(details: CompleteModeDetails) -> dict[str, object]:
    """Return the current comparison target surface.

    For complete-mode-detail alignment we intentionally ignore the page's
    distorted-superstructure and amplitude-state sections.  The useful static
    oracle surface is Parent / Subgroup / Undistorted superstructure /
    Displacive mode definitions.
    """
    return {
        "parent": details.parent,
        "parent_atom_count": len(details.parent_atoms),
        "subgroup": details.subgroup,
        "undistorted_atom_count": len(details.undistorted_atoms),
        "displacive_definition_count": len(details.displacive_definitions),
        "magnetic_definition_count": len(details.magnetic_definitions),
        "strain_definition_count": len(details.strain_definitions),
        "definition_labels": [mode.label for mode in details.displacive_definitions],
        "magnetic_definition_labels": [mode.label for mode in details.magnetic_definitions],
        "strain_definition_labels": [mode.label for mode in details.strain_definitions],
    }


def iter_paths(paths: Iterable[Path]) -> Iterable[Path]:
    for path in paths:
        if path.is_dir():
            yield from sorted([*path.rglob("*.html"), *path.rglob("*.txt.gz")])
        else:
            yield path


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("html", nargs="+", type=Path)
    parser.add_argument("--signature", action="store_true", help="emit kernel definition signature instead of a summary")
    parser.add_argument(
        "--target-signature",
        action="store_true",
        help="emit Parent/Subgroup/Undistorted/definition signature; ignore distorted/amplitudes",
    )
    args = parser.parse_args()

    payload: dict[str, object] = {}
    for path in iter_paths(args.html):
        if path.name.endswith(".txt.gz"):
            with gzip.open(path, "rt", encoding="utf-8", errors="strict") as handle:
                html = handle.read()
        else:
            html = path.read_text(encoding="utf-8", errors="replace")
        details = parse_complete_mode_details(html)
        if not details.displacive_definitions and not details.displacive_amplitudes:
            continue
        if args.target_signature:
            payload[str(path)] = target_signature(details)
        elif args.signature:
            payload[str(path)] = definition_signature(details)
        else:
            payload[str(path)] = summarize(details)
    print(json.dumps(payload, indent=2, ensure_ascii=False, default=_json_default))


if __name__ == "__main__":
    main()
