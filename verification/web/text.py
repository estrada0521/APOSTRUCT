"""Parse the structure sections of visible complete-mode-detail text for validation.

The parent, subgroup and undistorted-superstructure sections are read here; the
mode definition tables are read by the native definition comparison
(verification/compare/native/definitions.c).
"""

from __future__ import annotations

from dataclasses import dataclass
from html import unescape
import re


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
class CompleteModeDetails:
    parent: SpaceGroupHeader | None
    parent_lattice: Lattice | None
    parent_atoms: tuple[StructureAtom, ...]
    subgroup: SubgroupDetails | None
    undistorted_lattice: Lattice | None
    undistorted_atoms: tuple[StructureAtom, ...]


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
        "Occupancy mode definitions",
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


def parse_complete_mode_details(html: str) -> CompleteModeDetails:
    return parse_complete_mode_details_from_text(_plain_text(html))


def parse_complete_mode_details_from_text(text: str) -> CompleteModeDetails:
    """可視本文テキスト(observation の正本)から直接パースする。"""
    parent_section = _section(text, "Parent structure", "Subgroup details")
    undistorted_section = _section(text, "Undistorted superstructure", "Distorted superstructure")
    return CompleteModeDetails(
        parent=parse_parent_header(text),
        parent_lattice=parse_distorted_lattice(parent_section),
        parent_atoms=parse_structure_atoms(text, "Parent structure", "Subgroup details"),
        subgroup=parse_subgroup_details(text),
        undistorted_lattice=parse_distorted_lattice(undistorted_section),
        undistorted_atoms=parse_structure_atoms(text, "Undistorted superstructure", "Distorted superstructure"),
    )
