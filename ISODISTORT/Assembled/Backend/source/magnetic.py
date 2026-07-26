"""Reader for the magnetic ISOTROPY Source tables used by ISODISTORT.

The decoder exposes ``data_magnetic`` CSR records in their stored coordinate
conventions.  Subgroup basis and origin presentation remain downstream
responsibilities, so this module does not normalize values against rendered
output.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import sys
from typing import Any


ROOT = Path(__file__).resolve().parents[4]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from ISODISTORT.Assembled.Backend.source.iso_data import _load_table  # noqa: E402


SOURCE = ROOT / "Source"


@dataclass(frozen=True)
class MagneticIsotropyRow:
    row_id: int
    irrep_old_id: int
    parent_magnetic_group: int
    subgroup_magnetic_group: int
    subgroup_old: int
    subgroup_number_label: str
    subgroup_bns_label: str
    subgroup_unified_label: str
    basis: tuple[int, ...]
    origin: tuple[int, int, int, int]
    orderparam_dim: int
    orderparam_freeparam: int
    orderparam_label: str
    orderparam: tuple[int, ...]


@dataclass(frozen=True)
class SSGMagRecord:
    row_id: int
    ssg: int
    magnetic_group: int
    nlabel: str
    elements_count: int
    elements_pointer: int


class MagneticData:
    """Minimal decoder for ``Source/data_magnetic.txt``.

    The table stores magnetic isotropy rows as CSR arrays keyed by nonmagnetic
    ``data_irreps`` old-id.  Row ids are 1-based positions in
    ``mag_iso_irrep`` and match the leading id printed by ISODISTORT OPD pages
    for magnetic irreps.
    """

    def __init__(self, source: Path = SOURCE):
        self.source = Path(source)
        self.table = _load_table(self.source / "data_magnetic.txt")

    def irrep_has_magnetic_isotropy(self, old_id: int) -> bool:
        return bool(self.magnetic_isotropy_row_ids(old_id))

    def magnetic_isotropy_row_ids(self, old_id: int) -> range:
        pointers = self.table["mag_iso_irrep_pointer"]
        old_id = int(old_id)
        if old_id < 1 or old_id >= len(pointers):
            return range(0)
        start = int(pointers[old_id - 1])
        end = int(pointers[old_id])
        return range(start, end)

    def magnetic_isotropy_rows_for_irrep(self, old_id: int) -> list[MagneticIsotropyRow]:
        return [self.magnetic_isotropy_row(row_id) for row_id in self.magnetic_isotropy_row_ids(old_id)]

    def magnetic_isotropy_row(self, row_id: int) -> MagneticIsotropyRow:
        row_id = int(row_id)
        index = row_id - 1
        table: dict[str, list[Any]] = self.table
        subgroup = int(table["mag_iso_subgroup"][index])
        ptr = int(table["mag_iso_orderparam_pointer"][index])
        end = (
            int(table["mag_iso_orderparam_pointer"][index + 1])
            if index + 1 < len(table["mag_iso_orderparam_pointer"])
            else len(table["mag_iso_orderparam"]) + 1
        )
        return MagneticIsotropyRow(
            row_id=row_id,
            irrep_old_id=int(table["mag_iso_irrep"][index]),
            parent_magnetic_group=int(table["mag_iso_parent"][index]),
            subgroup_magnetic_group=subgroup,
            subgroup_old=int(table["mag_iso_subgroup_old"][index]),
            subgroup_number_label=str(table["mag_nlabel"][subgroup - 1]).strip(),
            subgroup_bns_label=str(table["mag_bns_label"][subgroup - 1]).strip(),
            subgroup_unified_label=str(table["mag_unified_label"][subgroup - 1]).strip(),
            basis=tuple(int(x) for x in table["mag_iso_basis"][index * 9:(index + 1) * 9]),
            origin=tuple(int(x) for x in table["mag_iso_origin"][index * 4:(index + 1) * 4]),  # type: ignore[arg-type]
            orderparam_dim=int(table["mag_iso_orderparam_dim"][index]),
            orderparam_freeparam=int(table["mag_iso_orderparam_freeparam"][index]),
            orderparam_label=str(table["mag_iso_orderparam_label"][index]).strip(),
            orderparam=tuple(int(x) for x in table["mag_iso_orderparam"][ptr - 1:end - 1]),
        )


_CACHE: MagneticData | None = None


def data() -> MagneticData:
    global _CACHE
    if _CACHE is None:
        _CACHE = MagneticData()
    return _CACHE


class SSGMagData:
    """Lazy decoder for ``Source/data_ssgmag.txt``.

    ``data_ssgmag`` is much larger than ``data_magnetic``.  Keep it out of the
    default magnetic isotropy path and instantiate this class only when the SSG
    mapping/operation records are explicitly needed.
    """

    def __init__(self, source: Path = SOURCE):
        self.source = Path(source)
        self.table = _load_table(self.source / "data_ssgmag.txt")

    def records_for_ssg(self, ssg: int) -> list[SSGMagRecord]:
        table: dict[str, list[Any]] = self.table
        pointers = table["issgmag_ssg_pointer"]
        ssg = int(ssg)
        if ssg < 1 or ssg >= len(pointers):
            return []
        start = int(pointers[ssg - 1])
        end = int(pointers[ssg])
        return [self.record(row_id) for row_id in range(start, end)]

    def record(self, row_id: int) -> SSGMagRecord:
        row_id = int(row_id)
        index = row_id - 1
        table: dict[str, list[Any]] = self.table
        return SSGMagRecord(
            row_id=row_id,
            ssg=int(table["issgmag_ssg"][index]),
            magnetic_group=int(table["issgmag_mag"][index]),
            nlabel=str(table["ssgmag_nlabel"][index]).strip(),
            elements_count=int(table["issgmag_elements_count"][index]),
            elements_pointer=int(table["issgmag_elements_pointer"][index]),
        )
