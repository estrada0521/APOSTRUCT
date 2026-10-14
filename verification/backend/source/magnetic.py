"""Reader for the magnetic ISOTROPY Source tables used by APOSTRUCT.

The decoder exposes ``data_magnetic`` CSR records in their stored coordinate
conventions.  Subgroup basis and origin presentation remain downstream
responsibilities, so this module does not normalize values against rendered
output.
"""

from __future__ import annotations

from pathlib import Path

from verification.backend.source.iso_data import ISOData


ROOT = Path(__file__).resolve().parents[3]
SOURCE = ROOT / "source"


class MagneticData:
    """Minimal decoder for ``Source/data_magnetic.txt``.

    The table stores magnetic isotropy rows as CSR arrays keyed by nonmagnetic
    ``data_irreps`` old-id.  Row ids are 1-based positions in
    ``mag_iso_irrep`` and match the leading id printed by magnetic OPD rows
    for magnetic irreps.
    """

    def __init__(self, source: Path = SOURCE):
        self.source = Path(source)
        if self.source.resolve() == SOURCE.resolve():
            from verification.backend.source.tables import source_tables

            self.table = source_tables().iso.magnetic
        else:
            self.table = ISOData(self.source).magnetic

_CACHE: MagneticData | None = None


def data() -> MagneticData:
    global _CACHE
    if _CACHE is None:
        _CACHE = MagneticData()
    return _CACHE
