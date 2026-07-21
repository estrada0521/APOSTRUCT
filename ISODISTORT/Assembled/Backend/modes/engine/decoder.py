"""Source-backed decoder facade for the mode projection runtime."""

from __future__ import annotations

from ISODISTORT.Assembled.Backend.modes.engine.symmetry.irreps import IrrepMixin
from ISODISTORT.Assembled.Backend.modes.engine.subgroup_structure import IsotropyMixin
from ISODISTORT.Assembled.Backend.modes.engine.project import ProjectMixin
from ISODISTORT.Assembled.Backend.modes.engine.source_catalog import ModeSourceCatalogMixin
from ISODISTORT.Assembled.Backend.modes.engine.symmetry.wyckoff_geometry import WyckoffGeometryMixin


class ModeDataDecoder(
    ModeSourceCatalogMixin,
    IrrepMixin,
    IsotropyMixin,
    WyckoffGeometryMixin,
    ProjectMixin,
):
    """Read Source tables once and evaluate the mode projection contracts."""

    # This is the small integer table addressed as DAT_0050da0c[irrep_type]
    # in the decompiled project_ routine.  SG205 R2-R2- fixes type 2 as 4:
    # frequency 2, image_type 2 doubles it to 4, and the mode kernel returns count 1.
    TYPE_DIVISOR = {1: 1, 2: 4, 3: 2}
