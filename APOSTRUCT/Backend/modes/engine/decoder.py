"""Source-backed decoder facade for the mode projection runtime."""

from __future__ import annotations

from APOSTRUCT.Backend.modes.engine.symmetry.irreps import IrrepMixin
from APOSTRUCT.Backend.modes.engine.subgroup_structure import IsotropyMixin
from APOSTRUCT.Backend.modes.engine.project import ProjectMixin
from APOSTRUCT.Backend.modes.engine.source_catalog import ModeSourceCatalogMixin
from APOSTRUCT.Backend.modes.engine.symmetry.wyckoff_geometry import WyckoffGeometryMixin


class ModeDataDecoder(
    ModeSourceCatalogMixin,
    IrrepMixin,
    IsotropyMixin,
    WyckoffGeometryMixin,
    ProjectMixin,
):
    """Evaluate one request against shared Source tables.

    Decoder instances are request-scoped. Caches keyed by ``decoder`` therefore
    only reuse work within that request; persistent Source storage must not
    extend the decoder lifetime.
    """

    # The project divisor follows the irrep real-form contract. SG205 R2-R2-
    # fixes type 2 as 4: frequency 2, image_type 2 doubles it to 4, and the
    # mode kernel returns count 1.
    TYPE_DIVISOR = {1: 1, 2: 4, 3: 2}
