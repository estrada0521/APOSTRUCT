"""Faithful mixin surface for subgroup-structure routines."""

from __future__ import annotations

from ISODISTORT.Assembled.Backend.modes.engine.subgroup_structure.atom_expansion import IsotropyAtomExpansionMixin
from ISODISTORT.Assembled.Backend.modes.engine.subgroup_structure.k_vectors import IsotropyKVectorMixin
from ISODISTORT.Assembled.Backend.modes.engine.subgroup_structure.newlat import NewlatMixin
from ISODISTORT.Assembled.Backend.modes.engine.subgroup_structure.setting_change import IsotropySettingChangeMixin


class IsotropyMixin(
    IsotropyKVectorMixin,
    IsotropyAtomExpansionMixin,
    IsotropySettingChangeMixin,
    NewlatMixin,
):
    """Runtime surface for ``find_isotropy_subgroup_``, ``newlat_``, and final atom setup."""
