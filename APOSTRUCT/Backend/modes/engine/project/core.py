from __future__ import annotations

from APOSTRUCT.Backend.modes.engine.project.rowops import RowOpsMixin
from APOSTRUCT.Backend.modes.engine.project.selection import ProjectSelectionMixin
from APOSTRUCT.Backend.modes.engine.project.setup import ProjectSetupMixin
from APOSTRUCT.Backend.modes.engine.project.transform_irrep2complex import TransformIrrep2ComplexMixin
from APOSTRUCT.Backend.modes.engine.project.vector import ProjectVectorMixin


class ProjectMixin(
    ProjectSetupMixin,
    ProjectSelectionMixin,
    TransformIrrep2ComplexMixin,
    ProjectVectorMixin,
    RowOpsMixin,
):
    """Combined site-irrep projection surface used by the mode engine."""
