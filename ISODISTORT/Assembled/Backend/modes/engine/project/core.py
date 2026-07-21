from __future__ import annotations

from ISODISTORT.Assembled.Backend.modes.engine.project.rowops import RowOpsMixin
from ISODISTORT.Assembled.Backend.modes.engine.project.selection import ProjectSelectionMixin
from ISODISTORT.Assembled.Backend.modes.engine.project.setup import ProjectSetupMixin
from ISODISTORT.Assembled.Backend.modes.engine.project.transform_irrep2complex import TransformIrrep2ComplexMixin
from ISODISTORT.Assembled.Backend.modes.engine.project.vector import ProjectVectorMixin


class ProjectMixin(
    ProjectSetupMixin,
    ProjectSelectionMixin,
    TransformIrrep2ComplexMixin,
    ProjectVectorMixin,
    RowOpsMixin,
):
    """Faithful port surface for binary ``project_`` and adjacent helpers."""
