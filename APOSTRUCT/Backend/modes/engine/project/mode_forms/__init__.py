"""Mode-vector shaping helpers for the mode-kernel pipeline."""

from APOSTRUCT.Backend.modes.engine.project.mode_forms.geometry import (
    _apply_fractional_vector_transform,
)
from APOSTRUCT.Backend.modes.engine.project.mode_forms.normalization import (
    _active_mode_indices,
    _add_mode_vectors,
    _dominant_mode_component,
    _normalize_mode_vectors,
    _same_mode,
)
from APOSTRUCT.Backend.modes.engine.project.mode_forms.print_order import (
    _regroup_repeated_component_modes,
    _type1_parametric_scalar_plane_print_order,
    _type3_parametric_component_print_order,
    _type3_parametric_component_print_order_by_block,
    _type3_parametric_kdim2_print_basis,
    _type3_real_print_modes,
)

__all__ = [
    "_active_mode_indices",
    "_add_mode_vectors",
    "_apply_fractional_vector_transform",
    "_dominant_mode_component",
    "_normalize_mode_vectors",
    "_regroup_repeated_component_modes",
    "_same_mode",
    "_type1_parametric_scalar_plane_print_order",
    "_type3_parametric_component_print_order",
    "_type3_parametric_component_print_order_by_block",
    "_type3_parametric_kdim2_print_basis",
    "_type3_real_print_modes",
]
