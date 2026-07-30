"""Mode-definition presentation helpers."""

from __future__ import annotations

import math
import re
from typing import Any
import gemmi
from ISODISTORT.Assembled.Backend.reciprocal import catalog as reciprocal_catalog
from ISODISTORT.Assembled.Backend.modes.engine.decoder import ModeDataDecoder
from ISODISTORT.Assembled.Backend.modes.engine.tensor import (
    strain_normfactor,
    strain_tensor_components,
    selected_rank2_tensor_row_groups,
    totally_symmetric_rank2_rows,
)
from ISODISTORT.Assembled.Backend.modes.presentation import (
    centering_translations,
    mode_normfactor,
)
from ISODISTORT.Assembled.Backend.modes.request_context import (
    _direction_matrix_text,
)
from ISODISTORT.Assembled.Backend.modes.common import (
    _k_label_from_irrep_label,
    _mode_decoder,
)


def _mode_label(
    parent_symbol: str,
    kvector: str,
    irrep_label: str,
    direction: str,
    site: dict[str, Any],
    site_irrep: str | None,
    mode_index: int,
    mode_kind: str = "dsp",
) -> str:
    atom_label = str(site.get("label") or site.get("type") or "Atom")
    wyckoff = str(site.get("wyckoff") or "?")
    suffix = "" if mode_index == 1 else f"_{mode_index}"
    tail = site_irrep or f"mode{suffix}"
    return f"{parent_symbol}[{kvector}]{irrep_label}{direction}[{atom_label}:{wyckoff}:{mode_kind}]{tail}"


def _strain_mode_definitions(
    sg: int,
    parent_symbol: str,
    decoder: ModeDataDecoder | None = None,
    mode_specs: list[dict[str, Any]] | None = None,
) -> list[dict[str, Any]]:
    """Build parent-cell strain modes fixed by the selected embedding."""

    decoder = decoder or _mode_decoder(None)
    selected_specs = [
        spec
        for spec in (mode_specs or [])
        if str(
            spec.get("k_label")
            or _k_label_from_irrep_label(str(spec.get("label") or ""))
        )
        == "GM"
        and spec.get("gid") is not None
        and spec.get("direction_matrix")
    ]
    if selected_specs:
        rows_by_spec = [
            (
                str(spec.get("label") or ""),
                _direction_matrix_text(spec["direction_matrix"]),
                selected_rank2_tensor_row_groups(
                    decoder,
                    int(sg),
                    int(spec["gid"]),
                    spec["direction_matrix"],
                ),
            )
            for spec in selected_specs
        ]
    else:
        rows_by_spec = []
    catalog = reciprocal_catalog.kpoints(int(sg))
    gamma = next(
        (
            item
            for item in catalog.get("kpoints") or []
            if str(item.get("label")) == "GM"
        ),
        None,
    )
    if gamma is None and not rows_by_spec:
        return []
    if not rows_by_spec:
        irrep = next(
            (
                item
                for item in gamma.get("irreps") or []
                if str(item.get("symbol")) in {"GM1+", "GM1"}
            ),
            None,
        )
        if irrep is None:
            return []
        rows_by_spec = [
            (
                str(irrep["symbol"]),
                "(a)",
                tuple(
                    (row,)
                    for row in totally_symmetric_rank2_rows(
                        decoder, int(sg), int(irrep["gid"])
                    )
                ),
            )
        ]
    definitions: list[dict[str, Any]] = []
    # faithful order is xx,xy,xz,yy,yz,zz; Web columns are
    # e1=xx,e2=yy,e3=zz,e4=yz,e5=xz,e6=xy.
    web_order = (0, 3, 5, 4, 2, 1)
    for symbol, direction_text, row_groups in rows_by_spec:
        presented_groups: list[list[tuple[float, ...]]] = []
        presented_rows: list[tuple[float, ...]] = []
        for source_group in row_groups:
            presented_group: list[tuple[float, ...]] = []
            for source in source_group:
                scale = max((abs(float(value)) for value in source), default=0.0)
                if scale <= 1e-15:
                    continue
                normalized = tuple(float(value) / scale for value in source)
                normalized = tuple(
                    float(round(value)) if abs(value - round(value)) <= 1e-12 else value
                    for value in normalized
                )
                pivot = next(
                    (
                        value
                        for value in normalized
                        if abs(abs(float(value)) - 1.0) <= 1e-10
                    ),
                    1.0,
                )
                if float(pivot) < 0.0:
                    normalized = tuple(-float(value) for value in normalized)
                if any(
                    all(
                        abs(float(left) - float(right)) <= 1e-10
                        for left, right in zip(normalized, seen, strict=True)
                    )
                    for seen in presented_rows
                ):
                    continue
                presented_rows.append(normalized)
                presented_group.append(normalized)
            if presented_group:
                presented_groups.append(presented_group)
        direction_parameters: list[str] = []
        for parameter in re.findall(
            r"(?<![A-Za-z])([a-z])(?![A-Za-z])", direction_text
        ):
            if parameter not in direction_parameters:
                direction_parameters.append(parameter)
        for group_index, presented_group in enumerate(presented_groups, start=1):
            suffix = f"_{group_index}" if len(presented_groups) > 1 else ""
            parameterized = len(direction_parameters) == len(presented_group)
            for row_index, normalized in enumerate(presented_group):
                parameter = direction_parameters[row_index] if parameterized else "a"
                definitions.append(
                    {
                        "label": f"{parent_symbol}[0,0,0]{symbol}{direction_text}strain{suffix}({parameter})",
                        "normfactor": strain_normfactor(normalized),
                        "components": tuple(
                            normalized[position] for position in web_order
                        ),
                        "tensor_components": strain_tensor_components(normalized),
                    }
                )
    return definitions


def _complete_mode_normfactor(
    rows: list[dict[str, Any]],
    child_cartesian: tuple[tuple[float, float, float], ...],
    child_sg: int | None,
) -> float | None:
    """Convert the conventional-row norm to Web's primitive-cell norm."""

    value = mode_normfactor(rows, child_cartesian)
    if value is None or child_sg is None:
        return value
    try:
        symbol = gemmi.find_spacegroup_by_number(int(child_sg)).hm
    except Exception:
        return value
    return value * math.sqrt(len(centering_translations(symbol)))


def _definition_site_irrep_key(label: str) -> tuple[str, str] | None:
    tail = str(label).rsplit("]", 1)[-1]
    match = re.match(r"^(?P<base>.*?)(?:_\d+)?\((?P<component>[^()]*)\)$", tail)
    if match is None:
        return None
    return match.group("base"), match.group("component")


def _orthogonalize_definition_modes(
    modes: list[dict[str, Any]],
    child_cartesian: list[list[float]] | tuple[tuple[float, ...], ...],
    child_sg: int | None,
) -> list[dict[str, Any]]:
    """Apply complete-mode Cartesian Gram-Schmidt within print families."""

    prior_by_key: dict[
        tuple[int, int, str, str],
        list[
            tuple[
                dict[str, Any],
                dict[int, tuple[float, float, float]],
                float,
            ]
        ],
    ] = {}
    out: list[dict[str, Any]] = []
    for source in modes:
        mode = {
            **source,
            "rows": [
                {**row, "dxyz": list(row["dxyz"])} for row in source.get("rows") or []
            ],
        }
        parsed = _definition_site_irrep_key(str(mode.get("label") or ""))
        if parsed is None:
            out.append(mode)
            continue
        key = (
            int(mode.get("_spec_order") or 0),
            int(mode.get("_site_order") or 0),
            parsed[0],
            parsed[1],
        )
        prior = prior_by_key.setdefault(key, [])
        mode_atom_ids = tuple(row.get("atom_id") for row in mode["rows"])
        if (
            not mode_atom_ids
            or not all(
                isinstance(atom_id, str) and atom_id for atom_id in mode_atom_ids
            )
            or len(set(mode_atom_ids)) != len(mode_atom_ids)
        ):
            raise ValueError("mode definition has no complete atom identity")
        for (
            orthogonal,
            reference_cartesian_by_row_id,
            exact_order_denominator,
        ) in prior:
            reference_rows = list(orthogonal.get("rows") or [])
            if mode_atom_ids != tuple(
                reference.get("atom_id") for reference in reference_rows
            ):
                raise ValueError("mode family atom identity/order changed")
            numerator = 0.0
            for row, reference in zip(mode["rows"], reference_rows, strict=True):
                reference_cart = reference_cartesian_by_row_id[id(reference)]
                current_cart = [
                    sum(
                        float(row["dxyz"][axis]) * float(child_cartesian[axis][column])
                        for axis in range(3)
                    )
                    for column in range(3)
                ]
                numerator += sum(
                    current_cart[axis] * reference_cart[axis] for axis in range(3)
                )
            if exact_order_denominator <= 1e-15:
                continue
            coefficient = numerator / exact_order_denominator
            for row, reference in zip(mode["rows"], reference_rows, strict=True):
                row["dxyz"] = [
                    float(row["dxyz"][axis])
                    - coefficient * float(reference["dxyz"][axis])
                    for axis in range(3)
                ]
        scale = max(
            (abs(float(value)) for row in mode["rows"] for value in row["dxyz"]),
            default=0.0,
        )
        if scale > 1e-15:
            for row in mode["rows"]:
                row["dxyz"] = [float(value) / scale for value in row["dxyz"]]
        mode["normfactor"] = _complete_mode_normfactor(
            mode["rows"], child_cartesian, child_sg
        )
        mode_cartesian_by_row_id = {
            id(row): tuple(
                sum(
                    float(row["dxyz"][axis]) * float(child_cartesian[axis][column])
                    for axis in range(3)
                )
                for column in range(3)
            )
            for row in mode["rows"]
        }
        mode_denominator = 0.0
        for row in mode["rows"]:
            vector = mode_cartesian_by_row_id[id(row)]
            mode_denominator += sum(value * value for value in vector)
        prior.append(
            (
                mode,
                mode_cartesian_by_row_id,
                mode_denominator,
            )
        )
        out.append(mode)
    return out


def _apply_dynamic_source_family_presentation(
    modes: list[dict[str, Any]],
) -> list[dict[str, Any]]:
    """Present type-1 real-carrier columns in component-major family order."""

    out = list(modes)
    grouped: dict[tuple[int, int, str, str], list[int]] = {}
    base_families: dict[tuple[int, int, str], set[int]] = {}
    for mode in modes:
        parsed = _definition_site_irrep_key(str(mode.get("label") or ""))
        if parsed is None or mode.get("_source_family") is None:
            continue
        base_key = (
            int(mode.get("_spec_order") or 0),
            int(mode.get("_site_order") or 0),
            parsed[0],
        )
        base_families.setdefault(base_key, set()).add(
            int(mode.get("_source_family") or 0)
        )
    for index, mode in enumerate(modes):
        if (
            mode.get("_source_family") is None
            or mode.get("_source_family_component") is None
        ):
            continue
        static_pair_candidate = (
            mode.get("_mode_kind") == "dsp"
            and not mode.get("_source_family_dynamic")
            and int(mode.get("_source_family_width") or 0) > 1
        )
        static_magnetic_pair_candidate = (
            mode.get("_mode_kind") == "mag"
            and not mode.get("_source_family_dynamic")
            and int(mode.get("_source_family_width") or 0) > 1
        )
        if (
            not mode.get("_source_family_phase")
            and not (
                mode.get("_mode_kind") == "dsp" and mode.get("_source_family_dynamic")
            )
            and not static_pair_candidate
            and not static_magnetic_pair_candidate
        ):
            continue
        parsed = _definition_site_irrep_key(str(mode.get("label") or ""))
        if parsed is None:
            continue
        base_key = (
            int(mode.get("_spec_order") or 0),
            int(mode.get("_site_order") or 0),
            parsed[0],
        )
        identity = mode.get("_source_print_identity") or {}
        multi_group = len(identity.get("opd_groups") or []) > 1
        little_type = int(identity.get("little_type") or 0)
        key = (
            *base_key,
            (
                ""
                if multi_group
                and (little_type == 3 or len(base_families.get(base_key) or ()) == 1)
                else parsed[1]
            ),
        )
        grouped.setdefault(key, []).append(index)

    for indices in grouped.values():
        sources = [modes[index] for index in indices]
        family_components = {
            (
                int(mode.get("_source_family") or 0),
                int(mode.get("_source_family_component") or 0),
            )
            for mode in sources
        }
        static_widths = {int(mode.get("_source_family_width") or 0) for mode in sources}
        source_families = {family for family, _component in family_components}
        static_rectangular_carrier = (
            all(mode.get("_mode_kind") == "dsp" for mode in sources)
            and not any(mode.get("_source_family_dynamic") for mode in sources)
            and len(static_widths) == 1
            and next(iter(static_widths), 0) > 1
            and source_families == set(range(len(source_families)))
            and family_components
            == {
                (family, component)
                for family in source_families
                for component in range(next(iter(static_widths), 0))
            }
        )
        static_rectangular_magnetic = (
            all(mode.get("_mode_kind") == "mag" for mode in sources)
            and not any(mode.get("_source_family_dynamic") for mode in sources)
            and len(static_widths) == 1
            and next(iter(static_widths), 0) > 1
            and source_families == set(range(len(source_families)))
            and family_components
            == {
                (family, component)
                for family in source_families
                for component in range(next(iter(static_widths), 0))
            }
        )
        static_width = next(iter(static_widths), 0)
        paired_static_carrier = (
            static_rectangular_carrier and static_width >= 4 and static_width % 2 == 0
        )
        single_family_multi_group = (
            all(mode.get("_mode_kind") == "dsp" for mode in sources)
            and not any(mode.get("_source_family_dynamic") for mode in sources)
            and all(
                int((mode.get("_source_print_identity") or {}).get("little_type") or 0)
                != 3
                for mode in sources
            )
            and len(source_families) == 1
            and len(static_widths) == 1
            and family_components
            == {
                (next(iter(source_families)), component)
                for component in range(static_width)
            }
            and all(
                len((mode.get("_source_print_identity") or {}).get("opd_groups") or [])
                > 1
                for mode in sources
            )
        )
        if paired_static_carrier:
            opd_row_counts = {
                sum(
                    len(group.get("rows") or [])
                    for group in (mode.get("_source_print_identity") or {}).get(
                        "opd_groups"
                    )
                    or []
                    if isinstance(group, dict)
                )
                for mode in sources
            }
            opd_row_count = (
                next(iter(opd_row_counts), 0) if len(opd_row_counts) == 1 else 0
            )
            source_little_types = {
                int((mode.get("_source_print_identity") or {}).get("little_type") or 0)
                for mode in sources
            }
            source_group_counts = {
                len((mode.get("_source_print_identity") or {}).get("opd_groups") or [])
                for mode in sources
            }
            if (
                source_little_types == {3}
                and opd_row_count > 0
                and static_width % opd_row_count == 0
            ):
                multi_group_type3 = any(count > 1 for count in source_group_counts)
                ordered = sorted(
                    sources,
                    key=lambda mode: (
                        int(
                            (
                                (mode.get("_source_print_identity") or {}).get(
                                    "component"
                                )
                                if multi_group_type3
                                else mode.get("_source_family_component")
                            )
                            or 0
                        )
                        // opd_row_count,
                        int(
                            (
                                (mode.get("_source_print_identity") or {}).get("family")
                                if multi_group_type3
                                else mode.get("_source_family")
                            )
                            or 0
                        ),
                        int(
                            (
                                (mode.get("_source_print_identity") or {}).get(
                                    "component"
                                )
                                if multi_group_type3
                                else mode.get("_source_family_component")
                            )
                            or 0
                        )
                        % opd_row_count,
                    ),
                )
            else:
                pair_width = static_width // 2
                ordered = sorted(
                    sources,
                    key=lambda mode: (
                        int(mode.get("_source_family_component") or 0) % pair_width,
                        int(mode.get("_source_family") or 0),
                        int(mode.get("_source_family_component") or 0) // pair_width,
                    ),
                )
        elif single_family_multi_group:
            ordered = sources
        elif static_rectangular_magnetic:
            ordered = sorted(
                sources,
                key=lambda mode: (
                    mode["_source_print_identity"]["print_component"],
                    int(mode.get("_source_family") or 0),
                ),
            )
        elif (
            static_rectangular_carrier
            or any(
                mode.get("_source_family_phase")
                or (
                    mode.get("_mode_kind") == "dsp"
                    and mode.get("_source_family_dynamic")
                )
                for mode in sources
            )
        ):
            ordered = sorted(
                sources,
                key=lambda mode: (
                    int(mode.get("_source_family_component") or 0),
                    int(mode.get("_source_family") or 0),
                ),
            )
        else:
            ordered = sources
        for target_index, source in zip(indices, ordered, strict=True):
            rows = [
                {
                    **row,
                    "dxyz": [float(value) for value in row.get("dxyz") or []],
                }
                for row in source.get("rows") or []
            ]
            if source.get("_source_family_positive_lead"):
                source_print_scalar = source.get("_source_print_scalar")
                if source_print_scalar is not None:
                    scalar = float(source_print_scalar)
                    rows = [
                        {
                            **row,
                            "dxyz": [
                                scalar * float(value) for value in row.get("dxyz") or []
                            ],
                        }
                        for row in rows
                    ]
            target = modes[target_index]
            label = str(target.get("label") or "")
            out[target_index] = {
                **target,
                "label": label,
                "rows": rows,
                "normfactor": source.get("normfactor"),
                "_source_family": source.get("_source_family"),
                "_source_family_component": source.get("_source_family_component"),
                "_source_family_width": source.get("_source_family_width"),
                "_source_print_scalar": source.get("_source_print_scalar"),
            }
    for index, mode in enumerate(out):
        display_site_irrep = mode.get("_display_site_irrep_label")
        if display_site_irrep:
            prefix, _tail = str(mode.get("label") or "").rsplit("]", 1)
            out[index] = {**mode, "label": f"{prefix}]{display_site_irrep}"}
    return out
