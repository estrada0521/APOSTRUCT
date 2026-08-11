"""Enumerate order-parameter directions and select an isotropy row."""

from __future__ import annotations

from fractions import Fraction
import re
from typing import Any, Mapping

from APOSTRUCT.Backend.isotropy.domains import (
    _fraction_mod1,
    _mat4_mul,
    _operation_matrix,
    _setting_matrix,
)

from APOSTRUCT.Backend.isotropy import catalog as isotropy_catalog
from APOSTRUCT.Backend.source import magnetic as magnetic_data
from APOSTRUCT.Backend.parent import parent_cell_tuple
from APOSTRUCT.Backend.reciprocal import catalog as reciprocal_catalog
from APOSTRUCT.Backend.reciprocal.k_points import normalize_k_params
from APOSTRUCT.Backend.reciprocal.irreps import select_irrep
from APOSTRUCT.Backend.isotropy.coupled import coupled_opd_rows
from APOSTRUCT.Backend.source.tables import source_tables


_COUPLED_OPD_PART_RE = re.compile(
    r"(?P<label>(?:P|C|S|(?:[4-9]|[1-9][0-9])D)[1-9][0-9]*)"
    r"\((?P<domain>[1-9][0-9]*)\)"
)


class OpdNotFoundError(ValueError):
    """The requested OPD is absent from the computed candidate rows."""


def _coupled_selected_labels(
    orderparam: int | str | None,
    *,
    slot_count: int,
    selected_opd_only: bool,
) -> tuple[str, ...] | None:
    """Parse exact per-slot Source labels from a complete OPD heading.

    The heading is a conditioned input.  It only bounds Source OPD discovery;
    no Web row values are accepted as computed output.
    """

    if not selected_opd_only or not isinstance(orderparam, str) or slot_count < 2:
        return None
    text = orderparam.strip()
    parts = tuple(_COUPLED_OPD_PART_RE.finditer(text))
    if len(parts) != int(slot_count):
        return None
    if "".join(part.group(0) for part in parts) != text:
        return None
    return tuple(part.group("label") for part in parts)


def _rows_for_selected_label(
    rows: list[dict[str, Any]],
    label: str | None,
) -> list[dict[str, Any]]:
    if label is None:
        return rows
    return [
        row
        for row in rows
        if label
        in {
            str((row.get("isotropy") or {}).get("opd_label") or "").strip(),
            str((row.get("direction") or {}).get("label") or "").strip(),
        }
    ]


def _select_orderparam(
    opd_rows: list[dict[str, Any]],
    *,
    orderparam: int | str | None = None,
    orderparam_index: int = 1,
) -> dict[str, Any] | None:
    requested = orderparam is not None
    if isinstance(orderparam, str) and orderparam.strip():
        label = orderparam.strip()
        selected = next(
            (
                row
                for row in opd_rows
                if label
                in {
                    str((row.get("isotropy") or {}).get("opd_label") or "").strip(),
                    str((row.get("direction") or {}).get("label") or "").strip(),
                }
            ),
            None,
        )
        if selected is not None:
            return selected

    orderparam_row_id = None
    if isinstance(orderparam, int) and not isinstance(orderparam, bool):
        orderparam_row_id = int(orderparam)
    elif isinstance(orderparam, str) and re.fullmatch(r"\s*\d+\s*", orderparam):
        orderparam_row_id = int(orderparam)
    if orderparam_row_id is not None:
        selected = next(
            (
                row
                for row in opd_rows
                if int(((row.get("isotropy") or {}).get("row_id") or 0)) == orderparam_row_id
            ),
            None,
        )
        if selected is not None:
            return selected
    if requested:
        available = list(
            dict.fromkeys(
                str((row.get("direction") or {}).get("label") or "").strip()
                for row in opd_rows
                if str((row.get("direction") or {}).get("label") or "").strip()
            )
        )
        raise OpdNotFoundError(
            f"OPD {orderparam!r} not found; available labels: {available}"
        )
    if 1 <= int(orderparam_index) <= len(opd_rows):
        return opd_rows[int(orderparam_index) - 1]
    return None


def _selected_opd_factors(
    parent_sg: int,
    slots: list[dict[str, Any]],
    selected_opd: dict[str, Any] | None,
) -> list[dict[str, Any]]:
    """Return the exact primary OPD/domain selection for downstream tools."""

    if not isinstance(selected_opd, dict):
        return []
    isotropy = selected_opd.get("isotropy")
    if not isinstance(isotropy, dict):
        return []
    if len(slots) == 1:
        labels = [str(isotropy.get("opd_label") or "")]
        domains = [int(isotropy.get("selected_domain") or 1)]
        parameter_counts = [int(isotropy.get("free") or 0)]
        display_opds = [
            str(isotropy.get("display_opd") or isotropy.get("source_opd") or "")
        ]
        domain_displays = [
            _domain_display(
                isotropy,
                label=labels[0],
                domain=domains[0],
                display_opd=display_opds[0],
            )
        ]
    else:
        components = (selected_opd.get("direction") or {}).get("components") or []
        if any(
            not isinstance(component, dict)
            or component.get("slot") != index
            for index, component in enumerate(components, 1)
        ):
            raise ValueError("selected OPD components do not match primary slots")
        labels = [str(component.get("opd") or "") for component in components]
        domains = [component.get("domain") for component in components]
        parameter_counts = []
        for slot, label in zip(slots, labels, strict=True):
            source_row = next(
                (
                    row
                    for row in slot.get("opd_rows") or []
                    if str(((row.get("isotropy") or {}).get("opd_label") or ""))
                    == label
                ),
                None,
            )
            parameter_counts.append(
                None
                if source_row is None
                else (source_row.get("isotropy") or {}).get("free")
            )
        display_opds = [
            str(component.get("display_opd") or "") for component in components
        ]
        domain_displays = [
            str(component.get("domain_display") or "") for component in components
        ]
    if (
        len(labels) != len(slots)
        or len(domains) != len(slots)
        or len(parameter_counts) != len(slots)
        or len(display_opds) != len(slots)
        or len(domain_displays) != len(slots)
        or any(not label for label in labels)
        or any(not display_opd for display_opd in display_opds)
        or any(not domain_display for domain_display in domain_displays)
        or any(type(domain) is not int or domain < 1 for domain in domains)
        or any(type(count) is not int or count < 1 for count in parameter_counts)
    ):
        raise ValueError(
            "selected OPD does not retain one exact displayed domain per primary slot"
        )

    def selected_k_parameters(slot: dict[str, Any]) -> list[str]:
        kpoint = slot.get("kpoint") or {}
        dimension = int(kpoint.get("dimension") or 0)
        if dimension <= 0:
            return []
        names = reciprocal_catalog.k_coordinate_map_for_slot(
            int(parent_sg), int(kpoint["kslot"])
        ).parameter_names
        if len(names) != dimension:
            raise ValueError("selected K slot does not retain its parameters")
        values = slot.get("display_k_params") or slot.get("k_params") or {}
        return [str(values.get(name) or "0").strip() or "0" for name in names]

    def selected_source_k_parameters(
        slot: dict[str, Any],
        label: str,
    ) -> list[str] | None:
        kpoint = slot.get("kpoint") or {}
        dimension = int(kpoint.get("dimension") or 0)
        if dimension <= 0:
            return []
        source = slot.get("source_kparam")
        if source is None and len(slots) == 1:
            source = isotropy.get("source_kparam")
        if source is None:
            source_row = next(
                (
                    row
                    for row in slot.get("opd_rows") or []
                    if str(((row.get("isotropy") or {}).get("opd_label") or ""))
                    == label
                ),
                None,
            )
            source = (source_row or {}).get("isotropy", {}).get("source_kparam")
        if not isinstance(source, (list, tuple)) or len(source) != 4:
            return None
        denominator = int(source[3])
        if denominator == 0:
            return None
        return [
            str(Fraction(int(source[index]), denominator))
            for index in range(dimension)
        ]

    factors: list[dict[str, Any]] = []
    for slot, label, domain, domain_display, parameter_count in zip(
        slots, labels, domains, domain_displays, parameter_counts, strict=True
    ):
        source_parameters = selected_source_k_parameters(slot, label)
        factor = {
            "slot": int(slot["slot"]),
            "gid": int(slot["irrep"]["gid"]),
            "label": str(slot["irrep"].get("symbol") or ""),
            "magnetic": bool(slot["irrep"].get("magnetic")),
            "k_parameters": selected_k_parameters(slot),
            "opd": label,
            "domain": int(domain),
            "domain_display": domain_display,
            "parameter_count": int(parameter_count),
        }
        if len(slots) == 1 and bool(isotropy.get("embedding_selected")):
            source_rows = list(isotropy.get("source_numeric_rows") or [])
            full_dimension = int(slot["irrep"]["full_dim"])
            if len(source_rows) != int(parameter_count) or any(
                len(row) != full_dimension for row in source_rows
            ):
                raise ValueError("selected embedding lost its exact direction matrix")
            factor["direction_matrix"] = [
                [source_rows[column][row] for column in range(len(source_rows))]
                for row in range(full_dimension)
            ]
        if source_parameters is not None:
            factor["source_k_parameters"] = source_parameters
        factors.append(factor)
    return factors


def select_embedding_direction(
    isotropy_state: dict[str, Any],
    directions_payload: Mapping[str, Any],
    direction_row: int,
) -> dict[str, Any]:
    """Replace a catalog OPD representative with one exact directions embedding."""

    if directions_payload.get("schema") != "APOSTRUCT.cli.directions":
        raise ValueError("--from-directions requires a compact directions payload")
    rows = directions_payload.get("directions")
    if not isinstance(rows, list) or not rows:
        raise ValueError("directions payload has no direction rows")
    if type(direction_row) is not int or not 1 <= direction_row <= len(rows):
        raise ValueError(f"direction row must be an integer from 1 through {len(rows)}")
    row = rows[direction_row - 1]
    if not isinstance(row, Mapping):
        raise ValueError(f"directions row {direction_row} is not an object")
    if row.get("role") != "primary":
        raise ValueError(f"directions row {direction_row} is not primary")

    selected = isotropy_state.get("selected") or {}
    slots = list(selected.get("selection_slots") or [])
    factors = list(selected.get("opd_factors") or [])
    selected_opd = selected.get("orderparam")
    if len(slots) != 1 or len(factors) != 1 or not isinstance(selected_opd, dict):
        raise ValueError("directions embedding reuse requires one selected irrep")
    factor = factors[0]
    if int(row.get("gid") or 0) != int(factor.get("gid") or 0):
        raise ValueError("directions row does not match the selected irrep")

    parent = int((isotropy_state.get("input") or {}).get("parent", {}).get("number") or 0)
    payload_parent = directions_payload.get("parent")
    if not isinstance(payload_parent, Mapping) or int(payload_parent.get("number") or 0) != parent:
        raise ValueError("directions parent differs from the mode parent")
    embedding = directions_payload.get("embedding")
    subgroup = directions_payload.get("subgroup")
    if not isinstance(embedding, Mapping) or not isinstance(subgroup, Mapping):
        raise ValueError("directions payload has no exact subgroup embedding")

    direction_matrix = row.get("direction_matrix")
    if not isinstance(direction_matrix, list) or not direction_matrix:
        raise ValueError("directions row has no direction matrix")
    try:
        source_rows = [
            [float(direction_matrix[coordinate][free]) for coordinate in range(len(direction_matrix))]
            for free in range(int(row["parameter_count"]))
        ]
    except (IndexError, KeyError, TypeError, ValueError) as exc:
        raise ValueError("directions row has an invalid direction matrix") from exc

    opd_label = str(row.get("opd") or "").strip()
    selected_domain = int(row.get("domain") or 0)
    if not opd_label or selected_domain < 1:
        raise ValueError("directions row lacks Source OPD/domain; regenerate it")

    basis_matrix = embedding.get("basis_matrix")
    origin_vector = embedding.get("origin_vector")
    if (
        not isinstance(basis_matrix, list)
        or len(basis_matrix) != 3
        or any(not isinstance(values, list) or len(values) != 3 for values in basis_matrix)
        or not isinstance(origin_vector, list)
        or len(origin_vector) != 3
    ):
        raise ValueError("directions payload lacks structured exact basis/origin values")
    try:
        public_basis = tuple(
            tuple(Fraction(str(value)) for value in values) for values in basis_matrix
        )
        public_origin = tuple(Fraction(str(value)) for value in origin_vector)
    except (ValueError, ZeroDivisionError) as exc:
        raise ValueError("directions embedding basis/origin is not exact rational data") from exc

    magnetic = bool(directions_payload.get("magnetic"))
    child = int(
        subgroup.get("ordinary_number") if magnetic else subgroup.get("number") or 0
    )
    parent_setting = payload_parent.get("setting")
    subgroup_setting = (
        subgroup.get("ordinary_setting") if magnetic else subgroup.get("setting")
    )
    if not isinstance(parent_setting, Mapping) or not isinstance(subgroup_setting, Mapping):
        raise ValueError("directions payload has no selected International settings")
    data = source_tables()
    source_basis, source_origin = data.subgroup_embedding_from_cinter(
        parent,
        child,
        public_basis,
        public_origin,
        parent_setting_id=int(parent_setting["id"]),
        subgroup_setting_id=int(subgroup_setting["id"]),
    )
    displayed = data.subgroup_change_setting_cinter(
        parent,
        child,
        source_basis,
        source_origin,
        parent_setting_id=int(
            (isotropy_state.get("_internal") or {}).get("parent_inter_setting_id")
            or parent_setting["id"]
        ),
        subgroup_setting_id=int(subgroup_setting["id"]),
    )
    denominator = int(displayed["basis_denominator"])
    basis_rows = [
        [Fraction(int(value), denominator) for value in values]
        for values in displayed["basis"]
    ]
    basis_text = isotropy_catalog._basis_text(basis_rows)  # noqa: SLF001
    origin_text = isotropy_catalog._display_origin_text(  # noqa: SLF001
        parent,
        child,
        basis_rows,
        tuple(int(value) for value in displayed["origin"]),
    )

    isotropy = dict(selected_opd.get("isotropy") or {})
    selected_subgroup = dict(isotropy.get("subgroup") or {})
    selected_subgroup["symbol"] = str(subgroup.get("symbol") or selected_subgroup.get("symbol") or "")
    isotropy.update(
        {
            "opd_label": opd_label,
            "direction_label": opd_label,
            "display_opd": str(row.get("direction") or ""),
            "source_opd": str(row.get("direction") or ""),
            "source_matrix": source_rows,
            "source_display_rows": source_rows,
            "source_numeric_rows": source_rows,
            "free": int(row["parameter_count"]),
            "s": int(embedding["cell_index"]),
            "i": int(embedding["index"]),
            "subgroup": selected_subgroup,
            "basis": isotropy_catalog._basis_values(basis_rows),  # noqa: SLF001
            "basis_text": basis_text,
            "basis_denominator": denominator,
            "origin": origin_text,
            "source_basis_values": list(source_basis),
            "source_origin_values": list(source_origin),
            "selected_domain": selected_domain,
            "embedding_selected": True,
            "embedding_direction_row": int(direction_row),
        }
    )
    direction = {
        **dict(selected_opd.get("direction") or {}),
        "label": opd_label,
        "opd": str(row.get("direction") or ""),
        "vectors": source_rows,
    }
    exact_opd = {**selected_opd, "direction": direction, "isotropy": isotropy}
    exact_selected = {
        **selected,
        "orderparam": exact_opd,
    }
    exact_selected["opd_factors"] = _selected_opd_factors(parent, slots, exact_opd)
    return {**isotropy_state, "selected": exact_selected}


def _domain_display(
    isotropy: dict[str, Any],
    *,
    label: str,
    domain: int,
    display_opd: str,
    basis: str | None = None,
    origin: str | None = None,
) -> str:
    display_opd = display_opd.replace(";", ",")
    subgroup = isotropy.get("subgroup") or {}
    subgroup_text = " ".join(
        part
        for part in (
            str(subgroup.get("display_label") or subgroup.get("number") or ""),
            str(subgroup.get("symbol") or ""),
        )
        if part
    )
    details = " ".join(
        part for part in (f"{label}({domain})", display_opd, subgroup_text) if part
    )
    basis_text = str(basis if basis is not None else isotropy.get("basis_text") or "")
    origin_text = str(origin if origin is not None else isotropy.get("origin") or "")
    if basis_text:
        details += f", basis={{{basis_text}}}"
    if origin_text:
        details += f", origin={origin_text}"
    return details


def _complete_coupled_domain_displays(
    *,
    data,
    parent_sg: int,
    parent_setting_id: int,
    slots: list[dict[str, Any]],
    selected_opd: dict[str, Any] | None,
) -> None:
    if not isinstance(selected_opd, dict):
        return
    components = (selected_opd.get("direction") or {}).get("components") or []
    if len(components) != len(slots):
        raise ValueError("selected OPD components do not match primary slots")
    source_data = isotropy_catalog._subgroup_core_data()  # noqa: SLF001
    for slot, component in zip(slots, components, strict=True):
        label = str(component.get("opd") or "")
        source_row = next(
            (
                row
                for row in slot.get("opd_rows") or []
                if str(((row.get("isotropy") or {}).get("opd_label") or "")) == label
            ),
            None,
        )
        if source_row is None:
            raise ValueError(f"selected Source OPD row not found: {label}")
        isotropy = source_row.get("isotropy") or {}
        domain = int(component["domain"])
        basis_text = str(isotropy.get("basis_text") or "")
        origin_text = str(isotropy.get("origin") or "")
        if domain > 1:
            raw_record = component.get("_domain_operation_record")
            if (
                not isinstance(raw_record, list)
                or len(raw_record) != 5
                or any(type(value) is not int for value in raw_record)
            ):
                raise ValueError("selected domain operation record is incomplete")
            record = tuple(raw_record)
            if component.get("_magnetic_domain_operation"):
                point_map = magnetic_data.data().table["mag_point_op_mag2nonmag"]
                record = (*record[:4], int(point_map[int(record[4]) - 1]))
            source_basis = tuple(int(value) for value in isotropy["source_basis_values"])
            source_origin = tuple(int(value) for value in isotropy["source_origin_values"])
            transform = _setting_matrix(source_basis, source_origin)  # type: ignore[arg-type]
            operation = _operation_matrix(source_data, int(parent_sg), record)
            basis_transform = _mat4_mul(operation, transform)
            origin_transform = _mat4_mul(transform, operation)
            basis_values = tuple(
                int(basis_transform[row][column])
                for row in range(3)
                for column in range(3)
            )
            origin = _fraction_mod1(
                tuple(origin_transform[3][column] for column in range(3))
            )
            subgroup = isotropy.get("subgroup") or {}
            ordinary_subgroup = int(
                subgroup.get("ordinary_number") or subgroup["number"]
            )
            displayed = data.subgroup_change_setting_cinter(
                int(parent_sg),
                ordinary_subgroup,
                basis_values,
                origin,
                parent_setting_id=int(parent_setting_id),
                subgroup_setting_id=isotropy_catalog._opd_subgroup_setting_id(  # noqa: SLF001
                    int(parent_sg), ordinary_subgroup, int(parent_setting_id)
                ),
            )
            basis_rows = displayed["basis"]
            basis_text = isotropy_catalog._basis_text(  # noqa: SLF001
                basis_rows, int(displayed["basis_denominator"])
            )
            origin_text = isotropy_catalog._display_origin_text(  # noqa: SLF001
                int(parent_sg), ordinary_subgroup, basis_rows, displayed["origin"]
            )
        component["domain_display"] = _domain_display(
            isotropy,
            label=label,
            domain=domain,
            display_opd=str(component.get("display_opd") or ""),
            basis=basis_text,
            origin=origin_text,
        )


def select_order_parameter_direction(
    reciprocal_state: dict[str, Any],
    *,
    k_params: dict[str, str] | None = None,
    display_k_params: dict[str, str] | None = None,
    irrep: str | None = None,
    irrep_index: int = 1,
    orderparam: int | str | None = None,
    orderparam_index: int = 1,
    selections: list[dict[str, Any]] | None = None,
    opd_row_limit: int | None = None,
    selected_opd_only: bool = False,
) -> dict[str, Any]:
    sg = int(reciprocal_state["input"]["parent"]["number"])
    if selections and len(selections) > 1:
        return _select_coupled_order_parameter(
            reciprocal_state,
            selections=selections,
            orderparam=orderparam,
            orderparam_index=orderparam_index,
            opd_row_limit=opd_row_limit,
            selected_opd_only=selected_opd_only,
        )
    k_params = normalize_k_params(k_params)
    display_k_params = normalize_k_params(display_k_params) or k_params
    local_data = source_tables()
    selection = (selections or [{}])[0]
    selected_irrep = select_irrep(
        reciprocal_state["selected"]["irreps"],
        selection.get("irrep", irrep),
        int(selection.get("irrep_index") or irrep_index),
    )
    old_id = int(selected_irrep["old_id"])
    if selected_irrep.get("magnetic"):
        opd_rows = isotropy_catalog.magnetic_opd_rows(
            old_id,
            parent_sg=sg,
            gid=int(selected_irrep["gid"]),
            k_params=k_params,
            display_k_params=display_k_params,
            selected_orderparam=orderparam if selected_opd_only else None,
            parent_setting_id=int(reciprocal_state["_internal"]["parent_inter_setting_id"]),
            parent_cell=parent_cell_tuple(reciprocal_state["input"]),
        )
    else:
        opd_rows = isotropy_catalog.opd_rows(
            old_id,
            parent_sg=sg,
            gid=int(selected_irrep["gid"]),
            k_params=k_params,
            display_k_params=display_k_params,
            selected_orderparam=orderparam if selected_opd_only else None,
            parent_setting_id=int(reciprocal_state["_internal"]["parent_inter_setting_id"]),
            parent_cell=parent_cell_tuple(reciprocal_state["input"]),
        )
    selected_opd = _select_orderparam(
        opd_rows,
        orderparam=orderparam,
        orderparam_index=orderparam_index,
    )
    selection_slots = [
        {
            "slot": 1,
            "kpoint": reciprocal_state["selected"].get("kpoint"),
            "irrep": selected_irrep,
            "k_params": k_params,
            "display_k_params": display_k_params,
        }
    ]
    return {
        **reciprocal_state,
        "schema": "APOSTRUCT.opd_direction.v1",
        "selected": {
            **reciprocal_state["selected"],
            "irrep": selected_irrep,
            "selection_slots": selection_slots,
            "opd_factors": _selected_opd_factors(sg, selection_slots, selected_opd),
            "image": None if old_id <= 0 else local_data.image_record(old_id),
            "opd_rows": opd_rows,
            "orderparam": selected_opd,
        },
    }


def _select_coupled_order_parameter(
    reciprocal_state: dict[str, Any],
    *,
    selections: list[dict[str, Any]],
    orderparam: int | str | None,
    orderparam_index: int,
    opd_row_limit: int | None,
    selected_opd_only: bool,
) -> dict[str, Any]:
    sg = int(reciprocal_state["input"]["parent"]["number"])
    data = source_tables()
    irrep_slots = list((reciprocal_state.get("selected") or {}).get("irrep_slots") or [])
    if len(irrep_slots) != len(selections):
        raise ValueError("K/IR slot count changed before coupled OPD generation")
    selected_labels = _coupled_selected_labels(
        orderparam,
        slot_count=len(selections),
        selected_opd_only=selected_opd_only,
    )
    slots: list[dict[str, Any]] = []
    for slot_index, (slot_state, selection) in enumerate(zip(irrep_slots, selections), start=1):
        selected_irrep = select_irrep(
            slot_state["irreps"],
            selection.get("irrep"),
            int(selection.get("irrep_index") or 1),
        )
        params = normalize_k_params(selection.get("k_params"))
        display_params = normalize_k_params(selection.get("display_k_params")) or params
        old_id = int(selected_irrep["old_id"])
        selected_label = (
            selected_labels[slot_index - 1]
            if selected_labels is not None
            else None
        )
        if selected_irrep.get("magnetic"):
            rows = isotropy_catalog.magnetic_opd_rows(
                old_id,
                parent_sg=sg,
                gid=int(selected_irrep["gid"]),
                k_params=params,
                display_k_params=display_params,
                selected_orderparam=selected_label,
                parent_setting_id=int(reciprocal_state["_internal"]["parent_inter_setting_id"]),
                parent_cell=parent_cell_tuple(reciprocal_state["input"]),
            )
        else:
            rows = isotropy_catalog.opd_rows(
                old_id,
                parent_sg=sg,
                gid=int(selected_irrep["gid"]),
                k_params=params,
                display_k_params=display_params,
                selected_orderparam=selected_label,
                parent_setting_id=int(reciprocal_state["_internal"]["parent_inter_setting_id"]),
                parent_cell=parent_cell_tuple(reciprocal_state["input"]),
            )
        rows = _rows_for_selected_label(rows, selected_label)
        source_kparam = isotropy_catalog._dynamic_kparam(  # noqa: SLF001
            data,
            gid=int(selected_irrep["gid"]),
            k_params=params,
        ) or (0, 0, 0, 1)
        slots.append({
            "slot": slot_index,
            "kpoint": slot_state["kpoint"],
            "irrep": selected_irrep,
            "k_params": params,
            "display_k_params": display_params,
            "source_kparam": source_kparam,
            "opd_rows": rows,
        })
    rows = coupled_opd_rows(
        isotropy_catalog._subgroup_core_data(),  # noqa: SLF001
        display_data=data,
        parent_sg=sg,
        slots=slots,
        parent_setting_id=int(reciprocal_state["_internal"]["parent_inter_setting_id"]),
        parent_cell=parent_cell_tuple(reciprocal_state["input"]),
        row_limit=opd_row_limit,
    )
    selected_opd = _select_orderparam(
        rows,
        orderparam=orderparam,
        orderparam_index=orderparam_index,
    )
    _complete_coupled_domain_displays(
        data=data,
        parent_sg=sg,
        parent_setting_id=int(reciprocal_state["_internal"]["parent_inter_setting_id"]),
        slots=slots,
        selected_opd=selected_opd,
    )
    for row in rows:
        for component in (row.get("direction") or {}).get("components") or []:
            component.pop("_domain_operation_record", None)
            component.pop("_magnetic_domain_operation", None)
    return {
        **reciprocal_state,
        "schema": "APOSTRUCT.opd_direction.multi.v1",
        "selected": {
            **reciprocal_state["selected"],
            "irrep": slots[0]["irrep"],
            "irreps_selected": [slot["irrep"] for slot in slots],
            "selection_slots": slots,
            "opd_factors": _selected_opd_factors(sg, slots, selected_opd),
            "opd_rows": rows,
            "orderparam": selected_opd,
        },
    }
