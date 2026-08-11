"""CIF-independent subgroup-to-order-parameter directions."""

from __future__ import annotations

from fractions import Fraction
from functools import lru_cache
from typing import Any, Sequence

import numpy as np

from APOSTRUCT.Backend.exactmath import integer_determinant3
from APOSTRUCT.Backend.isotropy.engine.id_subgroup import (
    id_subgroup_identify_with_generator_block,
)
from APOSTRUCT.Backend.isotropy.engine.id_subgroup_magnetic import (
    id_subgroup_magnetic_identify_with_generator_block,
    magnetic_nonmag_point_op,
)
from APOSTRUCT.Backend.isotropy.engine.dynamic_isotropy_file import (
    sort_dynamic_rows_for_file,
)
from APOSTRUCT.Backend.isotropy.engine.get_isotropy import (
    generate_dynamic_isotropy_rows_magnetic,
)
from APOSTRUCT.Backend.isotropy.engine.source_data import SourceData
from APOSTRUCT.Backend.isotropy.catalog import _magnetic_source_orderparam
from APOSTRUCT.Backend.invariants.authority import invariant_source
from APOSTRUCT.Backend.invariants.domains import (
    direction_domain_from_subgroup,
)
from APOSTRUCT.Backend.modes.engine.decoder import ModeDataDecoder
from APOSTRUCT.Backend.modes.engine.dynamic_subduction import (
    display_direction_magnetic_subduction_rows,
    display_direction_subduction_rows,
)
from APOSTRUCT.Backend.modes.engine.subgroup_structure.presentation_transport import (
    embedded_child_operation_records,
    embedded_magnetic_child_operation_records,
)
from APOSTRUCT.Backend.modes.request_context import (
    case_k_parameters_from_source_record,
)
from APOSTRUCT.Backend.parameter_names import parameter_name
from APOSTRUCT.Backend.reciprocal.catalog import k_coordinate_map_for_slot
from APOSTRUCT.Backend.settings import public_setting
from APOSTRUCT.Backend.source.iso_data import (
    reciprocal_to_cinter_matrix_from_table,
)
from APOSTRUCT.Backend.source.magnetic import data as magnetic_data
from APOSTRUCT.Backend.source.tables import source_tables


def _fraction_text(value: Fraction) -> str:
    return str(value.numerator) if value.denominator == 1 else f"{value.numerator}/{value.denominator}"


def _vector_text(values: Sequence[Fraction]) -> str:
    return "(" + ",".join(_fraction_text(Fraction(value)) for value in values) + ")"


def _basis_text(rows: Sequence[Sequence[Fraction]]) -> str:
    return ",".join(_vector_text(tuple(Fraction(value) for value in row)) for row in rows)


def _exact_vector_values(values: Sequence[Fraction]) -> list[str]:
    return [_fraction_text(Fraction(value)) for value in values]


def _exact_basis_values(rows: Sequence[Sequence[Fraction]]) -> list[list[str]]:
    return [_exact_vector_values(row) for row in rows]


def _display_coefficient(value: float) -> str:
    if abs(value) < 1e-10:
        return "0"
    if abs(value - round(value)) < 1e-10:
        return str(int(round(value)))
    return f"{value:.3f}"


def _direction_text(matrix: Sequence[Sequence[float]]) -> str:
    if not matrix:
        return "()"
    names = tuple(parameter_name(index) for index in range(len(matrix[0])))
    components: list[str] = []
    for row in matrix:
        expression = ""
        for value, name in zip(row, names, strict=True):
            coefficient = _display_coefficient(float(value))
            if coefficient == "0":
                continue
            negative = coefficient.startswith("-")
            magnitude = coefficient[1:] if negative else coefficient
            term = name if magnitude == "1" else f"{magnitude}{name}"
            expression += ("-" if negative else ("+" if expression else "")) + term
        components.append(expression or "0")
    return "(" + ",".join(components) + ")"


def _stride48_rows(matrix: Sequence[Sequence[float]]) -> tuple[float, ...]:
    if not matrix:
        return ()
    dimension = len(matrix)
    values: list[float] = []
    for parameter in range(len(matrix[0])):
        values.extend(float(matrix[coordinate][parameter]) for coordinate in range(dimension))
        values.extend(0.0 for _ in range(48 - dimension))
    return tuple(values)


def _reciprocal_cinter(
    decoder: ModeDataDecoder,
    parent_sg: int,
    vector: Sequence[Fraction],
    setting_id: int | None = None,
) -> tuple[Fraction, Fraction, Fraction]:
    matrix = reciprocal_to_cinter_matrix_from_table(
        decoder.iso.space,
        int(parent_sg),
        "pml",
        setting_id,
    )
    return tuple(
        sum(Fraction(vector[row]) * matrix[row][column] for row in range(3))
        for column in range(3)
    )  # type: ignore[return-value]


def _subgroup_index(
    decoder: ModeDataDecoder,
    parent_sg: int,
    child_sg: int,
    cell_index: int,
) -> int:
    numerator = int(cell_index) * decoder.tables.space_group_point_group_order(int(parent_sg))
    denominator = decoder.tables.space_group_point_group_order(int(child_sg))
    if denominator <= 0 or numerator % denominator:
        raise ValueError("subgroup embedding has a non-integral group index")
    return numerator // denominator


def _ordinary_stabilizer_symbol(
    tables: Any,
    stabilizer_sg: int,
    stabilizer_cell_index: int,
    *,
    parent_sg: int,
    parent_setting: dict[str, Any],
    requested_subgroup_sg: int,
    requested_cell_index: int,
    subgroup_setting: dict[str, Any],
) -> str:
    if (
        int(stabilizer_sg) == int(requested_subgroup_sg)
        and (
            int(requested_subgroup_sg) != int(parent_sg)
            or int(stabilizer_cell_index) == int(requested_cell_index)
        )
    ):
        return str(subgroup_setting["symbol"])
    if int(stabilizer_sg) == int(parent_sg):
        return str(parent_setting["symbol"])
    return tables.default_setting_space_symbol(int(stabilizer_sg))


def _validate_embedded_operations(
    decoder: ModeDataDecoder,
    parent_sg: int,
    operations: Sequence[Sequence[int]],
) -> None:
    parent_by_point = {
        int(record[4]): record for record in decoder.generate_space_group_records(int(parent_sg))
    }
    for operation in operations:
        parent = parent_by_point.get(int(operation[4]))
        if parent is None or any(
            (
                Fraction(int(operation[axis]), int(operation[3]))
                - Fraction(int(parent[axis]), int(parent[3]))
            ).denominator
            != 1
            for axis in range(3)
        ):
            raise ValueError("not all elements of the subgroup are elements of the parent group")


def _validate_embedded_magnetic_operations(
    decoder: ModeDataDecoder,
    parent_sg: int,
    operations: Sequence[Sequence[int]],
) -> None:
    parent_by_point = {
        int(record[4]): record for record in decoder.generate_space_group_records(int(parent_sg))
    }
    for operation in operations:
        parent = parent_by_point.get(magnetic_nonmag_point_op(int(operation[4])))
        if parent is None or any(
            (
                Fraction(int(operation[axis]), int(operation[3]))
                - Fraction(int(parent[axis]), int(parent[3]))
            ).denominator
            != 1
            for axis in range(3)
        ):
            raise ValueError(
                "not all elements of the magnetic subgroup are elements of the gray parent"
            )


def _magnetic_group_from_bns_number(number: str) -> int:
    wanted = str(number).strip()
    table = magnetic_data().table
    matches = [
        index
        for index, value in enumerate(table["mag_nlabel"], start=1)
        if str(value).strip() == wanted
    ]
    if len(matches) != 1:
        raise ValueError(f"unknown BNS magnetic group number: {number}")
    return int(matches[0])


def _magnetic_direction_label_domain(
    gid: int,
    matrix: Sequence[Sequence[float]],
    *,
    kparam: Sequence[int],
    include_domain: bool = True,
) -> tuple[str, int | None] | None:
    old_id = int(source_tables().little["little_irr_old"][int(gid) - 1])
    tolerance = 1e-7 if old_id > 0 else 1e-6
    authority = invariant_source()
    parent_sg = int(source_tables().little["little_irr_space_group"][int(gid) - 1])
    matches = []
    for label, candidate, child_sg, basis, origin in _magnetic_direction_records(
        int(gid), tuple(kparam)
    ):
        domain = direction_domain_from_subgroup(
            authority,
            sg=parent_sg,
            gid=int(gid),
            kparam=tuple(kparam),
            target=matrix,
            candidate=candidate,
            child_sg=child_sg,
            basis=basis,
            origin=origin,
            tolerance=tolerance,
        )
        if domain is not None:
            matches.append((label, domain))
    labels = {label for label, _domain in matches}
    if not matches:
        return None
    if len(labels) != 1:
        raise ValueError(f"carrier subspace matches multiple magnetic Source OPDs for gid={gid}")
    label = matches[0][0]
    domain = min(domain for matched, domain in matches if matched == label)
    return label, domain if include_domain else None


@lru_cache(maxsize=512)
def _magnetic_direction_records(
    gid: int,
    kparam: tuple[int, ...],
) -> tuple[tuple[str, np.ndarray, int, tuple[int, ...], tuple[int, ...]], ...]:
    tables = source_tables()
    old_id = int(tables.little["little_irr_old"][int(gid) - 1])
    magnetic = magnetic_data()
    if old_id > 0:
        return tuple(
            (
                str(row.orderparam_label),
                np.asarray(
                    (
                        source := _magnetic_source_orderparam(tables, row)
                    )["numeric_rows"][: int(source["free"])],
                    dtype=float,
                ).T,
                int(magnetic.table["mag_space_group"][row.subgroup_magnetic_group - 1]),
                tuple(int(value) for value in row.basis),
                tuple(int(value) for value in row.origin),
            )
            for row in magnetic.magnetic_isotropy_rows_for_irrep(old_id)
        )
    rows = sort_dynamic_rows_for_file(
        generate_dynamic_isotropy_rows_magnetic(
            invariant_source().source_data,
            gid=int(gid),
            kparam=kparam,
        )
    )
    return tuple(
        (
            str(row.direction).strip(),
            np.asarray(row.matrix, dtype=float).T,
            int(magnetic.table["mag_space_group"][int(row.subgroup_number) - 1]),
            tuple(int(value) for value in row.basis_values),
            tuple(int(value) for value in row.origin_values),
        )
        for row in rows
        if row.matrix
    )


def _magnetic_subgroup_index(
    parent_magnetic_group: int,
    subgroup_magnetic_group: int,
    cell_index: int,
) -> int:
    table = magnetic_data().table
    parent_point_group = int(table["mag_point_group"][int(parent_magnetic_group) - 1])
    child_point_group = int(table["mag_point_group"][int(subgroup_magnetic_group) - 1])
    parent_order = int(table["mag_point_group_order"][parent_point_group - 1])
    child_order = int(table["mag_point_group_order"][child_point_group - 1])
    numerator = int(cell_index) * parent_order
    if child_order <= 0 or numerator % child_order:
        raise ValueError("magnetic subgroup embedding has a non-integral group index")
    return numerator // child_order


def compatible_directions(
    parent_sg: int,
    subgroup_sg: int,
    basis: Sequence[Sequence[Fraction]],
    origin: Sequence[Fraction] = (Fraction(0), Fraction(0), Fraction(0)),
    *,
    parent_setting_id: int | None = None,
    subgroup_setting_id: int | None = None,
) -> dict[str, Any]:
    """Return ordinary OP directions compatible with one exact G:H embedding."""

    parent = int(parent_sg)
    subgroup = int(subgroup_sg)
    if not 1 <= parent <= 230 or not 1 <= subgroup <= 230:
        raise ValueError("parent and subgroup space groups must be integers from 1 through 230")

    tables = source_tables()
    parent_setting = public_setting(parent, parent_setting_id)
    subgroup_setting = public_setting(subgroup, subgroup_setting_id)
    decoder = ModeDataDecoder(tables=tables)
    source = SourceData(tables=decoder.iso)
    public_basis = tuple(tuple(Fraction(value) for value in row) for row in basis)
    public_origin = tuple(Fraction(value) for value in origin)
    source_basis, source_origin = tables.subgroup_embedding_from_cinter(
        parent,
        subgroup,
        public_basis,
        public_origin,
        parent_setting_id=int(parent_setting["id"]),
        subgroup_setting_id=int(subgroup_setting["id"]),
    )
    child_operations = embedded_child_operation_records(
        decoder,
        parent_sg=parent,
        child_sg=subgroup,
        subgroup_basis=source_basis,
        subgroup_origin=source_origin,
    )
    _validate_embedded_operations(decoder, parent, child_operations)

    input_cell_index = abs(integer_determinant3(source_basis))
    input_index = _subgroup_index(decoder, parent, subgroup, input_cell_index)
    rows = display_direction_subduction_rows(
        decoder,
        sg=parent,
        basis=source_basis,
        operations=child_operations,
        irrep_source=source,
    )
    direction_source = invariant_source()
    directions: list[dict[str, Any]] = []
    for row_number, row in enumerate(rows, start=1):
        matrix = tuple(tuple(float(value) for value in values) for values in row.direction_matrix)
        parameter_count = len(matrix[0]) if matrix else 0
        stabilizer_basis, stabilizer_operations = source.orderparam_to_subgroup(
            int(row.gid),
            tuple(int(value) for value in row.source_kparam),
            _stride48_rows(matrix),
            parameter_count,
        )
        identified = id_subgroup_identify_with_generator_block(
            source,
            parent,
            stabilizer_basis,
            stabilizer_operations,
            0,
        )
        if identified is None:
            raise ValueError(f"could not identify stabilizer for {row.irrep_label}")
        cell_index = abs(integer_determinant3(identified.basis))
        row_index = _subgroup_index(decoder, parent, int(identified.subgroup), cell_index)
        pml_k = tuple(Fraction(value) for value in row.reciprocal_vector_pml)
        cinter_k = _reciprocal_cinter(
            decoder,
            parent,
            pml_k,
            int(parent_setting["id"]),
        )
        kslot = int(decoder.iso.little["little_irr_k"][int(row.gid) - 1])
        lattice = int(decoder.iso.space["ispace_lattice"][parent - 1])
        k_dimension = int(decoder.iso.little["little_k_dim"][(lattice - 1) * 27 + kslot - 1])
        public_k_parameters: dict[str, str] = {}
        miller_love_parameters: dict[str, str] = {}
        if k_dimension > 0:
            public_values = case_k_parameters_from_source_record(
                decoder,
                gid=int(row.gid),
                source_kparam=row.source_kparam,
            )
            if public_values is None:
                raise ValueError(f"could not recover public k parameters for {row.irrep_label}")
            public_names = k_coordinate_map_for_slot(parent, kslot).parameter_names
            public_k_parameters = {
                name: _fraction_text(value)
                for name, value in zip(public_names, public_values, strict=True)
            }
            miller_love_parameters = {
                name: _fraction_text(
                    Fraction(int(row.source_kparam[index]), int(row.source_kparam[3]))
                )
                for index, name in enumerate(("alpha", "beta", "gamma")[:k_dimension])
            }
        primary = int(identified.subgroup) == subgroup and cell_index == input_cell_index
        if primary:
            source_direction = direction_source.direction_label_domain_for_matrix(
                parent,
                int(row.gid),
                matrix,
                kparam=row.source_kparam,
            )
        else:
            label = direction_source.direction_label_for_matrix(
                int(row.gid),
                matrix,
                kparam=row.source_kparam,
            )
            source_direction = (label, None) if label is not None else None
        directions.append(
            {
                "row": row_number,
                "k_label": str(row.k_label),
                "kvector": _vector_text(cinter_k),
                "miller_love_kvector": _vector_text(pml_k),
                "k_parameters": public_k_parameters,
                "miller_love_parameters": miller_love_parameters,
                "irrep": str(row.irrep_label),
                "gid": int(row.gid),
                "full_dim": len(matrix),
                "parameter_count": parameter_count,
                "parameters": [parameter_name(index) for index in range(parameter_count)],
                "opd": source_direction[0] if source_direction is not None else None,
                "domain": source_direction[1] if source_direction is not None else None,
                "direction": _direction_text(matrix),
                "direction_matrix": [list(values) for values in matrix],
                "role": "primary" if primary else "secondary",
                "subgroup": {
                    "number": int(identified.subgroup),
                    "symbol": _ordinary_stabilizer_symbol(
                        tables,
                        int(identified.subgroup),
                        cell_index,
                        parent_sg=parent,
                        parent_setting=parent_setting,
                        requested_subgroup_sg=subgroup,
                        requested_cell_index=input_cell_index,
                        subgroup_setting=subgroup_setting,
                    ),
                },
                "index": row_index,
                "cell_index": cell_index,
            }
        )

    return {
        "schema": "APOSTRUCT.cli.directions",
        "parent": {
            "number": parent,
            "symbol": str(parent_setting["symbol"]),
            "setting": parent_setting,
        },
        "subgroup": {
            "number": subgroup,
            "symbol": str(subgroup_setting["symbol"]),
            "setting": subgroup_setting,
        },
        "embedding": {
            "basis": _basis_text(public_basis),
            "origin": _vector_text(public_origin),
            "basis_matrix": _exact_basis_values(public_basis),
            "origin_vector": _exact_vector_values(public_origin),
            "index": input_index,
            "cell_index": input_cell_index,
        },
        "directions": directions,
    }


def compatible_magnetic_directions(
    parent_sg: int,
    subgroup_bns_number: str,
    basis: Sequence[Sequence[Fraction]],
    origin: Sequence[Fraction] = (Fraction(0), Fraction(0), Fraction(0)),
    *,
    parent_setting_id: int | None = None,
    subgroup_setting_id: int | None = None,
) -> dict[str, Any]:
    """Return time-odd OP directions compatible with one exact BNS embedding."""

    parent = int(parent_sg)
    if not 1 <= parent <= 230:
        raise ValueError("parent space group must be an integer from 1 through 230")

    tables = source_tables()
    decoder = ModeDataDecoder(tables=tables)
    source = SourceData(tables=decoder.iso)
    table = magnetic_data().table
    subgroup_magnetic_group = _magnetic_group_from_bns_number(subgroup_bns_number)
    subgroup_sg = int(table["mag_space_group"][subgroup_magnetic_group - 1])
    parent_setting = public_setting(parent, parent_setting_id)
    subgroup_setting = public_setting(subgroup_sg, subgroup_setting_id)
    public_basis = tuple(tuple(Fraction(value) for value in row) for row in basis)
    public_origin = tuple(Fraction(value) for value in origin)
    source_basis, source_origin = tables.subgroup_embedding_from_cinter(
        parent,
        subgroup_sg,
        public_basis,
        public_origin,
        parent_setting_id=int(parent_setting["id"]),
        subgroup_setting_id=int(subgroup_setting["id"]),
    )
    child_operations = embedded_magnetic_child_operation_records(
        decoder,
        parent_sg=parent,
        child_magnetic_group=subgroup_magnetic_group,
        subgroup_basis=source_basis,
        subgroup_origin=source_origin,
    )
    _validate_embedded_magnetic_operations(decoder, parent, child_operations)

    parent_magnetic_group = source.magnetic_parent_group_for_space_group(parent) + 1
    input_cell_index = abs(integer_determinant3(source_basis))
    input_index = _magnetic_subgroup_index(
        parent_magnetic_group,
        subgroup_magnetic_group,
        input_cell_index,
    )
    rows = display_direction_magnetic_subduction_rows(
        decoder,
        source,
        sg=parent,
        basis=source_basis,
        operations=child_operations,
    )
    directions: list[dict[str, Any]] = []
    for row_number, row in enumerate(rows, start=1):
        matrix = tuple(tuple(float(value) for value in values) for values in row.direction_matrix)
        parameter_count = len(matrix[0]) if matrix else 0
        stabilizer_basis, stabilizer_operations = source.orderparam_to_subgroup_magnetic(
            int(row.gid),
            tuple(int(value) for value in row.source_kparam),
            _stride48_rows(matrix),
            parameter_count,
        )
        identified = id_subgroup_magnetic_identify_with_generator_block(
            source,
            parent,
            stabilizer_basis,
            stabilizer_operations,
            0,
        )
        if identified is None:
            raise ValueError(f"could not identify magnetic stabilizer for m{row.irrep_label}")
        cell_index = abs(integer_determinant3(identified.basis))
        row_index = _magnetic_subgroup_index(
            parent_magnetic_group,
            int(identified.magnetic_group),
            cell_index,
        )
        pml_k = tuple(Fraction(value) for value in row.reciprocal_vector_pml)
        cinter_k = _reciprocal_cinter(
            decoder,
            parent,
            pml_k,
            int(parent_setting["id"]),
        )
        kslot = int(decoder.iso.little["little_irr_k"][int(row.gid) - 1])
        lattice = int(decoder.iso.space["ispace_lattice"][parent - 1])
        k_dimension = int(
            decoder.iso.little["little_k_dim"][(lattice - 1) * 27 + kslot - 1]
        )
        public_k_parameters: dict[str, str] = {}
        miller_love_parameters: dict[str, str] = {}
        if k_dimension > 0:
            public_values = case_k_parameters_from_source_record(
                decoder,
                gid=int(row.gid),
                source_kparam=row.source_kparam,
            )
            if public_values is None:
                raise ValueError(
                    f"could not recover public k parameters for m{row.irrep_label}"
                )
            public_names = k_coordinate_map_for_slot(parent, kslot).parameter_names
            public_k_parameters = {
                name: _fraction_text(value)
                for name, value in zip(public_names, public_values, strict=True)
            }
            miller_love_parameters = {
                name: _fraction_text(
                    Fraction(int(row.source_kparam[index]), int(row.source_kparam[3]))
                )
                for index, name in enumerate(("alpha", "beta", "gamma")[:k_dimension])
            }
        identified_group = int(identified.magnetic_group)
        primary = (
            identified_group == subgroup_magnetic_group
            and cell_index == input_cell_index
        )
        source_direction = _magnetic_direction_label_domain(
            int(row.gid),
            matrix,
            kparam=row.source_kparam,
            include_domain=primary,
        )
        directions.append(
            {
                "row": row_number,
                "k_label": str(row.k_label),
                "kvector": _vector_text(cinter_k),
                "miller_love_kvector": _vector_text(pml_k),
                "k_parameters": public_k_parameters,
                "miller_love_parameters": miller_love_parameters,
                "irrep": "m" + str(row.irrep_label),
                "gid": int(row.gid),
                "full_dim": len(matrix),
                "parameter_count": parameter_count,
                "parameters": [parameter_name(index) for index in range(parameter_count)],
                "opd": source_direction[0] if source_direction is not None else None,
                "domain": source_direction[1] if source_direction is not None else None,
                "direction": _direction_text(matrix),
                "direction_matrix": [list(values) for values in matrix],
                "role": "primary" if primary else "secondary",
                "subgroup": {
                    "number": identified_group,
                    "bns_number": str(table["mag_nlabel"][identified_group - 1]).strip(),
                    "symbol": str(table["mag_bns_label"][identified_group - 1]).strip(),
                    "ordinary_number": int(table["mag_space_group"][identified_group - 1]),
                },
                "index": row_index,
                "cell_index": cell_index,
            }
        )

    return {
        "schema": "APOSTRUCT.cli.directions",
        "magnetic": True,
        "parent": {
            "number": parent,
            "symbol": str(parent_setting["symbol"]),
            "setting": parent_setting,
            "gray_magnetic_group": parent_magnetic_group,
        },
        "subgroup": {
            "number": subgroup_magnetic_group,
            "bns_number": str(table["mag_nlabel"][subgroup_magnetic_group - 1]).strip(),
            "symbol": str(table["mag_bns_label"][subgroup_magnetic_group - 1]).strip(),
            "ordinary_number": subgroup_sg,
            "ordinary_setting": subgroup_setting,
        },
        "embedding": {
            "basis": _basis_text(public_basis),
            "origin": _vector_text(public_origin),
            "basis_matrix": _exact_basis_values(public_basis),
            "origin_vector": _exact_vector_values(public_origin),
            "index": input_index,
            "cell_index": input_cell_index,
        },
        "directions": directions,
    }


__all__ = ["compatible_directions", "compatible_magnetic_directions"]
