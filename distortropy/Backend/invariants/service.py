"""Validated service boundary for the standalone INVARIANTS interface."""

from __future__ import annotations

from fractions import Fraction
from functools import lru_cache
import math
import re
import threading
import time
from typing import Any, Mapping, Sequence

from distortropy.Backend.isotropy.catalog import (
    _basis_text,
    _display_origin_text,
    magnetic_opd_rows as catalog_magnetic_opd_rows,
    opd_rows as catalog_opd_rows,
)
from distortropy.Backend.isotropy.engine.get_isotropy import (
    _equation_key,
)
from distortropy.Backend.parameter_names import (
    parameter_name as _parameter_name,
)
from distortropy.Backend.reciprocal.catalog import (
    ensure_source_k_not_special,
    kpoints as reciprocal_kpoints,
)
from distortropy.Backend.source.magnetic import data as magnetic_source_data
from distortropy.Backend.source.tables import source_tables
from .domains import (
    _domain_transform_matrix,
    _fraction_record_mod1,
    _mat4_mul,
    _pml_operation_matrix,
    domain_count_from_isotropy_row,
    domain_count_from_subgroup,
    domain_operation_record_from_isotropy_row,
    domain_operation_record_from_subgroup,
)
from . import (
    coupled_fixed_gradient_action,
    coupled_fixed_irrep_matrices,
    coupled_parametric_gradient_action,
    coupled_parametric_irrep_matrices,
    fixed_irrep_dimensions,
    gradient_invariant_basis,
    gradient_polynomial_text,
    invariant_basis,
    polynomial_terms,
    polynomial_text,
    restricted_fixed_irrep_action,
    restricted_parametric_irrep_action,
)
from .source import (
    _exact_record_matrix,
    _quadratic3_record_matrix,
    _rational_entry,
    _source_display_entry,
)
from .subspace import (
    _exact_direction_matrix,
    _matrix_product,
    _project_action,
)
from .authority import InvariantSource, invariant_source


MAX_FACTORS = 12
MAX_DEGREE = 12
MAX_GRADIENT_ORDER = 6


def _integer(value: object, name: str, *, minimum: int, maximum: int) -> int:
    if type(value) is not int or not minimum <= value <= maximum:
        raise ValueError(f"{name} must be an integer from {minimum} through {maximum}")
    return value


def _fraction(value: object, name: str) -> Fraction:
    if type(value) is int:
        return Fraction(value)
    if not isinstance(value, str) or not value.strip():
        raise ValueError(f"{name} must be an exact rational string or integer")
    try:
        return Fraction(value.strip())
    except (ValueError, ZeroDivisionError) as exc:
        raise ValueError(f"{name} is not an exact rational value") from exc


def _fraction_text(value: Fraction) -> str:
    if value.denominator == 1:
        return str(value.numerator)
    return f"{value.numerator}/{value.denominator}"


def _canonical_k_expression(expression: str) -> tuple[str, tuple[str, ...]]:
    """Use the Web's a,b,g parameter names in first-occurrence order."""

    source_names: list[str] = []
    for name in re.findall(r"(?<![A-Za-z])([abg])(?![A-Za-z])", expression):
        if name not in source_names:
            source_names.append(name)
    public_names = ("a", "b", "g")[:len(source_names)]
    replacements = dict(zip(source_names, public_names, strict=True))
    display = re.sub(
        r"(?<![A-Za-z])([abg])(?![A-Za-z])",
        lambda match: replacements[match.group(1)],
        expression,
    )
    return display, tuple(source_names)


def _kovalev_k_label(irreps: Sequence[Mapping[str, Any]]) -> str:
    for irrep in irreps:
        match = re.match(r"(k\d+)", str(irrep.get("kov_label") or ""))
        if match:
            return match.group(1)
    return ""


def _display_coefficient(value: object) -> str:
    numeric = float(value)
    if abs(numeric) < 1e-10:
        return "0"
    if abs(numeric - round(numeric)) < 1e-10:
        return str(int(round(numeric)))
    return f"{numeric:.3f}"


def _display_direction_matrix(
    matrix: Sequence[Sequence[object]],
    *,
    parameter_offset: int = 0,
) -> str:
    """Format full-dimension x free-dimension columns as a Web OPD tuple."""

    if not matrix:
        return "()"
    width = len(matrix[0])
    names = tuple(_parameter_name(parameter_offset + index) for index in range(width))
    if any(len(row) != width for row in matrix):
        raise ValueError("OPD direction matrix cannot be displayed")
    components: list[str] = []
    for row in matrix:
        expression = ""
        for value, name in zip(row, names[:width], strict=True):
            coefficient = _display_coefficient(value)
            if coefficient == "0":
                continue
            negative = coefficient.startswith("-")
            magnitude = coefficient[1:] if negative else coefficient
            term = name if magnitude == "1" else f"{magnitude}{name}"
            if not expression:
                expression = f"-{term}" if negative else term
            else:
                expression += f"{'-' if negative else '+'}{term}"
        components.append(expression or "0")
    return "(" + ",".join(components) + ")"


def _parameter_polynomial_text(
    polynomial: Mapping[tuple[int, ...], object],
    parameter_names: Sequence[str] | None = None,
) -> str:
    internal = polynomial_text(polynomial)
    return re.sub(
        r"n([1-9]\d*)",
        lambda match: (
            parameter_names[int(match.group(1)) - 1]
            if parameter_names is not None
            else _parameter_name(int(match.group(1)) - 1)
        ),
        internal,
    )


def _time_reversal_action(
    matrices: Sequence[Sequence[Sequence[object]]],
    block_dimensions: Sequence[int],
    magnetic: Sequence[bool],
) -> tuple[tuple[tuple[object, ...], ...], ...]:
    """Extend a spatial action by the central time-reversal operation."""

    dimensions = tuple(int(value) for value in block_dimensions)
    flags = tuple(bool(value) for value in magnetic)
    if len(dimensions) != len(flags) or any(value < 1 for value in dimensions):
        raise ValueError("IR block dimensions and magnetic flags are inconsistent")
    original = tuple(tuple(tuple(row) for row in matrix) for matrix in matrices)
    if not any(flags):
        return original
    signs = tuple(
        Fraction(-1 if flag else 1)
        for width, flag in zip(dimensions, flags, strict=True)
        for _index in range(width)
    )
    if not original or len(original[0]) != len(signs):
        raise ValueError("time-reversal action does not cover the representation")
    reversed_action = tuple(
        tuple(
            tuple(signs[row] * value for value in values)
            for row, values in enumerate(matrix)
        )
        for matrix in original
    )
    return original + reversed_action


def _gradient_time_reversal_action(
    matrices: Sequence[Sequence[Sequence[object]]],
    variables: Sequence[tuple[int, tuple[int, ...]]],
    irrep_dimensions: Sequence[int],
    magnetic: Sequence[bool],
) -> tuple[tuple[tuple[object, ...], ...], ...]:
    dimensions = tuple(int(value) for value in irrep_dimensions)
    component_flags = tuple(
        flag
        for width, flag in zip(dimensions, magnetic, strict=True)
        for _index in range(width)
    )
    variable_flags = tuple(component_flags[int(component)] for component, _word in variables)
    return _time_reversal_action(matrices, (1,) * len(variable_flags), variable_flags)


class InvariantService:
    """Expose the exact DISPLAY INVARIANT kernel without changing its algebra."""

    def __init__(self, projection: InvariantSource | None = None):
        self.projection = projection or invariant_source()
        self.source_data = self.projection.source_data
        self._lock = threading.RLock()

    def bootstrap(self) -> dict[str, Any]:
        labels = self.projection.space["space_label_bc"]
        classes = self.projection.space["space_label"]
        return {
            "space_groups": [
                {
                    "number": number,
                    "symbol": str(labels[number - 1]).strip(),
                    "crystal_class": str(classes[number - 1]).strip(),
                    "display": (
                        f"{number} {str(labels[number - 1]).strip()} "
                        f"{str(classes[number - 1]).strip()}"
                    ),
                }
                for number in range(1, 231)
            ],
            "limits": {
                "factors": MAX_FACTORS,
                "degree": MAX_DEGREE,
                "gradient_order": MAX_GRADIENT_ORDER,
            },
        }

    def catalog(self, space_group: object) -> dict[str, Any]:
        sg = _integer(space_group, "space_group", minimum=1, maximum=230)
        reciprocal = reciprocal_kpoints(sg)
        magnetic = magnetic_source_data()
        k_groups: list[dict[str, Any]] = []
        for source_k in reciprocal["kpoints"]:
            irreps = source_k["irreps"]
            if not irreps:
                continue
            expression, source_parameter_names = _canonical_k_expression(
                str(source_k["display_kvector"])
            )
            kovalev = _kovalev_k_label(irreps)
            k_groups.append({
                "label": source_k["label"],
                "kovalev": kovalev,
                "expression": expression,
                "display": ", ".join(
                    part for part in (str(source_k["label"]), kovalev, expression) if part
                ),
                "kind": "fixed" if int(source_k["dimension"]) == 0 else "parametric",
                "parameter_count": int(source_k["dimension"]),
                "parameter_names": list(("a", "b", "g")[:int(source_k["dimension"])]),
                "source_parameter_names": list(source_parameter_names),
                "irreps": [
                    {
                        "gid": int(irrep["gid"]),
                        "old_id": int(irrep["old_id"]),
                        "label": str(irrep["symbol"]),
                        "kovalev": str(irrep["kov_label"]),
                        "display": ", ".join(
                            part for part in (str(irrep["symbol"]), str(irrep["kov_label"]))
                            if part
                        ),
                        "dimension": int(irrep["full_dim"]),
                        "type": int(irrep["type"]),
                        "magnetic_available": (
                            int(irrep["old_id"]) <= 0
                            or magnetic.irrep_has_magnetic_isotropy(int(irrep["old_id"]))
                        ),
                    }
                    for irrep in irreps
                ],
            })
        return {
            "space_group": sg,
            "symbol": str(self.projection.space["space_label_bc"][sg - 1]).strip(),
            "crystal_class": str(self.projection.space["space_label"][sg - 1]).strip(),
            "k_groups": k_groups,
        }

    def opds(self, payload: Mapping[str, Any]) -> dict[str, Any]:
        sg = _integer(payload.get("space_group"), "space_group", minimum=1, maximum=230)
        factor = self._resolve_factor(sg, payload)
        with self._lock:
            options = self._opd_options(sg, factor)
        return {"gid": factor["gid"], "label": factor["label"], "opds": options}

    def _opd_options(self, sg: int, factor: Mapping[str, Any]) -> list[dict[str, Any]]:
        gid = int(factor["gid"])
        old_id = int(factor["old_id"])
        is_magnetic = bool(factor["magnetic"])
        rich_by_label: dict[str, Mapping[str, Any]] = {}
        if old_id > 0:
            catalog = catalog_magnetic_opd_rows if is_magnetic else catalog_opd_rows
            labels = []
            for row in catalog(old_id, parent_sg=sg, gid=gid):
                isotropy = row.get("isotropy")
                if isinstance(isotropy, Mapping):
                    label = str(isotropy.get("opd_label") or "")
                    if label and label not in rich_by_label:
                        labels.append(label)
                        rich_by_label[label] = isotropy
            if not is_magnetic:
                labels = self.projection.list_direction_labels(gid)
        else:
            labels = []
            for row in self._dynamic_catalog_rows(
                sg,
                gid,
                tuple(factor["parameters"]),
                is_magnetic,
            ):
                isotropy = row.get("isotropy")
                if not isinstance(isotropy, Mapping):
                    continue
                label = str(isotropy.get("opd_label") or "")
                if label and label not in rich_by_label:
                    labels.append(label)
                    rich_by_label[label] = isotropy

        options = []
        for label in labels:
            rich = rich_by_label.get(label)
            if is_magnetic:
                magnetic_index = int((rich or {}).get("i") or 0)
                if magnetic_index < 2 or magnetic_index % 2:
                    raise ValueError(f"Source magnetic OPD {label!r} has invalid domain index")
                domain_count = magnetic_index // 2
            elif old_id <= 0 and rich is not None:
                domain_count = int(rich.get("i") or 0)
            else:
                domain_count = self._domain_count(sg, factor, label)
            if domain_count < 1:
                raise ValueError(f"Source OPD {label!r} has no domains")
            if rich is None:
                display_opd = _display_direction_matrix(
                    self._direction_matrix(sg, factor, label)
                )
                subgroup_text = ""
                basis_text = ""
                origin = ""
            else:
                display_opd = str(
                    rich.get("display_opd") or rich.get("source_opd") or ""
                ).replace(";", ",")
                subgroup = rich.get("subgroup") or {}
                subgroup_text = " ".join(
                    part
                    for part in (
                        str(subgroup.get("display_label") or subgroup.get("number") or ""),
                        str(subgroup.get("symbol") or ""),
                    )
                    if part
                )
                basis_text = str(rich.get("basis_text") or "")
                origin = str(rich.get("origin") or "")
            details = " ".join(part for part in (display_opd, subgroup_text) if part)
            if basis_text:
                details += f", basis={{{basis_text}}}"
            if origin:
                details += f", origin={origin}"
            options.append({
                "label": label,
                "display_opd": display_opd,
                "subgroup": subgroup_text,
                "basis": basis_text,
                "origin": origin,
                "display": f"{label} {details}".strip(),
                "domain_count": domain_count,
                "parameter_count": len(self._direction_matrix(sg, factor, label)[0]),
            })
        return options

    @lru_cache(maxsize=32)
    def _dynamic_catalog_rows(
        self,
        sg: int,
        gid: int,
        parameters: tuple[Fraction, ...],
        magnetic: bool = False,
    ) -> tuple[Mapping[str, Any], ...]:
        names = ("a", "b", "g")
        k_params = {
            names[index]: _fraction_text(value)
            for index, value in enumerate(parameters)
        }
        catalog = catalog_magnetic_opd_rows if magnetic else catalog_opd_rows
        return tuple(catalog(
            0,
            parent_sg=int(sg),
            gid=int(gid),
            k_params=k_params,
            display_k_params=k_params,
        ))

    @lru_cache(maxsize=256)
    def _static_catalog_rows(
        self,
        sg: int,
        old_id: int,
        gid: int,
        magnetic: bool = False,
    ) -> tuple[Mapping[str, Any], ...]:
        catalog = catalog_magnetic_opd_rows if magnetic else catalog_opd_rows
        return tuple(catalog(
            int(old_id),
            parent_sg=int(sg),
            gid=int(gid),
        ))

    def _dynamic_isotropy(
        self,
        sg: int,
        factor: Mapping[str, Any],
        direction: str,
    ) -> Mapping[str, Any]:
        if int(factor["old_id"]) > 0:
            rows = self._static_catalog_rows(
                int(sg),
                int(factor["old_id"]),
                int(factor["gid"]),
                bool(factor["magnetic"]),
            )
        else:
            rows = self._dynamic_catalog_rows(
                int(sg),
                int(factor["gid"]),
                tuple(factor["parameters"]),
                bool(factor["magnetic"]),
            )
        for row in rows:
            isotropy = row.get("isotropy")
            if (
                isinstance(isotropy, Mapping)
                and str(isotropy.get("opd_label") or "") == direction
            ):
                return isotropy
        raise ValueError(
            f"Source OPD row not found for gid={factor['gid']}: {direction}"
        )

    def domains(self, payload: Mapping[str, Any]) -> dict[str, Any]:
        sg = _integer(payload.get("space_group"), "space_group", minimum=1, maximum=230)
        factor = self._resolve_factor(sg, payload)
        direction = self._direction(factor)
        parameter_offset = _integer(
            payload.get("parameter_offset", 0),
            "parameter_offset",
            minimum=0,
            maximum=10_000,
        )
        with self._lock:
            options = {
                item["label"]: item["domain_count"]
                for item in self._opd_options(sg, factor)
            }
            if direction not in options:
                raise ValueError(f"unknown OPD {direction!r} for {factor['label']}")
            domains = self._domain_options(
                sg,
                factor,
                direction,
                options[direction],
                parameter_offset=parameter_offset,
            )
        return {
            "gid": factor["gid"],
            "label": factor["label"],
            "direction": direction,
            "domains": domains,
        }

    @staticmethod
    def _subspace_key(matrix: object) -> tuple[int, tuple[float, ...]]:
        if (
            isinstance(matrix, (str, bytes, bytearray))
            or not isinstance(matrix, Sequence)
            or not matrix
        ):
            raise ValueError("direction_matrix must be a non-empty matrix")
        rows: list[list[float]] = []
        for row in matrix:
            if (
                isinstance(row, (str, bytes, bytearray))
                or not isinstance(row, Sequence)
                or not row
            ):
                raise ValueError("direction_matrix must be rectangular")
            values = [float(value) for value in row]
            if any(not math.isfinite(value) for value in values):
                raise ValueError("direction_matrix must contain finite values")
            rows.append(values)
        width = len(rows[0])
        if any(len(row) != width for row in rows):
            raise ValueError("direction_matrix must be rectangular")
        stride = [0.0] * (width * 48)
        for free in range(width):
            for coordinate, row in enumerate(rows):
                stride[free * 48 + coordinate] = row[free]
        return _equation_key(len(rows), width, stride)

    def match_subspace(
        self, payload: Mapping[str, Any]
    ) -> dict[str, Any] | None:
        """Match a rendered fixed subspace to its first Source OPD domain."""

        sg = _integer(payload.get("space_group"), "space_group", minimum=1, maximum=230)
        factor = self._resolve_factor(sg, payload)
        direction_matrix = payload.get("direction_matrix")
        if not isinstance(direction_matrix, Sequence) or len(direction_matrix) != int(
            self.projection.little_record_by_gid(int(factor["gid"])).full_dim
        ):
            raise ValueError("direction_matrix dimension differs from its irrep")
        target = self._subspace_key(direction_matrix)
        matches: list[tuple[str, int, Mapping[str, Any]]] = []
        for option in self._opd_options(sg, factor):
            direction = str(option["label"])
            for item, matrix in self._domain_options_with_matrices(
                sg,
                factor,
                direction,
                int(option["domain_count"]),
            ):
                if self._subspace_key(matrix) == target:
                    matches.append((direction, int(item["number"]), item))
        directions = {direction for direction, _domain, _item in matches}
        if not matches:
            return None
        if len(directions) != 1:
            raise ValueError(f"rendered subspace matches multiple Source OPDs for {factor['label']}")
        direction = matches[0][0]
        domain, item = min(
            (domain, item)
            for label, domain, item in matches
            if label == direction
        )
        resolved = {
            "gid": int(factor["gid"]),
            "label": ("m" if factor["magnetic"] else "") + str(factor["label"]),
            "magnetic": bool(factor["magnetic"]),
            "k_parameters": [
                _fraction_text(value) for value in factor["parameters"]
            ],
            "opd": direction,
            "domain": domain,
        }
        if payload.get("include_domain_details") is True:
            resolved.update(
                {
                    "domain_display": str(item["display"]),
                    "parameter_count": int(item["parameter_count"]),
                    "ferroic_properties": list(item["ferroic_properties"]),
                    "domain_presentation": {
                        "display_opd": str(item["display_opd"]),
                        "subgroup": str(item["subgroup"]),
                        "basis": str(item["basis"]),
                        "origin": str(item["origin"]),
                    },
                }
            )
        return resolved

    def resolve_subspace(self, payload: Mapping[str, Any]) -> dict[str, Any]:
        """Require a rendered fixed subspace to have a Source OPD domain."""

        resolved = self.match_subspace(payload)
        if resolved is not None:
            return resolved
        sg = _integer(payload.get("space_group"), "space_group", minimum=1, maximum=230)
        factor = self._resolve_factor(sg, payload)
        raise ValueError(
            f"rendered subspace has no Source OPD domain for {factor['label']}"
        )

    def _direction_matrix(self, sg: int, factor: Mapping[str, Any], direction: str):
        if int(factor["old_id"]) > 0 and not factor["magnetic"]:
            return _exact_direction_matrix(self.projection, int(factor["gid"]), direction)[1]
        isotropy = self._dynamic_isotropy(sg, factor, direction)
        numeric_rows = isotropy.get("source_numeric_rows") or ()
        display_rows = isotropy.get("source_display_rows") or isotropy.get("source_matrix") or ()
        rows = numeric_rows or display_rows
        if not rows or any(not isinstance(row, Sequence) for row in rows):
            raise ValueError("dynamic OPD direction matrix is empty")
        width = len(rows[0])
        if width < 1 or any(len(row) != width for row in rows):
            raise ValueError("dynamic OPD direction matrix is ragged")
        def exact_entry(column: int, row: int):
            if numeric_rows:
                try:
                    return _rational_entry(complex(numeric_rows[column][row]))
                except ValueError:
                    pass
            return _source_display_entry(display_rows[column][row])
        return tuple(
            tuple(exact_entry(column, row) for column in range(len(rows)))
            for row in range(width)
        )

    def _domain_options(
        self,
        sg: int,
        factor: Mapping[str, Any],
        direction: str,
        count: int,
        *,
        parameter_offset: int = 0,
    ) -> list[dict[str, Any]]:
        return [
            item
            for item, _matrix in self._domain_options_with_matrices(
                sg,
                factor,
                direction,
                count,
                parameter_offset=parameter_offset,
            )
        ]

    def _domain_options_with_matrices(
        self,
        sg: int,
        factor: Mapping[str, Any],
        direction: str,
        count: int,
        *,
        parameter_offset: int = 0,
    ) -> list[tuple[dict[str, Any], Any]]:
        gid = int(factor["gid"])
        matrix = self._direction_matrix(sg, factor, direction)
        out = []
        for domain in range(1, count + 1):
            transformed = matrix
            record = None
            if domain > 1:
                if int(factor["old_id"]) > 0 and not factor["magnetic"]:
                    row_id = self.projection.isotropy_direction_row(gid, direction)
                    if row_id is None:
                        raise ValueError(f"Source OPD row not found for gid={gid}: {direction}")
                    record = domain_operation_record_from_isotropy_row(
                        self.projection,
                        self.source_data,
                        sg=sg,
                        row_id=row_id,
                        domain=domain,
                    )
                    action = _exact_record_matrix(self.source_data, gid, record)
                else:
                    isotropy = self._dynamic_isotropy(sg, factor, direction)
                    record = domain_operation_record_from_subgroup(
                        self.source_data,
                        sg=sg,
                        child_sg=int(
                            (isotropy.get("subgroup") or {}).get("ordinary_number")
                            or (isotropy.get("subgroup") or {})["number"]
                        ),
                        basis=tuple(int(value) for value in isotropy["source_basis_values"]),
                        origin=tuple(int(value) for value in isotropy["source_origin_values"]),
                        domain=domain,
                    )
                if record is None:
                    raise ValueError(f"Source domain {domain} has no operation record")
                if int(factor["old_id"]) <= 0:
                    source_values = tuple(int(value) for value in isotropy["source_kparam"])
                    action = _quadratic3_record_matrix(
                        self.projection, gid, record, source_values
                    )
                elif factor["magnetic"]:
                    action = _exact_record_matrix(self.source_data, gid, record)
                transformed = _matrix_product(action, matrix)
            display_opd = _display_direction_matrix(
                transformed,
                parameter_offset=parameter_offset,
            )
            presentation = self._domain_presentation(
                sg,
                factor,
                direction,
                record,
            )
            details = " ".join(
                part for part in (display_opd, presentation["subgroup"]) if part
            )
            if presentation["basis"]:
                details += f", basis={{{presentation['basis']}}}"
            if presentation["origin"]:
                details += f", origin={presentation['origin']}"
            out.append((
                {
                    "number": domain,
                    "display_opd": display_opd,
                    "parameter_count": len(matrix[0]),
                    **presentation,
                    "display": f"{direction}({domain}) {details}",
                },
                transformed,
            ))
        return out

    def _domain_presentation(
        self,
        sg: int,
        factor: Mapping[str, Any],
        direction: str,
        record: tuple[int, int, int, int, int] | None,
    ) -> dict[str, Any]:
        isotropy = self._dynamic_isotropy(sg, factor, direction)
        ferroic_properties = [
            str(value)
            for value in isotropy.get("ferroic_properties") or []
            if value
        ]
        subgroup = isotropy.get("subgroup") or {}
        subgroup_text = " ".join(
            part
            for part in (
                str(subgroup.get("display_label") or subgroup.get("number") or ""),
                str(subgroup.get("symbol") or ""),
            )
            if part
        )
        if record is None:
            return {
                "subgroup": subgroup_text,
                "basis": str(isotropy.get("basis_text") or ""),
                "origin": str(isotropy.get("origin") or ""),
                "ferroic_properties": ferroic_properties,
            }

        source_basis = tuple(int(value) for value in isotropy["source_basis_values"])
        source_origin = tuple(int(value) for value in isotropy["source_origin_values"])
        transform = _domain_transform_matrix(source_basis, source_origin)  # type: ignore[arg-type]
        operation = _pml_operation_matrix(self.source_data, int(sg), record)
        basis_transform = _mat4_mul(operation, transform)
        origin_transform = _mat4_mul(transform, operation)
        basis = tuple(
            int(basis_transform[row][column])
            for row in range(3)
            for column in range(3)
        )
        origin = _fraction_record_mod1(
            tuple(origin_transform[3][column] for column in range(3))
        )
        ordinary_subgroup = int(subgroup.get("ordinary_number") or subgroup["number"])
        displayed = source_tables().subgroup_change_setting_cinter(
            int(sg),
            ordinary_subgroup,
            basis,
            origin,
        )
        basis_rows = displayed["basis"]
        denominator = int(displayed["basis_denominator"])
        return {
            "subgroup": subgroup_text,
            "basis": _basis_text(basis_rows, denominator),
            "origin": _display_origin_text(
                int(sg),
                ordinary_subgroup,
                basis_rows,
                origin,
            ),
            "ferroic_properties": ferroic_properties,
        }

    def compute(self, payload: Mapping[str, Any]) -> dict[str, Any]:
        sg = _integer(payload.get("space_group"), "space_group", minimum=1, maximum=230)
        mode = payload.get("mode", "opd")
        if mode not in {"opd", "full", "gradient"}:
            raise ValueError("mode must be opd, full, or gradient")
        raw_factors = payload.get("factors")
        if not isinstance(raw_factors, list) or not 1 <= len(raw_factors) <= MAX_FACTORS:
            raise ValueError(f"factors must contain from 1 through {MAX_FACTORS} rows")
        factors = [self._resolve_factor(sg, row) for row in raw_factors]
        minimum = _integer(
            payload.get("minimum_degree", 2), "minimum_degree", minimum=1, maximum=MAX_DEGREE
        )
        maximum = _integer(
            payload.get("maximum_degree", minimum),
            "maximum_degree",
            minimum=minimum,
            maximum=MAX_DEGREE,
        )
        labels = tuple(factor["label"] for factor in factors)
        parameter_rows = tuple(factor["parameters"] for factor in factors)
        magnetic = tuple(bool(factor["magnetic"]) for factor in factors)
        has_parametric = any(factor["parameter_count"] for factor in factors)
        started = time.perf_counter()

        with self._lock:
            if mode == "gradient":
                if minimum != maximum:
                    raise ValueError("gradient mode requires one exact degree")
                gradient_order = _integer(
                    payload.get("gradient_order", 1),
                    "gradient_order",
                    minimum=1,
                    maximum=MAX_GRADIENT_ORDER,
                )
                if has_parametric:
                    matrices, variables = coupled_parametric_gradient_action(
                        sg,
                        labels,
                        parameter_rows,
                        gradient_order,
                        projection_source=self.projection,
                    )
                else:
                    matrices, variables = coupled_fixed_gradient_action(
                        sg,
                        labels,
                        gradient_order,
                        projection_source=self.projection,
                    )
                dimensions = fixed_irrep_dimensions(
                    sg, labels, source_data=self.source_data
                )
                matrices = _gradient_time_reversal_action(
                    matrices, variables, dimensions, magnetic
                )
                basis = gradient_invariant_basis(
                    matrices,
                    variables,
                    degree=minimum,
                    gradient_order=gradient_order,
                )
                action_count = len(matrices)
                result = {
                    "degrees": [
                        {
                            "degree": minimum,
                            "count": len(basis),
                            "invariants": [
                                gradient_polynomial_text(polynomial, variables)
                                for polynomial in basis
                            ],
                        }
                    ],
                    "variable_count": len(variables),
                    "action_count": action_count,
                    "gradient_order": gradient_order,
                }
            else:
                matrices, dimensions = self._ordinary_action(
                    sg, mode, factors, labels, parameter_rows, has_parametric
                )
                matrices = _time_reversal_action(matrices, dimensions, magnetic)
                parameter_names = self._assign_factor_variable_names(
                    factors, dimensions
                )
                bases = invariant_basis(
                    matrices,
                    minimum_degree=minimum,
                    maximum_degree=maximum,
                    block_dimensions=dimensions,
                )
                result = {
                    "degrees": [
                        {
                            "degree": degree,
                            "count": len(polynomials),
                            "invariants": [
                                _parameter_polynomial_text(value, parameter_names)
                                for value in polynomials
                            ],
                            "polynomials": [
                                {
                                    "terms": [
                                        {
                                            "coefficient": coefficient,
                                            "exponents": list(exponents),
                                        }
                                        for coefficient, exponents in polynomial_terms(
                                            value
                                        )
                                    ]
                                }
                                for value in polynomials
                            ],
                        }
                        for degree, polynomials in bases.items()
                    ],
                    "variable_count": sum(dimensions),
                    "action_count": len(matrices),
                }

        if mode == "gradient":
            parameter_names = self._assign_factor_variable_names(factors, dimensions)

        elapsed = round((time.perf_counter() - started) * 1000, 3)
        return {
            "space_group": sg,
            "symbol": str(self.projection.space["space_label_bc"][sg - 1]).strip(),
            "mode": mode,
            "minimum_degree": minimum,
            "maximum_degree": maximum,
            "factors": [self._public_factor(factor) for factor in factors],
            "variables": list(parameter_names),
            "elapsed_ms": elapsed,
            **result,
        }

    def _ordinary_action(
        self,
        sg: int,
        mode: str,
        factors: Sequence[dict[str, Any]],
        labels: tuple[str, ...],
        parameter_rows: tuple[tuple[Fraction, ...], ...],
        has_parametric: bool,
    ):
        if mode == "full":
            if has_parametric:
                matrices = coupled_parametric_irrep_matrices(
                    sg, labels, parameter_rows, projection_source=self.projection
                )
            else:
                matrices = coupled_fixed_irrep_matrices(
                    sg, labels, source_data=self.source_data
                )
            return matrices, fixed_irrep_dimensions(sg, labels, source_data=self.source_data)

        directions = tuple(self._direction(factor) for factor in factors)
        domains = tuple(self._domain(factor) for factor in factors)
        if any(factor["magnetic"] for factor in factors):
            for factor, direction, domain in zip(
                factors, directions, domains, strict=True
            ):
                options = {
                    item["label"]: item["domain_count"]
                    for item in self.opds({
                        "space_group": sg,
                        "gid": factor["gid"],
                        "k_parameters": [
                            _fraction_text(value) for value in factor["parameters"]
                        ],
                        "magnetic": factor["magnetic"],
                    })["opds"]
                }
                if direction not in options:
                    raise ValueError(
                        f"unknown OPD {direction!r} for {factor['label']}"
                    )
                if domain > options[direction]:
                    raise ValueError(
                        f"domain {domain} exceeds the upper bound "
                        f"{options[direction]} for {direction}"
                    )
            if has_parametric:
                full_action = coupled_parametric_irrep_matrices(
                    sg, labels, parameter_rows, projection_source=self.projection
                )
            else:
                full_action = coupled_fixed_irrep_matrices(
                    sg, labels, source_data=self.source_data
                )
            direction_blocks = [
                self._domain_direction_matrix(sg, factor, direction, domain)
                for factor, direction, domain in zip(
                    factors, directions, domains, strict=True
                )
            ]
            return _project_action(full_action, direction_blocks)
        if has_parametric:
            return restricted_parametric_irrep_action(
                sg,
                labels,
                parameter_rows,
                directions,
                domains,
                projection_source=self.projection,
                selected_direction_matrices=tuple(
                    factor.get("direction_matrix") for factor in factors
                ),
            )
        return restricted_fixed_irrep_action(
            sg,
            labels,
            directions,
            domains,
            projection_source=self.projection,
        )

    def _resolve_factor(self, sg: int, raw: object) -> dict[str, Any]:
        if not isinstance(raw, Mapping):
            raise ValueError("each factor must be an object")
        gid = _integer(raw.get("gid"), "gid", minimum=1, maximum=100_000)
        records = {item[0]: item for item in self.projection.list_irreps(sg)}
        if gid not in records:
            raise ValueError(f"gid {gid} is not an irrep of SG{sg}")
        _, label, k_label, kind = records[gid]
        parameter_count = self.projection.k_parameter_dimension_by_gid(gid)
        raw_parameters = raw.get("k_parameters", [])
        if not isinstance(raw_parameters, list):
            raise ValueError("k_parameters must be an array")
        if len(raw_parameters) != parameter_count:
            raise ValueError(
                f"{k_label} requires {parameter_count} exact k parameter values"
            )
        parameters = tuple(
            _fraction(value, f"k_parameters[{index}]")
            for index, value in enumerate(raw_parameters)
        )
        if parameter_count:
            kslot = int(self.projection.little["little_irr_k"][gid - 1])
            ensure_source_k_not_special(
                source_tables(),
                sg=sg,
                kslot=kslot,
                label=k_label,
                parameters=parameters,
            )
        magnetic = raw.get("magnetic", False)
        if type(magnetic) is not bool:
            raise ValueError("magnetic must be a boolean")
        old_id = int(self.projection.little_record_by_gid(gid).old_id)
        if (
            magnetic
            and old_id > 0
            and not magnetic_source_data().irrep_has_magnetic_isotropy(old_id)
        ):
            raise ValueError(f"Source magnetic OPD rows are unavailable for {label}")
        parameter_offset = raw.get("parameter_offset")
        if parameter_offset is not None and (
            type(parameter_offset) is not int or parameter_offset < 0
        ):
            raise ValueError("parameter_offset must be a nonnegative integer")
        return {
            "gid": gid,
            "old_id": old_id,
            "label": label,
            "k_label": k_label,
            "kind": kind,
            "parameter_count": parameter_count,
            "parameters": parameters,
            "direction_matrix": self._selected_direction_matrix(
                raw.get("selected_direction_matrix"),
                full_dimension=int(
                    self.projection.little["little_irr_full_dim"][gid - 1]
                ),
            ),
            "magnetic": magnetic,
            "parameter_offset": parameter_offset,
            "direction": raw.get("direction"),
            "domain": raw.get("domain"),
        }

    @staticmethod
    def _selected_direction_matrix(
        raw: object,
        *,
        full_dimension: int,
    ) -> tuple[tuple[object, ...], ...] | None:
        if raw is None:
            return None
        if not isinstance(raw, (list, tuple)) or len(raw) != full_dimension:
            raise ValueError("direction_matrix must have one row per irrep coordinate")
        rows = []
        width: int | None = None
        for raw_row in raw:
            if not isinstance(raw_row, (list, tuple)) or not raw_row:
                raise ValueError("direction_matrix rows must be nonempty arrays")
            row = tuple(_rational_entry(complex(value)) for value in raw_row)
            if width is None:
                width = len(row)
            elif len(row) != width:
                raise ValueError("direction_matrix rows must have equal length")
            rows.append(row)
        return tuple(rows)

    @staticmethod
    def _assign_factor_variable_names(
        factors: Sequence[dict[str, Any]],
        dimensions: Sequence[int],
    ) -> tuple[str, ...]:
        cursor = 0
        used: set[int] = set()
        names: list[str] = []
        for factor, width in zip(factors, dimensions, strict=True):
            offset = factor.get("parameter_offset")
            if offset is None:
                offset = cursor
            indices = range(int(offset), int(offset) + int(width))
            if any(index in used for index in indices):
                raise ValueError("factor parameter ranges overlap")
            assigned = tuple(_parameter_name(index) for index in indices)
            factor["variable_names"] = assigned
            names.extend(assigned)
            used.update(indices)
            cursor = max(cursor, int(offset) + int(width))
        return tuple(names)

    def _domain_direction_matrix(
        self,
        sg: int,
        factor: Mapping[str, Any],
        direction: str,
        domain: int,
    ):
        matrix = self._direction_matrix(sg, factor, direction)
        if domain == 1:
            return matrix
        gid = int(factor["gid"])
        if int(factor["old_id"]) > 0 and not factor["magnetic"]:
            row_id = self.projection.isotropy_direction_row(gid, direction)
            if row_id is None:
                raise ValueError(f"Source OPD row not found for gid={gid}: {direction}")
            record = domain_operation_record_from_isotropy_row(
                self.projection, self.source_data, sg=sg, row_id=row_id, domain=domain
            )
        else:
            isotropy = self._dynamic_isotropy(sg, factor, direction)
            subgroup = isotropy.get("subgroup") or {}
            record = domain_operation_record_from_subgroup(
                self.source_data,
                sg=sg,
                child_sg=int(subgroup.get("ordinary_number") or subgroup["number"]),
                basis=tuple(int(value) for value in isotropy["source_basis_values"]),
                origin=tuple(int(value) for value in isotropy["source_origin_values"]),
                domain=domain,
            )
        if record is None:
            raise ValueError(f"Source domain {domain} has no operation record")
        if int(factor["old_id"]) <= 0:
            source_values = tuple(int(value) for value in isotropy["source_kparam"])
            action = _quadratic3_record_matrix(self.projection, gid, record, source_values)
        else:
            action = _exact_record_matrix(self.source_data, gid, record)
        return _matrix_product(action, matrix)

    def _domain_count(self, sg: int, factor: Mapping[str, Any], direction: str) -> int:
        gid = int(factor["gid"])
        if factor["parameter_count"]:
            source_values = self.projection.source_kparam_for_gid(gid, factor["parameters"])
            row = self.projection._dynamic_row_for_gid_kparam_direction(
                gid, source_values, direction
            )
            if row is None:
                raise ValueError(f"dynamic OPD row not found for gid={gid}: {direction}")
            return domain_count_from_subgroup(
                self.source_data,
                sg=sg,
                child_sg=int(row.subgroup_number),
                basis=tuple(int(value) for value in row.basis_values),
            )
        row_id = self.projection.isotropy_direction_row(gid, direction)
        if row_id is None:
            raise ValueError(f"Source OPD row not found for gid={gid}: {direction}")
        return domain_count_from_isotropy_row(self.projection, self.source_data, sg, row_id)

    @staticmethod
    def _direction(factor: Mapping[str, Any]) -> str:
        value = factor.get("direction")
        if not isinstance(value, str) or not value.strip():
            raise ValueError(f"OPD is required for {factor['label']}")
        return value.strip()

    @staticmethod
    def _domain(factor: Mapping[str, Any]) -> int:
        return _integer(factor.get("domain"), "domain", minimum=1, maximum=100_000)

    @staticmethod
    def _public_factor(factor: Mapping[str, Any]) -> dict[str, Any]:
        out = {
            "gid": factor["gid"],
            "k": factor["k_label"],
            "irrep": ("m" if factor["magnetic"] else "") + factor["label"],
            "magnetic": factor["magnetic"],
            "k_parameters": [_fraction_text(value) for value in factor["parameters"]],
            "parameters": list(factor.get("variable_names") or ()),
        }
        if factor.get("direction") is not None:
            out["opd"] = factor["direction"]
        if factor.get("domain") is not None:
            out["domain"] = factor["domain"]
        return out


@lru_cache(maxsize=1)
def invariant_service() -> InvariantService:
    return InvariantService()


__all__ = [
    "InvariantService",
    "MAX_DEGREE",
    "MAX_FACTORS",
    "MAX_GRADIENT_ORDER",
    "invariant_service",
]
