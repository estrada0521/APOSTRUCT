"""Exact linear-independence certificates for displayed mode bases.

The check is intentionally a necessary condition, not a replacement for the
Web comparator.  For each displacive or magnetic parent-atom block, and for
the six-component parent-cell strain block, the displayed mode definitions
must be linearly independent over the rationals.
One full-rank modular image proves independence.  A suspected deficiency is
accepted only after an exact rational nullspace and an independent modular
rank lower bound agree.
"""

from __future__ import annotations

from collections import defaultdict
from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field
from decimal import Decimal, InvalidOperation
from functools import reduce
import hashlib
import math
import re
import struct
from typing import Any

from sympy import Rational
from sympy.polys.domains import GF
from sympy.polys.matrices import DomainMatrix


MODE_BASIS_THEOREM = "complete_mode_basis_linear_independence.v1"
POSITION_SCALE = 100_000
VECTOR_SCALE = 10_000
_MODULAR_PRIMES = (
    1_000_003,
    1_000_033,
    1_000_037,
    1_000_039,
    1_000_081,
    1_000_099,
    1_000_117,
    1_000_121,
)
_NUMBER = re.compile(r"[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?")
_MODE_SITE = re.compile(r"\[([^:\]]+):([^:\]]+):(dsp|mag)]")
_DEFINITION_HEADER = re.compile(
    rf"(?P<label>.+?)\s+normfactor\s*=\s*(?P<norm>\?|{_NUMBER.pattern})\s*"
)
_VECTOR_HEADER = re.compile(
    r"atom\s+x\s+y\s+z\s+(?:dx\s+dy\s+dz|dmx\s+dmy\s+dmz)"
)
_STRAIN_HEADER = re.compile(r"e1\s+e2\s+e3\s+e4\s+e5\s+e6")
_SITE_MULTIPLICITY = re.compile(r"([1-9]\d*)")
_DEFINITION_HEADINGS = {
    "Displacive mode definitions": "dsp",
    "Magnetic mode definitions": "mag",
}
_SECTION_STOPS = {
    "Displacive mode amplitudes",
    "Magnetic mode amplitudes",
    "Parent-cell strain mode definitions",
    "Parent-cell strain mode amplitudes",
}

Position = tuple[int, int, int]
Coordinate = tuple[int, int, int, int]
Vector = tuple[int, int, int]


@dataclass(frozen=True)
class QuantizedModeRow:
    position: Position
    vector: Vector
    atom_marker: str | None = field(default=None, compare=False)


@dataclass(frozen=True)
class QuantizedModeDefinition:
    ordinal: int
    kind: str
    atom: str
    label: str
    normfactor_zero: bool
    rows: tuple[QuantizedModeRow, ...]


@dataclass(frozen=True)
class QuantizedStrainDefinition:
    ordinal: int
    kind: str
    label: str
    normfactor_zero: bool
    components: tuple[int, int, int, int, int, int]


@dataclass(frozen=True)
class QuantizedStructureOrbit:
    label: str
    site: str
    multiplicity: int
    representative: Position
    expanded_positions: tuple[Position, ...]


def _decimal(token: str, *, field: str) -> Decimal:
    if _NUMBER.fullmatch(token) is None:
        raise ValueError(f"{field} must be one finite decimal token")
    try:
        value = Decimal(token)
    except InvalidOperation as exc:
        raise ValueError(f"{field} must be one finite decimal token") from exc
    if not value.is_finite():
        raise ValueError(f"{field} must be finite")
    return value


def _scaled_decimal(token: str, scale: int, *, field: str) -> int:
    value = _decimal(token, field=field) * scale
    integral = value.to_integral_value()
    if value != integral:
        raise ValueError(f"{field} exceeds the declared display precision")
    return int(integral)


def _finite_number(value: object, *, field: str) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise ValueError(f"{field} must be a finite JSON number")
    result = float(value)
    if not math.isfinite(result):
        raise ValueError(f"{field} must be finite")
    return result


def _quantized_float(value: object, digits: int, *, field: str) -> int:
    number = _finite_number(value, field=field)
    return _scaled_decimal(
        format(number, f".{digits}f"),
        10**digits,
        field=field,
    )


def _label_identity(label: str, *, expected_kind: str) -> tuple[str, str]:
    matches = list(_MODE_SITE.finditer(label))
    if len(matches) != 1:
        raise ValueError(f"mode label must contain one atom/site/kind block: {label!r}")
    atom, _site, kind = matches[0].groups()
    if kind != expected_kind:
        raise ValueError(
            f"mode label kind {kind!r} does not match {expected_kind!r}: {label!r}"
        )
    if not atom:
        raise ValueError(f"mode label has an empty parent atom: {label!r}")
    return atom, kind


def _definition(
    *,
    ordinal: int,
    kind: str,
    label: str,
    normfactor: Decimal,
    rows: Sequence[QuantizedModeRow],
) -> QuantizedModeDefinition:
    if not rows:
        raise ValueError(f"mode definition has no vector rows: {label!r}")
    atom, checked_kind = _label_identity(label, expected_kind=kind)
    return QuantizedModeDefinition(
        ordinal=ordinal,
        kind=checked_kind,
        atom=atom,
        label=label,
        normfactor_zero=normfactor == 0,
        rows=tuple(rows),
    )


def definitions_from_web_text(
    text: str,
) -> tuple[QuantizedModeDefinition | QuantizedStrainDefinition, ...]:
    """Parse the visible Web mode-definition tables without silent row skips."""

    if not isinstance(text, str):
        raise ValueError("Web complete-mode payload must be text")
    definitions: list[QuantizedModeDefinition | QuantizedStrainDefinition] = []
    current_kind: str | None = None
    current_label: str | None = None
    current_normfactor: Decimal | None = None
    current_rows: list[QuantizedModeRow] = []
    vector_header_seen = False
    seen_headings: set[str] = set()

    def flush() -> None:
        nonlocal current_label, current_normfactor, current_rows
        if current_label is None:
            return
        assert current_kind is not None
        assert current_normfactor is not None
        definitions.append(
            _definition(
                ordinal=len(definitions),
                kind=current_kind,
                label=current_label,
                normfactor=current_normfactor,
                rows=current_rows,
            )
        )
        current_label = None
        current_normfactor = None
        current_rows = []

    for line_number, raw_line in enumerate(text.splitlines(), 1):
        line = raw_line.strip()
        if line in _DEFINITION_HEADINGS:
            flush()
            if line in seen_headings:
                raise ValueError(f"duplicate mode-definition heading on line {line_number}")
            seen_headings.add(line)
            current_kind = _DEFINITION_HEADINGS[line]
            vector_header_seen = False
            continue
        if line in _SECTION_STOPS:
            flush()
            current_kind = None
            vector_header_seen = False
            continue
        if current_kind is None or not line:
            continue
        if _VECTOR_HEADER.fullmatch(line) is not None:
            if vector_header_seen or current_label is not None:
                raise ValueError(f"misplaced mode-vector heading on line {line_number}")
            vector_header_seen = True
            continue

        header = _DEFINITION_HEADER.fullmatch(line)
        if header is not None:
            if not vector_header_seen:
                raise ValueError(
                    f"mode definition precedes its vector heading on line {line_number}"
                )
            flush()
            label = header.group("label").strip()
            if not label:
                raise ValueError(f"empty mode label on line {line_number}")
            norm = header.group("norm")
            if norm == "?":
                raise ValueError(f"unknown mode normfactor on line {line_number}")
            current_label = label
            current_normfactor = _decimal(norm, field="mode normfactor")
            continue

        if current_label is None:
            raise ValueError(f"mode row precedes its definition on line {line_number}")
        parts = line.split()
        numeric: list[str]
        atom_marker: str | None
        if len(parts) == 6 and all(_NUMBER.fullmatch(part) for part in parts):
            numeric = parts
            atom_marker = None
        elif (
            len(parts) == 7
            and _NUMBER.fullmatch(parts[0]) is None
            and all(_NUMBER.fullmatch(part) for part in parts[1:])
        ):
            numeric = parts[1:]
            atom_marker = parts[0]
        else:
            raise ValueError(f"malformed mode row on line {line_number}")
        position = tuple(
            _scaled_decimal(token, POSITION_SCALE, field="mode position")
            % POSITION_SCALE
            for token in numeric[:3]
        )
        vector = tuple(
            _scaled_decimal(token, VECTOR_SCALE, field="mode vector")
            for token in numeric[3:]
        )
        current_rows.append(
            QuantizedModeRow(
                position=position,  # type: ignore[arg-type]
                vector=vector,  # type: ignore[arg-type]
                atom_marker=atom_marker,
            )
        )
    flush()
    strain_body = text.partition("Parent-cell strain mode definitions")[2].partition(
        "Parent-cell strain mode amplitudes"
    )[0]
    strain_header_seen = False
    strain_label: str | None = None
    strain_norm: Decimal | None = None
    for line_number, raw_line in enumerate(strain_body.splitlines(), 1):
        line = raw_line.strip()
        if not line:
            continue
        if _STRAIN_HEADER.fullmatch(line) is not None:
            if strain_header_seen or strain_label is not None:
                raise ValueError("misplaced strain-component heading")
            strain_header_seen = True
            continue
        header = _DEFINITION_HEADER.fullmatch(line)
        if header is not None:
            if not strain_header_seen or strain_label is not None:
                raise ValueError("misplaced strain definition")
            strain_label = header.group("label").strip()
            norm = header.group("norm")
            if not strain_label or norm == "?":
                raise ValueError("invalid strain definition header")
            strain_norm = _decimal(norm, field="strain normfactor")
            continue
        if strain_label is None or strain_norm is None:
            raise ValueError(f"strain row precedes its definition on line {line_number}")
        parts = line.split()
        if len(parts) != 6 or any(_NUMBER.fullmatch(part) is None for part in parts):
            raise ValueError(f"malformed strain row on line {line_number}")
        components = tuple(
            _scaled_decimal(token, VECTOR_SCALE, field="strain component")
            for token in parts
        )
        definitions.append(
            QuantizedStrainDefinition(
                ordinal=len(definitions),
                kind="strain",
                label=strain_label,
                normfactor_zero=strain_norm == 0,
                components=components,  # type: ignore[arg-type]
            )
        )
        strain_label = None
        strain_norm = None
    if strain_label is not None:
        raise ValueError("strain definition has no component row")
    if not seen_headings and not strain_header_seen:
        raise ValueError("complete-mode text has no mode-definition heading")
    return tuple(definitions)


def structure_from_web_text(text: str) -> tuple[QuantizedStructureOrbit, ...]:
    """Read the displayed undistorted asymmetric-unit rows."""

    body = text.partition("Undistorted superstructure")[2].partition(
        "Distorted superstructure"
    )[0]
    rows: list[QuantizedStructureOrbit] = []
    for raw_line in body.splitlines():
        line = raw_line.strip()
        if not line or line.startswith("a=") or line.startswith("atom site"):
            continue
        parts = line.split()
        if len(parts) < 5 or any(_NUMBER.fullmatch(token) is None for token in parts[2:5]):
            continue
        match = _SITE_MULTIPLICITY.match(parts[1])
        if match is None:
            raise ValueError(f"undistorted site has no multiplicity: {parts[1]!r}")
        position = tuple(
            _scaled_decimal(token, POSITION_SCALE, field="undistorted position")
            % POSITION_SCALE
            for token in parts[2:5]
        )
        rows.append(
            QuantizedStructureOrbit(
                label=parts[0],
                site=parts[1],
                multiplicity=int(match.group(1)),
                representative=position,  # type: ignore[arg-type]
                expanded_positions=(),
            )
        )
    if not rows:
        raise ValueError("complete-mode text has no undistorted structure rows")
    return tuple(rows)


def _required_mapping(value: object, *, field: str) -> Mapping[str, Any]:
    if not isinstance(value, Mapping):
        raise ValueError(f"{field} must be an object")
    return value


def _required_list(value: object, *, field: str) -> list[Any]:
    if not isinstance(value, list):
        raise ValueError(f"{field} must be an array")
    return value


def _local_definition(
    raw: object,
    *,
    ordinal: int,
    kind: str,
) -> QuantizedModeDefinition:
    value = _required_mapping(raw, field=f"{kind} definition")
    label = value.get("label")
    if not isinstance(label, str) or not label:
        raise ValueError(f"{kind} definition label must be nonempty text")
    if "normfactor" not in value:
        raise ValueError(f"{kind} definition is missing normfactor")
    normfactor = _finite_number(value["normfactor"], field="mode normfactor")
    rows = _required_list(value.get("rows"), field=f"{kind} definition rows")
    quantized_rows: list[QuantizedModeRow] = []
    for row_index, raw_row in enumerate(rows):
        row = _required_mapping(raw_row, field=f"{kind} row {row_index}")
        xyz = _required_list(row.get("xyz"), field=f"{kind} row {row_index} xyz")
        dxyz = _required_list(row.get("dxyz"), field=f"{kind} row {row_index} dxyz")
        if len(xyz) != 3 or len(dxyz) != 3:
            raise ValueError(f"{kind} row {row_index} needs vector3 xyz and dxyz")
        atom_marker = row.get("atom")
        if atom_marker is not None and not isinstance(atom_marker, str):
            raise ValueError(f"{kind} row {row_index} atom must be text or null")
        if atom_marker == "":
            raise ValueError(f"{kind} row {row_index} atom must not be empty")
        position = tuple(
            _quantized_float(component, 5, field="mode position") % POSITION_SCALE
            for component in xyz
        )
        vector = tuple(
            _quantized_float(component, 4, field="mode vector")
            for component in dxyz
        )
        quantized_rows.append(
            QuantizedModeRow(
                position=position,  # type: ignore[arg-type]
                vector=vector,  # type: ignore[arg-type]
                atom_marker=atom_marker,
            )
        )
    return _definition(
        ordinal=ordinal,
        kind=kind,
        label=label,
        normfactor=Decimal(format(normfactor, ".5f")),
        rows=quantized_rows,
    )


def _local_strain_definition(
    raw: object, *, ordinal: int
) -> QuantizedStrainDefinition:
    value = _required_mapping(raw, field="strain definition")
    label = value.get("label")
    if not isinstance(label, str) or not label:
        raise ValueError("strain definition label must be nonempty text")
    if "normfactor" not in value:
        raise ValueError("strain definition is missing normfactor")
    normfactor = _finite_number(value["normfactor"], field="strain normfactor")
    components = _required_list(
        value.get("components"), field="strain definition components"
    )
    if len(components) != 6:
        raise ValueError("strain definition needs six components")
    return QuantizedStrainDefinition(
        ordinal=ordinal,
        kind="strain",
        label=label,
        normfactor_zero=normfactor == 0,
        components=tuple(
            _quantized_float(component, 4, field="strain component")
            for component in components
        ),  # type: ignore[arg-type]
    )


def _local_structure_orbit(raw: object) -> QuantizedStructureOrbit:
    value = _required_mapping(raw, field="undistorted atom orbit")
    label = value.get("label")
    site = value.get("site")
    multiplicity = value.get("multiplicity")
    if not isinstance(label, str) or not label or not isinstance(site, str) or not site:
        raise ValueError("undistorted atom orbit needs label and site")
    if isinstance(multiplicity, bool) or not isinstance(multiplicity, int) or multiplicity <= 0:
        raise ValueError("undistorted atom orbit needs positive multiplicity")
    xyz = _required_list(value.get("xyz"), field="undistorted representative")
    if len(xyz) != 3:
        raise ValueError("undistorted representative needs vector3 xyz")
    representative = tuple(
        _quantized_float(component, 5, field="undistorted position") % POSITION_SCALE
        for component in xyz
    )
    raw_expanded = value.get("_presentation_orbit_points")
    expanded: list[Position] = []
    if raw_expanded is not None:
        for index, raw_position in enumerate(
            _required_list(raw_expanded, field="undistorted orbit points")
        ):
            position = _required_list(
                raw_position, field=f"undistorted orbit point {index}"
            )
            if len(position) != 3:
                raise ValueError("undistorted orbit point needs vector3 xyz")
            expanded.append(
                tuple(
                    _quantized_float(component, 5, field="undistorted position")
                    % POSITION_SCALE
                    for component in position
                )  # type: ignore[arg-type]
            )
    return QuantizedStructureOrbit(
        label=label,
        site=site,
        multiplicity=multiplicity,
        representative=representative,  # type: ignore[arg-type]
        expanded_positions=tuple(expanded),
    )


def definitions_from_local_payload(
    payload: Mapping[str, Any],
) -> tuple[QuantizedModeDefinition | QuantizedStrainDefinition, ...]:
    """Read the Local native mode definitions at the public display precision."""

    root = _required_mapping(payload, field="Local payload")
    preview = _required_mapping(root.get("preview"), field="Local preview")
    selected = _required_mapping(preview.get("selected"), field="Local selected state")
    details = _required_mapping(
        selected.get("mode_details"), field="Local complete-mode details"
    )
    definitions: list[QuantizedModeDefinition | QuantizedStrainDefinition] = []
    for kind, key in (
        ("dsp", "displacive_definitions"),
        ("mag", "magnetic_definitions"),
    ):
        rows = _required_list(details.get(key), field=f"Local {key}")
        for raw in rows:
            definitions.append(
                _local_definition(
                    raw,
                    ordinal=len(definitions),
                    kind=kind,
                )
            )
    strain_rows = _required_list(
        details.get("strain_definitions"), field="Local strain_definitions"
    )
    for raw in strain_rows:
        definitions.append(
            _local_strain_definition(raw, ordinal=len(definitions))
        )
    return tuple(definitions)


def structure_from_local_payload(
    payload: Mapping[str, Any],
) -> tuple[QuantizedStructureOrbit, ...]:
    root = _required_mapping(payload, field="Local payload")
    preview = _required_mapping(root.get("preview"), field="Local preview")
    selected = _required_mapping(preview.get("selected"), field="Local selected state")
    details = _required_mapping(
        selected.get("mode_details"), field="Local complete-mode details"
    )
    rows = _required_list(
        details.get("undistorted_atoms"), field="Local undistorted_atoms"
    )
    result = tuple(_local_structure_orbit(row) for row in rows)
    if not result:
        raise ValueError("Local output has no undistorted atom orbits")
    return result


def _columns_for_group(
    definitions: Sequence[QuantizedModeDefinition],
) -> tuple[list[dict[Coordinate, int]], tuple[Position, ...], str]:
    columns: list[dict[Coordinate, int]] = []
    positions: set[Position] = set()
    digest = hashlib.sha256()
    digest.update(b"isodistort.mode-basis.integer-matrix.v1\0")
    digest.update(struct.pack("<I", len(definitions)))
    for definition in definitions:
        by_position: dict[Position, Vector] = {}
        for row in definition.rows:
            if row.position in by_position:
                raise ValueError(
                    f"mode definition repeats a periodic position: {definition.label!r}"
                )
            by_position[row.position] = row.vector
            positions.add(row.position)
        digest.update(struct.pack("<I?", len(by_position), definition.normfactor_zero))
        column: dict[Coordinate, int] = {}
        for position, vector in sorted(by_position.items()):
            digest.update(struct.pack("<6q", *position, *vector))
            if definition.normfactor_zero:
                continue
            for axis, component in enumerate(vector):
                if component:
                    column[(*position, axis)] = component
        columns.append(column)
    return columns, tuple(sorted(positions)), digest.hexdigest()


def _rank_mod_prime(columns: Sequence[Mapping[Coordinate, int]], prime: int) -> int:
    coordinates = sorted({coordinate for column in columns for coordinate in column})
    row_index = {coordinate: index for index, coordinate in enumerate(coordinates)}
    sparse_rows: dict[int, dict[int, int]] = {}
    for column_index, column in enumerate(columns):
        for coordinate, value in column.items():
            reduced = int(value) % prime
            if reduced:
                sparse_rows.setdefault(row_index[coordinate], {})[
                    column_index
                ] = reduced
    matrix = DomainMatrix.from_dict_sympy(
        len(coordinates), len(columns), sparse_rows
    ).convert_to(GF(prime))
    return matrix.rank()


def _primitive_integer_relation(values: Sequence[object]) -> tuple[int, ...]:
    rationals = [Rational(value) for value in values]
    denominator = math.lcm(*(int(value.q) for value in rationals))
    integers = [int(value * denominator) for value in rationals]
    divisor = reduce(math.gcd, (abs(value) for value in integers if value), 0)
    if divisor == 0:
        raise ValueError("exact nullspace returned a zero relation")
    integers = [value // divisor for value in integers]
    first = next(value for value in integers if value)
    if first < 0:
        integers = [-value for value in integers]
    return tuple(integers)


def _verify_relation(
    columns: Sequence[Mapping[Coordinate, int]], coefficients: Sequence[int]
) -> None:
    total: dict[Coordinate, int] = defaultdict(int)
    for coefficient, column in zip(coefficients, columns, strict=True):
        if not coefficient:
            continue
        for coordinate, value in column.items():
            total[coordinate] += coefficient * value
    if any(total.values()):
        raise ValueError("exact nullspace relation failed integer substitution")


def _exact_null_relations(
    columns: Sequence[Mapping[Coordinate, int]],
    coordinates: Sequence[Coordinate],
) -> tuple[tuple[int, ...], ...]:
    row_index = {coordinate: index for index, coordinate in enumerate(coordinates)}
    sparse_rows: dict[int, dict[int, int]] = {}
    for column_index, column in enumerate(columns):
        for coordinate, value in column.items():
            index = row_index[coordinate]
            sparse_rows.setdefault(index, {})[column_index] = value
    matrix = DomainMatrix.from_dict_sympy(
        len(coordinates), len(columns), sparse_rows
    ).to_field()
    nullspace = matrix.nullspace().to_Matrix()
    relations = tuple(
        _primitive_integer_relation(tuple(nullspace[row, column] for column in range(nullspace.cols)))
        for row in range(nullspace.rows)
    )
    for relation in relations:
        _verify_relation(columns, relation)
    return relations


def _relation_independence_prime(
    relations: Sequence[Sequence[int]],
) -> int | None:
    relation_columns = [
        {
            (definition_index, 0, 0, 0): coefficient
            for definition_index, coefficient in enumerate(relation)
            if coefficient
        }
        for relation in relations
    ]
    for prime in _MODULAR_PRIMES:
        if _rank_mod_prime(relation_columns, prime) == len(relations):
            return prime
    return None


def _group_certificate(
    kind: str,
    atom: str,
    definitions: Sequence[QuantizedModeDefinition],
) -> dict[str, Any]:
    columns, positions, matrix_sha256 = _columns_for_group(definitions)
    column_count = len(columns)
    modular_ranks: dict[str, int] = {}
    for prime in _MODULAR_PRIMES:
        rank = _rank_mod_prime(columns, prime)
        modular_ranks[str(prime)] = rank
        if rank == column_count:
            return {
                "kind": kind,
                "atom": atom,
                "definition_count": column_count,
                "matrix_rows": 3 * len(positions),
                "matrix_sha256": matrix_sha256,
                "exact_rank_over_q": column_count,
                "deficiency": 0,
                "proof": "full_column_rank_mod_prime",
                "modular_ranks": modular_ranks,
            }

    coordinates = tuple(
        (*position, axis) for position in positions for axis in range(3)
    )
    relations = _exact_null_relations(columns, coordinates)
    exact_rank = column_count - len(relations)
    if max(modular_ranks.values(), default=0) != exact_rank:
        raise ValueError(
            "fixed modular primes did not certify the exact rational rank lower bound"
        )
    relation_prime = _relation_independence_prime(relations)
    if relations and relation_prime is None:
        raise ValueError(
            "fixed modular primes did not certify null-relation independence"
        )
    relation_rows: list[dict[str, Any]] = []
    for relation in relations:
        terms = [
            {
                "definition_index": definition.ordinal,
                "coefficient": coefficient,
                "label": definition.label,
            }
            for definition, coefficient in zip(definitions, relation, strict=True)
            if coefficient
        ]
        relation_rows.append({"terms": terms})
    return {
        "kind": kind,
        "atom": atom,
        "definition_count": column_count,
        "matrix_rows": 3 * len(positions),
        "matrix_sha256": matrix_sha256,
        "exact_rank_over_q": exact_rank,
        "deficiency": len(relations),
        "proof": "exact_nullspace_plus_modular_lower_bound",
        "modular_ranks": modular_ranks,
        "null_relation_independence_prime": relation_prime,
        "exact_integer_null_relations": relation_rows,
    }


def _strain_group_certificate(
    definitions: Sequence[QuantizedStrainDefinition],
) -> dict[str, Any]:
    coordinates: tuple[Coordinate, ...] = tuple(
        (axis, 0, 0, 0) for axis in range(6)
    )
    digest = hashlib.sha256()
    digest.update(b"isodistort.strain-basis.integer-matrix.v1\0")
    digest.update(struct.pack("<I", len(definitions)))
    columns: list[dict[Coordinate, int]] = []
    for definition in definitions:
        digest.update(
            struct.pack(
                "<6q?", *definition.components, definition.normfactor_zero
            )
        )
        columns.append(
            {
                coordinates[axis]: component
                for axis, component in enumerate(definition.components)
                if component and not definition.normfactor_zero
            }
        )
    column_count = len(columns)
    modular_ranks: dict[str, int] = {}
    for prime in _MODULAR_PRIMES:
        rank = _rank_mod_prime(columns, prime)
        modular_ranks[str(prime)] = rank
        if rank == column_count:
            return {
                "kind": "strain",
                "atom": "parent_cell",
                "definition_count": column_count,
                "matrix_rows": 6,
                "matrix_sha256": digest.hexdigest(),
                "exact_rank_over_q": column_count,
                "deficiency": 0,
                "proof": "full_column_rank_mod_prime",
                "modular_ranks": modular_ranks,
            }
    relations = _exact_null_relations(columns, coordinates)
    exact_rank = column_count - len(relations)
    if max(modular_ranks.values(), default=0) != exact_rank:
        raise ValueError(
            "fixed modular primes did not certify the strain rank lower bound"
        )
    relation_prime = _relation_independence_prime(relations)
    if relations and relation_prime is None:
        raise ValueError(
            "fixed modular primes did not certify strain-relation independence"
        )
    relation_rows = []
    for relation in relations:
        relation_rows.append(
            {
                "terms": [
                    {
                        "definition_index": definition.ordinal,
                        "coefficient": coefficient,
                        "label": definition.label,
                    }
                    for definition, coefficient in zip(
                        definitions, relation, strict=True
                    )
                    if coefficient
                ]
            }
        )
    return {
        "kind": "strain",
        "atom": "parent_cell",
        "definition_count": column_count,
        "matrix_rows": 6,
        "matrix_sha256": digest.hexdigest(),
        "exact_rank_over_q": exact_rank,
        "deficiency": len(relations),
        "proof": "exact_nullspace_plus_modular_lower_bound",
        "modular_ranks": modular_ranks,
        "null_relation_independence_prime": relation_prime,
        "exact_integer_null_relations": relation_rows,
    }


def assess_mode_basis(
    definitions: Sequence[QuantizedModeDefinition | QuantizedStrainDefinition],
) -> dict[str, Any]:
    """Return exact certificates for the mode-basis independence theorem."""

    if not isinstance(definitions, Sequence):
        raise ValueError("mode definitions must be an ordered sequence")
    grouped: dict[tuple[str, str], list[QuantizedModeDefinition]] = {}
    strain: list[QuantizedStrainDefinition] = []
    row_count = 0
    for expected_ordinal, definition in enumerate(definitions):
        if not isinstance(
            definition, (QuantizedModeDefinition, QuantizedStrainDefinition)
        ):
            raise ValueError("mode definitions use an unsupported definition type")
        if definition.ordinal != expected_ordinal:
            raise ValueError("mode definition ordinals must be contiguous and ordered")
        if isinstance(definition, QuantizedStrainDefinition):
            if definition.kind != "strain":
                raise ValueError("strain definition has an invalid kind")
            strain.append(definition)
            row_count += 1
        else:
            grouped.setdefault((definition.kind, definition.atom), []).append(definition)
            row_count += len(definition.rows)

    for (kind, atom), group_definitions in grouped.items():
        zero_definition = next(
            (
                definition
                for definition in group_definitions
                if definition.normfactor_zero
                or not any(component for row in definition.rows for component in row.vector)
            ),
            None,
        )
        if zero_definition is not None:
            _, positions, matrix_sha256 = _columns_for_group(group_definitions)
            return {
                "theorem": MODE_BASIS_THEOREM,
                "status": "refuted",
                "proof_scope": "exact_zero_column",
                "definition_count": len(definitions),
                "definition_row_count": row_count,
                "strain_definition_count": len(strain),
                "group_count": len(grouped) + bool(strain),
                "refuted_group_count": 1,
                "groups": [],
                "zero_column": {
                    "kind": kind,
                    "atom": atom,
                    "group_definition_count": len(group_definitions),
                    "matrix_rows": 3 * len(positions),
                    "matrix_sha256": matrix_sha256,
                    "definition_index": zero_definition.ordinal,
                    "label": zero_definition.label,
                },
            }

    certificates = [
        _group_certificate(kind, atom, rows)
        for (kind, atom), rows in grouped.items()
    ]
    if strain:
        certificates.append(_strain_group_certificate(strain))
    refuted = [row for row in certificates if int(row["deficiency"]) > 0]
    return {
        "theorem": MODE_BASIS_THEOREM,
        "status": "refuted" if refuted else "satisfied",
        "definition_count": len(definitions),
        "definition_row_count": row_count,
        "strain_definition_count": len(strain),
        "group_count": len(certificates),
        "refuted_group_count": len(refuted),
        "groups": certificates,
    }


__all__ = [
    "MODE_BASIS_THEOREM",
    "QuantizedModeDefinition",
    "QuantizedModeRow",
    "QuantizedStrainDefinition",
    "QuantizedStructureOrbit",
    "assess_mode_basis",
    "definitions_from_local_payload",
    "definitions_from_web_text",
    "structure_from_local_payload",
    "structure_from_web_text",
]
