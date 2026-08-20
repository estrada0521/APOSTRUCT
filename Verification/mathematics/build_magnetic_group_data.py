#!/usr/bin/env python3
"""Rebuild and audit the exact magnetic child-group authority artifact."""

from __future__ import annotations

import argparse
from collections import deque
from collections.abc import Iterable, Mapping, Sequence
import csv
from fractions import Fraction
import gzip
import hashlib
import io
import json
from pathlib import Path
import re
import shlex
from typing import Any, TextIO

import spglib


ROOT = Path(__file__).resolve().parents[2]
DATA_DIR = Path(__file__).with_name("data")
UPSTREAM_DIR = DATA_DIR / "upstream/spglib-2.7.0"
HALL_PATH = UPSTREAM_DIR / "magnetic_hall_symbols.yaml"
NUMBERS_PATH = UPSTREAM_DIR / "msg_numbers.csv"
LICENSE_PATH = DATA_DIR / "SPGLIB-LICENSE.txt"
OUTPUT_PATH = DATA_DIR / "magnetic_group_generators.json.gz"
SOURCE_SPACE_PATH = ROOT / "Source/data_space.txt"

SPGLIB_VERSION = "2.7.0"
SPGLIB_COMMIT = "12355c77fb7c505a55f52cae36341d73b781a065"
HALL_SHA256 = "6e3d7e2de1540c20a3ef13ad3554fac3e68e399c10c435e2feb86c933229d088"
NUMBERS_SHA256 = "7f1c92556be68f20999971a6c18631ba460c367f3b2d0e2f892643e811d0f539"
LICENSE_SHA256 = "9c4e602ce15bf48206dad105b666ddaeedba747eeda6465807106bc11f4f587f"
OFFICIAL_ISO_MAG_SHA256 = (
    "2b11217ae10687b0836d8151846db2fad57b0d77211da90bb138ccd123a7b1fc"
)
SPGLIB_SETTING_VERSION = "2.7.0"
SOURCE_SPACE_SHA256 = "47d8de79c643ad2f6c8e57424e179cd5702ddad669f02032296c0e3f8c2e594a"

# These prime placements are the unique Hall variants whose exact closures
# equal the pinned ISO-MAG operation table.  ``audit_official_iso_mag`` checks
# all 1,651 groups rather than treating this list as sufficient evidence.
OFFICIAL_PARITY_OVERRIDES = {
    "16.3": "P 2 2'",
    "27.80": "P 2' -2c'",
    "32.137": "P 2' -2ab'",
    "34.158": "P 2' -2n'",
    "35.167": "C 2' -2'",
    "37.182": "C 2' -2c'",
    "42.221": "F 2' -2'",
    "43.226": "F 2' -2d'",
    "44.231": "I 2' -2'",
}

Matrix = tuple[tuple[int, int, int], tuple[int, int, int], tuple[int, int, int]]
Vector = tuple[Fraction, Fraction, Fraction]
Operation = tuple[Matrix, Vector, int]

ZERO = Fraction(0)
HALF = Fraction(1, 2)
QUARTER = Fraction(1, 4)
IDENTITY_MATRIX: Matrix = ((1, 0, 0), (0, 1, 0), (0, 0, 1))
ZERO_VECTOR: Vector = (ZERO, ZERO, ZERO)
IDENTITY: Operation = (IDENTITY_MATRIX, ZERO_VECTOR, 1)

LATTICE_SYMBOLS: dict[str, tuple[Vector, ...]] = {
    "P": (ZERO_VECTOR,),
    "A": (ZERO_VECTOR, (ZERO, HALF, HALF)),
    "B": (ZERO_VECTOR, (HALF, ZERO, HALF)),
    "C": (ZERO_VECTOR, (HALF, HALF, ZERO)),
    "I": (ZERO_VECTOR, (HALF, HALF, HALF)),
    "R": (
        ZERO_VECTOR,
        (Fraction(2, 3), Fraction(1, 3), Fraction(1, 3)),
        (Fraction(1, 3), Fraction(2, 3), Fraction(2, 3)),
    ),
    "H": (
        ZERO_VECTOR,
        (Fraction(2, 3), Fraction(1, 3), ZERO),
        (Fraction(1, 3), Fraction(2, 3), ZERO),
    ),
    "F": (
        ZERO_VECTOR,
        (ZERO, HALF, HALF),
        (HALF, ZERO, HALF),
        (HALF, HALF, ZERO),
    ),
    "p": (ZERO_VECTOR,),
    "c": (ZERO_VECTOR, (HALF, HALF, ZERO)),
}

ROTATION_MATRICES: dict[str, Matrix] = {
    "1x": IDENTITY_MATRIX,
    "1y": IDENTITY_MATRIX,
    "1z": IDENTITY_MATRIX,
    "2x": ((1, 0, 0), (0, -1, 0), (0, 0, -1)),
    "2y": ((-1, 0, 0), (0, 1, 0), (0, 0, -1)),
    "2z": ((-1, 0, 0), (0, -1, 0), (0, 0, 1)),
    "3x": ((1, 0, 0), (0, 0, -1), (0, 1, -1)),
    "3y": ((-1, 0, 1), (0, 1, 0), (-1, 0, 0)),
    "3z": ((0, -1, 0), (1, -1, 0), (0, 0, 1)),
    "4x": ((1, 0, 0), (0, 0, -1), (0, 1, 0)),
    "4y": ((0, 0, 1), (0, 1, 0), (-1, 0, 0)),
    "4z": ((0, -1, 0), (1, 0, 0), (0, 0, 1)),
    "6x": ((1, 0, 0), (0, 1, -1), (0, 1, 0)),
    "6y": ((0, 0, 1), (0, 1, 0), (-1, 0, 1)),
    "6z": ((1, -1, 0), (1, 0, 0), (0, 0, 1)),
    "2px": ((-1, 0, 0), (0, 0, -1), (0, -1, 0)),
    "2ppx": ((-1, 0, 0), (0, 0, 1), (0, 1, 0)),
    "2py": ((0, 0, -1), (0, -1, 0), (-1, 0, 0)),
    "2ppy": ((0, 0, 1), (0, -1, 0), (1, 0, 0)),
    "2pz": ((0, -1, 0), (-1, 0, 0), (0, 0, -1)),
    "2ppz": ((0, 1, 0), (1, 0, 0), (0, 0, -1)),
    "3*": ((0, 0, 1), (1, 0, 0), (0, 1, 0)),
}

TRANSLATIONS: dict[str, Vector] = {
    "a": (HALF, ZERO, ZERO),
    "b": (ZERO, HALF, ZERO),
    "c": (ZERO, ZERO, HALF),
    "n": (HALF, HALF, HALF),
    "u": (QUARTER, ZERO, ZERO),
    "v": (ZERO, QUARTER, ZERO),
    "w": (ZERO, ZERO, QUARTER),
    "d": (QUARTER, QUARTER, QUARTER),
}


def sha256_bytes(value: bytes) -> str:
    return hashlib.sha256(value).hexdigest()


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def require_source(path: Path, expected_sha256: str) -> None:
    if not path.is_file() or sha256_file(path) != expected_sha256:
        raise ValueError(f"pinned authority source differs: {path}")


def matrix_vector(matrix: Matrix, vector: Vector) -> Vector:
    return tuple(
        sum(Fraction(matrix[row][axis]) * vector[axis] for axis in range(3))
        for row in range(3)
    )  # type: ignore[return-value]


def matrix_multiply(left: Matrix, right: Matrix) -> Matrix:
    return tuple(
        tuple(
            sum(left[row][axis] * right[axis][column] for axis in range(3))
            for column in range(3)
        )
        for row in range(3)
    )  # type: ignore[return-value]


def matrix_determinant(matrix: Matrix) -> int:
    first, second, third = matrix
    return (
        first[0] * (second[1] * third[2] - second[2] * third[1])
        - first[1] * (second[0] * third[2] - second[2] * third[0])
        + first[2] * (second[0] * third[1] - second[1] * third[0])
    )


def vector_add(left: Vector, right: Vector) -> Vector:
    return tuple(left[axis] + right[axis] for axis in range(3))  # type: ignore[return-value]


def vector_subtract(left: Vector, right: Vector) -> Vector:
    return tuple(left[axis] - right[axis] for axis in range(3))  # type: ignore[return-value]


def normalize_vector(vector: Iterable[Fraction | int]) -> Vector:
    values = tuple(Fraction(value) % 1 for value in vector)
    if len(values) != 3:
        raise ValueError("operation translation must contain three components")
    return values  # type: ignore[return-value]


def negate_matrix(matrix: Matrix) -> Matrix:
    return tuple(tuple(-matrix[row][column] for column in range(3)) for row in range(3))  # type: ignore[return-value]


def compose(left: Operation, right: Operation) -> Operation:
    left_rotation, left_translation, left_time = left
    right_rotation, right_translation, right_time = right
    rotated = matrix_vector(left_rotation, right_translation)
    return (
        matrix_multiply(left_rotation, right_rotation),
        normalize_vector(left_translation[axis] + rotated[axis] for axis in range(3)),
        left_time * right_time,
    )


def closure(
    generators: Iterable[Operation], *, limit: int = 4096
) -> frozenset[Operation]:
    generators = tuple(generators)
    seen = {IDENTITY}
    queue = deque([IDENTITY])
    while queue:
        current = queue.popleft()
        for generator in generators:
            product = compose(current, generator)
            if product in seen:
                continue
            seen.add(product)
            queue.append(product)
            if len(seen) > limit:
                raise ValueError(f"magnetic Hall closure exceeds {limit} operations")
    return frozenset(seen)


def _parse_origin(tokens: Sequence[str]) -> Vector:
    if len(tokens) != 3 or not tokens[0].startswith("(") or not tokens[2].endswith(")"):
        raise ValueError(f"invalid Hall origin shift: {tokens!r}")
    return (
        Fraction(int(tokens[0][1:]), 12),
        Fraction(int(tokens[1]), 12),
        Fraction(int(tokens[2][:-1]), 12),
    )


def _parse_rotation_token(
    token: str,
    count: int,
    previous_nfold: str | None,
    previous_axis: str | None,
) -> tuple[Matrix, Vector, str, str, int]:
    position = 0
    improper = token[position] == "-"
    if improper:
        position += 1
    nfold = token[position]
    if nfold not in "12346":
        raise ValueError(f"invalid Hall rotation order in {token!r}")
    position += 1

    axis = ""
    if position < len(token):
        if token[position] == "^":
            axis = "p"
            position += 1
        elif token[position] in ('"', "="):
            axis = "pp"
            position += 1
    if position < len(token) and token[position] in "xyz*":
        axis = token[position]
        position += 1
    if axis in {"p", "pp"} and previous_axis in {"x", "y", "z"}:
        axis += previous_axis

    if nfold == "1":
        axis = "z"
    if axis in {"", "p", "pp"}:
        if count == 0:
            axis += "z"
        elif count == 1 and previous_nfold in {"2", "4"}:
            axis += "x"
        elif count == 1 and previous_nfold in {"3", "6"}:
            axis = "pz"
        elif count == 2 and nfold == "3":
            axis += "*"
        else:
            raise ValueError(f"unreachable default Hall axis in {token!r}")

    key = nfold + axis
    try:
        rotation = ROTATION_MATRICES[key]
    except KeyError as exc:
        raise ValueError(f"unsupported Hall rotation key {key!r}") from exc
    if improper:
        rotation = negate_matrix(rotation)

    translation = ZERO_VECTOR
    while position < len(token):
        component = token[position]
        if component in "123456":
            translation = (ZERO, ZERO, Fraction(int(component), int(nfold)))
        elif component in TRANSLATIONS:
            translation = vector_add(translation, TRANSLATIONS[component])
        else:
            break
        position += 1

    time_reversal = 1
    if position < len(token) and token[position] == "'":
        time_reversal = -1
        position += 1
    if position != len(token):
        raise ValueError(f"unparsed Hall token suffix in {token!r}")
    return rotation, translation, nfold, axis, time_reversal


def hall_generators(hall_symbol: str) -> tuple[Operation, ...]:
    tokens = hall_symbol.split()
    if not tokens:
        raise ValueError("empty Hall symbol")
    lattice_token = tokens[0]
    inversion = lattice_token.startswith("-")
    lattice = lattice_token[1:] if inversion else lattice_token
    if lattice not in LATTICE_SYMBOLS:
        raise ValueError(f"unknown Hall lattice symbol {lattice!r}")

    parsed: list[Operation] = []
    shift = ZERO_VECTOR
    previous_nfold = None
    previous_axis = None
    count = 0
    index = 1
    while index < len(tokens):
        if tokens[index].startswith("("):
            shift = _parse_origin(tokens[index:])
            index = len(tokens)
            break
        rotation, translation, nfold, axis, time = _parse_rotation_token(
            tokens[index], count, previous_nfold, previous_axis
        )
        parsed.append((rotation, translation, time))
        previous_nfold = nfold
        previous_axis = axis
        count += 1
        index += 1

    generators: list[Operation] = []
    for translation in LATTICE_SYMBOLS[lattice]:
        if translation != ZERO_VECTOR:
            generators.append((IDENTITY_MATRIX, translation, 1))
    if inversion:
        rotation = negate_matrix(IDENTITY_MATRIX)
        generators.append(
            (
                rotation,
                normalize_vector(
                    vector_subtract(shift, matrix_vector(rotation, shift))
                ),
                1,
            )
        )
    for rotation, translation, time in parsed:
        generators.append(
            (
                rotation,
                normalize_vector(
                    vector_add(
                        translation,
                        vector_subtract(shift, matrix_vector(rotation, shift)),
                    )
                ),
                time,
            )
        )
    return tuple(generators)


def load_hall_symbols(path: Path) -> dict[str, tuple[int, str]]:
    result: dict[str, tuple[int, str]] = {}
    hall_number: int | None = None
    for line_number, raw in enumerate(path.read_text(encoding="utf-8").splitlines(), 1):
        content = raw.split("#", 1)[0].rstrip()
        if not content:
            continue
        parent = re.fullmatch(r"(\d+):", content)
        if parent is not None:
            hall_number = int(parent.group(1))
            continue
        child = re.fullmatch(r'  "([0-9]+\.[0-9]+)":\s*(\S(?:.*\S)?)', content)
        if child is None or hall_number is None:
            raise ValueError(f"unsupported Hall table line {line_number}: {raw!r}")
        bns_number, hall_symbol = child.groups()
        if bns_number in result:
            raise ValueError(f"duplicate BNS number {bns_number}")
        result[bns_number] = hall_number, hall_symbol
    if len(result) != 1651:
        raise ValueError("magnetic Hall table must contain 1,651 entries")
    return result


def operation_json(operation: Operation) -> dict[str, object]:
    rotation, translation, time_reversal = operation
    return {
        "rotation": [list(row) for row in rotation],
        "translation": [str(value) for value in translation],
        "time_reversal": time_reversal,
    }


def operation_digest(operations: Iterable[Operation]) -> str:
    rows = [operation_json(operation) for operation in sorted(operations)]
    encoded = json.dumps(rows, sort_keys=True, separators=(",", ":")).encode()
    return sha256_bytes(encoded)


def _generator_order(operation: Operation) -> tuple[object, ...]:
    rotation, translation, time_reversal = operation
    return (
        0 if rotation == IDENTITY_MATRIX else 1,
        0 if time_reversal == -1 else 1,
        sum(abs(value) for value in translation),
        sum(
            abs(rotation[row][column] - IDENTITY_MATRIX[row][column])
            for row in range(3)
            for column in range(3)
        ),
        rotation,
        translation,
        time_reversal,
    )


def minimal_generators(operations: frozenset[Operation]) -> tuple[Operation, ...]:
    generators: list[Operation] = []
    reached = frozenset({IDENTITY})
    for operation in sorted(operations, key=_generator_order):
        if operation in reached:
            continue
        generators.append(operation)
        reached = closure(generators)
        if reached == operations:
            break
    if reached != operations:
        raise ValueError("generator selection did not reach the target group")
    changed = True
    while changed:
        changed = False
        for index in range(len(generators) - 1, -1, -1):
            trial = generators[:index] + generators[index + 1 :]
            if closure(trial) == operations:
                generators = trial
                changed = True
                break
    return tuple(generators)


def number_metadata() -> tuple[dict[str, str], dict[str, int]]:
    ordinary: dict[str, str] = {}
    uni_numbers: dict[str, int] = {}
    with NUMBERS_PATH.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            bns = row["BNS"]
            og = row["OG"]
            uni = int(row["UNI"])
            if bns in uni_numbers:
                raise ValueError(f"duplicate BNS number {bns}")
            uni_numbers[bns] = uni
            number = bns.split(".", 1)[0]
            if og.split(".", 2)[1] != "1":
                continue
            if number in ordinary:
                raise ValueError(f"duplicate type-I group for ordinary group {number}")
            ordinary[number] = bns
    if set(ordinary) != {str(number) for number in range(1, 231)}:
        raise ValueError("type-I BNS mapping must cover ordinary groups 1..230")
    if set(uni_numbers.values()) != set(range(1, 1652)):
        raise ValueError("UNI mapping must cover 1..1651")
    return ordinary, uni_numbers


def _spglib_fraction(value: float) -> Fraction:
    observed = float(value)
    exact = Fraction(observed).limit_denominator(6)
    if abs(float(exact) - observed) > 1e-12:
        raise ValueError(f"spglib setting translation is not sixth-rational: {value!r}")
    return exact % 1


def _spglib_setting_operations(
    uni_number: int,
    hall_number: int,
) -> frozenset[Operation]:
    raw = spglib.get_magnetic_symmetry_from_database(int(uni_number), int(hall_number))
    if raw is None:
        raise ValueError(
            f"spglib has no magnetic setting UNI{uni_number}/Hall{hall_number}"
        )
    rotations = raw["rotations"]
    translations = raw["translations"]
    time_reversals = raw["time_reversals"]
    if not (len(rotations) == len(translations) == len(time_reversals)):
        raise ValueError("spglib magnetic setting arrays have different lengths")
    operations = frozenset(
        (
            tuple(tuple(int(component) for component in row) for row in rotation),
            tuple(_spglib_fraction(component) for component in translation),
            -1 if int(time_reversal) else 1,
        )
        for rotation, translation, time_reversal in zip(
            rotations, translations, time_reversals, strict=True
        )
    )
    if not operations or any(
        abs(matrix_determinant(row[0])) != 1 for row in operations
    ):
        raise ValueError("spglib magnetic setting contains an invalid operation")
    return operations


def _ordinary_hall_settings() -> dict[int, tuple[int, ...]]:
    settings: dict[int, list[int]] = {}
    for hall_number in range(1, 531):
        setting = spglib.get_spacegroup_type(hall_number)
        if setting is None or int(setting.hall_number) != hall_number:
            raise ValueError(f"spglib has no ordinary Hall setting {hall_number}")
        settings.setdefault(int(setting.number), []).append(hall_number)
    if set(settings) != set(range(1, 231)):
        raise ValueError("spglib Hall settings do not cover ordinary groups 1..230")
    return {number: tuple(values) for number, values in settings.items()}


def _setting_json(hall_number: int, operations: frozenset[Operation]) -> dict[str, Any]:
    generators = minimal_generators(operations)
    if closure(generators) != operations:
        raise AssertionError("setting generators changed the target group")
    setting_type = spglib.get_spacegroup_type(int(hall_number))
    if setting_type is None:
        raise ValueError(f"spglib has no Hall metadata {hall_number}")
    return {
        "setting_key": f"hall:{int(hall_number)}",
        "authority": "hall",
        "hall_number": int(hall_number),
        "source_setting_ids": [],
        "choice": str(setting_type.choice),
        "operation_count": len(operations),
        "operation_sha256": operation_digest(operations),
        "generators": [operation_json(operation) for operation in generators],
    }


def _source_setting_rows() -> dict[int, tuple[dict[str, Any], ...]]:
    """Return the Source-only presentation sets absent from Hall 1..530.

    Source contains 742 current (non-version-2) setting rows, but signed-axis
    aliases collapse to 542 exact operation sets.  Hall already supplies 527
    of those sets.  Keeping only the 15 exact differences prevents presentation
    coverage from turning into a 742-way runtime scan.
    """

    from APOSTRUCT.Backend.exactmath import fraction_matrix_inverse3
    from APOSTRUCT.Backend.source.magnetic_operations import (
        _pml_point_operation_matrix,
        fraction_values,
    )
    from APOSTRUCT.Backend.source.tables import source_tables

    require_source(SOURCE_SPACE_PATH, SOURCE_SPACE_SHA256)
    tables = source_tables()
    centering = {
        "P": (ZERO_VECTOR,),
        "A": (ZERO_VECTOR, (ZERO, HALF, HALF)),
        "B": (ZERO_VECTOR, (HALF, ZERO, HALF)),
        "C": (ZERO_VECTOR, (HALF, HALF, ZERO)),
        "I": (ZERO_VECTOR, (HALF, HALF, HALF)),
        "F": (
            ZERO_VECTOR,
            (ZERO, HALF, HALF),
            (HALF, ZERO, HALF),
            (HALF, HALF, ZERO),
        ),
        "R": (
            ZERO_VECTOR,
            (Fraction(2, 3), Fraction(1, 3), Fraction(1, 3)),
            (Fraction(1, 3), Fraction(2, 3), Fraction(2, 3)),
        ),
    }

    def source_operations(sg: int, setting_id: int) -> frozenset[Operation]:
        matrix = tables.pml_to_cinter_matrix(sg, setting_id)
        inverse = fraction_matrix_inverse3(matrix)
        origin = tables.cml_to_cinter_origin(sg, setting_id)
        record = tables.inter_setting_record(setting_id)
        lattice = str(record["label_short"])[0].upper()
        centers = (
            (ZERO_VECTOR,)
            if lattice == "R" and str(record["axis"]) == "R"
            else centering[lattice]
        )
        operations: set[Operation] = set()
        for raw in tables.generate_space_group_records(sg):
            pml_rotation = _pml_point_operation_matrix(sg, int(raw[4]))
            row_rotation = matrix_multiply(
                matrix_multiply(inverse, pml_rotation), matrix
            )
            translated = matrix_vector(
                tuple(tuple(value for value in row) for row in zip(*matrix)),
                fraction_values(raw[:4]),
            )
            rotated_origin = tuple(
                sum(origin[axis] * row_rotation[axis][column] for axis in range(3))
                for column in range(3)
            )
            base_translation = tuple(
                (translated[axis] - rotated_origin[axis] + origin[axis]) % 1
                for axis in range(3)
            )
            column_rotation = tuple(
                tuple(int(row_rotation[column][row]) for column in range(3))
                for row in range(3)
            )
            for center in centers:
                operations.add(
                    (
                        column_rotation,  # type: ignore[arg-type]
                        tuple(
                            (base_translation[axis] + center[axis]) % 1
                            for axis in range(3)
                        ),  # type: ignore[arg-type]
                        1,
                    )
                )
        return frozenset(operations)

    # Ordinary operations do not use a UNI number.  Build the exact Hall index
    # directly from the ordinary database instead.
    hall_by_operations = {}
    hall_operations_by_number = {}
    for hall_number in range(1, 531):
        raw = spglib.get_symmetry_from_database(hall_number)
        if raw is None:
            raise ValueError(f"spglib has no ordinary Hall setting {hall_number}")
        operations = frozenset(
            (
                tuple(tuple(int(value) for value in row) for row in rotation),
                tuple(_spglib_fraction(value) for value in translation),
                1,
            )
            for rotation, translation in zip(
                raw["rotations"], raw["translations"], strict=True
            )
        )
        hall_operations_by_number[hall_number] = operations
        hall_by_operations.setdefault(operations, hall_number)

    source_sets: dict[
        frozenset[Operation], tuple[int, list[int]]
    ] = {}
    for setting_id, raw_sg in enumerate(tables.space["ispace_inter_number"], 1):
        sg = int(raw_sg)
        if not sg:
            continue
        record = tables.inter_setting_record(setting_id)
        if int(record["version"]) == 2:
            continue
        operations = source_operations(sg, setting_id)
        prior = source_sets.get(operations)
        if prior is None:
            source_sets[operations] = (sg, [setting_id])
        else:
            if prior[0] != sg:
                raise ValueError("Source operation set aliases different space groups")
            prior[1].append(setting_id)

    rows: dict[int, list[dict[str, Any]]] = {}
    for operations, (sg, setting_ids) in source_sets.items():
        if operations in hall_by_operations:
            continue
        target_id = min(setting_ids)
        target_matrix = tables.pml_to_cinter_matrix(sg, target_id)
        target_origin = tables.cml_to_cinter_origin(sg, target_id)
        base_candidates: list[tuple[int, int]] = []
        for base_id in tables.inter_setting_ids_for_space_group(sg):
            base_record = tables.inter_setting_record(base_id)
            if int(base_record["version"]) == 2:
                continue
            hall_number = hall_by_operations.get(source_operations(sg, base_id))
            if hall_number is not None:
                base_candidates.append((base_id, hall_number))
        if not base_candidates:
            raise ValueError(f"Source-only SG{sg} setting has no Hall-overlap base")
        base_id, base_hall = min(base_candidates)
        base_matrix = tables.pml_to_cinter_matrix(sg, base_id)
        base_origin = tables.cml_to_cinter_origin(sg, base_id)
        row_matrix = matrix_multiply(
            fraction_matrix_inverse3(base_matrix), target_matrix
        )
        row_shift = tuple(
            target_origin[column]
            - sum(base_origin[axis] * row_matrix[axis][column] for axis in range(3))
            for column in range(3)
        )
        column_matrix = tuple(
            tuple(row_matrix[column][row] for column in range(3))
            for row in range(3)
        )
        row = {
            "setting_key": f"source:{target_id}",
            "authority": "source",
            "hall_number": base_hall,
            "source_setting_ids": sorted(setting_ids),
            "choice": str(tables.inter_setting_record(target_id)["label_short"]),
            "coordinate_matrix": [
                [str(value) for value in line] for line in column_matrix
            ],
            "coordinate_shift": [str(value) for value in row_shift],
            "ordinary_operation_sha256": operation_digest(operations),
        }
        transported = _transform_setting_operations(
            hall_operations_by_number[base_hall], row
        )
        if operation_digest(transported) != row["ordinary_operation_sha256"]:
            raise ValueError("Source affine transport does not reproduce its setting")
        rows.setdefault(sg, []).append(row)
    if sum(map(len, rows.values())) != 15:
        raise ValueError("Source-only presentation coverage must contain 15 sets")
    return {
        sg: tuple(sorted(values, key=lambda row: str(row["setting_key"])))
        for sg, values in rows.items()
    }


def _fraction_token(token: str) -> Fraction:
    value = Fraction(token)
    if str(value) != token:
        raise ValueError(f"noncanonical Source affine fraction {token!r}")
    return value


def _transform_setting_operations(
    operations: frozenset[Operation], row: Mapping[str, Any]
) -> frozenset[Operation]:
    from APOSTRUCT.Backend.exactmath import fraction_matrix_inverse3

    matrix = tuple(
        tuple(_fraction_token(str(value)) for value in raw_row)
        for raw_row in row["coordinate_matrix"]
    )
    inverse = fraction_matrix_inverse3(matrix)
    shift = tuple(_fraction_token(str(value)) for value in row["coordinate_shift"])
    transformed: set[Operation] = set()
    for rotation, translation, time_reversal in operations:
        changed_rotation = matrix_multiply(
            matrix_multiply(matrix, rotation), inverse
        )
        changed_translation = tuple(
            (
                sum(matrix[axis][column] * translation[column] for column in range(3))
                + shift[axis]
                - sum(changed_rotation[axis][column] * shift[column] for column in range(3))
            )
            % 1
            for axis in range(3)
        )
        if any(value.denominator != 1 for line in changed_rotation for value in line):
            raise ValueError("Source presentation produced a nonintegral rotation")
        transformed.add(
            (
                tuple(tuple(int(value) for value in line) for line in changed_rotation),  # type: ignore[arg-type]
                changed_translation,  # type: ignore[arg-type]
                time_reversal,
            )
        )
    return frozenset(transformed)


def _source_setting_json(
    row: Mapping[str, Any], operations: frozenset[Operation]
) -> dict[str, Any]:
    generators = minimal_generators(operations)
    return {
        "setting_key": row["setting_key"],
        "authority": "source",
        "hall_number": row["hall_number"],
        "source_setting_ids": row["source_setting_ids"],
        "choice": row["choice"],
        "operation_count": len(operations),
        "operation_sha256": operation_digest(operations),
        "generators": [operation_json(operation) for operation in generators],
    }


def build_document() -> tuple[dict[str, Any], dict[str, frozenset[Operation]]]:
    if str(spglib.__version__) != SPGLIB_SETTING_VERSION:
        raise ValueError(
            f"setting authority requires spglib {SPGLIB_SETTING_VERSION}, "
            f"found {spglib.__version__}"
        )
    require_source(HALL_PATH, HALL_SHA256)
    require_source(NUMBERS_PATH, NUMBERS_SHA256)
    require_source(LICENSE_PATH, LICENSE_SHA256)
    symbols = load_hall_symbols(HALL_PATH)
    ordinary, uni_numbers = number_metadata()
    entries: dict[str, dict[str, object]] = {}
    complete: dict[str, frozenset[Operation]] = {}
    total_operations = 0
    setting_count = 0
    setting_operation_count = 0
    incomplete_setting_groups: list[str] = []
    hall_settings = _ordinary_hall_settings()
    source_settings = _source_setting_rows()
    for label in sorted(symbols, key=lambda value: tuple(map(int, value.split(".")))):
        hall_number, raw_symbol = symbols[label]
        hall_symbol = OFFICIAL_PARITY_OVERRIDES.get(label, raw_symbol)
        operations = closure(hall_generators(hall_symbol))
        generators = minimal_generators(operations)
        if closure(generators) != operations:
            raise AssertionError("minimal generators changed the target group")
        complete[label] = operations
        total_operations += len(operations)
        ordinary_number = int(label.split(".", 1)[0])
        api_default = _spglib_setting_operations(uni_numbers[label], hall_number)
        setting_coverage_complete = api_default == operations
        if setting_coverage_complete:
            setting_rows = [
                _setting_json(
                    setting_hall,
                    _spglib_setting_operations(uni_numbers[label], setting_hall),
                )
                for setting_hall in hall_settings[ordinary_number]
            ]
            for source_row in source_settings.get(ordinary_number, ()):
                base_operations = _spglib_setting_operations(
                    uni_numbers[label], int(source_row["hall_number"])
                )
                transformed = _transform_setting_operations(
                    base_operations, source_row
                )
                setting_rows.append(_source_setting_json(source_row, transformed))
        else:
            incomplete_setting_groups.append(label)
            setting_rows = [_setting_json(hall_number, operations)]
        setting_count += len(setting_rows)
        setting_operation_count += sum(
            int(setting["operation_count"]) for setting in setting_rows
        )
        entries[label] = {
            "uni_number": uni_numbers[label],
            "hall_number": hall_number,
            "hall_symbol": hall_symbol,
            "operation_count": len(operations),
            "operation_sha256": operation_digest(operations),
            "generators": [operation_json(operation) for operation in generators],
            "setting_coverage_complete": setting_coverage_complete,
            "settings": setting_rows,
        }
    if len(entries) != 1651 or total_operations != 38_307:
        raise ValueError("magnetic group authority coverage changed")
    return (
        {
            "schema": "isodistort.validation.magnetic-group-generators.v3",
            "sources": {
                "license": "BSD-3-Clause",
                "spglib_version": SPGLIB_VERSION,
                "spglib_commit": SPGLIB_COMMIT,
                "spglib_magnetic_hall_symbols_sha256": HALL_SHA256,
                "spglib_msg_numbers_sha256": NUMBERS_SHA256,
                "spglib_license_sha256": LICENSE_SHA256,
                "spglib_setting_version": SPGLIB_SETTING_VERSION,
                "official_iso_mag_sha256": OFFICIAL_ISO_MAG_SHA256,
                "official_parity_overrides": OFFICIAL_PARITY_OVERRIDES,
                "source_data_space_sha256": SOURCE_SPACE_SHA256,
            },
            "group_count": len(entries),
            "total_operation_count": total_operations,
            "setting_count": setting_count,
            "setting_operation_count": setting_operation_count,
            "incomplete_setting_groups": incomplete_setting_groups,
            "ordinary_groups": ordinary,
            "entries": entries,
        },
        complete,
    )


def encode_document(document: dict[str, Any]) -> bytes:
    payload = (
        json.dumps(document, sort_keys=True, separators=(",", ":")) + "\n"
    ).encode()
    raw = io.BytesIO()
    with gzip.GzipFile(
        filename="", mode="wb", compresslevel=9, fileobj=raw, mtime=0
    ) as handle:
        handle.write(payload)
    return raw.getvalue()


class _TokenStream:
    def __init__(self, handle: TextIO):
        self._handle = handle
        self._buffer: list[str] = []

    def next(self) -> str:
        while not self._buffer:
            line = self._handle.readline()
            if not line:
                raise EOFError("unexpected end of ISO-MAG table")
            self._buffer = shlex.split(line, comments=False, posix=True)
        return self._buffer.pop(0)

    def integer(self) -> int:
        return int(self.next())

    def require_end(self) -> None:
        if self._buffer or any(line.split() for line in self._handle):
            raise ValueError("ISO-MAG table contains trailing records")


def _matrix3(values: Sequence[int]) -> Matrix:
    if len(values) != 9:
        raise ValueError("ISO-MAG rotation must contain nine values")
    return tuple(
        tuple(values[3 * row + column] for column in range(3)) for row in range(3)
    )  # type: ignore[return-value]


def _translation_cosets(generators: Iterable[Vector]) -> frozenset[Vector]:
    generators = tuple(normalize_vector(vector) for vector in generators)
    seen = {ZERO_VECTOR}
    queue = deque([ZERO_VECTOR])
    while queue:
        current = queue.popleft()
        for generator in generators:
            translated = normalize_vector(
                current[axis] + generator[axis] for axis in range(3)
            )
            if translated in seen:
                continue
            seen.add(translated)
            queue.append(translated)
            if len(seen) > 128:
                raise ValueError("ISO-MAG translation quotient exceeds 128")
    return frozenset(seen)


def _skip_wyckoff_rows(tokens: _TokenStream, count: int) -> None:
    for _index in range(count):
        position_count = tokens.integer()
        tokens.integer()
        tokens.next()
        for _position in range(position_count):
            for _value in range(22):
                tokens.integer()


def official_iso_mag_operations(path: Path) -> dict[str, frozenset[Operation]]:
    if sha256_file(path) != OFFICIAL_ISO_MAG_SHA256:
        raise ValueError("official ISO-MAG table differs from the pinned source")
    _ordinary, uni_numbers = number_metadata()
    with path.open("r", encoding="utf-8") as handle:
        tokens = _TokenStream(handle)
        nonhexagonal: list[Matrix] = []
        for expected in range(1, 49):
            if tokens.integer() != expected:
                raise ValueError("ISO-MAG nonhexagonal point-operation order changed")
            tokens.next()
            tokens.next()
            nonhexagonal.append(_matrix3([tokens.integer() for _value in range(9)]))
        hexagonal: list[Matrix] = []
        for expected in range(1, 25):
            if tokens.integer() != expected:
                raise ValueError("ISO-MAG hexagonal point-operation order changed")
            tokens.next()
            tokens.next()
            hexagonal.append(_matrix3([tokens.integer() for _value in range(9)]))

        groups: dict[str, frozenset[Operation]] = {}
        for expected_uni in range(1, 1652):
            bns_first = tokens.integer()
            bns_second = tokens.integer()
            bns_label = tokens.next()
            if bns_label != f"{bns_first}.{bns_second}":
                raise ValueError("ISO-MAG BNS components and label disagree")
            if uni_numbers.get(bns_label) != expected_uni:
                raise ValueError("ISO-MAG and spglib UNI numbering disagree")
            tokens.next()
            tokens.next()
            tokens.integer()
            tokens.integer()
            tokens.integer()
            tokens.next()
            tokens.next()
            magnetic_type = tokens.integer()
            if magnetic_type == 4:
                for _value in range(13):
                    tokens.integer()

            operation_rows = []
            for _operation in range(tokens.integer()):
                point_index = tokens.integer()
                numerators = [tokens.integer() for _axis in range(3)]
                denominator = tokens.integer()
                time_reversal = tokens.integer()
                operation_rows.append(
                    (point_index, numerators, denominator, time_reversal)
                )
            lattice_vectors = []
            for _vector in range(tokens.integer()):
                numerators = [tokens.integer() for _axis in range(3)]
                denominator = tokens.integer()
                lattice_vectors.append(
                    tuple(Fraction(value, denominator) for value in numerators)
                )
            _skip_wyckoff_rows(tokens, tokens.integer())
            if magnetic_type == 4:
                for _operation in range(tokens.integer()):
                    for _value in range(6):
                        tokens.integer()
                for _vector in range(tokens.integer()):
                    for _value in range(4):
                        tokens.integer()
                _skip_wyckoff_rows(tokens, tokens.integer())

            point_table = hexagonal if 143 <= bns_first <= 194 else nonhexagonal
            translations = _translation_cosets(lattice_vectors)
            operations = {
                (
                    point_table[point_index - 1],
                    normalize_vector(
                        Fraction(numerators[axis], denominator) + centering[axis]
                        for axis in range(3)
                    ),
                    time_reversal,
                )
                for point_index, numerators, denominator, time_reversal in operation_rows
                for centering in translations
            }
            if bns_label in groups:
                raise ValueError(f"duplicate ISO-MAG group {bns_label}")
            groups[bns_label] = frozenset(operations)
            if not operations:
                raise ValueError(f"ISO-MAG group {bns_label} has no operations")
        tokens.require_end()
    if len(groups) != 1651:
        raise ValueError("ISO-MAG table must contain 1,651 groups")
    return groups


def audit_official_iso_mag(
    expected: dict[str, frozenset[Operation]], path: Path
) -> None:
    official = official_iso_mag_operations(path)
    if set(official) != set(expected):
        raise ValueError("ISO-MAG and Hall tables have different BNS populations")
    mismatches = [
        label
        for label in sorted(
            expected, key=lambda value: tuple(map(int, value.split(".")))
        )
        if official[label] != expected[label]
    ]
    if mismatches:
        raise ValueError(
            "ISO-MAG parity differs for corrected Hall groups: "
            + ", ".join(mismatches[:20])
        )
    if sum(map(len, official.values())) != 38_307:
        raise ValueError("ISO-MAG operation total differs from 38,307")


def _write_atomic(path: Path, value: bytes) -> None:
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.unlink(missing_ok=True)
    try:
        temporary.write_bytes(value)
        temporary.replace(path)
    finally:
        temporary.unlink(missing_ok=True)


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--write",
        action="store_true",
        help="replace the committed artifact after a complete ISO-MAG audit",
    )
    parser.add_argument(
        "--official-iso-mag",
        type=Path,
        help="pinned magnetic_data.txt used for the exact 1,651-group parity audit",
    )
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    args = _parser().parse_args(argv)
    document, operations = build_document()
    encoded = encode_document(document)
    if args.official_iso_mag is not None:
        audit_official_iso_mag(operations, args.official_iso_mag.resolve())
    if args.write:
        if args.official_iso_mag is None:
            raise ValueError("--write requires --official-iso-mag")
        _write_atomic(OUTPUT_PATH, encoded)
    elif not OUTPUT_PATH.is_file() or OUTPUT_PATH.read_bytes() != encoded:
        raise ValueError("committed magnetic group authority is not reproducible")
    print(
        json.dumps(
            {
                "artifact_sha256": sha256_bytes(encoded),
                "groups": document["group_count"],
                "official_audit": args.official_iso_mag is not None,
                "operations": document["total_operation_count"],
                "wrote": bool(args.write),
            },
            sort_keys=True,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
