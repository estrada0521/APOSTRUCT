"""Combine published compact mode definitions without rerunning the pipeline."""

from __future__ import annotations

from fractions import Fraction
import json
from pathlib import Path
import sys
from typing import Any, Sequence


_DEFINITION_FIELDS = (
    ("displacive", "displacive_definitions"),
    ("magnetic", "magnetic_definitions"),
)


def load_compact_modes(path: Path) -> dict[str, Any]:
    try:
        if str(path) == "-":
            value = json.load(sys.stdin)
        else:
            with path.expanduser().resolve().open(encoding="utf-8") as stream:
                value = json.load(stream)
    except json.JSONDecodeError as exc:
        raise ValueError(f"invalid compact modes JSON: {exc}") from exc
    if not isinstance(value, dict):
        raise ValueError("compact modes JSON must be an object")
    if value.get("schema") != "distortropy.cli.modes":
        raise ValueError("expected a distortropy.cli.modes payload")
    return value


def _weight_rows(raw_weights: Sequence[str]) -> list[tuple[str, Fraction]]:
    if not raw_weights:
        raise ValueError("combine-modes requires at least one weight")
    rows: list[tuple[str, Fraction]] = []
    seen: set[str] = set()
    for raw in raw_weights:
        definition_id, separator, raw_value = str(raw).partition("=")
        definition_id = definition_id.strip()
        raw_value = raw_value.strip()
        if not separator or not definition_id or not raw_value:
            raise ValueError("weights must use DEFINITION_ID=VALUE")
        if definition_id in seen:
            raise ValueError(f"duplicate weight for definition {definition_id}")
        seen.add(definition_id)
        try:
            value = Fraction(raw_value)
        except (ValueError, ZeroDivisionError) as exc:
            raise ValueError(
                f"weight for definition {definition_id} must be an exact number"
            ) from exc
        rows.append((definition_id, value))
    return rows


def _vector(raw: Sequence[float]) -> list[float]:
    return [float(value) for value in raw]


def _definitions(mode_details: dict[str, Any]) -> dict[str, tuple[str, dict[str, Any]]]:
    indexed: dict[str, tuple[str, dict[str, Any]]] = {}
    for kind, field in _DEFINITION_FIELDS:
        for definition in mode_details.get(field) or []:
            definition_id = definition["definition_id"]
            indexed[definition_id] = (kind, definition)
    return indexed


def _viewer_atoms(
    mode_details: dict[str, Any],
) -> tuple[list[dict[str, Any]], dict[str, int]]:
    atoms: list[dict[str, Any]] = []
    positions: dict[str, int] = {}
    for raw in mode_details.get("viewer_atoms") or []:
        atom_id = raw["atom_id"]
        positions[atom_id] = len(atoms)
        atoms.append(
            {
                "atom_id": atom_id,
                "child_site": raw["child_site"],
                "element": raw.get("element"),
                "site_order": raw["site_order"],
                "xyz": _vector(raw["xyz"]),
            }
        )
    return atoms, positions


def combine_mode_payload(
    payload: dict[str, Any],
    *,
    raw_weights: Sequence[str],
    apply_normfactor: bool = False,
) -> dict[str, Any]:
    """Apply explicit weights to definition rows in one compact modes payload."""

    if payload.get("schema") != "distortropy.cli.modes":
        raise ValueError("combine-modes requires a distortropy.cli.modes payload")
    mode_details = payload.get("mode_details")
    if not isinstance(mode_details, dict) or mode_details.get("status") != "ok":
        raise ValueError("compact modes payload has no completed mode details")
    indexed = _definitions(mode_details)
    weights = _weight_rows(raw_weights)
    atoms, atom_positions = _viewer_atoms(mode_details)
    displacements = [[0.0, 0.0, 0.0] for _ in atoms]
    magnetic_vectors = [[0.0, 0.0, 0.0] for _ in atoms]
    touched_displacive: set[int] = set()
    touched_magnetic: set[int] = set()
    applied: list[dict[str, Any]] = []

    for definition_id, exact_weight in weights:
        try:
            kind, definition = indexed[definition_id]
        except KeyError as exc:
            available = ", ".join(indexed)
            raise ValueError(
                f"definition {definition_id} not found; available IDs: [{available}]"
            ) from exc
        input_weight = float(exact_weight)
        normfactor: float | None = None
        if apply_normfactor:
            normfactor = float(definition["normfactor"])
        weight = input_weight * (normfactor if normfactor is not None else 1.0)
        target = displacements if kind == "displacive" else magnetic_vectors
        touched = touched_displacive if kind == "displacive" else touched_magnetic
        for row in definition["rows"]:
            position = atom_positions[row["atom_id"]]
            vector = _vector(row["dxyz"])
            touched.add(position)
            for axis in range(3):
                target[position][axis] += weight * vector[axis]
        applied_row = {
            "definition_id": definition_id,
            "kind": kind,
            "weight": str(exact_weight),
            "effective_row_weight": weight,
        }
        if normfactor is not None:
            applied_row["normfactor"] = normfactor
        applied.append(applied_row)

    combined_atoms: list[dict[str, Any]] = []
    for index, atom in enumerate(atoms):
        if index not in touched_displacive and index not in touched_magnetic:
            continue
        combined = dict(atom)
        if index in touched_displacive:
            combined["displacement"] = displacements[index]
        if index in touched_magnetic:
            combined["magnetic_vector"] = magnetic_vectors[index]
        combined_atoms.append(combined)

    result: dict[str, Any] = {
        "schema": "distortropy.cli.mode_combination",
        "source_schema": "distortropy.cli.modes",
        "coordinate_system": "child_crystallographic_dxyz",
        "weight_convention": (
            "amplitudes use each definition's published normfactor; no common physical unit is assigned"
            if apply_normfactor
            else "weights multiply published definition rows; normfactor is not applied"
        ),
        "weights": applied,
        "lattice": mode_details.get("lattice"),
        "atoms": combined_atoms,
    }
    if touched_magnetic:
        net_magnetic_vector: list[float] = []
        for axis in range(3):
            total = 0.0
            for vector in magnetic_vectors:
                total += vector[axis]
            net_magnetic_vector.append(total)
        result["net_magnetic_vector"] = net_magnetic_vector
        result["net_magnetic_vector_scope"] = "returned_conventional_child_cell"
    return result


__all__ = ["combine_mode_payload", "load_compact_modes"]
