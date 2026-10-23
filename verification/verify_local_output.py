"""Independently check one `apo modes --full-state` JSON output.

This is the means of independent verification, not a report of results: it
runs the same two theorems Validation/math_validation.py runs internally
(exact rational linear independence of the mode basis; group-invariance of
the displayed mode fields under the declared subgroup), but takes a plain
`apo` CLI JSON file directly and needs nothing else — no content-addressed
store, no Records/, no network.

The `apo modes ... --full-state -o output.json` CLI writes mode-detail data
at `selected.mode_details`. Its vector definitions are expanded into per-atom
rows (verification.local_output, the same expansion the private Validation
harness uses) and wrapped as `{"preview": {"selected": ...}}`, the shape the
theorem parsers expect; a non-"ok" mode_details.status is rejected.

Usage:
    apo modes structure.cif --k X --irrep X3- --displacive Li --opd P1 \
        --full-state -o output.json
    python3 -m verification.verify_local_output output.json
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from verification.local_output import with_definition_rows  # noqa: E402
from verification.mathematics.group_invariance import (  # noqa: E402
    GROUP_INVARIANCE_THEOREM,
    assess_group_invariance,
    subgroup_label_from_local_payload,
)
from verification.mathematics.mode_basis import (  # noqa: E402
    MODE_BASIS_THEOREM,
    _local_definition,
    _local_structure_orbit,
    assess_mode_basis,
)


def load_definitions_and_structure(mode_details: dict) -> tuple[list, list]:
    definitions = []
    for raw in mode_details.get("displacive_definitions", []):
        definitions.append(_local_definition(raw, ordinal=len(definitions), kind="dsp"))
    for raw in mode_details.get("magnetic_definitions", []):
        definitions.append(_local_definition(raw, ordinal=len(definitions), kind="mag"))
    for raw in mode_details.get("occupancy_definitions", []):
        definitions.append(_local_definition(raw, ordinal=len(definitions), kind="occ"))
    structure = [_local_structure_orbit(raw) for raw in mode_details.get("undistorted_atoms", [])]
    return definitions, structure


def verify(apo_output_path: Path) -> dict[str, dict]:
    raw = json.loads(apo_output_path.read_text())
    if "selected" not in raw or "mode_details" not in raw.get("selected", {}):
        raise ValueError(
            "not an `apo modes --full-state` output "
            "(missing selected.mode_details; rerun with --full-state)"
        )
    mode_details = with_definition_rows(raw["selected"]["mode_details"])
    status = mode_details.get("status")
    if status != "ok":
        reason = str(mode_details.get("reason") or "").strip()
        detail = f": {reason}" if reason else ""
        raise ValueError(f"mode details are not complete (status={status!r}){detail}")
    payload = {"preview": {"selected": {"mode_details": mode_details}}}

    definitions, structure = load_definitions_and_structure(mode_details)
    subgroup_label = subgroup_label_from_local_payload(payload)

    return {
        MODE_BASIS_THEOREM: assess_mode_basis(definitions),
        GROUP_INVARIANCE_THEOREM: assess_group_invariance(definitions, subgroup_label, structure),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("apo_output", type=Path, help="JSON file from `apo modes --full-state -o ...`")
    args = parser.parse_args()

    results = verify(args.apo_output)
    ok = True
    for theorem, result in results.items():
        status = result.get("status", "unknown")
        print(f"{theorem}: {status}")
        if status not in ("satisfied", "not_applicable"):
            ok = False
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
