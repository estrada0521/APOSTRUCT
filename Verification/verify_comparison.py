"""Independently reproduce a strict_pass/physical_pass/physical_fail verdict.

This is the actual comparator this project's Strict-pass and Physical-pass
figures come from — not the independent math check verify_local_output.py
and verify_web_text.py run, but the Local-vs-Web equivalence judgment
itself. Without this, a third party can only see the published percentages,
never what "physical pass" concretely means for one branch; this makes that
opaque.

This calls the same `Verification.compare.compare_local_to_web` and
`Verification.result.result_for_comparison` functions used by the private
monorepo Validation operation. Nothing about either is reimplemented here.
The only thing this file adds is building their plain arguments from an `apo`
CLI output and a saved Web text export directly, instead of from private
content-addressed stores and provenance checks.

Usage:
    apo modes structure.cif --k X --irrep X3- --displacive Li --opd P1 \
        --full-state -o local.json
    # Separately obtain web_text.txt: the "Complete Mode Details" text for
    # the same structure/k/irrep/OPD from ISODISTORT's public Web version.
    python3 Verification/verify_comparison.py local.json web_text.txt \
        --input case.json --opd P1
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from Verification.compare import compare_local_to_web  # noqa: E402
from Verification.result import (  # noqa: E402
    DEFAULT_TOLERANCE,
    result_for_comparison,
)
from Verification.version import comparator_version  # noqa: E402


def _request_selections(input_payload: dict[str, Any]) -> list[dict[str, Any]]:
    return [
        {
            "k_label": slot["label"],
            "k_params": dict(slot.get("params") or {}),
            "irrep": slot["ir"],
        }
        for slot in input_payload["k"]
    ]


def verify(
    local_apo_output: dict[str, Any],
    web_text: str,
    input_payload: dict[str, Any],
    *,
    input_id: str,
    opd_label: str,
    tol: float,
) -> dict[str, Any]:
    mode_details = local_apo_output["selected"]["mode_details"]
    if mode_details.get("status") != "ok":
        raise ValueError(
            f"local mode details are not complete (status={mode_details.get('status')!r})"
        )
    preview = {
        "space_group": local_apo_output.get("space_group") or {},
        "input": local_apo_output.get("input") or {},
        "selected": local_apo_output["selected"],
    }

    opd, comparison = compare_local_to_web(
        opd_exists=True,
        preview=preview,
        mode_text=web_text,
        opd_label=opd_label,
        k_selections=_request_selections(input_payload),
        case_label=input_id,
        cif_label=str(input_payload.get("cif", "unknown.cif")),
        comparison_source="verify_comparison",
        tol=tol,
    )

    result = result_for_comparison(opd, comparison, tol=tol)
    result["comparator_version"] = comparator_version(ROOT)
    return result


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("local_output", type=Path, help="JSON from `apo modes --full-state -o ...`")
    parser.add_argument("web_text", type=Path, help="Saved 'Complete Mode Details' text from the Web app")
    parser.add_argument("--input", type=Path, required=True, help="Canonical input JSON (case.json)")
    parser.add_argument("--opd", required=True, help="Web-selected OPD label, e.g. P1")
    parser.add_argument("--input-id", default="unknown", help="Optional label for the result row")
    parser.add_argument("--tol", type=float, default=DEFAULT_TOLERANCE)
    args = parser.parse_args()

    local_apo_output = json.loads(args.local_output.read_text())
    web_text = args.web_text.read_text(encoding="utf-8")
    input_payload = json.loads(args.input.read_text())

    result = verify(
        local_apo_output,
        web_text,
        input_payload,
        input_id=args.input_id,
        opd_label=args.opd,
        tol=args.tol,
    )
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0 if result.get("verdict") in ("strict_pass", "physical_pass") else 1


if __name__ == "__main__":
    raise SystemExit(main())
