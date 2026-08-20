"""Independently check one Web "Complete Mode Details" text export.

Companion to verify_local_output.py, for the Web side of a comparison: paste
in the "Complete Mode Details" text ISODISTORT's public Web version renders
for the same CIF/k/irrep/OPD, and this runs the same two theorems against
it. Needs nothing else — no content-addressed store, no Records/, no network
access of its own; fetching the Web page is left to whoever runs this.

Usage:
    python3 Verification/verify_web_text.py complete_mode_details.txt
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from Verification.mathematics.group_invariance import (  # noqa: E402
    GROUP_INVARIANCE_THEOREM,
    assess_group_invariance,
    subgroup_label_from_web_text,
)
from Verification.mathematics.mode_basis import (  # noqa: E402
    MODE_BASIS_THEOREM,
    assess_mode_basis,
    definitions_from_web_text,
    structure_from_web_text,
)


def verify(web_text_path: Path) -> dict[str, dict]:
    text = web_text_path.read_text(encoding="utf-8")
    definitions = definitions_from_web_text(text)
    structure = structure_from_web_text(text)
    subgroup_label = subgroup_label_from_web_text(text)

    return {
        MODE_BASIS_THEOREM: assess_mode_basis(definitions),
        GROUP_INVARIANCE_THEOREM: assess_group_invariance(definitions, subgroup_label, structure),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("web_text", type=Path, help="Saved 'Complete Mode Details' text from the Web app")
    args = parser.parse_args()

    results = verify(args.web_text)
    ok = True
    for theorem, result in results.items():
        status = result.get("status", "unknown")
        print(f"{theorem}: {status}")
        if status not in ("satisfied", "not_applicable"):
            ok = False
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
