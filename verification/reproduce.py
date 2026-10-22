"""Recompute Verdicts/Population tables from branches/results.jsonl.gz.

This script depends only on the small published branch result file
and the standard library. It does not read the internal `log/` trees, does
not need APOSTRUCT, and does not run any computation of its own: it applies
the documented population-eligibility rule and verdict counting to
already-judged per-branch records.

It prints numbers; it does not tell you whether they are the right numbers.
Comparing the output to any particular published document is left to
whoever runs it.

Usage:
    python3 -m verification.reproduce
"""

from __future__ import annotations

import argparse
import gzip
import json
from collections import Counter, defaultdict
from pathlib import Path

from verification.result import canonical_verdict


def default_results_path() -> Path:
    checkout_path = Path.cwd() / "branches" / "results.jsonl.gz"
    if checkout_path.is_file():
        return checkout_path
    return Path(__file__).resolve().parents[1] / "branches" / "results.jsonl.gz"

THEOREMS = (
    "complete_mode_basis_linear_independence.v1",
    "complete_mode_field_group_invariance.v1",
)


def load_records(path: Path) -> list[dict]:
    with gzip.open(path, "rt") as f:
        return [json.loads(line) for line in f]


def is_eligible(record: dict) -> bool:
    """Web output not refuted by either independent mathematics check."""

    web_math = record.get("web_math")
    if not web_math or not record.get("local_output"):
        return False
    for theorem in THEOREMS:
        status = web_math.get(theorem)
        if status not in ("satisfied", "not_applicable"):
            return False
    return True


def summarize(records: list[dict]) -> dict[str, dict]:
    by_pool: dict[str, list[dict]] = defaultdict(list)
    for record in records:
        if record["pool"] in ("mp", "magndata") and is_eligible(record):
            by_pool[record["pool"]].append(record)

    summary: dict[str, dict] = {}
    for pool, rows in by_pool.items():
        verdicts = Counter(canonical_verdict(r["verdict"]) for r in rows)
        identical = verdicts.get("identical", 0)
        equivalent = identical + verdicts.get("equivalent", 0)
        population = len(rows)

        by_k: dict[str, tuple[int, int, int]] = {}
        rows_by_k: dict[str, list[dict]] = defaultdict(list)
        for r in rows:
            rows_by_k[r["k_signature"]].append(r)
        for k_sig, k_rows in rows_by_k.items():
            k_verdicts = Counter(canonical_verdict(r["verdict"]) for r in k_rows)
            k_identical = k_verdicts.get("identical", 0)
            k_equivalent = k_identical + k_verdicts.get("equivalent", 0)
            by_k[k_sig] = (len(k_rows), k_equivalent, k_identical)

        summary[pool] = {
            "population": population,
            "equivalent": equivalent,
            "identical": identical,
            "by_k": by_k,
        }
    return summary


def _fmt_pct(numerator: int, denominator: int) -> str:
    if denominator == 0:
        return "n/a"
    return f"{100.0 * numerator / denominator:.2f}%"


def print_report(summary: dict[str, dict]) -> None:
    for pool in ("mp", "magndata"):
        if pool not in summary:
            continue
        s = summary[pool]
        pop = s["population"]
        print(f"\n=== {pool} pool ===")
        print(f"Valid Validation population: {pop}")
        print(f"Equivalent (incl. identical): {s['equivalent']} / {pop} = {_fmt_pct(s['equivalent'], pop)}")
        print(f"Identical:                    {s['identical']} / {pop} = {_fmt_pct(s['identical'], pop)}")
        print("\nBy K-signature:")
        for k_sig in sorted(s["by_k"], key=lambda k: (len(k), k)):
            total, equivalent, identical = s["by_k"][k_sig]
            print(
                f"  {k_sig:<6} population={total:<5} "
                f"equivalent={equivalent} ({_fmt_pct(equivalent, total)})  "
                f"identical={identical} ({_fmt_pct(identical, total)})"
            )


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("results", nargs="?", type=Path, default=default_results_path())
    args = parser.parse_args()
    records = load_records(args.results)
    summary = summarize(records)
    print_report(summary)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
