"""Recompute every derived field from the stored measurements; check or repair.

The one-command reproducibility check for the database. Every field that is
a function of the measurements (see `rel_quest/derived.py` for the list and
the function that owns each) is recomputed with the canonical code and
compared with what is stored.

  --check    exit 1 if anything differs (CI runs this)
  --dry-run  report what would change
  (default)  rewrite the differing blocks in place, back up the old file,
             stamp metadata.last_updated, write an audit

Only blocks that differ are replaced, so an unchanged record is written back
byte-for-byte.

Runs that repaired something (each appended to the audit):
  2026-09-23  `performance_summary` made reproducible — empty on 62 records,
              stale on 34, returned verbatim by the REST API.
  2026-09-23  `source_quality`/`sigma_log_D` on 8 SAFE conditions, after
              classify_measurement_quality moved to reference.Locator.
  2026-10-23  `sigma_log_D` on 3,792 figure-prior cells: the figure prior is now
              measured (uncertainty.SIGMA_PRIORS, 0.033 -> 0.0715 dex).

Inputs / outputs:
  data/sx_ligands_merged_v14_flat.jsonl            (rewritten in place)
  data/archive/<stem>.derivedbak_<date>.jsonl      (backup, gitignored)
  outputs/derived_fields_recompute_audit.txt

Run:
  python scripts/recompute_derived_fields.py [--check | --dry-run] [--database PATH]
"""

from __future__ import annotations

import argparse
import collections
import json
import shutil
import sys
from datetime import date
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from rel_quest.derived import derived_differences, recompute_record  # noqa: E402

DB_PATH = ROOT / "data" / "sx_ligands_merged_v14_flat.jsonl"
AUDIT_PATH = ROOT / "outputs" / "derived_fields_recompute_audit.txt"
TODAY = date.today().isoformat()


def apply_fresh(stored: dict, fresh: dict, diffs: list) -> None:
    """Replace, in place, only the stored blocks that `diffs` names."""
    blocks = {d.split(".", 1)[0] if d.split(".", 1)[0] in (
        "selectivity", "data_quality") else d for d in diffs}
    for key in ("elements_measured", "performance_summary", "selectivity", "data_quality"):
        if key in blocks:
            stored[key] = fresh[key]
    if any("." in d and d.split(".", 1)[0] not in ("selectivity", "data_quality") for d in diffs):
        stored["experimental_data"] = fresh["experimental_data"]


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument("--check", action="store_true", help="exit 1 on any difference")
    mode.add_argument("--dry-run", action="store_true", help="report without writing")
    parser.add_argument("--database", type=Path, default=DB_PATH)
    args = parser.parse_args()

    lines = args.database.read_text(encoding="utf-8").splitlines()
    records = [json.loads(line) for line in lines if line.strip()]

    changed = []
    by_field = collections.Counter()
    for record in records:
        fresh = recompute_record(record)
        diffs = derived_differences(record, fresh)
        if not diffs:
            continue
        for diff in diffs:
            by_field[diff.split(".")[-1] if "." in diff else diff] += 1
        changed.append((record["ligand"]["Name"], diffs, record, fresh))

    report = [
        "Derived-field recompute",
        f"Date: {TODAY}",
        f"Database: {args.database}",
        "Code path: rel_quest.derived.recompute_record",
        "",
        f"Records: {len(records)}",
        f"Records with a derived field that does not reproduce: {len(changed)}",
    ]
    if by_field:
        report.append("Differences by field:")
        for field, n in by_field.most_common():
            report.append(f"  {field:<28} {n:>4}")
        report.append("")
        report.append("Per record:")
        for name, diffs, record, fresh in changed:
            report.append(f"  {name}")
            for diff in diffs[:8]:
                if diff == "performance_summary":
                    report.append(f"      performance_summary  {record.get('performance_summary')!r}")
                    report.append(f"                        -> {fresh['performance_summary']!r}")
                else:
                    report.append(f"      {diff}")
            if len(diffs) > 8:
                report.append(f"      ... and {len(diffs) - 8} more")

    if args.check:
        print("\n".join(report[:12 + len(by_field)]))
        if changed:
            print(f"\nFAIL: {len(changed)} record(s) do not reproduce "
                  f"(run without --check to see and repair)")
            return 1
        print("\nOK: every derived field reproduces from the stored measurements")
        return 0

    if args.dry_run or not changed:
        print("\n".join(report))
        if not changed:
            print("\nNothing to do.")
        return 0

    for name, diffs, record, fresh in changed:
        apply_fresh(record, fresh, diffs)
        record["metadata"]["last_updated"] = TODAY

    backup = ROOT / "data" / "archive" / f"{args.database.stem}.derivedbak_{TODAY}.jsonl"
    run = 2
    while backup.exists():  # a second run the same day keeps the first backup
        backup = backup.with_name(f"{args.database.stem}.derivedbak_{TODAY}_{run}.jsonl")
        run += 1
    backup.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(args.database, backup)
    with args.database.open("w", encoding="utf-8", newline="\n") as handle:
        for record in records:
            handle.write(json.dumps(record) + "\n")

    # Appended, not overwritten: each run that repairs something is a dated
    # entry in the same history.
    AUDIT_PATH.parent.mkdir(parents=True, exist_ok=True)
    has_history = AUDIT_PATH.exists() and AUDIT_PATH.stat().st_size > 0
    with AUDIT_PATH.open("a", encoding="utf-8") as handle:
        if has_history:
            handle.write("\n" + "=" * 72 + "\n\n")
        handle.write("\n".join(report) + "\n")
    print("\n".join(report[:12 + len(by_field)]))
    print(f"\nRewrote {len(changed)} records.")
    print(f"Backup:  {backup.relative_to(ROOT)}")
    print(f"Audit:   {AUDIT_PATH.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
