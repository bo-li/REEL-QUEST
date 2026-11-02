"""Check contribution files before they are applied (the CI job on a data pull request).

  python scripts/check_contribution.py [contributions/<paper>.json ...] [--database PATH]
         [--contributions-dir contributions] [--report-dir outputs/contributions]
         [--offline] [--summary]

Without files: every file in the contributions directory not yet applied. For
each file: the schema (schema/contribution.schema.json), the DOI (a bare DOI,
not a compilation's, resolving at doi.org unless --offline), the new-ligand
declarations against the database (a name or structure already there), each
condition's grammar, duplicates of stored conditions, then a dry-run apply
whose records must pass the validators and the JSON Schema. Nothing is written
to the database or the file. Writes <report-dir>/<stem>_check.md per file and
exits 1 on any problem; --summary also appends the reports to
$GITHUB_STEP_SUMMARY (the job summary on GitHub).
"""

from __future__ import annotations

import argparse
import os
import sys
from datetime import date
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from rel_quest import contributions as lib  # noqa: E402

DB = ROOT / "data" / "sx_ligands_merged_v14_flat.jsonl"
CONTRIBUTIONS = ROOT / "contributions"
REPORTS = ROOT / "outputs" / "contributions"


def doi_exists(doi: str):
    """True when doi.org redirects the DOI, False on 404, None when it cannot be told."""
    try:
        import requests
        r = requests.head(f"https://doi.org/{doi}", allow_redirects=False, timeout=10,
                          headers={"User-Agent": "reel-quest contribution check"})
    except Exception:          # noqa: BLE001 - any network failure means "could not tell"
        return None
    if r.status_code in (301, 302, 303, 307, 308):
        return True
    if r.status_code == 404:
        return False
    return None


def check_files(files, database, report_dir, offline=False, today=None):
    """[(path, CheckResult, report text)], reports written to report_dir, and the records every
    file would leave behind. Files are checked in order on the running state: a file that passes
    is applied in memory before the next is checked, so a later file may add conditions to a
    ligand an earlier one declares. Nothing is written to the database or the files."""
    records = lib.load_records(database)
    today = today or date.today().isoformat()
    report_dir = Path(report_dir)
    report_dir.mkdir(parents=True, exist_ok=True)
    out = []
    for path in files:
        path = Path(path)
        try:
            paper = lib.load_contribution(path)
        except (OSError, ValueError) as e:
            paper = None
            result = lib.CheckResult(path=str(path), problems=[f"cannot read the file: {e}"])
        else:
            result = lib.check_contribution(paper, records, today, path=str(path),
                                            doi_exists=None if offline else doi_exists)
            if result.ok:
                records = lib.apply_contribution(paper, records, today).records
        text = "\n".join(lib.report_lines(result, paper or {})) + "\n"
        (report_dir / f"{path.stem}_check.md").write_text(text, encoding="utf-8")
        out.append((path, result, text))
    return out, records


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("files", nargs="*", help="contribution files (default: every pending file)")
    parser.add_argument("--database", type=Path, default=DB)
    parser.add_argument("--contributions-dir", type=Path, default=CONTRIBUTIONS)
    parser.add_argument("--report-dir", type=Path, default=REPORTS)
    parser.add_argument("--offline", action="store_true", help="do not ask doi.org whether the DOI exists")
    parser.add_argument("--summary", action="store_true", help="append the reports to $GITHUB_STEP_SUMMARY")
    args = parser.parse_args(argv)
    files = [Path(f) for f in args.files] or lib.pending_files(args.contributions_dir)
    if not files:
        print(f"nothing to check: no pending file in {args.contributions_dir}")
        return 0
    results, _ = check_files(files, args.database, args.report_dir, offline=args.offline)
    failed = 0
    for path, result, text in results:
        verdict = "OK" if result.ok else f"{len(result.problems)} problem(s)"
        print(f"{path.name}: {verdict}" + (f", {len(result.warnings)} warning(s)" if result.warnings else ""))
        for p in result.problems:
            print(f"  - {p}")
        for w in result.warnings:
            print(f"  ? {w}")
        failed += not result.ok
        if args.summary and os.environ.get("GITHUB_STEP_SUMMARY"):
            with open(os.environ["GITHUB_STEP_SUMMARY"], "a", encoding="utf-8") as h:
                h.write(text + "\n")
    print(f"reports: {args.report_dir}")
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
