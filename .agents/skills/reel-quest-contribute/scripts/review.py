#!/usr/bin/env python3
"""Check a second read of a contribution file, and write its summary.

    python review.py check outputs/contributions/<doi-slug>_review.json

Refuses (exit 1, the problems on stdout) a review that misses a condition or
a new ligand of the file, names one the file has not got, uses a verdict or
confidence outside the vocabulary, gives no evidence, gives no correction for
a wrong value or condition, or was made on another version of the file (its
contribution_sha256). Otherwise writes <review>.md beside it: the summary to
post in the pull request.
"""

import argparse
import hashlib
import json
import sys
from pathlib import Path

VERDICTS = ("Correct", "Wrong value", "Wrong molecule", "Wrong condition", "Wrong acid",
            "Not measured here", "Not in paper", "Duplicate", "Unsure")
LIGAND_VERDICTS = ("Correct", "Wrong molecule", "Already in the database", "Unsure")
CONFIDENCE = ("high", "medium", "low")
NEEDS_CORRECTION = ("Wrong value", "Wrong condition", "Wrong acid", "Wrong molecule")


def repository_root() -> Path:
    for parent in Path(__file__).resolve().parents:
        if (parent / "rel_quest" / "__init__.py").is_file():
            return parent
    return Path.cwd()


def locate(path: str, review_path: Path) -> Path:
    for candidate in (Path(path), Path.cwd() / path, repository_root() / path):
        if candidate.is_file():
            return candidate
    raise FileNotFoundError(path)


def problems_of(review: dict, review_path: Path):
    problems = []
    try:
        path = locate(review.get("contribution", ""), review_path)
    except FileNotFoundError:
        return [f"contribution {review.get('contribution')!r} not found"], None
    paper = json.loads(path.read_text(encoding="utf-8"))
    if review.get("contribution_sha256") != hashlib.sha256(path.read_bytes()).hexdigest():
        problems.append("contribution_sha256 is not the file's: the review was made on another version")
    if (review.get("doi") or "").lower() != (paper.get("DOI") or "").lower():
        problems.append(f"doi {review.get('doi')!r} is not the file's {paper.get('DOI')!r}")
    if not review.get("pages_viewed"):
        problems.append("pages_viewed is empty: say which pages were read")

    def entries(key, ids, verdicts, label):
        seen = {}
        for item in review.get(key) or []:
            ident = item.get(label)
            if ident in seen:
                problems.append(f"{key}: {ident} reviewed twice")
            seen[ident] = item
            if ident not in ids:
                problems.append(f"{key}: {ident!r} is not in the file")
            if item.get("verdict") not in verdicts:
                problems.append(f"{key} {ident}: verdict {item.get('verdict')!r} is not one of {', '.join(verdicts)}")
            if item.get("confidence") not in CONFIDENCE:
                problems.append(f"{key} {ident}: confidence {item.get('confidence')!r} is not high, medium or low")
            if item.get("verdict") == "Unsure" and item.get("confidence") != "low":
                problems.append(f"{key} {ident}: an Unsure verdict has low confidence")
            if not str(item.get("evidence") or "").strip():
                problems.append(f"{key} {ident}: no evidence (page, figure or table, what was read)")
            if item.get("verdict") in NEEDS_CORRECTION and not (item.get("correction") or item.get("note")):
                problems.append(f"{key} {ident}: {item['verdict']} without a correction or note saying what the paper has")
        for ident in ids:
            if ident not in seen:
                problems.append(f"{key}: {ident} is not reviewed")
        return seen

    entries("conditions", [c["id"] for c in paper.get("conditions", [])], VERDICTS, "id")
    entries("new_ligands", [lig["Name"] for lig in paper.get("new_ligands", [])], LIGAND_VERDICTS, "Name")
    return problems, paper


def summary(review: dict, paper: dict) -> str:
    rows = review.get("conditions", [])
    counts = {}
    for row in rows:
        counts[row["verdict"]] = counts.get(row["verdict"], 0) + 1
    lines = [f"## Second read of `{Path(review['contribution']).name}`", "",
             f"Paper: https://doi.org/{paper['DOI']} - pages read: {', '.join(map(str, review['pages_viewed']))}", "",
             "Verdicts: " + ", ".join(f"{v} {n}" for v, n in sorted(counts.items())), ""]
    if review.get("paper_notes"):
        lines += [review["paper_notes"], ""]
    if paper.get("new_ligands"):
        lines += ["| New ligand | Verdict | Confidence | Evidence |", "|---|---|---|---|"]
        for item in review.get("new_ligands", []):
            lines.append(f"| {item['Name']} | {item['verdict']} | {item['confidence']} | {item['evidence']} |")
        lines.append("")
    lines += ["| Condition | Ligand | Verdict | Confidence | Evidence | Correction |", "|---|---|---|---|---|---|"]
    ligand = {c["id"]: c.get("ligand", "") for c in paper.get("conditions", [])}
    for row in rows:
        correction = json.dumps(row.get("correction") or {}, ensure_ascii=False) if row.get("correction") else ""
        note = f" {row['note']}" if row.get("note") else ""
        lines.append(f"| {row['id']} | {ligand.get(row['id'], '')} | {row['verdict']} | {row['confidence']} | "
                     f"{row['evidence']}{note} | {correction} |")
    lines += ["", "This is an AI second read. Before the file is applied, a person spot-checks it against the "
                  "paper (GOVERNANCE.md, the review rule) and adds `reviewed_by` to the file."]
    return "\n".join(lines) + "\n"


def cmd_check(args) -> int:
    path = Path(args.review)
    review = json.loads(path.read_text(encoding="utf-8"))
    problems, paper = problems_of(review, path)
    if problems:
        print(f"{path.name}: {len(problems)} problem(s)")
        for p in problems:
            print(f"  - {p}")
        return 1
    out = path.with_suffix(".md")
    out.write_text(summary(review, paper), encoding="utf-8")
    print(f"{path.name}: OK; summary at {out}")
    return 0


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0],
                                     formatter_class=argparse.RawDescriptionHelpFormatter, epilog=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    p = sub.add_parser("check", help="check a review file and write its summary")
    p.add_argument("review")
    p.set_defaults(func=cmd_check)
    args = parser.parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    raise SystemExit(main())
