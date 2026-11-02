"""Technical validation of the source-transcribed values.

Derived fields are checked by recomputation (rel_quest/derived.py). The
transcribed values — each D read from a paper's figure or table — can only be
checked against the papers. This module supports the two checks reported in
the paper's Technical Validation section:

1. **Re-digitisation audit.** A seeded, stratified sample of conditions is
   drawn (`draw_sample`); a curator re-reads every cell from the source PDF;
   `score_sample` reports the error by stratum and compares it with the
   `sigma_log_D` prior the database assigns.

2. **Cross-source agreement.** Where two compilations transcribed the same
   ligand under the same conditions (`cross_source_pairs`), their log D
   values are compared. Pairs citing the same paper measure transcription
   agreement; pairs citing different papers measure inter-laboratory spread.

Scripts: scripts/validation/draw_revalidation_sample.py,
scripts/validation/score_revalidation.py,
scripts/validation/cross_source_agreement.py.
"""

from __future__ import annotations

import math
import os
import re
from pathlib import Path
from typing import Callable, Dict, Iterable, List, Optional

import numpy as np

ROOT = Path(__file__).resolve().parents[1]

# Conditions per stratum. Fixed rather than proportional so the small strata
# (tabulated SAFE values, the Pramanik control) are not starved: 38
# conditions, about 5.7 % of the database, roughly 250 cells.
DEFAULT_ALLOCATION = {
    ("Liu2022", "figure"): 12,
    ("SAFE", "figure"): 10,
    ("SAFE", "table"): 6,
    ("IDEAL", "figure"): 8,
    ("Pramanik2024", "table"): 2,   # already verified cell by cell: a control
}
DEFAULT_SEED = 20260922

SAMPLE_COLUMNS = [
    "sample_id", "ligand", "condition_id", "stratum", "doi", "pdf_path", "locator",
    "element", "stored_value", "stored_log_value", "sigma_log_D",
    "redigitised_value", "redigitised_log_value", "checker", "date", "note",
]


def stratum(condition: Dict) -> tuple:
    """(condition_id prefix, source_quality) — the sampling stratum."""
    prefix = condition["condition_id"].split("_", 1)[0]
    cells = condition.get("distribution_ratios") or {}
    quality = next(iter(cells.values()), {}).get("source_quality", "unknown")
    return prefix, quality


def doi_file_key(text: str) -> str:
    """The letters and digits of a DOI or file name, lower-cased.

    A PDF saved by hand keeps whatever the browser or the person made of the
    DOI's punctuation: `10.1016_s1002-0721(14)60437-3.pdf` keeps the brackets
    the downloader's naming turns into `_`, and `10.1088_1757-899x91012068.pdf`
    drops the slashes of 10.1088/1757-899x/9/1/012068. Only the letters and
    digits survive every one of those, so files are matched on them.
    """
    return re.sub(r"[^a-z0-9]", "", text.lower())


def downloaded_pdf_finder(pdf_dir: Path) -> Callable[[Optional[str]], Optional[str]]:
    """Return a function mapping a DOI to the name of its PDF in `pdf_dir` (or None).

    A file matches when its name, less `.pdf`, has the DOI's letters and
    digits (`doi_file_key`); failing that, when it is the DOI in the
    downloader's naming (`<DOI, runs of other characters as _>`) plus a copy
    suffix such as `_2` or ` (1)` — not a supporting-information file, and not
    a longer DOI that merely starts with this one.
    """
    names = sorted(os.listdir(pdf_dir)) if Path(pdf_dir).is_dir() else []
    by_key = {}
    for name in names:
        if name.lower().endswith(".pdf"):
            by_key.setdefault(doi_file_key(name[:-4]), name)

    def find(doi):
        if not doi:
            return None
        if doi_file_key(doi) in by_key:
            return by_key[doi_file_key(doi)]
        stem = re.sub(r"[^A-Za-z0-9._-]+", "_", doi).lower()
        for name in names:
            rest = name.lower()[len(stem):]
            if name.lower().startswith(stem) and (
                    rest == ".pdf" or (rest[:1] in ("_", " ") and "_si" not in rest)):
                return name
        return None

    return find


def pdf_lookup(root: Path = ROOT) -> Callable[[Optional[str]], Optional[str]]:
    """Return a function mapping a DOI to a local PDF path (or None).

    Uses data/doi_to_pdf_index.json, then papers/downloaded/ through
    `downloaded_pdf_finder`, which also finds PDFs saved by hand.
    """
    import json
    from .references import normalize_doi

    index_path = root / "data" / "doi_to_pdf_index.json"
    index = {}
    if index_path.exists():
        for doi, paths in json.loads(index_path.read_text(encoding="utf-8")).items():
            key = normalize_doi(doi)
            if key:
                index[key] = paths
    find_downloaded = downloaded_pdf_finder(root / "papers" / "downloaded")

    def lookup(doi):
        if not doi:
            return None
        for path in index.get(doi, []):
            if (root / path).exists():
                return path
        name = find_downloaded(doi)
        return f"papers/downloaded/{name}" if name else None

    return lookup


def draw_sample(records: List[Dict], allocation: Dict = None, seed: int = DEFAULT_SEED,
                lookup: Callable = None) -> Dict:
    """Seeded stratified sample of conditions that have a local PDF.

    Returns {"rows": [...one per cell...], "eligible": {stratum: n},
    "population": {stratum: n}, "drawn": {stratum: n}}.
    """
    allocation = allocation or DEFAULT_ALLOCATION
    lookup = lookup or (lambda doi: None)
    population, eligible = {}, {}
    for record in records:
        for condition in record["experimental_data"]:
            key = stratum(condition)
            population[key] = population.get(key, 0) + 1
            pdf = lookup(condition["reference"].get("DOI"))
            # a bound has no value to re-digitise; a condition of bounds only gives no row
            if pdf and any(not cell.get("censor") for cell in condition["distribution_ratios"].values()):
                eligible.setdefault(key, []).append((record, condition, pdf))

    rng = np.random.default_rng(seed)
    rows, drawn = [], {}
    for key in sorted(allocation):
        pool = sorted(eligible.get(key, []), key=lambda item: item[1]["condition_id"])
        n = min(allocation[key], len(pool))
        picks = sorted(rng.choice(len(pool), size=n, replace=False)) if n else []
        drawn[key] = n
        for pick in picks:
            record, condition, pdf = pool[pick]
            reference = condition["reference"]
            for element, cell in condition["distribution_ratios"].items():
                if cell.get("censor"):
                    continue
                rows.append({
                    "sample_id": f"S{len(rows) + 1:04d}",
                    "ligand": record["ligand"]["Name"],
                    "condition_id": condition["condition_id"],
                    "stratum": f"{key[0]}/{key[1]}",
                    "doi": reference.get("DOI"),
                    "pdf_path": pdf,
                    "locator": reference.get("Locator") or "",
                    "element": element,
                    "stored_value": cell.get("value"),
                    "stored_log_value": cell.get("log_value"),
                    "sigma_log_D": cell.get("sigma_log_D"),
                    "redigitised_value": "",
                    "redigitised_log_value": "",
                    "checker": "",
                    "date": "",
                    "note": "",
                })
    return {"rows": rows,
            "population": {f"{k[0]}/{k[1]}": v for k, v in sorted(population.items())},
            "eligible": {f"{k[0]}/{k[1]}": len(v) for k, v in sorted(eligible.items())},
            "drawn": {f"{k[0]}/{k[1]}": v for k, v in sorted(drawn.items())}}


def _as_float(text) -> Optional[float]:
    try:
        value = float(text)
    except (TypeError, ValueError):
        return None
    return value if math.isfinite(value) else None


def score_sample(rows: Iterable[Dict]) -> Dict:
    """Error statistics of a filled re-digitisation sheet, by stratum.

    A row counts when it has a stored log D and a re-digitised value (either
    `redigitised_log_value` or a positive `redigitised_value`).
    """
    by_stratum: Dict[str, List[tuple]] = {}
    for row in rows:
        stored = _as_float(row.get("stored_log_value"))
        fresh = _as_float(row.get("redigitised_log_value"))
        if fresh is None:
            value = _as_float(row.get("redigitised_value"))
            fresh = math.log10(value) if value and value > 0 else None
        if stored is None or fresh is None:
            continue
        prior = _as_float(row.get("sigma_log_D"))
        by_stratum.setdefault(row["stratum"], []).append((fresh - stored, prior))

    def summarise(pairs):
        errors = np.array([abs(d) for d, _ in pairs])
        signed = np.array([d for d, _ in pairs])
        priors = [p for _, p in pairs if p is not None]
        return {
            "n_cells": int(len(errors)),
            "median_abs_dex": round(float(np.median(errors)), 4),
            "p90_abs_dex": round(float(np.percentile(errors, 90)), 4),
            "rms_dex": round(float(np.sqrt(np.mean(signed ** 2))), 4),
            "mean_signed_dex": round(float(np.mean(signed)), 4),
            "within_0.05_dex": round(float(np.mean(errors <= 0.05)), 3),
            "within_0.10_dex": round(float(np.mean(errors <= 0.10)), 3),
            "gross_errors_gt_0.3_dex": int(np.sum(errors > 0.3)),
            "prior_sigma_dex": round(float(np.median(priors)), 4) if priors else None,
        }

    summary = {name: summarise(pairs) for name, pairs in sorted(by_stratum.items())}
    everything = [pair for pairs in by_stratum.values() for pair in pairs]
    if everything:
        summary["ALL"] = summarise(everything)
    return summary


def _close(a, b, rel=0.02, abs_tol=1e-9):
    if a is None or b is None:
        return False
    return math.isclose(float(a), float(b), rel_tol=rel, abs_tol=abs_tol)


def _norm(text) -> str:
    return re.sub(r"\s+", " ", str(text or "")).strip().lower()


def _norm_solvent(text) -> str:
    """Diluent name without parenthetical aliases: 'phenyl trifluoromethyl
    sulfone (FS-13)' and 'phenyl trifluoromethyl sulfone' are one diluent."""
    return _norm(re.sub(r"\([^)]*\)", "", str(text or "")))


def condition_key_matches(a: Dict, b: Dict, rel: float = 0.02, temp_tol: float = 2.0) -> bool:
    """Strict same-condition test: acid, [acid], [extractant], diluent, T, modifier,
    aqueous additives; a mixture (`Co_Extractants`) never matches.

    Both concentrations must be present: two conditions with an unknown
    concentration are not known to match. A condition with a co-extractant or
    an aqueous complexant measures a different system (v14's T2EHDGA + HDEHP
    with HEDTA paired with T2EHDGA alone before 2026-10-19). (All SAFE conditions currently lack
    them — see scripts/audit_safe_condition_grouping.py — so SAFE rows cannot
    pair until that is repaired.)
    """
    return (_norm(a.get("Acid_Name")) == _norm(b.get("Acid_Name"))
            and _close(a.get("Acid_Concentration_M"), b.get("Acid_Concentration_M"), rel)
            and _close(a.get("Extractant_Concentration_M"), b.get("Extractant_Concentration_M"), rel)
            and _norm_solvent(a.get("Solvent_Name")) == _norm_solvent(b.get("Solvent_Name"))
            and a.get("Temperature_C") is not None and b.get("Temperature_C") is not None
            and abs(float(a["Temperature_C"]) - float(b["Temperature_C"])) <= temp_tol
            and _norm(a.get("Phase_Modifier_Name")) == _norm(b.get("Phase_Modifier_Name"))
            and _norm(a.get("Additives")) == _norm(b.get("Additives"))
            and not a.get("Co_Extractants") and not b.get("Co_Extractants"))


def cross_source_pairs(records: List[Dict], temp_tol: float = 2.0) -> List[Dict]:
    """Cells measured for the same ligand, element and condition by two sources.

    Returns one dict per paired cell with the two log D values, their
    difference, and whether both cite the same primary paper. `temp_tol`
    (degC) is loosened for the relaxed tier because the Liu SI often records
    298 K where the paper reports room temperature.
    """
    pairs = []
    for record in records:
        conditions = record["experimental_data"]
        for i, a in enumerate(conditions):
            source_a = a["condition_id"].split("_", 1)[0]
            for b in conditions[i + 1:]:
                source_b = b["condition_id"].split("_", 1)[0]
                if source_a == source_b or not condition_key_matches(a, b, temp_tol=temp_tol):
                    continue
                doi_a, doi_b = a["reference"].get("DOI"), b["reference"].get("DOI")
                for element in sorted(set(a["distribution_ratios"]) & set(b["distribution_ratios"])):
                    log_a = a["distribution_ratios"][element].get("log_value")
                    log_b = b["distribution_ratios"][element].get("log_value")
                    if log_a is None or log_b is None:
                        continue
                    pairs.append({
                        "ligand": record["ligand"]["Name"],
                        "element": element,
                        "condition_a": a["condition_id"],
                        "condition_b": b["condition_id"],
                        "source_a": source_a,
                        "source_b": source_b,
                        "doi_a": doi_a,
                        "doi_b": doi_b,
                        "same_paper": bool(doi_a) and doi_a == doi_b,
                        "log_D_a": log_a,
                        "log_D_b": log_b,
                        "delta_dex": round(log_b - log_a, 4),
                    })
    return pairs


def summarise_pairs(pairs: List[Dict]) -> Dict:
    """Agreement statistics, split by same-paper vs different-paper pairs."""
    out = {}
    for label, subset in (("same_paper", [p for p in pairs if p["same_paper"]]),
                          ("different_paper", [p for p in pairs if not p["same_paper"]]),
                          ("all", pairs)):
        if not subset:
            out[label] = {"n_cells": 0}
            continue
        deltas = np.array([p["delta_dex"] for p in subset])
        errors = np.abs(deltas)
        out[label] = {
            "n_cells": int(len(deltas)),
            "n_ligands": len({p["ligand"] for p in subset}),
            "median_abs_dex": round(float(np.median(errors)), 4),
            "p90_abs_dex": round(float(np.percentile(errors, 90)), 4),
            "max_abs_dex": round(float(np.max(errors)), 4),
            "bland_altman_bias_dex": round(float(np.mean(deltas)), 4),
            "bland_altman_loa_dex": round(float(1.96 * np.std(deltas, ddof=1)), 4)
            if len(deltas) > 1 else None,
        }
    return out
