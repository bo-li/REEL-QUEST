"""Every number the manuscript quotes, computed from the database.

The manuscript must not carry a hand-typed count: this repository has a long
history of numbers that went stale in prose (166 vs 170 ligands, Solvation
116 vs 120). `compute_stats()` derives each one from
data/sx_ligands_merged_v14_flat.jsonl, and from the validation outputs when
they exist; `paper/stats.json` is the committed copy the manuscript cites,
and `tests/test_paper_stats.py` fails when it no longer matches the data.

Run:
  python paper/stats.py            print and write paper/stats.json
  python paper/stats.py --check    exit 1 if paper/stats.json is stale
"""

from __future__ import annotations

import argparse
import collections
import contextlib
import io
import json
import sys
import zlib
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from rel_quest.schema import SCHEMA_VERSION, TREND_CLASSIFICATIONS  # noqa: E402
from rel_quest.trends import classify_selectivity_trend  # noqa: E402
from rel_quest.uncertainty import (  # noqa: E402
    QUALITY_PRIORS, SIGMA_PRIORS, TABLE_PRIOR_ON_D, condition_trend_points, figure_prior_evidence, is_censored,
    is_mixture)

DB_PATH = ROOT / "data" / "sx_ligands_merged_v14_flat.jsonl"
STATS_PATH = ROOT / "paper" / "stats.json"
VALIDATION = ROOT / "outputs" / "validation"
LANTHANIDES = ["La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy",
               "Ho", "Er", "Tm", "Yb", "Lu"]
PAIRS = [("Eu", "La"), ("Nd", "Pr"), ("Tb", "Dy"), ("Er", "Ho"), ("Yb", "Lu"), ("Dy", "Nd")]


def load_records():
    return [json.loads(line) for line in DB_PATH.open(encoding="utf-8") if line.strip()]


def counter(values):
    return dict(sorted(collections.Counter(values).items(), key=lambda kv: (-kv[1], str(kv[0]))))


def quartiles(values):
    if not values:
        return None
    q = np.percentile(values, [25, 50, 75])
    return {"q1": round(float(q[0]), 4), "median": round(float(q[1]), 4),
            "q3": round(float(q[2]), 4)}


def extremum_sharpness(record):
    """Largest sharpness of a Peaked/U_Shaped record's extremum, in dex.

    Sharpness = the extremum's log D minus the mean of its two neighbours
    (sign-flipped for U_Shaped), in the condition where it is largest. Log D
    moves by well under 0.3 dex between adjacent lanthanides in a smooth
    curve, so an extremum standing > 0.5 dex above its neighbours is one
    element off the curve — more likely a digitisation or transcription
    outlier than a peak. None for other trend classes.
    """
    trend = record["selectivity"]["trend_classification"]
    if trend not in ("Peaked", "U_Shaped"):
        return None
    sign = 1 if trend == "Peaked" else -1
    sharpest = None
    for condition in record["experimental_data"]:
        ys = [y * sign for _, y in sorted(condition_trend_points(condition))]
        if len(ys) < 4:  # the peak test's own minimum
            continue
        i = ys.index(max(ys))
        if 0 < i < len(ys) - 1:
            value = ys[i] - (ys[i - 1] + ys[i + 1]) / 2
            sharpest = value if sharpest is None else max(sharpest, value)
    return sharpest


SHARP_EXTREMUM_DEX = 0.5


def sharp_extremum(record, threshold=SHARP_EXTREMUM_DEX):
    """True if a curved label rests on one element far off its neighbours."""
    sharpness = extremum_sharpness(record)
    return sharpness is not None and sharpness > threshold


def duplicated_conditions(records):
    """Conditions whose measurements are identical to another ligand's.

    Groups conditions by (acid, [acid], [extractant], solvent, T, every
    element's log D); a group spanning two or more ligands is one set of
    measurements attributed to several molecules.
    """
    groups = collections.defaultdict(list)
    for r in records:
        for c in r["experimental_data"]:
            key = (c.get("Acid_Name"), c.get("Acid_Concentration_M"),
                   c.get("Extractant_Concentration_M"), c.get("Solvent_Name"),
                   c.get("Temperature_C"),
                   tuple(sorted((el, v.get("log_value")) for el, v in c["distribution_ratios"].items())))
            groups[key].append((r["ligand"]["Name"], c["condition_id"]))
    return [g for g in groups.values() if len({name for name, _ in g}) > 1]


def distinct_conditions(records):
    """(subtype, condition), counting each repeated measurement set once.

    A set of measurements stored on several molecules (one paper's values
    attributed to each of its structural variants) is one experiment. The key
    is duplicated_conditions' key within a subtype, so it also collapses a
    ligand's own repeated rows. Figures that count conditions use this.
    """
    seen, rows = set(), []
    for r in records:
        subtype = r["ligand"]["Subtype"]
        for c in r["experimental_data"]:
            key = (subtype, c.get("Acid_Name"), c.get("Acid_Concentration_M"),
                   c.get("Extractant_Concentration_M"), c.get("Solvent_Name"),
                   c.get("Temperature_C"),
                   tuple(sorted((el, v.get("log_value"))
                                for el, v in c["distribution_ratios"].items())))
            if key not in seen:
                seen.add(key)
                rows.append((subtype, c))
    return rows


def pooled_conditions(records):
    """distinct_conditions without the mixtures, for every number pooled across
    conditions: a mixed-extractant D (uncertainty.is_mixture) is not the ligand's alone."""
    return [(subtype, c) for subtype, c in distinct_conditions(records) if not is_mixture(c)]


MC_DRAWS = 200
MC_SEED = 20260929


def perturbed_draws(record, draws=MC_DRAWS, seed=MC_SEED):
    """Yield (conditions, label) with every measured log D moved by its sigma.

    Each draw adds N(0, sigma_log_D) to every measured cell (bounds and cells
    without a sigma stay put) and re-runs the canonical classifier. Seeded
    per ligand name, so every consumer sees the same draws.
    """
    rng = np.random.default_rng([seed, zlib.crc32(record["ligand"]["Name"].encode())])
    cells = [(i, el, v) for i, c in enumerate(record["experimental_data"])
             for el, v in c["distribution_ratios"].items()
             if isinstance(v, dict) and v.get("log_value") is not None and not is_censored(v)]
    sigma = np.array([v.get("sigma_log_D") or 0.0 for _, _, v in cells])
    base = np.array([v["log_value"] for _, _, v in cells])
    for noise in rng.normal(0.0, 1.0, size=(draws, len(cells))):
        conditions = [{**c, "distribution_ratios": dict(c["distribution_ratios"])}
                      for c in record["experimental_data"]]
        for (i, el, v), value in zip(cells, base + noise * sigma):
            conditions[i]["distribution_ratios"][el] = {**v, "log_value": float(value)}
        with contextlib.redirect_stdout(io.StringIO()):
            yield conditions, classify_selectivity_trend(conditions)["trend_classification"]


def perturbed_label_counts(record, draws=MC_DRAWS, seed=MC_SEED):
    """Counter of trend labels over the perturbed draws; None if nothing is fitted."""
    if record["selectivity"].get("slope") is None:
        return None
    return collections.Counter(label for _, label in perturbed_draws(record, draws, seed))


def label_stability(records):
    """{name: (stored label, share of draws that keep it, Counter)} for fitted ligands."""
    out = {}
    for r in records:
        counts = perturbed_label_counts(r)
        if counts:
            stored = r["selectivity"]["trend_classification"]
            out[r["ligand"]["Name"]] = (stored, counts[stored] / sum(counts.values()), counts)
    return out


def stability_summary(records, stability=None):
    stability = stability if stability is not None else label_stability(records)
    kept = [share for _, share, _ in stability.values()]
    return {"draws": MC_DRAWS, "n_fitted": len(kept),
            "label_kept_in_at_least_90pct": sum(1 for k in kept if k >= 0.9),
            "label_kept_in_50_to_90pct": sum(1 for k in kept if 0.5 <= k < 0.9),
            "label_kept_in_under_50pct": sum(1 for k in kept if k < 0.5)}


TABLE_PRIOR_FACTOR = 2  # the sensitivity check: the assumed table prior, doubled


def reported_sigma(records):
    """The experimental sigma papers give for the values they report.

    A measured sigma on a reported value is either the paper's own (a stated
    uncertainty, or the spread of its replicates) or the table prior, written in
    when verification replaced a cell with its printed value; only the former
    says anything about the prior.
    """
    found = [(c["reference"]["DOI"], cell["sigma_log_D"]) for r in records for c in r["experimental_data"]
             for cell in c["distribution_ratios"].values()
             if cell.get("sigma_source") == "measured" and cell.get("value_origin") == "reported"
             and not is_censored(cell) and cell.get("sigma_log_D") != SIGMA_PRIORS["table"]]
    return {"n_values": len(found), "dois": sorted({doi for doi, _ in found}),
            "median_pct_on_D": round(100 * float(np.median([s for _, s in found])) * np.log(10), 1)
            if found else None}


def has_table_prior(cell):
    """A measured cell whose sigma is the assumed table prior."""
    return (cell.get("source_quality") == "table" and cell.get("sigma_source") != "measured"
            and cell.get("log_value") is not None and not is_censored(cell))


def table_prior_sensitivity(records, stability, factor=TABLE_PRIOR_FACTOR):
    """The sigma-dependent numbers again, with the table prior multiplied by `factor`.

    Too few papers give an experimental sigma to estimate the table prior
    (reported_sigma), so the manuscript shows instead that its conclusions do not
    rest on it. Trend classes and slope intervals do not read
    sigma; the Monte Carlo label stability and the adjacent-pair 2-sigma counts
    do. The draws are seeded per ligand, so both runs see the same noise.
    """
    scale = lambda cell: {**cell, "sigma_log_D": cell["sigma_log_D"] * factor} if has_table_prior(cell) else cell  # noqa: E731
    scaled = [{**r, "experimental_data": [
        {**c, "distribution_ratios": {el: scale(cell) for el, cell in c["distribution_ratios"].items()}}
        for c in r["experimental_data"]]} for r in records]
    rescaled = label_stability(scaled)
    at_90 = lambda s, name: s[name][1] >= 0.9  # noqa: E731
    return {"factor": factor, "table_pct_on_D": round(100 * factor * TABLE_PRIOR_ON_D, 1),
            "n_values": sum(1 for r in records for c in r["experimental_data"]
                            for cell in c["distribution_ratios"].values() if has_table_prior(cell)),
            "trend_label_stability": stability_summary(scaled, rescaled),
            "ligands_crossing_90pct": sum(1 for name in stability if at_90(stability, name) != at_90(rescaled, name)),
            "adjacent_pair_log_sf": adjacent_sf_summary(scaled)}


ADJACENT = [(LANTHANIDES[i + 1], LANTHANIDES[i]) for i in range(len(LANTHANIDES) - 1)]


def measured_log(condition, element):
    """log D of one cell, or None when it is missing or only a bound."""
    cell = (condition.get("distribution_ratios") or {}).get(element)
    if not isinstance(cell, dict) or cell.get("log_value") is None or is_censored(cell):
        return None
    return cell["log_value"]


def showcase_groups(records, a="Eu", b="La", min_n=5):
    """Fig. 5: (subtypes with at least min_n conditions measuring both, ordered by median
    log SF; {subtype: [log D(a)]}; {subtype: [log SF(a/b)]}) over pooled conditions."""
    log_d, log_sf = collections.defaultdict(list), collections.defaultdict(list)
    for subtype, c in pooled_conditions(records):
        da, db = measured_log(c, a), measured_log(c, b)
        if da is not None:
            log_d[subtype].append(da)
            if db is not None:
                log_sf[subtype].append(da - db)
    kept = sorted((s for s in log_sf if len(log_sf[s]) >= min_n),
                  key=lambda s: (-np.median(log_sf[s]), s))
    return kept, log_d, log_sf


def log_d_by_subtype(records, element):
    """{subtype: [log D of element]} over distinct conditions measuring it, mixtures left out."""
    out = collections.defaultdict(list)
    for subtype, c in pooled_conditions(records):
        value = measured_log(c, element)
        if value is not None:
            out[subtype].append(value)
    return out


def eu_la_by_subtype(records):
    """Fig. 5 in numbers, in its row order: per subtype, conditions and medians of log D(Eu),
    log D(La) and log SF(Eu/La), and the median SF to two significant figures."""
    kept, log_eu, log_sf = showcase_groups(records)
    log_la = log_d_by_subtype(records, "La")
    return {s: {"n_sf": len(log_sf[s]), "median_log_sf": round(float(np.median(log_sf[s])), 2),
                "median_sf": float(f"{10 ** np.median(log_sf[s]):.2g}"),
                "n_log_d_eu": len(log_eu[s]), "median_log_d_eu": round(float(np.median(log_eu[s])), 2),
                "n_log_d_la": len(log_la[s]), "median_log_d_la": round(float(np.median(log_la[s])), 2)}
            for s in kept}


def adjacent_sf(condition):
    """[(pair, log SF, sigma of log SF)] for measured adjacent lanthanide pairs."""
    out = []
    ratios = condition["distribution_ratios"]
    for heavy, light in ADJACENT:
        a, b = ratios.get(heavy), ratios.get(light)
        if not (isinstance(a, dict) and isinstance(b, dict)) or is_censored(a) or is_censored(b):
            continue
        if a.get("log_value") is None or b.get("log_value") is None:
            continue
        sigma = float(np.hypot(a.get("sigma_log_D") or 0.0, b.get("sigma_log_D") or 0.0))
        out.append((f"{heavy}/{light}", a["log_value"] - b["log_value"], sigma))
    return out


def pair_sf_summary(records):
    """Lanthanide pairs within 2 sigma of unity by how far apart they are (1 = neighbors), and for Eu/La."""
    by_gap, eu_la = collections.defaultdict(lambda: [0, 0]), [0, 0]
    for _, c in pooled_conditions(records):
        ratios = c["distribution_ratios"]
        measured = [el for el in LANTHANIDES if isinstance(ratios.get(el), dict)
                    and ratios[el].get("log_value") is not None and not is_censored(ratios[el])]
        for i, light in enumerate(measured):
            for heavy in measured[i + 1:]:
                a, b = ratios[heavy], ratios[light]
                within = abs(a["log_value"] - b["log_value"]) < 2 * float(
                    np.hypot(a.get("sigma_log_D") or 0.0, b.get("sigma_log_D") or 0.0))
                for tally in [by_gap[LANTHANIDES.index(heavy) - LANTHANIDES.index(light)]] + \
                        ([eu_la] if {light, heavy} == {"Eu", "La"} else []):
                    tally[0] += 1
                    tally[1] += within
    as_dict = lambda t: {"n_values": t[0], "within_2_sigma": t[1]}  # noqa: E731
    return {"by_separation": {str(gap): as_dict(by_gap[gap]) for gap in sorted(by_gap)},
            "Eu/La": as_dict(eu_la)}


def adjacent_sf_summary(records):
    rows = [sf for _, c in pooled_conditions(records) for sf in adjacent_sf(c)]
    return {"n_values": len(rows),
            "sign_undetermined_at_2_sigma": sum(1 for _, v, s in rows if abs(v) < 2 * s)}


def lanthanides_measured(record):
    """Lanthanides with at least one measured (not censored) D; Pm never is, so at most 14."""
    return {el for c in record["experimental_data"] for el, v in c["distribution_ratios"].items()
            if el in LANTHANIDES and isinstance(v, dict) and v.get("log_value") is not None
            and not is_censored(v)}


def subtype_gaps(records, stability):
    """Per subtype, where more experiments would change what the database can say.

    Counts, not shares, so every denominator is explicit: ligands with no trend
    class (no condition with 4 lanthanides); the median number of lanthanides
    measured per ligand (of 14) and of conditions per ligand; fitted ligands
    whose class survives fewer than 90 % of the sigma draws (label_stability);
    adjacent-pair log SF values within 2 sigma of zero, over distinct
    conditions without mixtures (pooled_conditions; a value can sit there because the pair is barely separated or
    because sigma is large: either way its direction is not settled). Ordered
    as `subtypes` (most ligands first).
    """
    by_subtype = collections.defaultdict(list)
    for r in records:
        by_subtype[r["ligand"]["Subtype"]].append(r)
    sf_rows = collections.defaultdict(list)
    for subtype, c in pooled_conditions(records):
        sf_rows[subtype].extend(adjacent_sf(c))
    out = {}
    for subtype in counter(r["ligand"]["Subtype"] for r in records):
        rs = by_subtype[subtype]
        kept = [stability[r["ligand"]["Name"]][1] for r in rs if r["ligand"]["Name"] in stability]
        rows = sf_rows[subtype]
        out[subtype] = {
            "n_ligands": len(rs),
            "trend_classes": {t: n for t in TREND_CLASSIFICATIONS
                              if (n := sum(1 for r in rs
                                           if r["selectivity"]["trend_classification"] == t))},
            "no_trend_class": sum(1 for r in rs
                                  if r["selectivity"]["trend_classification"] == "Insufficient_Data"),
            "median_lanthanides_measured": float(np.median([len(lanthanides_measured(r))
                                                            for r in rs])),
            "median_conditions": float(np.median([len(r["experimental_data"]) for r in rs])),
            "n_fitted": len(kept),
            "fitted_class_kept_under_90pct": sum(1 for k in kept if k < 0.9),
            "n_ligands_with_adjacent_sf": sum(1 for r in rs if any(
                adjacent_sf(c) for c in r["experimental_data"] if not is_mixture(c))),
            "n_adjacent_sf": len(rows),
            "adjacent_sf_within_2_sigma": sum(1 for _, v, s in rows if abs(v) < 2 * s),
        }
    return out


def cell_provenance(conditions):
    """How each stored D was obtained and how its sigma was set (schema v0.3 fields)."""
    cells = [cell for _, c in conditions for cell in c["distribution_ratios"].values()
             if isinstance(cell, dict)]
    return {
        "value_origin": counter(cell.get("value_origin") or "unset" for cell in cells),
        "n_bounds": sum(1 for cell in cells if is_censored(cell)),
        "n_measured_sigma": sum(1 for cell in cells if cell.get("sigma_source") == "measured"),
        "n_prior_sigma": sum(1 for cell in cells
                             if cell.get("sigma_source") != "measured" and not is_censored(cell)),
        "n_mixture_conditions": sum(1 for _, c in conditions
                                    if (c.get("uncertainty") or {}).get("is_mixture")
                                    or c.get("Co_Extractants")),
    }


LEDGER = ROOT / "outputs" / "v14_apply_ledger.csv"
SPOT_CHECK = ROOT / "outputs" / "verification" / "spot_check_answers.csv"


def verification_summary():
    """The v14 check of every entry against its paper, from the tracked ledger and spot-check.

    One ledger row per verification-form entry (plus additions), with the verdict,
    who reached it and what was done. The spot-check settles the curator-estimated
    cells ('Estimates') and a random audit of the AI-check 'Correct' verdicts that
    decision 6 accepted ('Audit').
    """
    import csv
    from scipy.stats import beta

    if not LEDGER.exists():
        return None
    rows = list(csv.DictReader(LEDGER.open(encoding="utf-8")))

    def checker(value):
        if value.startswith("Curator decision"):
            return "curator decision"
        if value.startswith("Claude"):
            return "transcribed for curator review"
        return value or "unset"

    def action(value):
        for prefix in ("merged into", "recovered", "resolved"):
            if value.startswith(prefix):
                return prefix.split()[0]
        return value or "unset"

    summary = {
        "n_entries": len({row["entry_id"] for row in rows}),
        "verdicts": counter(row["verdict"] or "unset" for row in rows),
        "checked_by": counter(checker(row["checked_by"]) for row in rows),
        "actions": counter(action(row["action"]) for row in rows),
        "accepted_by_decision_6": sum(1 for row in rows if row["checked_by"] == "Curator decision 6"),
    }
    if SPOT_CHECK.exists():
        answers = list(csv.DictReader(SPOT_CHECK.open(encoding="utf-8-sig")))
        estimates = [a for a in answers if a["sheet"] == "Estimates" and a["answer"].strip()]
        audit = [a for a in answers if a["sheet"] == "Audit" and a["answer"].strip()]
        wrong = sum(1 for a in audit if a["answer"].strip() == "Wrong")
        summary["estimates_answered"] = len(estimates)
        summary["estimates_accepted"] = sum(1 for a in estimates if a["answer"].strip() == "accept")
        summary["audit_answered"] = len(audit)
        summary["audit_wrong"] = wrong
        summary["audit_error_upper_95"] = (round(float(beta.ppf(0.95, wrong + 1, len(audit) - wrong)), 4)
                                           if audit else None)
    return summary


# Keys that need the developer repository: the source labels (metadata.data_source,
# reference.Provenance_Level, condition-id prefixes, merge lineage) or the v14
# verification ledger. The public release (rel_quest.schema.PUBLIC_FIELDS) has
# neither, so its stats.json keeps these as computed here; `complete` carries them
# over and every other key is recomputed and compared.
DEV_ONLY_KEYS = ("provenance_level_conditions", "ligands_by_data_source", "conditions_by_origin",
                 "d_values_by_origin", "n_merged_records", "safe_conditions_without_concentrations",
                 "verification")


def has_source_labels(records):
    return any((r.get("metadata") or {}).get("data_source") for r in records)


def complete(fresh, stored):
    """`fresh` with the DEV_ONLY_KEYS it could not compute taken from `stored`, in stored order."""
    merged = {}
    for key, value in (stored or {}).items():
        if key in fresh:
            merged[key] = fresh[key]
        elif key in DEV_ONLY_KEYS:
            merged[key] = value
    for key, value in fresh.items():
        merged.setdefault(key, value)
    return merged


def load_stored():
    return json.loads(STATS_PATH.read_text(encoding="utf-8")) if STATS_PATH.exists() else None


def compute_stats(records=None, stability=None):
    """stability: label_stability(records), when the caller has it (it is the slow part)."""
    records = records if records is not None else load_records()
    stability = stability if stability is not None else label_stability(records)
    conditions = [(r, c) for r in records for c in r["experimental_data"]]
    cells = [(r, c, el, cell) for r, c in conditions
             for el, cell in c["distribution_ratios"].items()]
    origin = lambda c: c["condition_id"].split("_", 1)[0]  # noqa: E731
    labelled = has_source_labels(records)
    fitted = [r for r in records if r["selectivity"].get("slope") is not None]
    ci_width = [r["selectivity"]["slope_CI_upper"] - r["selectivity"]["slope_CI_lower"]
                for r in fitted if r["selectivity"].get("slope_CI_lower") is not None]
    references = [c["reference"] for _, c in conditions]
    n_cond = collections.Counter(len(r["experimental_data"]) for r in records)

    def both_measured(a, b):
        return sum(1 for _, c in conditions
                   if (c["distribution_ratios"].get(a) or {}).get("log_value") is not None
                   and (c["distribution_ratios"].get(b) or {}).get("log_value") is not None)

    stats = {
        "schema_version": SCHEMA_VERSION,
        "n_ligands": len(records),
        "n_conditions": len(conditions),
        "n_d_values": len(cells),
        "n_distinct_primary_dois": len({ref["DOI"] for ref in references if ref.get("DOI")}),
        "n_conditions_url_only": sum(1 for ref in references if not ref.get("DOI")),
        "n_conditions_with_locator": sum(1 for ref in references if ref.get("Locator")),
        **({"provenance_level_conditions": counter(ref["Provenance_Level"] for ref in references),
            "ligands_by_data_source": counter(r["metadata"]["data_source"] for r in records),
            "conditions_by_origin": counter(origin(c) for _, c in conditions),
            "d_values_by_origin": counter(origin(c) for _, c, _, _ in cells)} if labelled else {}),
        "d_values_by_source_quality": counter(cell["source_quality"] for *_, cell in cells),
        "d_values_by_element": {el: sum(1 for *_, e, _ in cells if e == el)
                                for el in LANTHANIDES + ["Y", "Sc"]},
        "ligands_measuring_element": {el: sum(1 for r in records if el in r["elements_measured"])
                                      for el in LANTHANIDES + ["Y", "Sc"]},
        "mechanisms": counter(r["ligand"]["Mechanism"] for r in records),
        "subtypes": counter(r["ligand"]["Subtype"] for r in records),
        "n_subtypes_present": len({r["ligand"]["Subtype"] for r in records}),
        "trend_classes": {t: sum(1 for r in records
                                 if r["selectivity"]["trend_classification"] == t)
                          for t in TREND_CLASSIFICATIONS},
        "trend_confidence": counter(r["selectivity"]["trend_confidence"] for r in records),
        "n_ligands_with_fitted_slope": len(fitted),
        "slope_ci_width_dex_per_Z": quartiles(ci_width),
        "ligands_by_condition_count": {"1": n_cond.get(1, 0), "2-4": sum(v for k, v in n_cond.items() if 2 <= k <= 4),
                                       "5+": sum(v for k, v in n_cond.items() if k >= 5)},
        "conditions_with_both": {f"{a}/{b}": both_measured(a, b) for a, b in PAIRS},
        "n_distinct_conditions": len(distinct_conditions(records)),
        "trend_label_stability": stability_summary(records, stability),
        "adjacent_pair_log_sf": adjacent_sf_summary(records),
        "pair_log_sf": pair_sf_summary(records),
        **({"n_merged_records": sum(1 for r in records if r["metadata"].get("merged_from"))}
           if labelled else {}),
        "identifier_coverage": {
            "SMILES": sum(1 for r in records if r["ligand"].get("SMILES")),
            "InChI_Key": sum(1 for r in records if r["ligand"].get("InChI_Key")),
            "IUPAC_Name": sum(1 for r in records if r["ligand"].get("IUPAC_Name")),
            "PubChem_CID": sum(1 for r in records if r["ligand"].get("PubChem_CID")),
            "CAS_Number": sum(1 for r in records if r["ligand"].get("CAS_Number")),
        },
        "quality_grade": counter(r["data_quality"]["quality_grade"] for r in records),
        "n_outlier_flags": sum(len(r["data_quality"]["outlier_flags"]) for r in records),
        "curved_trends_on_a_sharp_extremum": {
            "threshold_dex": SHARP_EXTREMUM_DEX,
            "Peaked": sum(1 for r in records if sharp_extremum(r)
                          and r["selectivity"]["trend_classification"] == "Peaked"),
            "U_Shaped": sum(1 for r in records if sharp_extremum(r)
                            and r["selectivity"]["trend_classification"] == "U_Shaped"),
        },
        "duplicated_condition_groups": len(duplicated_conditions(records)),
        "conditions_in_duplicated_groups": sum(len(g) for g in duplicated_conditions(records)),
        "ligands_in_duplicated_groups": len({n for g in duplicated_conditions(records) for n, _ in g}),
        **({"safe_conditions_without_concentrations": sum(
            1 for _, c in conditions if origin(c) == "SAFE"
            and c.get("Acid_Concentration_M") is None and c.get("Extractant_Concentration_M") is None)}
           if labelled else {}),
        "by_subtype": subtype_gaps(records, stability),
        "eu_la_by_subtype": eu_la_by_subtype(records),
        "cell_provenance": cell_provenance(conditions),
        "sigma_priors": {"figure_dex": SIGMA_PRIORS["figure"], "table_dex": SIGMA_PRIORS["table"],
                         "figure_pct_on_D": round(100 * QUALITY_PRIORS["figure"], 2),
                         "table_pct_on_D": round(100 * TABLE_PRIOR_ON_D, 1), "figure_evidence": figure_prior_evidence(records)},
        "reported_sigma": reported_sigma(records),
        "table_prior_sensitivity": table_prior_sensitivity(records, stability),
    }
    verification = verification_summary()
    if verification is not None:
        stats["verification"] = verification
    for name in ("revalidation_summary", "cross_source_agreement"):
        path = VALIDATION / f"{name}.json"
        if path.exists():
            stats[name] = json.loads(path.read_text(encoding="utf-8"))
    return stats


def render(stats) -> str:
    return json.dumps(stats, indent=2, ensure_ascii=False) + "\n"


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    text = render(complete(compute_stats(), load_stored()))
    if args.check:
        if not STATS_PATH.exists() or STATS_PATH.read_text(encoding="utf-8") != text:
            print("paper/stats.json is stale — run python paper/stats.py")
            return 1
        print("paper/stats.json is current")
        return 0
    STATS_PATH.write_text(text, encoding="utf-8")
    print(text)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
