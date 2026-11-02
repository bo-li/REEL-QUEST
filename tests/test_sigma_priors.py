"""The sigma(log D) prior for a figure value is measured, not assumed.

A value without a measured sigma gets a prior by how it was read
(rel_quest.uncertainty.SIGMA_PRIORS). The figure prior is the median, over
calibrated figures, of the reading uncertainty measured on them
(uncertainty.figure_prior_evidence); the table prior stays an assumption,
since too few papers report an experimental sigma to measure one.

Ways this goes wrong, each pinned below:
  1. the prior in the code drifts from the evidence in the canonical data;
  2. the evidence admits a sigma that is not a figure reading: a prior, a
     tabulated or interpolated value, a bound, or a reading clamped to the
     vector-figure floor, which says "at least the floor";
  3. one heavily read figure outweighs the rest (the prior stands in for an
     uncalibrated figure, so one figure is one observation);
  4. the code that writes new cells uses another constant than the one
     recompute_derived_fields applies;
  5. the manuscript types the numbers instead of citing paper/stats.json.
"""

import json
import unittest
from pathlib import Path

from rel_quest import contributions
from rel_quest.uncertainty import SIGMA_PRIORS, VECTOR_FIGURE_FLOOR, figure_prior_evidence

ROOT = Path(__file__).resolve().parent.parent
FLOOR = VECTOR_FIGURE_FLOOR


def records():
    path = ROOT / "data" / "sx_ligands_merged_v14_flat.jsonl"
    return [json.loads(line) for line in path.open(encoding="utf-8") if line.strip()]


def cell(sigma, origin="digitised", source="measured", **extra):
    return {"value": 1.0, "log_value": 0.0, "sigma_log_D": sigma, "value_origin": origin,
            "sigma_source": source, **extra}


def record(*conditions):
    return {"experimental_data": [
        {"reference": {"DOI": doi, "Locator": locator}, "distribution_ratios": cells}
        for doi, locator, cells in conditions]}


class TestFigurePriorEvidence(unittest.TestCase):

    def test_prior_matches_the_canonical_data(self):
        evidence = figure_prior_evidence(records())
        self.assertGreaterEqual(evidence["n_figures"], 30)
        self.assertAlmostEqual(SIGMA_PRIORS["figure"], evidence["median_dex"], delta=0.005,
                               msg="the figure prior no longer matches the calibrated readings: "
                                   "set SIGMA_PRIORS['figure'] and run recompute_derived_fields.py")

    def test_only_unclamped_figure_readings_count(self):
        ev = figure_prior_evidence([record(("10.1/a", "Figure 1", {
            "La": cell(0.05),
            "Ce": cell(0.9, source="prior"),
            "Pr": cell(0.9, origin="reported"),
            "Nd": cell(0.9, origin="interpolated"),
            "Sm": cell(0.9, censor="<"),
            "Eu": cell(FLOOR),
        }))])
        self.assertEqual((ev["n_cells"], ev["n_figures"], ev["n_floor_excluded"]), (1, 1, 1))
        self.assertEqual(ev["median_dex"], 0.05)

    def test_each_figure_counts_once(self):
        heavy = {el: cell(0.2) for el in ("La", "Ce", "Pr", "Nd", "Sm", "Eu", "Gd", "Tb", "Dy")}
        ev = figure_prior_evidence([record(("10.1/a", "Figure 1", heavy),
                                           ("10.1/a", "Figure 2", {"La": cell(0.05)}),
                                           ("10.1/b", "Figure 1", {"La": cell(0.06)}))])
        self.assertEqual(ev["n_figures"], 3)
        self.assertEqual(ev["median_dex"], 0.06)

    def test_figure_without_readings_is_not_counted(self):
        ev = figure_prior_evidence([record(("10.1/a", "Figure 1", {"La": cell(FLOOR)}),
                                           ("10.1/b", "Figure 1", {"La": cell(0.07)}))])
        self.assertEqual((ev["n_figures"], ev["n_figures_calibrated"], ev["median_dex"]), (1, 2, 0.07))


class TestPriorsAreOneConstant(unittest.TestCase):

    def test_writers_use_the_derived_priors(self):
        self.assertEqual(contributions.FIGURE_PRIOR, SIGMA_PRIORS["figure"])
        self.assertEqual(contributions.TABLE_PRIOR, SIGMA_PRIORS["table"])

    def test_stored_prior_cells_carry_the_prior(self):
        seen = set()
        for r in records():
            for c in r["experimental_data"]:
                for v in c["distribution_ratios"].values():
                    if isinstance(v, dict) and v.get("sigma_source") != "measured" and v.get("source_quality"):
                        seen.add((v["source_quality"], v["sigma_log_D"]))
        self.assertEqual(seen, {(q, SIGMA_PRIORS[q]) for q, _ in seen})

    def test_manuscript_cites_stats_not_typed_numbers(self):
        path = ROOT / "paper" / "manuscript_draft.md"
        if not path.exists():
            self.skipTest("public release: no manuscript template")
        text = path.read_text(encoding="utf-8")
        self.assertIn("{{sigma_priors.figure_dex}}", text)
        self.assertNotIn("7.5%", text)
