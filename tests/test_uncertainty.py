"""Pin the slope-CI estimator to the slope it describes.

`selectivity.slope` is the mean of per-condition regression slopes
(`trends.classify_selectivity_trend`). Until 2026-08-21
`uncertainty.bootstrap_trend_ci` pooled every point from every condition
into one regression, a different estimator, and eight records carried a
slope outside its own interval. These tests hold the two together.
"""

import json
import unittest
from pathlib import Path

from rel_quest.trends import classify_selectivity_trend
from rel_quest.uncertainty import (
    MIN_POINTS_PER_CONDITION,
    bootstrap_trend_ci,
    detect_outliers,
    per_condition_slopes,
)

ROOT = Path(__file__).resolve().parents[1]
DB_PATH = ROOT / "data" / "sx_ligands_merged_v14_flat.jsonl"

LN = ["La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu"]


def condition(slope, intercept, elements=LN, condition_id="c", noise=0.0):
    """A synthetic condition whose log D is a line in Z, optionally with a
    fixed alternating +/-noise so the fit is not exact."""
    ratios = {}
    for i, el in enumerate(elements):
        log_d = intercept + slope * (57 + LN.index(el)) + noise * (-1) ** i
        ratios[el] = {"value": 10 ** log_d, "log_value": log_d}
    return {"condition_id": condition_id, "distribution_ratios": ratios}


class SyntheticCases(unittest.TestCase):

    def test_two_conditions_interval_spans_their_slopes_and_holds_the_mean(self):
        data = [condition(0.10, 0.0, condition_id="a"),
                condition(0.20, 3.0, condition_id="b")]  # same slopes, shifted
        ci = bootstrap_trend_ci(data)
        self.assertIsNotNone(ci)
        self.assertAlmostEqual(ci["slope_CI_lower"], 0.10, places=3)
        self.assertAlmostEqual(ci["slope_CI_upper"], 0.20, places=3)
        self.assertTrue(ci["slope_CI_lower"] <= 0.15 <= ci["slope_CI_upper"])
        # The pooled estimator would have been dragged by the 3-dex offset
        # between the two conditions; the condition bootstrap is not.
        self.assertGreater(ci["slope_std"], 0)

    def test_identical_conditions_give_a_zero_width_interval(self):
        data = [condition(0.10, 0.0, condition_id=str(i)) for i in range(5)]
        ci = bootstrap_trend_ci(data)
        self.assertAlmostEqual(ci["slope_CI_lower"], 0.10, places=6)
        self.assertAlmostEqual(ci["slope_CI_upper"], 0.10, places=6)

    def test_single_condition_reduces_to_its_own_case_bootstrap(self):
        data = [condition(0.05, 0.0, noise=0.05)]
        fits = per_condition_slopes(data)
        self.assertEqual(len(fits), 1)
        ci = bootstrap_trend_ci(data)
        self.assertIsNotNone(ci)
        self.assertGreater(ci["slope_CI_upper"], ci["slope_CI_lower"])
        self.assertTrue(ci["slope_CI_lower"] <= fits[0][0] <= ci["slope_CI_upper"])

    def test_conditions_below_the_minimum_are_ignored_like_the_classifier(self):
        too_few = condition(0.30, 0.0, elements=LN[:MIN_POINTS_PER_CONDITION - 1],
                            condition_id="short")
        self.assertEqual(per_condition_slopes([too_few]), [])
        self.assertIsNone(bootstrap_trend_ci([too_few]))
        # ...and adding it to a real condition does not move the interval.
        full = condition(0.10, 0.0, condition_id="full")
        self.assertEqual(bootstrap_trend_ci([full]), bootstrap_trend_ci([full, too_few]))

    def test_y_and_sc_are_not_trend_points(self):
        data = [condition(0.10, 0.0)]
        data[0]["distribution_ratios"]["Y"] = {"value": 1e6, "log_value": 6.0}
        data[0]["distribution_ratios"]["Sc"] = {"value": 1e-6, "log_value": -6.0}
        self.assertEqual(len(per_condition_slopes(data)[0][1]), len(LN))

    def test_seeded_and_reproducible(self):
        data = [condition(0.10, 0.0, condition_id="a"), condition(0.12, 1.0, condition_id="b"),
                condition(0.08, -1.0, condition_id="c")]
        self.assertEqual(bootstrap_trend_ci(data), bootstrap_trend_ci(data))
        self.assertNotEqual(bootstrap_trend_ci(data, seed=1), bootstrap_trend_ci(data, seed=2))


class DetectOutliersPerCondition(unittest.TestCase):
    """A flag must mean 'off its own condition's trend', never 'off the
    pooled record trend'. Until 2026-08-25 the detector fitted one pooled
    regression, so between-condition offsets and slope disagreement were
    read as outliers."""

    def test_a_clean_steep_condition_among_flat_ones_is_not_flagged(self):
        # Five flat conditions plus one steep one: the pooled fit ran
        # between them, so the steep condition's La and Lu carried the
        # largest residuals of the record despite lying exactly on their
        # own line. Per-condition fitting leaves nothing to flag.
        data = [condition(0.0, 0.0, condition_id=f"flat{i}", noise=0.05)
                for i in range(5)]
        data.append(condition(0.30, -19.2, condition_id="steep", noise=0.05))
        self.assertEqual(detect_outliers(data), [])

    def test_a_within_condition_spike_is_flagged_against_its_own_fit(self):
        data = [condition(0.10, 0.0, condition_id="clean", noise=0.05),
                condition(0.10, 2.0, condition_id="spiked", noise=0.05)]
        data[1]["distribution_ratios"]["Eu"]["log_value"] += 3.0
        flags = detect_outliers(data)
        self.assertEqual([(f["element"], f["condition_id"]) for f in flags],
                         [("Eu", "spiked")])
        # The prediction is the spiked condition's own line, so the
        # residual is the spike itself, not a pooled-fit artefact.
        self.assertGreater(flags[0]["residual"], 2.0)
        self.assertIn(flags[0]["severity"], ("Medium", "High"))

    def test_conditions_below_the_minimum_contribute_no_flags(self):
        # A wild point in a condition too small to define a trend is
        # un-judgeable: there is no fit of its own to deviate from.
        small = condition(0.0, 0.0, elements=LN[:MIN_POINTS_PER_CONDITION - 1],
                          condition_id="small")
        small["distribution_ratios"]["La"]["log_value"] = -9.0
        clean = condition(0.10, 0.0, condition_id="clean", noise=0.05)
        self.assertEqual(detect_outliers([small]), [])
        self.assertEqual(detect_outliers([clean, small]),
                         detect_outliers([clean]))

    def test_y_and_sc_are_never_flagged(self):
        data = [condition(0.10, 0.0, condition_id=str(i), noise=0.05)
                for i in range(3)]
        data[0]["distribution_ratios"]["Y"] = {"value": 1e9, "log_value": 9.0}
        self.assertEqual(detect_outliers(data), [])

    def test_flat_residuals_produce_no_flags(self):
        data = [condition(0.10, 0.0, condition_id="exact")]  # zero residuals
        self.assertEqual(detect_outliers(data), [])


class AgainstTheDatabase(unittest.TestCase):
    """The bootstrap's units must be the classifier's units."""

    @classmethod
    def setUpClass(cls):
        with DB_PATH.open(encoding="utf-8") as handle:
            cls.records = [json.loads(line) for line in handle if line.strip()]

    def test_mean_of_per_condition_slopes_is_the_stored_slope(self):
        off = []
        for record in self.records:
            stored = (record.get("selectivity") or {}).get("slope")
            fits = per_condition_slopes(record.get("experimental_data", []))
            if stored is None:
                if fits:
                    off.append(f"{record['ligand']['Name']}: slope None but "
                               f"{len(fits)} eligible conditions")
                continue
            mean_slope = round(sum(s for s, _ in fits) / len(fits), 4)
            if abs(mean_slope - stored) > 1e-4:
                off.append(f"{record['ligand']['Name']}: {mean_slope} != {stored}")
        self.assertEqual(off, [], off[:5])

    def test_eligible_condition_count_matches_the_classifier_notes(self):
        off = []
        for record in self.records:
            fresh = classify_selectivity_trend(record.get("experimental_data", []))
            notes = fresh.get("notes") or ""
            n_fits = len(per_condition_slopes(record.get("experimental_data", [])))
            if notes.startswith("Based on "):
                n_notes = int(notes.split()[2])
                if n_notes != n_fits:
                    off.append(f"{record['ligand']['Name']}: notes {n_notes}, bootstrap {n_fits}")
            elif n_fits:
                off.append(f"{record['ligand']['Name']}: Insufficient_Data but {n_fits} eligible")
        self.assertEqual(off, [], off[:5])


if __name__ == "__main__":
    unittest.main()
