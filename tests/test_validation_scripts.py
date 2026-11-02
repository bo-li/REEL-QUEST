"""The technical-validation helpers (rel_quest/validation.py) do what the
paper's Technical Validation section says they do."""

import json
import math
import unittest
from pathlib import Path
from rel_quest.schema import PROFILE

from rel_quest.validation import (
    DEFAULT_ALLOCATION, condition_key_matches, cross_source_pairs, draw_sample,
    score_sample, summarise_pairs)

ROOT = Path(__file__).resolve().parent.parent


def records():
    path = ROOT / "data" / "sx_ligands_merged_v14_flat.jsonl"
    return [json.loads(line) for line in path.open(encoding="utf-8") if line.strip()]


@unittest.skipIf(PROFILE == "public", "source labels are in the developer repository only")
class TestSample(unittest.TestCase):

    def setUp(self):
        self.records = records()
        self.every_pdf = lambda doi: "x.pdf" if doi else None

    def test_seeded_and_reproducible(self):
        a = draw_sample(self.records, lookup=self.every_pdf)
        b = draw_sample(self.records, lookup=self.every_pdf)
        self.assertEqual(a["rows"], b["rows"])
        c = draw_sample(self.records, seed=1, lookup=self.every_pdf)
        self.assertNotEqual([r["condition_id"] for r in a["rows"]],
                            [r["condition_id"] for r in c["rows"]])

    def test_allocation_is_honoured(self):
        sample = draw_sample(self.records, lookup=self.every_pdf)
        for (prefix, quality), n in DEFAULT_ALLOCATION.items():
            conditions = {r["condition_id"] for r in sample["rows"]
                          if r["stratum"] == f"{prefix}/{quality}"}
            with self.subTest(stratum=(prefix, quality)):
                self.assertEqual(len(conditions), n)

    def test_bounds_are_not_sampled(self):
        # A bound (D < x, D > x) has no value to re-digitise; a condition holding only bounds
        # would be drawn and give no row, so it is not eligible either.
        db = {c["condition_id"]: c for r in self.records for c in r["experimental_data"]}
        sample = draw_sample(self.records, lookup=self.every_pdf)
        sampled = [db[row["condition_id"]]["distribution_ratios"][row["element"]] for row in sample["rows"]]
        self.assertFalse([c for c in sampled if c.get("censor")])
        drawn = {row["condition_id"] for row in sample["rows"]}
        self.assertEqual(sum(sample["drawn"].values()), len(drawn))

    def test_only_conditions_with_a_pdf_are_eligible(self):
        sample = draw_sample(self.records, lookup=lambda doi: None)
        self.assertEqual(sample["rows"], [])


class TestScoring(unittest.TestCase):

    def test_errors_by_stratum(self):
        rows = [
            {"stratum": "A", "stored_log_value": "0.0", "redigitised_log_value": "0.1",
             "sigma_log_D": "0.0326"},
            {"stratum": "A", "stored_log_value": "1.0", "redigitised_value": str(10 ** 0.8),
             "sigma_log_D": "0.0326"},
            {"stratum": "B", "stored_log_value": "-1.0", "redigitised_log_value": "",
             "redigitised_value": ""},   # not filled: ignored
        ]
        summary = score_sample(rows)
        self.assertEqual(set(summary), {"A", "ALL"})
        self.assertEqual(summary["A"]["n_cells"], 2)
        self.assertAlmostEqual(summary["A"]["median_abs_dex"], 0.15, places=3)
        self.assertAlmostEqual(summary["A"]["rms_dex"], math.sqrt((0.01 + 0.04) / 2), places=3)
        self.assertEqual(summary["A"]["within_0.10_dex"], 0.5)

    def test_empty_sheet(self):
        self.assertEqual(score_sample([]), {})


class TestCrossSource(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.records = records()
        cls.strict = cross_source_pairs(cls.records, temp_tol=2.0)
        cls.relaxed = cross_source_pairs(cls.records, temp_tol=5.0)

    def test_strict_pairs_are_same_paper_transcriptions(self):
        # v14 keeps one copy of each duplicated transcription (the verification's
        # 'Duplicate - drop'), so v13's DMDBMA and HONTA pairs are gone; TODGA's and
        # T2EHDGA's went when aqueous additives and mixtures entered the key (they
        # paired different oxalic acid levels, and a T2EHDGA + HDEHP mixture).
        ligands = {p["ligand"] for p in self.strict}
        self.assertTrue({"C8-DOODA", "DMDDdDGA"} <= ligands, ligands)
        self.assertTrue(all(p["source_a"] != p["source_b"] for p in self.strict))

    def test_relaxed_tier_adds_inter_laboratory_pairs(self):
        different = {p["ligand"] for p in self.relaxed if not p["same_paper"]}
        self.assertTrue(different, "expected different-paper pairs in the relaxed tier")
        self.assertGreaterEqual(len(self.relaxed), len(self.strict))
        summary = summarise_pairs(self.relaxed)
        self.assertEqual(summary["all"]["n_cells"],
                         summary["same_paper"]["n_cells"] + summary["different_paper"]["n_cells"])

    def test_unknown_concentrations_never_match(self):
        a = {"Acid_Name": "HNO3", "Acid_Concentration_M": None, "Extractant_Concentration_M": 0.1,
             "Solvent_Name": "n-dodecane", "Temperature_C": 25, "Phase_Modifier_Name": ""}
        self.assertFalse(condition_key_matches(a, dict(a)))

    def test_solvent_aliases_in_parentheses_are_ignored(self):
        a = {"Acid_Name": "HNO3", "Acid_Concentration_M": 1.0, "Extractant_Concentration_M": 0.1,
             "Solvent_Name": "phenyl trifluoromethyl sulfone", "Temperature_C": 25,
             "Phase_Modifier_Name": ""}
        b = dict(a, Solvent_Name="phenyl trifluoromethyl sulfone (FS-13)", Phase_Modifier_Name=None)
        self.assertTrue(condition_key_matches(a, b))


if __name__ == "__main__":
    unittest.main()
