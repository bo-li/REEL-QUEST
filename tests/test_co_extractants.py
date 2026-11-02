"""Mixed-extractant conditions (`Co_Extractants`): kept on the main ligand, left out of pooling.

A condition measured with a second extractant in the organic phase (e.g.
CyMe4-BTBP + TODGA, TODGA + DMDOHEMA) records it in `Co_Extractants`
([{Name, SMILES, Concentration_M}]). The D belongs to the mixture, not to
the ligand alone, so the condition stays on the record (it is a real
measurement, shown in the ligand's own condition table) but must not shape
anything pooled across a ligand's conditions. Ways that could fail:

1. A mixture condition votes in the selectivity trend, enters the slope or
   its bootstrap interval, or is scored by the outlier check; or its
   elements appear in `element_coverage` (the points the fit consumed);
   or the notes do not say a mixture was left out, so a mixture-only record
   reads as Insufficient_Data with nothing measured.
2. Its elements vanish from `elements_measured` (they were measured).
3. It enters a ligand's aggregated (median) D, or a point in the
   single-element or logSF viewers (they compare ligands with each other).
4. The JSON Schema rejects a well-formed `Co_Extractants`, or accepts one
   without a name.
5. The ligand's condition table does not say which conditions are mixtures.
"""

import copy
import json
import unittest
from pathlib import Path

from rel_quest import LigandDatabase
from rel_quest.derived import recompute_record
from rel_quest.plotting import extract_element_data, extract_separation_factor_data
from rel_quest.uncertainty import is_mixture
from rel_quest.validators import validate_against_json_schema

ROOT = Path(__file__).resolve().parent.parent
CANONICAL = ROOT / "data" / "sx_ligands_merged_v14_flat.jsonl"
LN = ["La", "Ce", "Pr", "Nd", "Sm", "Eu", "Gd"]
TODGA = {"Name": "TODGA", "SMILES": "CCCCCCCCN(CCCCCCCC)C(=O)COCC(=O)N(CCCCCCCC)CCCCCCCC",
         "Concentration_M": 0.005}


def cells(slope, start=0.0, elements=LN):
    return {el: {"value": round(10 ** (start + slope * i), 6), "log_value": round(start + slope * i, 4),
                 "n_replicates": 1, "sigma_log_D": 0.0326, "source_quality": "figure"}
            for i, el in enumerate(elements)}


class TestCoExtractants(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        base = next(json.loads(l) for l in CANONICAL.open(encoding="utf-8") if json.loads(l)["ligand"]["Name"] == "HONTA")
        t = base["experimental_data"][0]
        alone = [dict(copy.deepcopy(t), condition_id=f"A{i}", distribution_ratios=cells(0.3, i * 0.1))
                 for i in range(3)]                                  # ascending, without a co-extractant
        mixed = dict(copy.deepcopy(t), condition_id="M1", Co_Extractants=[TODGA],
                     distribution_ratios=cells(-0.4, 3.0, LN + ["Tb", "Dy", "Y"]))  # reverse, with TODGA
        cls.alone = recompute_record(dict(copy.deepcopy(base), experimental_data=alone))
        cls.both = recompute_record(dict(copy.deepcopy(base), experimental_data=alone + [mixed]))

    def test_is_mixture(self):
        self.assertTrue(is_mixture({"Co_Extractants": [TODGA]}))
        self.assertFalse(is_mixture({"Co_Extractants": []}))
        self.assertFalse(is_mixture({"Co_Extractants": None}))
        self.assertFalse(is_mixture({}))

    def test_1_trend_slope_interval_and_outliers_ignore_the_mixture(self):
        strip = lambda sel: {k: v for k, v in sel.items() if k != "notes"}   # noqa: E731
        self.assertEqual(strip(self.both["selectivity"]), strip(self.alone["selectivity"]))
        self.assertIn("1 condition(s) measured with a co-extractant", self.both["selectivity"]["notes"])
        self.assertEqual(self.both["selectivity"]["trend_classification"], "Ascending")
        self.assertNotIn("Tb", self.both["selectivity"]["element_coverage"])
        self.assertEqual(self.both["data_quality"]["outlier_flags"], self.alone["data_quality"]["outlier_flags"])

    def test_1_a_mixture_only_record_says_why_it_has_no_trend(self):
        only = recompute_record(dict(copy.deepcopy(self.both), experimental_data=self.both["experimental_data"][-1:]))
        self.assertEqual(only["selectivity"]["trend_classification"], "Insufficient_Data")
        self.assertIn("All 1 condition(s) were measured with a co-extractant", only["selectivity"]["notes"])

    def test_2_its_elements_are_still_measured(self):
        for el in ("Tb", "Dy", "Y"):
            self.assertIn(el, self.both["elements_measured"])
        self.assertEqual(self.both["performance_summary"]["conditions_tested"], 4)

    def test_3_aggregates_and_viewers_ignore_the_mixture(self):
        import tempfile
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "db.jsonl"
            path.write_text(json.dumps(self.both) + "\n", encoding="utf-8")
            db = LigandDatabase(str(path))
            self.assertEqual(db.load_errors, [])
            lig = db.get_all_ligands()[0]
            agg = lig.get_aggregated_d_values("median")
            self.assertNotIn("Tb", agg)
            self.assertEqual(agg["La"]["n_conditions"], 3)
            eu = extract_element_data(db, "Eu")
            self.assertNotIn("M1", set(eu["condition_id"]))
            self.assertEqual(len(eu), 3)
            sf = extract_separation_factor_data(db, "Eu", "La")
            self.assertNotIn("M1", set(sf["condition_id"]))

    def test_4_schema(self):
        self.assertEqual(validate_against_json_schema(self.both), [])
        bad = copy.deepcopy(self.both)
        bad["experimental_data"][-1]["Co_Extractants"] = [{"SMILES": "C", "Concentration_M": 0.1}]
        self.assertNotEqual(validate_against_json_schema(bad), [])

    def test_5_condition_table_names_the_co_extractant(self):
        import app
        self.assertEqual(app.co_extractants_text(self.both["experimental_data"][-1]), "TODGA 0.005 M")
        self.assertEqual(app.co_extractants_text(self.both["experimental_data"][0]), "")


if __name__ == "__main__":
    unittest.main()
