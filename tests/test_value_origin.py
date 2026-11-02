"""End to end: how each stored D was obtained, and how well it is known.

The record schema (still v0.3: the change is additive) gains three optional
fields on a distribution-ratio cell:

  value_origin   reported | digitised | interpolated | derived | censored
  censor         ">" or "<": `value`/`log_value` hold the bound (D > 100)
  sigma_source   "prior" (sigma_log_D is the flat prior by source type) or
                 "measured" (sigma_log_D was computed for this value, e.g.
                 from the figure calibration's marker size)

Absent fields mean what v13 always meant: origin from `source_quality`, no
censoring, prior sigma. Ways this could go wrong, each pinned below:

1. A measured sigma is overwritten by the prior when derived fields are
   recomputed, or `recompute_derived_fields.py --check` calls it stale.
2. A censored bound enters the trend fit (the classifier or the bootstrap
   interval), outlier detection or the condition statistics as if it were
   a measurement — Sypula 2012's ">100" rows would read as a flat trend.
3. A separation factor is computed from a censored bound.
4. The schema rejects a record carrying the new fields, or accepts a bad
   one: an origin outside the vocabulary, a censor other than > or <, a
   censored origin without a bound direction, a measured source without a
   sigma.
5. The v13 database (no new fields) changes on recompute, or fails the
   schema.
6. A condition whose every value is a bound (Han 2018's Pr rows: all below
   the figure's resolution) fails the schema, because its statistics block
   (measured values only) is empty.

The recomputed records and schema errors are produced in memory; the check
against the real database runs the CI script itself.
"""

import copy
import math
import json
import subprocess
import sys
import unittest
from pathlib import Path

from rel_quest.derived import derived_differences, recompute_record
from rel_quest.uncertainty import condition_trend_points, detect_outliers, is_censored
from rel_quest.utils import condition_statistics
from rel_quest.validators import validate_against_json_schema

ROOT = Path(__file__).resolve().parent.parent
DB = ROOT / "data" / "sx_ligands_merged_v14_flat.jsonl"   # fixture records are copied from it
LN = ["La", "Ce", "Pr", "Nd", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu"]


def record(name="TODGA"):
    for line in DB.open(encoding="utf-8"):
        r = json.loads(line)
        if r["ligand"]["Name"] == name:
            return r
    raise LookupError(name)


def ascending_condition(censor_top=False):
    """A clean Ascending condition (log D = 0.1 * i); optionally the four heaviest are '>100' bounds."""
    base = copy.deepcopy(record()["experimental_data"][0])
    cells = {}
    for i, el in enumerate(LN):
        log_d = round(0.1 * i, 3)
        cell = {"value": round(10 ** log_d, 6), "log_value": log_d, "n_replicates": 1,
                "sigma_log_D": 0.0326, "source_quality": "figure"}
        if censor_top and i >= 10:
            cell.update(value=100.0, log_value=2.0, censor=">", value_origin="censored")
        cells[el] = cell
    base["distribution_ratios"] = cells
    return base


class TestValueOrigin(unittest.TestCase):

    def test_1_a_measured_sigma_survives_recompute(self):
        r = record()
        cell = next(iter(r["experimental_data"][0]["distribution_ratios"].values()))
        cell.update(sigma_log_D=0.41, sigma_source="measured", value_origin="digitised")
        fresh = recompute_record(r)
        fresh_cell = next(iter(fresh["experimental_data"][0]["distribution_ratios"].values()))
        self.assertEqual(fresh_cell["sigma_log_D"], 0.41)
        self.assertNotIn("sigma_log_D", " ".join(derived_differences(r, fresh)))

    def test_2_censored_bounds_stay_out_of_fits_and_statistics(self):
        clean, censored = ascending_condition(), ascending_condition(censor_top=True)
        self.assertTrue(is_censored(censored["distribution_ratios"]["Lu"]))
        self.assertFalse(is_censored(clean["distribution_ratios"]["Lu"]))
        self.assertEqual(len(condition_trend_points(censored)), 10)
        self.assertEqual(condition_statistics(censored["distribution_ratios"])["n_elements"], 10)
        self.assertEqual(condition_statistics(censored["distribution_ratios"])["max_log_d"], 0.9)
        r = record()
        r["experimental_data"] = [censored]
        sel = recompute_record(r)["selectivity"]
        self.assertEqual(sel["trend_classification"], "Ascending")
        # The fit regresses on atomic number over the ten uncensored points only.
        from scipy.stats import linregress
        pts = condition_trend_points(censored)
        self.assertAlmostEqual(sel["slope"], round(linregress(*zip(*pts)).slope, 4), places=4)
        self.assertNotIn("Lu", sel["element_coverage"])
        self.assertEqual(detect_outliers([censored]), [])

    def test_3_no_separation_factor_from_a_bound(self):
        from rel_quest.core import LigandDatabase
        r = record()
        r["experimental_data"] = [ascending_condition(censor_top=True)]
        q = LigandDatabase()._record_to_query(r)
        self.assertIsNone(q.get_separation_factor("Lu", "La"))
        self.assertIsNotNone(q.get_separation_factor("Nd", "La"))

    def test_4_schema_accepts_the_fields_and_rejects_bad_values(self):
        r = recompute_record(record())
        cells = r["experimental_data"][0]["distribution_ratios"]
        first = next(iter(cells))
        cells[first].update(value_origin="censored", censor=">", sigma_source="prior")
        self.assertEqual(validate_against_json_schema(r), [])
        for bad in ({"value_origin": "guessed"}, {"censor": "~"}, {"sigma_source": "vibes"},
                    {"value_origin": "censored", "censor": None}):
            broken = copy.deepcopy(r)
            broken["experimental_data"][0]["distribution_ratios"][first].update(bad)
            self.assertTrue(validate_against_json_schema(broken), bad)

    def test_6_a_condition_of_bounds_only_is_valid(self):
        r = record()
        cond = ascending_condition()
        for c in cond["distribution_ratios"].values():
            c.update(value=0.5, log_value=round(math.log10(0.5), 4), censor="<", value_origin="censored")
        r["experimental_data"] = [cond] + r["experimental_data"][:1]
        fresh = recompute_record(r)
        self.assertEqual(fresh["experimental_data"][0]["statistics"], {})
        self.assertEqual(validate_against_json_schema(fresh), [])

    def test_5_v13_is_unchanged_and_valid(self):
        result = subprocess.run([sys.executable, "scripts/recompute_derived_fields.py", "--check"],
                                cwd=ROOT, capture_output=True, text=True)
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)


if __name__ == "__main__":
    unittest.main()
