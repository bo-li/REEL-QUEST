"""A mixed-extractant condition stays out of every number pooled across conditions.

A condition measured with other extractants in the organic phase (Co_Extractants,
uncertainty.is_mixture) gives the D of the mixture, not of the record's ligand.
The trend fit and the dashboard leave it out; the manuscript's pooled numbers did
not, so Fig. 5c ranked the multidentate N-donors first on 18 BTBP + co-extractant
conditions (median SF ~100 instead of ~10) while the dashboard (Fig. 2a) did not.

Ways this goes wrong, each pinned below:
  1. a mixture's log D enters Fig. 5 (log D(Eu), log D(La) or SF(Eu/La) by subtype);
  2. a mixture's adjacent-pair SF enters the manuscript's 2-sigma count;
  3. ... or the per-subtype counts behind Fig. 8;
  4. the test is vacuous: the canonical data must hold mixtures these would catch.
Each check deletes every mixture condition and requires the number not to move.
"""

import copy
import importlib.util
import json
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
GAP_KEYS = ("n_ligands_with_adjacent_sf", "n_adjacent_sf", "adjacent_sf_within_2_sigma")


def load(path, name):
    spec = importlib.util.spec_from_file_location(name, ROOT / path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def records():
    path = ROOT / "data" / "sx_ligands_merged_v14_flat.jsonl"
    return [json.loads(line) for line in path.open(encoding="utf-8") if line.strip()]


def without_mixtures(recs):
    from rel_quest.uncertainty import is_mixture
    out = copy.deepcopy(recs)
    for r in out:
        r["experimental_data"] = [c for c in r["experimental_data"] if not is_mixture(c)]
    return out


class TestPooledNumbersExcludeMixtures(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.stats = load("paper/stats.py", "paper_stats_mixtures")
        cls.records = records()
        cls.stripped = without_mixtures(cls.records)

    def test_the_data_holds_mixtures_that_would_count(self):              # 4
        from rel_quest.uncertainty import is_mixture
        mixtures = [c for r in self.records for c in r["experimental_data"] if is_mixture(c)]
        self.assertTrue(any(self.stats.adjacent_sf(c) for c in mixtures))
        self.assertTrue(any({"Eu", "La"} <= set(c["distribution_ratios"]) for c in mixtures))

    def test_figure_5_ignores_mixtures(self):                              # 1
        try:
            import matplotlib  # noqa: F401
        except ImportError:  # pragma: no cover - optional `paper` extra
            self.skipTest("matplotlib not installed (the `paper` extra)")
        figures = load("paper/figures.py", "paper_figures_mixtures")
        self.assertEqual(figures.showcase_groups(self.records), figures.showcase_groups(self.stripped))
        for element in ("Eu", "La"):
            self.assertEqual(figures.log_d_by_subtype(self.records, element),
                             figures.log_d_by_subtype(self.stripped, element))

    def test_adjacent_sf_count_ignores_mixtures(self):                    # 2
        self.assertEqual(self.stats.adjacent_sf_summary(self.records),
                         self.stats.adjacent_sf_summary(self.stripped))

    def test_subtype_gaps_ignore_mixtures(self):                          # 3
        full = self.stats.subtype_gaps(self.records, {})
        bare = self.stats.subtype_gaps(self.stripped, {})
        for subtype in full:
            with self.subTest(subtype=subtype):
                self.assertEqual({k: full[subtype][k] for k in GAP_KEYS},
                                 {k: bare[subtype][k] for k in GAP_KEYS})


if __name__ == "__main__":
    unittest.main()
