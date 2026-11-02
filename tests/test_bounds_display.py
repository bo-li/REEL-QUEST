"""End to end: censored cells (bounds, D > x or D < x) as the dashboard shows them.

v14 stores a bound in a cell's `value`/`log_value` with `censor` '>' or '<'
(rel_quest.uncertainty.is_censored). A bound is a measurement of that element
but not a value, so the dashboard must neither drop it nor treat it as one.
Ways it could mislead, each pinned below:

1. A bound enters a statistic as if it were a value: the box/violin of the
   single-element and logSF viewers, their category medians (which also set
   the category order), the aggregated (median) D of a ligand behind the
   Eu/La scatter and the aggregated trend plot, a per-condition trend fit,
   or the R² column of the condition table.
2. A bound disappears: no marker for it, or a category whose every point is
   a bound vanishes from the plot.
3. A bound is drawn but reads as a value: no direction in its symbol, or a
   hover that says "log D = 0.75" instead of "D > 5.62".
4. The separation factor of a pair with a bound gets the wrong direction.
   logSF = log D1 - log D2: D1 '>' or D2 '<' make it a lower bound, D1 '<'
   or D2 '>' an upper bound; two bounds that point the same way give no
   information and must be left out (and counted), never plotted.
5. The sidebar's (or the REST API's) log D range filter counts a bound as
   "a value in range".
6. v13, which has no bounds, renders differently afterwards. Pinned here:
   v13 has no bounds, so the new paths are never taken. The stronger check
   (frames and aggregates equal to a pickle taken before the change) was run
   once, 2026-10-05, and runs again when RELQUEST_V13_BASELINE names one.

A small fixture database (two ligands copied from v13, cells replaced) is
written to the test's temporary directory; the v13 check reads the canonical
file.
"""

import copy
import json
import math
import tempfile
import unittest
from pathlib import Path

from rel_quest import LigandDatabase
from rel_quest.plotting import (
    extract_element_data,
    extract_separation_factor_data,
    separation_factor_censor,
)

ROOT = Path(__file__).resolve().parent.parent
V13 = ROOT / "data" / "archive" / "sx_ligands_merged_v13_flat.jsonl"
# Fixture records are copied from the canonical file, which every checkout has
CANONICAL = ROOT / "data" / "sx_ligands_merged_v14_flat.jsonl"


def cell(log_d, censor=None):
    c = {"value": round(10 ** log_d, 6), "log_value": log_d, "n_replicates": 1, "sigma_log_D": 0.05,
         "source_quality": "figure"}
    if censor:
        c.update(censor=censor, value_origin="censored")
    return c


def condition(template, cid, cells, acid=3.0):
    c = copy.deepcopy(template)
    c.update(condition_id=cid, Acid_Concentration_M=acid, distribution_ratios=cells)
    return c


class Fixture(unittest.TestCase):
    """FixA (a DGA): values and bounds mixed. FixB (Malonamide): Eu measured only as bounds."""

    @classmethod
    def setUpClass(cls):
        records = {json.loads(l)["ligand"]["Name"]: json.loads(l) for l in CANONICAL.open(encoding="utf-8")}
        a, b = copy.deepcopy(records["TODGA"]), copy.deepcopy(records["DMDOHEMA"])
        t = a["experimental_data"][0]
        a["ligand"]["Name"], b["ligand"]["Name"] = "FixA", "FixB"
        a["experimental_data"] = [
            # a trend condition: Nd is a lower bound, the rest values
            condition(t, "A1", {"La": cell(0.0), "Ce": cell(0.3), "Pr": cell(0.6), "Nd": cell(2.5, ">"),
                                "Sm": cell(1.5), "Eu": cell(1.8)}),
            condition(t, "A2", {"Eu": cell(3.0, ">"), "La": cell(0.2)}),    # SF Eu/La > 2.8
            condition(t, "A3", {"Eu": cell(0.5), "La": cell(-2.0, "<")}),   # SF Eu/La > 2.5
            condition(t, "A4", {"Eu": cell(2.0, ">"), "La": cell(1.0, ">")}),  # no information
            condition(t, "A5", {"Eu": cell(-1.0, "<"), "La": cell(0.0, ">")}),  # SF Eu/La < -1.0
        ]
        b["experimental_data"] = [condition(b["experimental_data"][0], "B1", {"Eu": cell(-1.5, "<"),
                                                                              "La": cell(-0.5)})]
        cls.tmp = tempfile.TemporaryDirectory()
        path = Path(cls.tmp.name) / "fixture.jsonl"
        path.write_text("".join(json.dumps(r) + "\n" for r in (a, b)), encoding="utf-8")
        cls.db = LigandDatabase(str(path))
        assert cls.db.load_errors == [], cls.db.load_errors
        cls.ligands = cls.db.get_all_ligands()
        cls.eu = extract_element_data(cls.db, "Eu", cls.ligands)
        cls.sf = extract_separation_factor_data(cls.db, "Eu", "La", cls.ligands)
        # pandas 3 stores a missing censor as NaN, pandas 2 as None: the plots must take both
        cls.sf_nan = cls.sf.copy()
        for c in ("censor", "censor_elem1", "censor_elem2"):
            cls.sf_nan[c] = cls.sf_nan[c].where(cls.sf_nan[c].notna(), float("nan"))
        import app
        cls.app = app

    @classmethod
    def tearDownClass(cls):
        cls.tmp.cleanup()

    @staticmethod
    def traces(fig, kind):
        return [t for t in fig.data if t.type == kind]

    @staticmethod
    def bound_traces(fig):
        return [t for t in fig.data if t.type == "scatter" and t.marker is not None
                and t.marker.symbol is not None and "triangle" in str(t.marker.symbol)]


class TestSeparationFactorAlgebra(unittest.TestCase):

    def test_4_every_combination(self):
        self.assertEqual(separation_factor_censor(None, None), "=")
        self.assertEqual(separation_factor_censor(">", None), ">")
        self.assertEqual(separation_factor_censor("<", None), "<")
        self.assertEqual(separation_factor_censor(None, ">"), "<")
        self.assertEqual(separation_factor_censor(None, "<"), ">")
        self.assertEqual(separation_factor_censor(">", "<"), ">")
        self.assertEqual(separation_factor_censor("<", ">"), "<")
        self.assertIsNone(separation_factor_censor(">", ">"))
        self.assertIsNone(separation_factor_censor("<", "<"))


class TestFrames(Fixture):

    def test_2_bounds_stay_in_the_frames_marked(self):
        eu = self.eu.set_index("condition_id")
        self.assertEqual(eu.loc["A2", "censor"], ">")
        self.assertEqual(eu.loc["B1", "censor"], "<")
        self.assertTrue(eu.loc[["A1", "A3"], "censor"].isna().all())

    def test_4_separation_factor_rows(self):
        sf = self.sf.set_index("condition_id")
        self.assertEqual((sf.loc["A2", "censor"], round(sf.loc["A2", "log_sf"], 3)), (">", 2.8))
        self.assertEqual((sf.loc["A3", "censor"], round(sf.loc["A3", "log_sf"], 3)), (">", 2.5))
        self.assertEqual((sf.loc["A5", "censor"], round(sf.loc["A5", "log_sf"], 3)), ("<", -1.0))
        self.assertNotIn("A4", sf.index)
        self.assertEqual(self.sf.attrs["no_information"], 1)

    def test_1_aggregates_use_values_only(self):
        agg = next(l for l in self.ligands if l.name == "FixA").get_aggregated_d_values("median", True)
        self.assertAlmostEqual(agg["Eu"]["log_value"], 1.15)  # median of 1.8, 0.5 (not 3.0, 2.0, -1.0)
        self.assertEqual((agg["Eu"]["n_conditions"], agg["Eu"]["n_bounds"]), (2, 3))
        self.assertNotIn("Nd", agg)                           # measured only as a bound
        b = next(l for l in self.ligands if l.name == "FixB").get_aggregated_d_values("median", True)
        self.assertNotIn("Eu", b)


class TestPlots(Fixture):

    def test_1_and_2_element_plot(self):
        for mode in ("Ligand Types", "Solvents", "Ligand Acidity"):
            with self.subTest(mode=mode):
                fig = self.app.plot_element_distribution(self.eu, "Eu", mode)
                plotted = [y for t in self.traces(fig, "box") + self.traces(fig, "violin") for y in t.y]
                for bound in (3.0, 2.0, -1.0, -1.5):
                    self.assertNotIn(bound, plotted)
                ys = sorted(y for t in self.bound_traces(fig) for y in t.y)
                self.assertEqual(ys, [-1.5, -1.0, 2.0, 3.0])
        fig = self.app.plot_element_distribution(self.eu, "Eu", "Ligand Types")
        # FixB's family (all bounds) still appears on the axis
        self.assertIn("Malonamide", list(fig.layout.xaxis.categoryarray))

    def test_3_bound_markers_read_as_bounds(self):
        fig = self.app.plot_element_distribution(self.eu, "Eu", "Ligand Types")
        symbols = {s for t in self.bound_traces(fig) for s in
                   (t.marker.symbol if isinstance(t.marker.symbol, (list, tuple)) else [t.marker.symbol])}
        self.assertEqual(symbols, {"triangle-up-open", "triangle-down-open"})
        hovers = " ".join(str(h) for t in self.bound_traces(fig) for h in (t.hovertext or t.text or []))
        self.assertIn("D[Eu] > 100 (log 2.00)", hovers)
        self.assertIn("D[Eu] < 0.1 (log -1.00)", hovers)

    def test_1_and_2_separation_factor_plot(self):
        for mode in ("Ligand Types", "Acid Concentration"):
            with self.subTest(mode=mode):
                fig = self.app.plot_separation_factor_distribution(self.sf, "Eu", "La", mode)
                self.app.plot_separation_factor_distribution(self.sf_nan, "Eu", "La", mode)
                plotted = [y for t in self.traces(fig, "box") + self.traces(fig, "violin") for y in t.y]
                if mode == "Acid Concentration":
                    plotted = [y for t in fig.data if t.type == "scatter" and t not in self.bound_traces(fig)
                               for y in t.y]
                for bound in (2.8, 2.5, -1.0):   # no value in the fixture sits on a bound
                    self.assertFalse(any(math.isclose(y, bound) for y in plotted), (mode, bound))
                ys = sorted(round(y, 3) for t in self.bound_traces(fig) for y in t.y)
                self.assertEqual(ys, [-1.0, -1.0, 2.5, 2.8])   # A5, B1 (Eu < -1.5, La -0.5), A3, A2

    def test_1_trend_plot_and_r2_leave_the_bound_out(self):
        fixa = next(l for l in self.ligands if l.name == "FixA")
        fig = self.app.plot_selectivity_trend(fixa, "A1")
        data = next(t for t in fig.data if t.name == "Data")
        self.assertNotIn(60, list(data.x))                     # Nd (Z 60) is not a data point
        (bound,) = self.bound_traces(fig)
        self.assertEqual((list(bound.x), list(bound.y)), ([60], [2.5]))
        r2, slope = self.app.calculate_condition_r2(fixa.experimental_conditions[0]["distribution_ratios"])
        # OLS of La..Eu without Nd: points (57,0) (58,.3) (59,.6) (62,1.5) (63,1.8) lie on slope 0.3
        self.assertAlmostEqual(slope, 0.3, places=6)
        self.assertAlmostEqual(r2, 1.0, places=6)

    def test_5_range_filter_counts_values_only(self):
        fixb = next(l for l in self.ligands if l.name == "FixB")
        self.assertFalse(self.app.has_value_in_range(fixb, -2.0, -1.0))   # only the Eu < -1.5 bound
        self.assertTrue(self.app.has_value_in_range(fixb, -1.0, 0.0))     # La -0.5
        from rel_quest import api
        names = lambda lo, hi: {l.name for l in api.apply_filters(self.ligands, logd_min=lo, logd_max=hi)}
        self.assertNotIn("FixB", names(-2.0, -1.0))
        self.assertIn("FixB", names(-1.0, 0.0))

    def test_1_and_2_element_plot_scatter_mode(self):
        fig = self.app.plot_element_distribution(self.eu, "Eu", "Acid Concentration")
        (bounds,) = self.bound_traces(fig)
        self.assertEqual(sorted(bounds.y), [-1.5, -1.0, 2.0, 3.0])
        self.assertEqual(set(bounds.x), {3.0})                 # the conditions' acid_conc
        values = [y for t in fig.data if t is not bounds for y in t.y]
        self.assertEqual(sorted(values), [0.5, 1.8])


@unittest.skipUnless(V13.exists(), "the v13 archive is in the developer repository only")
class TestV13Unchanged(unittest.TestCase):

    def test_6_v13_has_no_bounds(self):
        db = LigandDatabase(str(V13))
        ligands = db.get_all_ligands()
        for el in ("La", "Nd", "Eu", "Lu", "Y"):
            self.assertTrue(extract_element_data(db, el, ligands)["censor"].isna().all(), el)
        sf = extract_separation_factor_data(db, "Pr", "Nd", ligands)
        self.assertTrue(sf["censor"].isna().all())
        self.assertEqual(sf.attrs["no_information"], 0)
        self.assertTrue(all(d["n_bounds"] == 0 for l in ligands
                            for d in l.get_aggregated_d_values("median").values()))

    def test_6_frames_and_aggregates(self):
        import pickle
        import os
        path = os.environ.get("RELQUEST_V13_BASELINE")
        if not path or not Path(path).exists():
            self.skipTest("set RELQUEST_V13_BASELINE to a pickle of the pre-change frames")
        base = pickle.load(open(path, "rb"))
        db = LigandDatabase(str(V13))
        ligands = db.get_all_ligands()
        for el, frame in base["el"].items():
            now = extract_element_data(db, el, ligands)
            self.assertTrue(now["censor"].isna().all())
            self.assertTrue(now.drop(columns="censor").equals(frame), el)
        for pair, frame in base["sf"].items():
            now = extract_separation_factor_data(db, *pair, ligands)
            self.assertTrue(now[["censor", "censor_elem1", "censor_elem2"]].isna().all().all())
            self.assertTrue(now.drop(columns=["censor", "censor_elem1", "censor_elem2"]).equals(frame), pair)
        for lig in ligands:
            now = lig.get_aggregated_d_values("median", True)
            self.assertEqual({e: {k: v for k, v in d.items() if k != "n_bounds"} for e, d in now.items()},
                             base["agg"][lig.name], lig.name)


if __name__ == "__main__":
    unittest.main()
