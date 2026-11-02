"""Tests for the Ligand Browser filter logic.

app._filter_ligands() is pure logic with no Streamlit calls, so the sidebar
filters can be tested directly instead of through the UI.
"""

import unittest

from rel_quest import LigandDatabase
from rel_quest.schema import PROFILE

NO_FILTERS = dict(
    search_text="",
    extractant_types_filter=[],
    mechanisms=[],
    trends=[],
    show_liu_only=False,
    sources=[],
    required_elements=[],
    logd_min=-4.0,
    logd_max=3.0,
)


class TestFilterLigands(unittest.TestCase):
    """Each sidebar filter narrows the result set as advertised."""

    @classmethod
    def setUpClass(cls):
        import app

        # Held as a module reference: assigning the function to a class
        # attribute would bind it as a method and shift the arguments.
        cls.app = app
        cls.db = LigandDatabase()
        cls.total = len(cls.db.get_all_ligands())

    def run_filter(self, **overrides):
        kwargs = {**NO_FILTERS, **overrides}
        return self.app._filter_ligands(self.db, **kwargs)

    def test_no_filters_returns_everything(self):
        """The default slider range must not drop any ligand."""
        self.assertEqual(len(self.run_filter()), self.total)

    def test_search_is_case_insensitive_substring(self):
        upper = self.run_filter(search_text="TODGA")
        lower = self.run_filter(search_text="todga")
        self.assertEqual({l.name for l in upper}, {l.name for l in lower})
        self.assertTrue(upper, "expected at least one TODGA match")
        for lig in upper:
            self.assertIn("todga", lig.name.lower())

    def test_mechanism_filter(self):
        solvation = self.run_filter(mechanisms=["Solvation"])
        self.assertTrue(solvation)
        self.assertTrue(all(l.mechanism == "Solvation" for l in solvation))
        self.assertLess(len(solvation), self.total)

    def test_trend_filter_accepts_quadratic_classes(self):
        """Peaked/U_Shaped are real trend values and must be filterable."""
        for trend in ("Peaked", "U_Shaped"):
            with self.subTest(trend=trend):
                hits = self.run_filter(trends=[trend])
                self.assertTrue(hits, f"no ligands with trend {trend}")
                self.assertTrue(all(l.selectivity_trend == trend for l in hits))

    @unittest.skipIf(PROFILE == "public", "source labels are in the developer repository only")
    def test_liu_only_filter(self):
        liu = self.run_filter(show_liu_only=True)
        self.assertTrue(liu)
        self.assertTrue(all("Liu_2022_JACSAu" in l.data_sources for l in liu))

    def test_required_elements_requires_all(self):
        both = self.run_filter(required_elements=["Eu", "La"])
        self.assertTrue(both)
        for lig in both:
            for element in ("Eu", "La"):
                self.assertTrue(
                    any(element in c.get("distribution_ratios", {})
                        for c in lig.experimental_conditions),
                    f"{lig.name} missing {element}",
                )

    def test_filters_compose(self):
        """Combining filters is at least as restrictive as either alone."""
        a = self.run_filter(mechanisms=["Solvation"])
        b = self.run_filter(show_liu_only=True)
        both = self.run_filter(mechanisms=["Solvation"], show_liu_only=True)
        self.assertLessEqual(len(both), min(len(a), len(b)))

    def test_narrow_logd_range_excludes_ligands(self):
        """A narrow logD window must drop ligands with no value inside it."""
        narrow = self.run_filter(logd_min=2.5, logd_max=3.0)
        self.assertLess(len(narrow), self.total)
        for lig in narrow:
            values = [
                data.get("log_value")
                for cond in lig.experimental_conditions
                for data in cond.get("distribution_ratios", {}).values()
            ]
            self.assertTrue(
                any(v is not None and 2.5 <= v <= 3.0 for v in values),
                f"{lig.name} has no logD in [2.5, 3.0]",
            )


if __name__ == "__main__":
    unittest.main()
