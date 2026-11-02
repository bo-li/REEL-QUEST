"""Smoke test for the Streamlit dashboard.

app.py is 3,000+ lines and had no automated coverage at all — the plotting
consolidation in 2026-07 was verified by hand. This runs the whole script
through Streamlit's AppTest harness (no browser required) and fails if any
callback raises or if the charts stop rendering.

It is deliberately shallow: it asserts the app boots, loads the canonical
database and produces plots. It does not assert plot contents.
"""

import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
APP = ROOT / "app.py"

# A full AppTest run loads the database and renders every tab; allow headroom
# for slower CI machines.
TIMEOUT_SECONDS = 300


class TestDashboardSmoke(unittest.TestCase):
    """End-to-end: does app.py actually run?"""

    @classmethod
    def setUpClass(cls):
        try:
            from streamlit.testing.v1 import AppTest
        except ImportError:  # pragma: no cover - depends on environment
            raise unittest.SkipTest("streamlit not installed")

        cls.app = AppTest.from_file(str(APP), default_timeout=TIMEOUT_SECONDS)
        cls.app.run()

    def test_runs_without_exceptions(self):
        """No uncaught exception anywhere in the script run."""
        exceptions = [
            getattr(exc, "value", exc) for exc in self.app.exception
        ]
        self.assertEqual(exceptions, [], f"Dashboard raised: {exceptions}")
        print("✓ Dashboard ran with no exceptions")

    def test_renders_plots(self):
        """The overview/analysis tabs produced Plotly charts."""
        charts = self.app.get("plotly_chart")
        self.assertGreater(len(charts), 0, "Dashboard rendered no plots")
        print(f"✓ Dashboard rendered {len(charts)} plots")

    def test_reports_no_errors(self):
        """No st.error() surfaced (e.g. the database-not-found path)."""
        errors = [e.value for e in self.app.error]
        self.assertEqual(errors, [], f"Dashboard surfaced st.error: {errors}")
        print("✓ Dashboard surfaced no errors")

    def test_element_selectors_render(self):
        """All three periodic-table selectors render their 15 buttons.

        They share render_element_selector(); this catches a bad state_key or
        a duplicated widget-key prefix, which Streamlit would reject at runtime.
        """
        keys = {b.key for b in self.app.get("button") if b.key}
        for prefix in ("elem", "ln1", "ln2"):
            found = {k for k in keys if k.startswith(f"{prefix}_")}
            self.assertEqual(
                len(found), 15,
                f"selector '{prefix}' rendered {len(found)} buttons, expected 15",
            )
        print("✓ All 3 element selectors rendered 15 buttons each")

    def test_ligand_detail_view_renders(self):
        """The detail view renders for ligands that hit different branches.

        It is ~440 lines and unreachable from a plain app run (it needs a table
        row selection, which AppTest cannot set), so it is driven through
        _detail_probe_app.py. Rows are chosen to cover records with and without
        PubChem IDs, alternative names, and multiple data sources.
        """
        from streamlit.testing.v1 import AppTest

        probe = Path(__file__).resolve().parent / "_detail_probe_app.py"
        if not probe.exists():  # pragma: no cover
            self.skipTest("detail probe app missing")

        for row in (0, 5, 8, 18, 165):
            with self.subTest(row=row):
                at = AppTest.from_file(str(probe), default_timeout=TIMEOUT_SECONDS)
                at.session_state["probe_row"] = row
                at.run()

                exceptions = [str(getattr(e, "value", e)) for e in at.exception]
                self.assertEqual(exceptions, [], f"detail view row {row} raised")
                self.assertGreater(
                    len(at.get("plotly_chart")), 0,
                    f"detail view row {row} rendered no trend plot",
                )
        print("✓ Ligand detail view renders for 5 branch-diverse ligands")

    def test_detail_view_without_selection_shows_placeholder(self):
        """With no row selected the detail view prompts instead of rendering."""
        from streamlit.testing.v1 import AppTest

        probe = Path(__file__).resolve().parent / "_detail_probe_app.py"
        at = AppTest.from_file(str(probe), default_timeout=TIMEOUT_SECONDS)
        at.session_state["probe_row"] = -1          # no selection
        at.run()

        self.assertEqual(
            [str(getattr(e, "value", e)) for e in at.exception], [],
            "empty-selection path raised",
        )
        self.assertEqual(
            len(at.get("plotly_chart")), 0,
            "no ligand is selected, so nothing should be plotted",
        )
        messages = [i.value for i in at.get("info")]
        self.assertTrue(
            any("Select a ligand" in m for m in messages),
            f"expected a 'select a ligand' prompt, got {messages}",
        )
        print("✓ Detail view shows the placeholder when nothing is selected")

    def test_ai_assistant_page_short_circuits(self):
        """Switching to the AI Assistant page bypasses the database content."""
        from streamlit.testing.v1 import AppTest

        at = AppTest.from_file(str(APP), default_timeout=TIMEOUT_SECONDS)
        at.run()
        at.sidebar.radio[0].set_value("🤖 AI Assistant").run()

        self.assertEqual(
            [str(getattr(e, "value", e)) for e in at.exception], [],
            "AI Assistant page raised",
        )
        self.assertEqual(
            len(at.get("plotly_chart")), 0,
            "AI Assistant page should not render database charts",
        )
        print("✓ AI Assistant page short-circuits cleanly")


if __name__ == "__main__":
    unittest.main()
