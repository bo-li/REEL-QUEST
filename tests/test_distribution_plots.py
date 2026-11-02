"""Tests for the two distribution plot builders.

plot_element_distribution and plot_separation_factor_distribution drive the
dashboard's two main viewers and had no coverage. These call them directly
(they return plotly figures and need no Streamlit context) across every x-axis
mode, asserting structure rather than pinned pixel values so the tests survive
data updates.
"""

import unittest

import pandas as pd

from rel_quest import LigandDatabase
from rel_quest.plotting import (
    extract_element_data,
    extract_separation_factor_data,
)

X_AXIS_MODES = [
    "Ligand Types", "Solvents", "Solvent Categories", "Ligand Acidity",
    "Acid Concentration", "Extractant Concentration", "Temperature",
]
CATEGORICAL_MODES = X_AXIS_MODES[:4]


class TestDistributionPlots(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        import app

        cls.app = app
        cls.db = LigandDatabase()
        ligands = cls.db.get_all_ligands()
        cls.element_df = extract_element_data(cls.db, "Eu", ligands)
        cls.sf_df = extract_separation_factor_data(cls.db, "Eu", "La", ligands)

    def test_element_plot_renders_every_mode(self):
        for mode in X_AXIS_MODES:
            with self.subTest(mode=mode):
                fig = self.app.plot_element_distribution(self.element_df, "Eu", mode)
                self.assertGreater(len(fig.data), 0, f"{mode} produced no traces")

    def test_separation_factor_plot_renders_every_mode(self):
        for mode in X_AXIS_MODES:
            with self.subTest(mode=mode):
                fig = self.app.plot_separation_factor_distribution(
                    self.sf_df, "Eu", "La", mode)
                self.assertGreater(len(fig.data), 0, f"{mode} produced no traces")

    def test_categorical_modes_use_box_and_violin(self):
        for mode in CATEGORICAL_MODES:
            with self.subTest(mode=mode):
                types = {t.type for t in
                         self.app.plot_element_distribution(self.element_df, "Eu", mode).data}
                self.assertIn("box", types, f"{mode} drew no box traces")
                self.assertIn("violin", types, f"{mode} drew no violin traces")

    def test_separation_factor_has_zero_reference_line(self):
        """logSF plots mark the no-selectivity baseline at y=0."""
        fig = self.app.plot_separation_factor_distribution(
            self.sf_df, "Eu", "La", "Ligand Types")
        self.assertGreaterEqual(
            len(fig.layout.shapes or []), 1,
            "expected a y=0 reference line on the separation-factor plot",
        )

    def test_category_color_map_returns_the_real_palettes(self):
        """_category_color_map must hand back the module palettes themselves.

        Asserted against the palette constants rather than against the
        function's own output, which would make the test self-fulfilling.
        """
        ligand_types = self.app._category_color_map("Ligand Types")
        self.assertEqual(ligand_types["DGA"], self.app.CHEMICAL_FAMILY_COLORS["DGA"])
        self.assertEqual(
            ligand_types["Phosphine_oxide"],
            self.app.CHEMICAL_FAMILY_COLORS["Phosphine_oxide"],
        )
        self.assertEqual(
            self.app._category_color_map("Ligand Acidity"),
            self.app.LIGAND_ACIDITY_COLORS,
        )
        self.assertEqual(
            self.app._category_color_map("Solvent Categories"),
            self.app.SOLVENT_CATEGORY_COLORS,
        )
        self.assertEqual(
            self.app._category_color_map("Solvents"), self.app.SOLVENT_COLORS)

    def test_ligand_type_colors_come_from_the_shared_palette(self):
        """Both plots colour Ligand Types from the real palette constants."""
        palette = (set(self.app.CHEMICAL_FAMILY_COLORS.values())
                   | set(self.app.EXTRACTANT_TYPE_COLORS.values())
                   | {self.app.BOUND_COLOR})   # bounds (D > x, D < x) are drawn neutral on purpose
        fallbacks = {"#999999", "#CCCCCC"}

        for name, fig in (
            ("element", self.app.plot_element_distribution(
                self.element_df, "Eu", "Ligand Types")),
            ("logSF", self.app.plot_separation_factor_distribution(
                self.sf_df, "Eu", "La", "Ligand Types")),
        ):
            with self.subTest(plot=name):
                colors = [
                    str(t.marker.color) for t in fig.data
                    if getattr(t, "marker", None) and t.marker.color
                ]
                self.assertTrue(colors, f"{name} plot set no marker colours")

                unknown = set(colors) - palette - fallbacks
                self.assertEqual(
                    unknown, set(),
                    f"{name} used colours outside the palette: {unknown}",
                )
                # Guard against a broken palette silently falling back for
                # everything, which a subset check alone would accept.
                real = [c for c in colors if c not in fallbacks]
                self.assertGreater(
                    len(real), len(colors) // 2,
                    f"{name}: most categories fell back to a default colour "
                    f"({len(real)}/{len(colors)} resolved)",
                )

    def test_empty_dataframes_return_placeholder_figures(self):
        """Both plots degrade to an annotated placeholder, not an exception."""
        empty = pd.DataFrame()

        fig = self.app.plot_element_distribution(empty, "Eu", "Ligand Types")
        self.assertEqual(len(fig.data), 0)
        self.assertTrue(fig.layout.annotations, "no 'no data' message shown")
        self.assertIn("Eu", fig.layout.annotations[0].text)

        fig = self.app.plot_separation_factor_distribution(
            empty, "Eu", "La", "Ligand Types")
        self.assertEqual(len(fig.data), 0)
        self.assertTrue(fig.layout.annotations, "no 'no data' message shown")


if __name__ == "__main__":
    unittest.main()
