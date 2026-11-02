"""Manuscript numbers, figures and benchmark splits are derived, not typed.

paper/stats.json is what the manuscript cites; it must equal a fresh
computation from the database. data/splits/benchmark_splits.json must never
put identical measurements on both sides of a split.
"""

import importlib.util
import json
import re
import unittest
from pathlib import Path

from rel_quest._release import PROFILE

ROOT = Path(__file__).resolve().parent.parent


def load(path, name):
    spec = importlib.util.spec_from_file_location(name, ROOT / path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def records():
    path = ROOT / "data" / "sx_ligands_merged_v14_flat.jsonl"
    return [json.loads(line) for line in path.open(encoding="utf-8") if line.strip()]


class TestPaperStats(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.stats_module = load("paper/stats.py", "paper_stats")
        cls.fresh = cls.stats_module.compute_stats()

    def test_committed_stats_are_current(self):
        # Without source labels (the public release) the DEV_ONLY_KEYS are carried over
        # from the committed file and everything else must recompute.
        committed = (ROOT / "paper" / "stats.json").read_text(encoding="utf-8")
        fresh = self.stats_module.complete(self.fresh, json.loads(committed))
        self.assertEqual(committed, self.stats_module.render(fresh),
                         "paper/stats.json is stale — run python paper/stats.py")

    def test_counts_are_internally_consistent(self):
        s = self.stats_module.complete(self.fresh, self.stats_module.load_stored())
        self.assertEqual(sum(s["ligands_by_data_source"].values()), s["n_ligands"])
        self.assertEqual(sum(s["conditions_by_origin"].values()), s["n_conditions"])
        self.assertEqual(sum(s["d_values_by_origin"].values()), s["n_d_values"])
        self.assertEqual(sum(s["trend_classes"].values()), s["n_ligands"])
        self.assertEqual(sum(s["mechanisms"].values()), s["n_ligands"])
        self.assertEqual(sum(s["provenance_level_conditions"].values()), s["n_conditions"])

    def test_dev_only_keys_are_computed_where_labels_exist(self):
        if not self.stats_module.has_source_labels(records()):
            self.skipTest("public release: no source labels")
        for key in self.stats_module.DEV_ONLY_KEYS:
            with self.subTest(key=key):
                self.assertIn(key, self.fresh)

    def test_readme_headline_numbers_match(self):
        if PROFILE == "public":
            self.skipTest("the public README carries no headline numbers")
        readme = (ROOT / "README.md").read_text(encoding="utf-8")
        s = self.fresh
        for label, value in (("unique ligands", s["n_ligands"]),
                             ("experimental conditions", s["n_conditions"]),
                             ("D-values", s["n_d_values"])):
            match = re.search(rf"\*\*([\d,]+) {re.escape(label)}\*\*", readme)
            with self.subTest(label=label):
                self.assertIsNotNone(match, f"README headline lacks '**N {label}**'")
                self.assertEqual(int(match.group(1).replace(",", "")), value)
        for mechanism, n in s["mechanisms"].items():
            with self.subTest(mechanism=mechanism):
                self.assertIn(f"{mechanism} {n}", readme)


class TestBenchmarkSplits(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.module = load("scripts/make_benchmark_splits.py", "make_benchmark_splits")
        # The splits are checked against the database they were drawn from (the script's DB_PATH)
        with cls.module.DB_PATH.open(encoding="utf-8") as h:
            cls.records = [json.loads(line) for line in h if line.strip()]
        cls.committed = json.loads((ROOT / "data" / "splits" / "benchmark_splits.json").read_text())

    def test_committed_splits_are_current(self):
        self.assertEqual(self.committed, json.loads(json.dumps(self.module.build(self.records))))

    def test_every_ligand_in_exactly_one_fold(self):
        names = [e["name"] for fold in self.committed["grouped_5fold"]["folds"] for e in fold]
        self.assertEqual(len(names), len(set(names)))
        self.assertEqual(set(names), {r["ligand"]["Name"] for r in self.records})

    def test_identical_measurements_never_straddle_folds(self):
        fold_of = {e["name"]: i for i, fold in enumerate(self.committed["grouped_5fold"]["folds"])
                   for e in fold}
        owners = {}
        for r in self.records:
            for c in r["experimental_data"]:
                owners.setdefault(self.module.condition_key(c), set()).add(r["ligand"]["Name"])
        straddling = [sorted(n) for n in owners.values() if len({fold_of[x] for x in n}) > 1]
        self.assertEqual(straddling, [])


class TestFigures(unittest.TestCase):

    def test_every_figure_builds(self):
        try:
            import matplotlib  # noqa: F401
        except ImportError:  # pragma: no cover - optional `paper` extra
            self.skipTest("matplotlib not installed (the `paper` extra)")
        figures = load("paper/figures.py", "paper_figures")
        made = figures.build(export=False)
        self.assertEqual(set(made), set(figures.FIGURES))


if __name__ == "__main__":
    unittest.main()
