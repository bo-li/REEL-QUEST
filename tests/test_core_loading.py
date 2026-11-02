"""Loading must be all-or-nothing per record.

PROJECT_REVIEW_2026-07-31.md finding 6. `_build_index()` appended a record to
`self._records` and *then* read its required fields, both inside one broad
`except`. A record missing `ligand.Mechanism` was therefore warned about,
left in `_records`, and never indexed — so `get_statistics()` counted a
ligand no query could return, and `get_all_ligands()` raised `KeyError`
on it later, far from the cause.

The loader now validates first and appends-and-indexes together, so a record
is either fully present or fully absent and recorded in `load_errors`.

These build tiny database files rather than mutating the canonical one.
"""

import json
import tempfile
import unittest
import warnings
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
DB_PATH = ROOT / "data" / "sx_ligands_merged_v14_flat.jsonl"


def a_valid_record():
    with DB_PATH.open() as handle:
        return json.loads(next(line for line in handle if line.strip()))


class TestPartialLoadIsImpossible(unittest.TestCase):

    def setUp(self):
        self.good = a_valid_record()
        self._tmp = tempfile.TemporaryDirectory()
        self.tmp = Path(self._tmp.name)

    def tearDown(self):
        self._tmp.cleanup()

    def _load(self, *records, **kwargs):
        from rel_quest import LigandDatabase

        path = self.tmp / "probe.jsonl"
        path.write_text("".join(json.dumps(r) + "\n" for r in records))
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            db = LigandDatabase(str(path), **kwargs)
        return db, caught

    # -- the malformed shapes that used to slip through -------------------

    def _malformed_cases(self):
        return {
            "missing ligand.Mechanism": {
                **self.good,
                "ligand": {k: v for k, v in self.good["ligand"].items()
                           if k != "Mechanism"}},
            "missing ligand.Name": {
                **self.good,
                "ligand": {k: v for k, v in self.good["ligand"].items()
                           if k != "Name"}},
            "missing selectivity block": {
                k: v for k, v in self.good.items() if k != "selectivity"},
            "selectivity is a string": {**self.good, "selectivity": "oops"},
            "ligand is a list": {**self.good, "ligand": []},
            "record is a list": [1, 2, 3],
            "record is a bare string": "not a record",
            "missing trend_classification": {
                **self.good,
                "selectivity": {k: v for k, v in self.good["selectivity"].items()
                                if k != "trend_classification"}},
        }

    def test_a_malformed_record_is_absent_everywhere(self):
        for label, bad in self._malformed_cases().items():
            with self.subTest(case=label):
                db, _ = self._load(bad, self.good)
                self.assertEqual(len(db._records), 1, "malformed record was kept")
                self.assertEqual(db.get_statistics()["total_ligands"], 1)
                self.assertEqual(len(db.get_all_ligands()), 1)
                self.assertEqual(len(db.load_errors), 1)

    def test_get_all_ligands_no_longer_raises(self):
        """The reported symptom: statistics counted it, conversion then failed."""
        for label, bad in self._malformed_cases().items():
            with self.subTest(case=label):
                db, _ = self._load(bad, self.good)
                try:
                    ligands = db.get_all_ligands()
                except Exception as exc:  # pragma: no cover - the regression
                    self.fail(f"get_all_ligands() raised {type(exc).__name__}: {exc}")
                self.assertEqual(len(ligands), 1)

    def test_statistics_and_records_always_agree(self):
        for label, bad in self._malformed_cases().items():
            with self.subTest(case=label):
                db, _ = self._load(bad, self.good, self.good)
                self.assertEqual(db.get_statistics()["total_ligands"],
                                 len(db.get_all_ligands()))

    # -- how the problem is surfaced --------------------------------------

    def test_skipping_warns_with_line_and_reason(self):
        bad = {**self.good, "selectivity": "oops"}
        db, caught = self._load(bad, self.good)
        self.assertTrue(caught, "skipping a record must warn")
        message = str(caught[-1].message)
        self.assertIn("line 1", message)
        self.assertIn("could not", message)
        self.assertEqual(db.load_errors[0]["line"], 1)
        self.assertIn("selectivity", db.load_errors[0]["reason"])

    def test_strict_mode_raises_instead_of_skipping(self):
        from rel_quest import DatabaseLoadError, LigandDatabase

        path = self.tmp / "strict.jsonl"
        bad = {**self.good, "ligand": {}}
        path.write_text(json.dumps(bad) + "\n" + json.dumps(self.good) + "\n")
        with self.assertRaises(DatabaseLoadError) as ctx:
            LigandDatabase(str(path), strict=True)
        self.assertIn("line 1", str(ctx.exception))

    def test_invalid_json_is_reported_with_its_line_number(self):
        from rel_quest import LigandDatabase

        path = self.tmp / "badjson.jsonl"
        path.write_text(json.dumps(self.good) + "\n{ not json\n"
                        + json.dumps(self.good) + "\n")
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            db = LigandDatabase(str(path))
        self.assertEqual(len(db.get_all_ligands()), 2)
        self.assertEqual([e["line"] for e in db.load_errors], [2])
        self.assertIn("invalid JSON", db.load_errors[0]["reason"])

    # -- a clean load must stay clean and silent --------------------------

    def test_a_clean_file_produces_no_errors_and_no_warning(self):
        db, caught = self._load(self.good)
        self.assertEqual(db.load_errors, [])
        self.assertEqual(
            [w for w in caught if issubclass(w.category, RuntimeWarning)], [])

    def test_blank_lines_are_not_errors(self):
        from rel_quest import LigandDatabase

        path = self.tmp / "blanks.jsonl"
        path.write_text("\n" + json.dumps(self.good) + "\n\n")
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            db = LigandDatabase(str(path))
        self.assertEqual(len(db.get_all_ligands()), 1)
        self.assertEqual(db.load_errors, [])
        self.assertEqual(
            [w for w in caught if issubclass(w.category, RuntimeWarning)], [])

    def test_the_canonical_database_loads_without_errors(self):
        from rel_quest import LigandDatabase

        db = LigandDatabase()
        self.assertEqual(db.load_errors, [])
        self.assertEqual(db.get_statistics()["total_ligands"],
                         len(db.get_all_ligands()))


if __name__ == "__main__":
    unittest.main()
