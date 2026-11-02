"""The data-only release bundle is complete, checksummed and reproducible.

scripts/export_release_bundle.py builds what gets deposited on Zenodo and
cited in the paper's Data Records, including the SAFE-excluded fallback
described in docs/DATA_PROVENANCE_AND_TERMS.md.
"""

import csv
import hashlib
import importlib.util
import json
import tempfile
import unittest
from pathlib import Path
from rel_quest.schema import PROFILE

ROOT = Path(__file__).resolve().parent.parent


def exporter():
    spec = importlib.util.spec_from_file_location(
        "export_release_bundle", ROOT / "scripts" / "export_release_bundle.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class TestReleaseBundle(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.tmp = tempfile.TemporaryDirectory()
        cls.module = exporter()
        cls.bundle, cls.archive, cls.records, cls.n_cells, _ = cls.module.build(Path(cls.tmp.name))

    @classmethod
    def tearDownClass(cls):
        cls.tmp.cleanup()

    def test_contents(self):
        for name in ("sx_ligands_merged_v14_flat.jsonl", "reel_quest_measurements.csv",
                     "reel_quest_ligands.csv", "schema/ligand_record.schema.json",
                     "schema/data_dictionary.csv", "LICENSE-DATA", "CITATION.cff",
                     "README.md", "SHA256SUMS", "zenodo_dataset.json",
                     "DATA_PROVENANCE_AND_TERMS.md", "references.bib"):
            with self.subTest(name=name):
                self.assertTrue((self.bundle / name).is_file())

    def test_database_is_copied_unchanged(self):
        self.assertEqual((self.bundle / "sx_ligands_merged_v14_flat.jsonl").read_bytes(),
                         (ROOT / "data" / "sx_ligands_merged_v14_flat.jsonl").read_bytes())

    def test_measurement_csv_has_one_row_per_d_value(self):
        stored = sum(len(c["distribution_ratios"]) for r in self.records
                     for c in r["experimental_data"])
        with (self.bundle / "reel_quest_measurements.csv").open(encoding="utf-8") as handle:
            rows = list(csv.DictReader(handle))
        self.assertEqual(len(rows), stored)
        self.assertEqual(self.n_cells, stored)
        self.assertTrue(all(row["doi"] or row["url"] for row in rows))

    def test_checksums_verify(self):
        for line in (self.bundle / "SHA256SUMS").read_text().splitlines():
            digest, name = line.split("  ", 1)
            with self.subTest(name=name):
                self.assertEqual(hashlib.sha256((self.bundle / name).read_bytes()).hexdigest(),
                                 digest)

    def test_zip_is_byte_reproducible(self):
        first = self.archive.read_bytes()
        with tempfile.TemporaryDirectory() as other:
            _, again, _, _, _ = self.module.build(Path(other))
            self.assertEqual(again.read_bytes(), first)

    def test_dataset_metadata_is_cc_by(self):
        zenodo = json.loads((self.bundle / "zenodo_dataset.json").read_text())
        self.assertEqual((zenodo["upload_type"], zenodo["license"]), ("dataset", "cc-by-4.0"))


@unittest.skipIf(PROFILE == "public", "source labels are in the developer repository only")
class TestExcludeAggregator(unittest.TestCase):

    def test_without_safe(self):
        module = exporter()
        with tempfile.TemporaryDirectory() as tmp:
            bundle, _, records, n_cells, dropped = module.build(Path(tmp), ["SAFE_01262026"])
            names = {(c["reference"].get("Aggregator") or {}).get("name")
                     for r in records for c in r["experimental_data"]}
            self.assertNotIn("SAFE_01262026", names)
            with (ROOT / "data" / "sx_ligands_merged_v14_flat.jsonl").open(encoding="utf-8") as h:
                total = sum(1 for line in h if line.strip())
            self.assertEqual(len(dropped) + len(records), total)
            self.assertTrue(all(r["metadata"]["data_source"] != "SAFE_01262026" for r in records))
            # Derived fields describe what remains, not what was dropped.
            from rel_quest.derived import derived_differences, recompute_record
            stale = [r["ligand"]["Name"] for r in records
                     if derived_differences(r, recompute_record(r))]
            self.assertEqual(stale, [])
            self.assertIn("Excluded sources", (bundle / "README.md").read_text())


if __name__ == "__main__":
    unittest.main()
