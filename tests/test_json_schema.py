"""The published JSON Schema describes the data and cannot go stale.

schema/ligand_record.schema.json and schema/data_dictionary.csv are generated
by scripts/generate_json_schema.py from the code's own vocabularies. These
tests check that the files on disk equal a fresh render, that every stored
record and the example template validate, and that the schema rejects the
defects it exists to catch.
"""

import copy
import csv
import importlib.util
import json
import unittest
from pathlib import Path

from jsonschema import Draft202012Validator
from rel_quest.schema import PROFILE

ROOT = Path(__file__).resolve().parent.parent
SCHEMA_PATH = ROOT / "schema" / "ligand_record.schema.json"
DICTIONARY_PATH = ROOT / "schema" / "data_dictionary.csv"
DB_PATH = ROOT / "data" / "sx_ligands_merged_v14_flat.jsonl"


def load_generator():
    spec = importlib.util.spec_from_file_location(
        "generate_json_schema", ROOT / "scripts" / "generate_json_schema.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def records():
    return [json.loads(line) for line in DB_PATH.open(encoding="utf-8") if line.strip()]


class TestSchemaFiles(unittest.TestCase):

    def test_files_on_disk_match_a_fresh_render(self):
        schema_text, dictionary_text = load_generator().render()
        self.assertEqual(SCHEMA_PATH.read_text(encoding="utf-8"), schema_text,
                         "run scripts/generate_json_schema.py")
        self.assertEqual(DICTIONARY_PATH.read_text(encoding="utf-8"), dictionary_text,
                         "run scripts/generate_json_schema.py")

    def test_schema_is_valid_draft_2020_12(self):
        Draft202012Validator.check_schema(json.loads(SCHEMA_PATH.read_text()))

    def test_enumerations_come_from_the_code(self):
        from rel_quest.schema import TREND_CLASSIFICATIONS, canonical_mechanisms
        schema = json.loads(SCHEMA_PATH.read_text())
        props = schema["properties"]
        self.assertEqual(props["selectivity"]["properties"]["trend_classification"]["enum"],
                         list(TREND_CLASSIFICATIONS))
        self.assertEqual(props["ligand"]["properties"]["Mechanism"]["enum"],
                         list(canonical_mechanisms()))
        subtypes = json.loads((ROOT / "vocab" / "chemical_families.json").read_text())["subtypes"]
        self.assertEqual(set(props["ligand"]["properties"]["Subtype"]["enum"]), set(subtypes))

    def test_dictionary_documents_units_for_every_concentration(self):
        with DICTIONARY_PATH.open(encoding="utf-8") as handle:
            rows = {row["path"]: row for row in csv.DictReader(handle)}
        for path in ("experimental_data[].Acid_Concentration_M",
                     "experimental_data[].Extractant_Concentration_M",
                     "experimental_data[].Temperature_C",
                     "experimental_data[].distribution_ratios.<element>.sigma_log_D",
                     "selectivity.slope"):
            with self.subTest(path=path):
                self.assertTrue(rows[path]["units"], f"{path} has no units")
        undocumented = [path for path, row in rows.items() if not row["description"]]
        self.assertEqual(undocumented, [])


class TestStoredDataValidates(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.validator = Draft202012Validator(json.loads(SCHEMA_PATH.read_text()))

    def test_every_record(self):
        failures = {}
        for record in records():
            errors = [e.message for e in self.validator.iter_errors(record)]
            if errors:
                failures[record["ligand"]["Name"]] = errors[:2]
        self.assertEqual(failures, {}, list(failures.items())[:3])

    def test_template_is_a_valid_current_record(self):
        template = json.loads((ROOT / "data" / "template.json").read_text())
        self.assertEqual([e.message for e in self.validator.iter_errors(template)], [])

    def test_template_is_the_stored_record_it_claims_to_be(self):
        template = json.loads((ROOT / "data" / "template.json").read_text())
        stored = {r["ligand"]["Name"]: r for r in records()}
        self.assertEqual(template, stored[template["ligand"]["Name"]],
                         "data/template.json has drifted from the database; regenerate it")


class TestSchemaRejectsDefects(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.validator = Draft202012Validator(json.loads(SCHEMA_PATH.read_text()))
        cls.base = json.loads((ROOT / "data" / "template.json").read_text())

    def assertRejected(self, mutate):
        record = copy.deepcopy(self.base)
        mutate(record)
        self.assertTrue(list(self.validator.iter_errors(record)), "schema accepted a defect")

    def test_comma_joined_doi(self):
        def mutate(r):
            r["experimental_data"][0]["reference"]["DOI"] = (
                "10.1000/a, https://doi.org/10.1021/jacs.5c19738")
        self.assertRejected(mutate)

    def test_unknown_top_level_field(self):
        self.assertRejected(lambda r: r.update({"extra_field": 1}))

    def test_unknown_element_symbol(self):
        def mutate(r):
            ratios = r["experimental_data"][0]["distribution_ratios"]
            ratios["Am"] = dict(next(iter(ratios.values())))
        self.assertRejected(mutate)

    def test_non_positive_distribution_ratio(self):
        def mutate(r):
            next(iter(r["experimental_data"][0]["distribution_ratios"].values()))["value"] = 0
        self.assertRejected(mutate)

    @unittest.skipIf(PROFILE == "public", "source labels are in the developer repository only")
    def test_primary_row_with_an_aggregator(self):
        def mutate(r):
            r["experimental_data"][0]["reference"]["Aggregator"] = {
                "name": "SAFE_01262026", "DOI": None, "URL": None, "entry": None,
                "curation_DOI": None}
        # Pramanik rows are primary, so an Aggregator contradicts the level.
        self.assertEqual(self.base["experimental_data"][0]["reference"]["Provenance_Level"],
                         "primary")
        self.assertRejected(mutate)

    def test_trend_outside_the_vocabulary(self):
        self.assertRejected(lambda r: r["selectivity"].update({"trend_classification": "Wavy"}))


class TestValidatorUsesTheSchema(unittest.TestCase):

    def test_validate_database_file_reports_schema_errors(self):
        import tempfile
        from rel_quest.validators import validate_database_file

        record = json.loads((ROOT / "data" / "template.json").read_text())
        record["extra_field"] = 1
        with tempfile.NamedTemporaryFile("w", suffix=".jsonl", delete=False) as handle:
            handle.write(json.dumps(record) + "\n")
        try:
            valid, invalid, errors = validate_database_file(handle.name)
            self.assertEqual((valid, invalid), (0, 1))
            self.assertIn("schema:", errors[0])
            valid, invalid, _ = validate_database_file(handle.name, use_json_schema=False)
            self.assertEqual((valid, invalid), (1, 0))
        finally:
            Path(handle.name).unlink()


if __name__ == "__main__":
    unittest.main()
