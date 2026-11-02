"""End to end: the live contribution route.

One file per paper (contributions/<doi-slug>.json, grammar in
rel_quest/contributions.py and schema/contribution.schema.json) is checked by
scripts/check_contribution.py (the CI job on a data pull request) and applied to
the released database by scripts/apply_contributions.py once a maintainer has
reviewed it. Ways this could corrupt the database or mislead a contributor, each
pinned below:

 1. A DOI that is not a DOI, or a compilation's DOI (SAFE, SAFE-MolGen), passes
    the check.
 2. A new ligand whose structure the database already holds (same InChIKey,
    under another name) creates a second record for one molecule; one that
    shares its connectivity with a stored record (a stereoisomer) passes
    without a warning.
 3. A cell given as both D and log_D that disagree is written as one of them;
    a D of zero or less is written at all.
 4. A condition whose cells equal a stored condition of the same ligand is
    added again (double counting); one equal to another ligand's goes
    unremarked.
 5. A condition that cannot be applied as written (no locator, a ligand that
    no record or new_ligands entry names, a condition field the schema has
    not got, a cell without value_origin) is skipped silently or applied with
    a guess.
 6. The file violates the contribution schema (no proposed_by, an id that is
    not A + four digits, a Subtype outside the vocabulary) and is applied.
 7. A new ligand with an unparseable SMILES, or one that no condition in the
    file uses, creates a record.
 8. The check changes the database or the contribution file.
 9. Applying stamps last_updated on records the contribution did not touch,
    or leaves it unstamped on the records it did.
10. The applied records fail the validators or the JSON Schema, or a derived
    field no longer reproduces (rel_quest.derived.recompute_record is not a
    fixed point), or every derived field is not recomputed.
11. A `Lit_` condition arrives without primary provenance: a bare lower-case
    DOI, the locator, Provenance_Level 'primary', no aggregator, and a
    per-condition data_source label on a record from another source.
12. A trend label changes on an existing record and the report does not say
    so; a new record's class is not reported.
13. The same file is applied twice: the file is marked applied, its
    condition_ids exist, and the second apply must refuse.
14. Nothing is written when any file in a run has a problem (all or
    nothing); the audit file does not list what was done.
17. The review rule (GOVERNANCE.md §4) has no enforcement point: a file that
    names no reviewer is applied without a trace, or the reviewer's name does
    not reach the data (the new record's curator, each condition's notes).
15. The contribution schema on disk is stale (differs from the module's
    render), or the literature additions that built v14
    (outputs/verification/additions/) no longer fit the grammar (the two
    routes must share one format).
18. A new ligand whose every condition holds a single element (Eu in an
    Am/Eu study, say) passes without the reviewer being told: GOVERNANCE.md
    §1 admits single-element series only for ligands already in the
    database. Single-element conditions on a stored ligand, or a new ligand
    with one multi-element series beside single-element ones, are not flagged.
16. The refactor of scripts/apply_verification.py (shared helpers moved to
    rel_quest.contributions) changes what the v14 build does: pinned by
    tests/test_apply_verification.py, which must keep passing unchanged.

The repeatable artifact of a run is outputs/contributions/apply_<date>.txt
(here written into a temporary directory and read back).
"""

import copy
import importlib.util
import json
import tempfile
import unittest
from datetime import date
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
DB = ROOT / "data" / "sx_ligands_merged_v14_flat.jsonl"
ADDITIONS = ROOT / "outputs" / "verification" / "additions"
FIXTURE_RECORDS = ("HONTA", "DMDOHEMA", "BLPhen")
DOI = "10.1234/fixture"
SAFE_DOI = "10.1080/07366299.2025.2564381"
# A diglycolamide no compilation carries: tetrakis(cyclopropylmethyl)-DGA.
NEW_SMILES = "C1CC1CN(CC1CC1)C(=O)COCC(=O)N(CC1CC1)CC1CC1"


def load_script(name):
    spec = importlib.util.spec_from_file_location(name, ROOT / "scripts" / f"{name}.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def fixture_paper():
    return {
        "DOI": DOI,
        "citation": "Fixture, Example, J. Fixtures 2026, 1, 1",
        "proposed_by": "tests/test_contributions.py, from Table 2 and Figure 3",
        "reviewed_by": "a second reader, 2026-10-20, against the fixture paper",
        "new_ligands": [{
            "Name": "FixtureDGA", "SMILES": NEW_SMILES, "Mechanism": "Solvation",
            "Extractant_Type": "Diglycolamide", "Subtype": "DGA", "Chemical_Family": "DGA",
            "reason": "the paper's ligand 1 (fixture)"}],
        "conditions": [
            {"id": "A0001", "ligand": "HONTA",
             "conditions": {"Acid_Name": "HNO3", "Acid_Concentration_M": 1.0, "Extractant_Concentration_M": 0.1,
                            "Solvent_Name": "n-dodecane", "Temperature_C": 25.0, "Contact_Time": "30 min"},
             "additives": "0.5 M NaNO3",
             "reference": {"Locator": "Table 2"},
             "cells": {"La": {"log_D": -0.5, "value_origin": "reported"},
                       "Eu": {"D": 2.0, "value_origin": "reported"},
                       "Nd": {"log_D": 0.1, "sigma_log_D": 0.05, "value_origin": "digitised"},
                       "Lu": {"censor": ">", "log_bound": 2.0}},
             "why": "Table 2, the 1 M HNO3 row (fixture)", "note": None},
            {"id": "A0002", "ligand": "FixtureDGA",
             "conditions": {"Acid_Name": "HNO3", "Acid_Concentration_M": 3.0, "Extractant_Concentration_M": 0.2,
                            "Solvent_Name": "n-dodecane", "Temperature_C": 25.0},
             "reference": {"Locator": "Figure 3"},
             "cells": {el: {"log_D": round(-1.0 + 0.25 * i, 3), "value_origin": "digitised"}
                       for i, el in enumerate(["La", "Nd", "Eu", "Tb", "Er", "Lu"])},
             "why": "Figure 3, the 3 M HNO3 series (fixture)", "note": None},
        ],
    }


class TestContributionRoute(unittest.TestCase):

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        root = Path(self.tmp.name)
        self.db = root / "db.jsonl"
        self.records = [r for r in map(json.loads, DB.open(encoding="utf-8")) if r["ligand"]["Name"] in FIXTURE_RECORDS]
        self.assertEqual(len(self.records), 3)
        self.db.write_text("".join(json.dumps(r) + "\n" for r in self.records), encoding="utf-8")
        self.dir, self.reports, self.audits, self.archive = (root / d for d in ("contributions", "reports", "audits", "archive"))
        self.dir.mkdir()
        self.check = load_script("check_contribution")
        # Applying is a maintainer's step in the developer repository; the public one only checks
        self.apply = (load_script("apply_contributions")
                      if (ROOT / "scripts" / "apply_contributions.py").exists() else None)
        from rel_quest import contributions
        self.lib = contributions

    def tearDown(self):
        self.tmp.cleanup()

    def write(self, paper, name="10.1234_fixture.json"):
        path = self.dir / name
        path.write_text(json.dumps(paper, indent=2) + "\n", encoding="utf-8")
        return path

    def run_check(self, *files):
        return self.check.main([*map(str, files), "--database", str(self.db), "--report-dir", str(self.reports),
                                "--offline"])

    def run_apply(self, *files, dry_run=False):
        if self.apply is None:
            self.skipTest("scripts/apply_contributions.py is in the developer repository only")
        return self.apply.main([*map(str, files), "--database", str(self.db), "--report-dir", str(self.reports),
                                "--audit-dir", str(self.audits), "--archive-dir", str(self.archive),
                                "--offline", *(["--dry-run"] if dry_run else [])])

    def report(self, path):
        return (self.reports / f"{Path(path).stem}_check.md").read_text(encoding="utf-8")

    def stored(self):
        return {r["ligand"]["Name"]: r for r in map(json.loads, self.db.open(encoding="utf-8"))}

    # --- the route end to end -------------------------------------------------

    def test_check_then_apply(self):
        from rel_quest.derived import derived_differences, recompute_record
        from rel_quest.validators import validate_against_json_schema, validate_ligand_record
        path = self.write(fixture_paper())
        db_before, file_before = self.db.read_bytes(), path.read_bytes()
        before = self.stored()

        self.assertEqual(self.run_check(path), 0)
        report = self.report(path)
        self.assertIn("No problems", report)
        self.assertIn("HONTA", report)
        self.assertEqual(self.db.read_bytes(), db_before, "the check changed the database")      # 8
        self.assertEqual(path.read_bytes(), file_before, "the check changed the file")

        self.assertEqual(self.run_apply(path, dry_run=True), 0)
        self.assertEqual(self.db.read_bytes(), db_before, "a dry run wrote the database")
        self.assertEqual(path.read_bytes(), file_before, "a dry run marked the file")

        self.assertEqual(self.run_apply(path), 0)
        after = self.stored()
        today = date.today().isoformat()
        self.assertEqual(set(after), set(before) | {"FixtureDGA"})

        # 11: the HONTA condition carries primary provenance and a source label
        honta = after["HONTA"]
        cid = f"Lit_HONTA_10.1234_fixture_A0001"
        cond = next(c for c in honta["experimental_data"] if c["condition_id"] == cid)
        ref = cond["reference"]
        self.assertEqual((ref["DOI"], ref["Locator"], ref["Provenance_Level"], ref.get("Aggregator")),
                         (DOI, "Table 2", "primary", None))
        self.assertEqual(cond["data_source"], "Literature_curated")
        self.assertEqual(cond["Additives"], "0.5 M NaNO3")
        cells = cond["distribution_ratios"]
        self.assertAlmostEqual(cells["La"]["log_value"], -0.5)
        self.assertAlmostEqual(cells["La"]["value"], 10 ** -0.5, places=5)
        self.assertAlmostEqual(cells["Eu"]["log_value"], 0.301, places=3)          # given as D
        self.assertEqual((cells["Lu"]["censor"], cells["Lu"]["value_origin"]), (">", "censored"))
        self.assertEqual((cells["Nd"]["sigma_log_D"], cells["Nd"]["sigma_source"]), (0.05, "measured"))
        self.assertEqual(cells["La"]["source_quality"], "table")
        self.assertIn("mean_log_d", cond["statistics"])
        self.assertIn(cid, [c["condition_id"] for c in honta["experimental_data"]])
        self.assertEqual(len(honta["experimental_data"]), len(before["HONTA"]["experimental_data"]) + 1)

        # 9: only the touched records are re-dated
        self.assertEqual(honta["metadata"]["last_updated"], today)
        for name in ("DMDOHEMA", "BLPhen"):
            self.assertEqual(after[name], before[name], f"{name} changed although the file never names it")

        # 7 / 11: the new record
        new = after["FixtureDGA"]
        self.assertTrue(new["ligand"]["InChI_Key"])
        self.assertEqual(new["ligand"]["Molecular_Formula"], "C20H32N2O3")
        self.assertEqual(new["metadata"]["data_source"], "Literature_curated")
        self.assertEqual((new["metadata"]["added_date"], new["metadata"]["last_updated"]), (today, today))
        self.assertEqual(new["metadata"]["curation_status"], "Validated")
        self.assertIn("tests/test_contributions.py", new["metadata"]["curator"])
        self.assertIn("reviewed by a second reader", new["metadata"]["curator"])             # 17
        self.assertIn("reviewed by a second reader", ref["Notes"])
        self.assertEqual([c["condition_id"] for c in new["experimental_data"]],
                         ["Lit_FixtureDGA_10.1234_fixture_A0002"])
        self.assertNotIn("data_source", new["experimental_data"][0])    # the record's own source
        self.assertEqual(new["selectivity"]["trend_classification"], "Ascending")

        # 10: every record validates and every derived field reproduces
        for r in after.values():
            ok, errors = validate_ligand_record(r)
            errors += validate_against_json_schema(r)
            self.assertEqual(errors, [], r["ligand"]["Name"])
            self.assertEqual(derived_differences(r, recompute_record(r)), [], r["ligand"]["Name"])

        # 13: the file is marked and a second apply refuses
        marked = json.loads(path.read_text(encoding="utf-8"))
        self.assertEqual(marked["status"], "applied")
        self.assertEqual(marked["applied_on"], today)
        self.assertEqual(marked["condition_ids"], [cid, "Lit_FixtureDGA_10.1234_fixture_A0002"])
        self.assertEqual(self.lib.pending_files(self.dir), [])
        db_applied = self.db.read_bytes()
        self.assertNotEqual(self.run_apply(path), 0)
        self.assertEqual(self.db.read_bytes(), db_applied)
        self.assertIn("already", self.report(path))

        # 14: the audit lists what was done; 12: labels before and after
        audits = sorted(self.audits.glob("apply_*.txt"))
        self.assertEqual(len(audits), 1)
        audit = audits[0].read_text(encoding="utf-8")
        for text in (cid, "Lit_FixtureDGA_10.1234_fixture_A0002", "FixtureDGA: new record (Ascending)",
                     f"HONTA: {before['HONTA']['selectivity']['trend_classification']} -> "
                     f"{honta['selectivity']['trend_classification']}", "CHANGELOG"):
            self.assertIn(text, audit)
        backups = list(self.archive.glob("*.jsonl"))
        self.assertEqual(len(backups), 1)
        self.assertEqual(backups[0].read_bytes(), db_before)

    # --- the ways a file is refused -------------------------------------------

    def refused(self, mutate, expected, name="10.1234_fixture.json"):
        paper = fixture_paper()
        mutate(paper)
        path = self.write(paper, name)
        db_before = self.db.read_bytes()
        code = self.run_check(path)
        report = self.report(path)
        self.assertEqual(code, 1, f"accepted a file that should be refused: {expected!r}\n{report}")
        self.assertIn(expected, report)
        self.assertEqual(self.db.read_bytes(), db_before)
        self.assertNotEqual(self.run_apply(path), 0)
        self.assertEqual(self.db.read_bytes(), db_before)
        self.assertNotIn("status", json.loads(path.read_text(encoding="utf-8")))
        return report

    def test_refusals(self):
        honta = self.records[[r["ligand"]["Name"] for r in self.records].index("HONTA")]
        stored = next(c for c in honta["experimental_data"] if len(c["distribution_ratios"]) >= 3
                      and not any(v.get("censor") for v in c["distribution_ratios"].values()))
        same_cells = {el: {"log_D": v["log_value"], "value_origin": "reported"}
                      for el, v in stored["distribution_ratios"].items()}

        def set_(key, value, where="paper"):
            def f(p):
                target = p if where == "paper" else p["conditions"][0] if where == "cond" else p["new_ligands"][0]
                if value is KeyError:
                    del target[key]
                else:
                    target[key] = value
            return f

        cases = [
            (set_("DOI", "not-a-doi"), "DOI"),                                                      # 1
            (set_("DOI", SAFE_DOI), "compilation"),                                                 # 1
            (set_("SMILES", honta["ligand"]["SMILES"], "lig"), "same structure as HONTA"),          # 2
            (lambda p: p["conditions"][0]["cells"].update(Eu={"D": 2.0, "log_D": 0.9, "value_origin": "reported"}),
             "disagree"),                                                                            # 3
            (lambda p: p["conditions"][0]["cells"].update(Eu={"D": 0.0, "value_origin": "reported"}), "positive"),
            (lambda p: p["conditions"][0].update(cells=same_cells), f"duplicate of {stored['condition_id']}"),  # 4
            (set_("reference", {}, "cond"), "Locator"),                                             # 5
            (set_("ligand", "Nobody", "cond"), "no record"),
            (set_("conditions", {"Acid": "HNO3"}, "cond"), "cannot be set"),
            (lambda p: p["conditions"][0]["cells"].update(La={"log_D": -0.5}), "value_origin"),
            (set_("proposed_by", KeyError), "proposed_by"),                                          # 6
            (set_("id", "X1", "cond"), "id"),
            (set_("Subtype", "Foo", "lig"), "Subtype"),
            (set_("SMILES", "C1CC", "lig"), "SMILES"),                                              # 7
            (lambda p: p["conditions"].pop(1), "no condition"),
        ]
        for i, (mutate, expected) in enumerate(cases):
            with self.subTest(expected=expected):
                self.refused(mutate, expected, name=f"10.1234_case{i}.json")

    def test_duplicate_of_another_ligands_condition_is_a_warning(self):
        dmdohema = self.records[[r["ligand"]["Name"] for r in self.records].index("DMDOHEMA")]
        stored = next(c for c in dmdohema["experimental_data"] if len(c["distribution_ratios"]) >= 3
                      and not any(v.get("censor") for v in c["distribution_ratios"].values()))
        paper = fixture_paper()
        paper["conditions"][0]["cells"] = {el: {"log_D": v["log_value"], "value_origin": "reported"}
                                           for el, v in stored["distribution_ratios"].items()}
        path = self.write(paper)
        self.assertEqual(self.run_check(path), 0)
        self.assertIn(f"same cells as DMDOHEMA's {stored['condition_id']}", self.report(path))

    def test_new_ligand_with_only_single_element_series_is_a_warning(self):        # 18
        def eu_only(i, ligand, acid):
            return {"id": f"A000{i}", "ligand": ligand,
                    "conditions": {"Acid_Name": "HNO3", "Acid_Concentration_M": acid,
                                   "Extractant_Concentration_M": 0.2, "Solvent_Name": "n-dodecane"},
                    "reference": {"Locator": "Figure 4"},
                    "cells": {"Eu": {"log_D": round(0.1 * acid, 3), "value_origin": "digitised"}},
                    "why": "Figure 4, Eu against [HNO3] (fixture)", "note": None}

        flag = "FixtureDGA: every condition holds a single element"
        sparse = fixture_paper()
        sparse["conditions"] = [sparse["conditions"][0], eu_only(2, "FixtureDGA", 1.0), eu_only(3, "FixtureDGA", 3.0)]
        path = self.write(sparse, name="10.1234_sparse.json")
        self.assertEqual(self.run_check(path), 0)          # a warning, not a refusal: the maintainer decides
        self.assertIn(flag, self.report(path))

        mixed = fixture_paper()                            # one multi-element series beside an Eu-only one
        mixed["conditions"].append(eu_only(3, "FixtureDGA", 1.0))
        mixed["conditions"].append(eu_only(4, "HONTA", 2.0))   # single-element data on a stored ligand
        path = self.write(mixed, name="10.1234_mixed.json")
        self.assertEqual(self.run_check(path), 0)
        self.assertNotIn("holds a single element", self.report(path))

    def test_stereoisomer_of_a_stored_record_is_a_warning(self):
        from rdkit import Chem
        for r in self.records:
            mol = Chem.MolFromSmiles(r["ligand"]["SMILES"])
            centres = Chem.FindMolChiralCenters(mol, includeUnassigned=True, useLegacyImplementation=False)
            if centres:
                mol.GetAtomWithIdx(centres[0][0]).SetChiralTag(Chem.ChiralType.CHI_TETRAHEDRAL_CW)
                smiles, name = Chem.MolToSmiles(mol), r["ligand"]["Name"]
                break
        else:
            self.skipTest("no fixture record has a stereocentre")
        paper = fixture_paper()
        paper["new_ligands"][0]["SMILES"] = smiles
        path = self.write(paper)
        self.assertEqual(self.run_check(path), 0)
        self.assertIn(f"same connectivity as {name}", self.report(path))

    def test_unreviewed_file_is_refused_unless_waived(self):
        paper = fixture_paper()
        del paper["reviewed_by"]
        path = self.write(paper)
        self.assertEqual(self.run_check(path), 0)            # the check does not need a reviewer yet
        db_before = self.db.read_bytes()
        self.assertNotEqual(self.run_apply(path), 0)
        self.assertEqual(self.db.read_bytes(), db_before)
        self.assertIn("reviewed_by", self.report(path))
        self.assertEqual(self.apply.main([str(path), "--database", str(self.db), "--report-dir", str(self.reports),
                                          "--audit-dir", str(self.audits), "--archive-dir", str(self.archive),
                                          "--offline", "--unreviewed"]), 0)
        audit = next(self.audits.glob("apply_*.txt")).read_text(encoding="utf-8")
        self.assertIn("applied --unreviewed", audit)

    def test_all_or_nothing(self):
        good = self.write(fixture_paper(), "10.1234_good.json")
        bad_paper = fixture_paper()
        bad_paper["DOI"] = "10.1234/other"
        bad_paper["new_ligands"] = []
        bad_paper["conditions"] = [dict(bad_paper["conditions"][0], ligand="Nobody")]
        bad = self.write(bad_paper, "10.1234_bad.json")
        db_before = self.db.read_bytes()
        self.assertNotEqual(self.run_apply(good, bad), 0)
        self.assertEqual(self.db.read_bytes(), db_before, "a run with a bad file wrote the database")
        for path in (good, bad):
            self.assertNotIn("status", json.loads(path.read_text(encoding="utf-8")))
        self.assertEqual(list(self.audits.glob("*")), [])
        # without arguments the scripts take every pending file of the directory
        code = self.check.main(["--contributions-dir", str(self.dir), "--database", str(self.db),
                                "--report-dir", str(self.reports), "--offline"])
        self.assertEqual(code, 1)
        self.assertTrue((self.reports / "10.1234_good_check.md").exists())
        self.assertTrue((self.reports / "10.1234_bad_check.md").exists())

    def test_two_files_in_one_run_see_each_other(self):
        first = self.write(fixture_paper(), "10.1234_first.json")
        second_paper = fixture_paper()
        second_paper["DOI"] = "10.1234/second"
        second_paper["new_ligands"] = []                      # FixtureDGA now comes from the first file
        shifted = {el: dict(c, log_D=round(c["log_D"] + 0.3, 3)) for el, c in second_paper["conditions"][1]["cells"].items()}
        second_paper["conditions"] = [dict(second_paper["conditions"][1], id="A0007", cells=shifted,
                                           reference={"Locator": "Figure 4"})]
        second = self.write(second_paper, "10.1234_second.json")
        repeat_paper = fixture_paper()
        repeat_paper["DOI"] = "10.1234/repeat"
        repeat_paper["new_ligands"] = []
        repeat_paper["conditions"] = [dict(repeat_paper["conditions"][1], id="A0009")]
        repeat = self.write(repeat_paper, "10.1234_repeat.json")
        db_before = self.db.read_bytes()
        self.assertNotEqual(self.run_apply(first, second, repeat), 0)
        self.assertEqual(self.db.read_bytes(), db_before)
        self.assertIn("duplicate of Lit_FixtureDGA_10.1234_fixture_A0002", self.report(repeat))
        repeat.unlink()
        self.assertEqual(self.run_apply(first, second), 0)
        after = self.stored()
        self.assertEqual([c["condition_id"] for c in after["FixtureDGA"]["experimental_data"]],
                         ["Lit_FixtureDGA_10.1234_fixture_A0002", "Lit_FixtureDGA_10.1234_second_A0007"])

    # --- the format is one format ----------------------------------------------

    def test_schema_on_disk_is_current(self):
        generator = load_script("generate_json_schema")
        on_disk = (ROOT / "schema" / "contribution.schema.json").read_text(encoding="utf-8")
        self.assertEqual(on_disk, generator.render_contribution(), "run scripts/generate_json_schema.py")
        from jsonschema import Draft202012Validator
        Draft202012Validator.check_schema(json.loads(on_disk))

    def test_v14_additions_fit_the_grammar(self):
        # The literature additions that built v14 declared their new ligands in
        # data/new_ligands_v14.json, so only their conditions are the shared format.
        if not ADDITIONS.exists():
            self.skipTest("the v14 additions are in the developer repository only")
        files = sorted(ADDITIONS.glob("*.json"))
        self.assertTrue(files)
        for path in files:
            paper = self.lib.load_contribution(path)
            self.assertEqual(self.lib.normalize_doi(paper["DOI"]), paper["DOI"], path.name)
            for a in paper["conditions"]:
                self.assertEqual(self.lib.addition_problems(paper, a), [], f"{path.name} {a.get('id')}")


if __name__ == "__main__":
    unittest.main()
