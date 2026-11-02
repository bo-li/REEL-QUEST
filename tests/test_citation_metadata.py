"""Citation, licence and archival metadata agree with each other and the code.

A data paper's "Code/Data Availability" statement points at these files, and
Zenodo archives whatever they say at release time. They repeat the version
and the licence split, so these tests stop them drifting — the same failure
`tests/test_packaging.py` guards against for dependencies.
"""

import json
import re
import subprocess
import sys
import unittest
from pathlib import Path

from rel_quest.schema import PROFILE

ROOT = Path(__file__).resolve().parent.parent


def cff_field(name):
    match = re.search(rf"^{name}:\s*(.+)$", (ROOT / "CITATION.cff").read_text(), re.M)
    return match.group(1).strip().strip('"') if match else None


class TestLicenceFiles(unittest.TestCase):

    def test_code_licence_is_mit(self):
        text = (ROOT / "LICENSE").read_text()
        self.assertTrue(text.startswith("MIT License"))
        self.assertIn("LICENSE-DATA", text)

    def test_data_licence_carries_the_cc_by_legal_code(self):
        text = (ROOT / "LICENSE-DATA").read_text()
        self.assertIn("Attribution 4.0 International", text)
        self.assertIn("Section 1 -- Definitions.", text)
        self.assertIn("data/sx_ligands_merged_v14_flat.jsonl", text)
        self.assertIn("paper/references.bib", text)

    def test_safe_paper_is_listed_among_the_references(self):
        # The SAFE database paper is credited in CITATION.cff's references (as a
        # source of the data), not in the licence.
        self.assertIn("10.1080/07366299.2025.2564381", (ROOT / "CITATION.cff").read_text())


class TestBibliography(unittest.TestCase):
    """paper/references.bib: every paper the manuscript and ESI cite, as BibTeX.

    Ways it could go wrong:
    1. a paper the database's values come from is missing, so it cannot be credited;
    2. it is stale against paper/references.json and the drafts (developer repo only);
    3. an entry carries a curation field (data_source) the public data has dropped;
    4. two entries share a key, or a key breaks BibTeX (a DOI's '/' or ':');
    5. the drafts' markup (^sup^, [[TODO: ...]]) or a bare '&' or '%' reaches a field.
    """

    @classmethod
    def setUpClass(cls):
        cls.text = (ROOT / "paper" / "references.bib").read_text(encoding="utf-8")
        cls.keys = re.findall(r"^@\w+\{([^,\s]+),$", cls.text, re.M)
        cls.dois = set(re.findall(r"^\s*doi = \{([^}]+)\},?$", cls.text, re.M))

    def test_every_measuring_paper_is_listed(self):                      # 1
        database = ROOT / "data" / "sx_ligands_merged_v14_flat.jsonl"
        measured = {c["reference"]["DOI"] for line in database.open(encoding="utf-8")
                    for c in json.loads(line)["experimental_data"] if c["reference"].get("DOI")}
        self.assertEqual(sorted(measured - self.dois), [])
        self.assertIn("10.1080/07366299.2025.2564381", self.dois)

    @unittest.skipIf(PROFILE == "public", "built from the drafts in the developer repository")
    def test_is_current(self):                                           # 2
        result = subprocess.run([sys.executable, "paper/references_bib.py", "--check"], cwd=ROOT,
                                capture_output=True, text=True)
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)

    def test_carries_no_curation_fields(self):                           # 3
        self.assertNotRegex(self.text, r"data_source|Aggregator|Provenance_Level")

    def test_keys_are_unique_and_plain(self):                            # 4
        self.assertGreater(len(self.keys), 100)
        self.assertEqual(len(self.keys), len(set(self.keys)))
        self.assertEqual([k for k in self.keys if not re.fullmatch(r"[a-z][a-z0-9]*", k)], [])

    def test_no_draft_markup_or_bare_specials(self):                     # 5
        self.assertNotIn("[[TODO", self.text)
        self.assertNotRegex(self.text, r"\^[^{]")
        entries = "\n".join(line for line in self.text.splitlines() if not line.startswith("%"))
        self.assertNotRegex(entries, r"(?<!\\)[&%]")


class TestVersionAgreement(unittest.TestCase):

    def setUp(self):
        import rel_quest
        self.version = rel_quest.__version__

    def test_citation_cff(self):
        self.assertEqual(cff_field("version"), self.version)
        self.assertEqual(cff_field("cff-version"), "1.2.0")

    def test_zenodo(self):
        zenodo = json.loads((ROOT / ".zenodo.json").read_text())
        self.assertEqual(zenodo["version"], self.version)
        self.assertEqual(zenodo["license"], "MIT")
        self.assertIn("CC BY 4.0", zenodo["notes"])

    def test_codemeta(self):
        codemeta = json.loads((ROOT / "codemeta.json").read_text())
        self.assertEqual(codemeta["version"], self.version)
        self.assertEqual(set(codemeta["license"]),
                         {"https://spdx.org/licenses/MIT",
                          "https://spdx.org/licenses/CC-BY-4.0"})

    def test_changelog_leads_with_this_version(self):
        first = re.search(r"^## \[([^\]]+)\]", (ROOT / "CHANGELOG.md").read_text(), re.M)
        self.assertEqual(first.group(1), self.version)

    def test_codemeta_requirements_match_requirements_txt(self):
        codemeta = json.loads((ROOT / "codemeta.json").read_text())
        requirements = [line.strip() for line in (ROOT / "requirements.txt").read_text().splitlines()
                        if line.strip() and not line.startswith("#")]
        self.assertEqual(codemeta["softwareRequirements"], requirements)


class TestRelatedIdentifiers(unittest.TestCase):

    def test_source_dois_are_well_formed_and_consistent(self):
        from rel_quest.references import DOI_PATTERN, LIU_2022_DOI, SAFE_DATABASE_DOI

        zenodo = json.loads((ROOT / ".zenodo.json").read_text())
        dois = {r["identifier"] for r in zenodo["related_identifiers"] if r["scheme"] == "doi"}
        for doi in dois:
            with self.subTest(doi=doi):
                self.assertRegex(doi, DOI_PATTERN)
        self.assertIn(LIU_2022_DOI, dois)
        self.assertIn(SAFE_DATABASE_DOI, dois)
        cff = (ROOT / "CITATION.cff").read_text()
        for doi in dois:
            self.assertIn(doi, cff)

    def test_cff_parses_when_yaml_is_available(self):
        try:
            import yaml
        except ImportError:  # pragma: no cover - optional
            self.skipTest("PyYAML not installed")
        data = yaml.safe_load((ROOT / "CITATION.cff").read_text())
        self.assertEqual(data["license"], ["MIT", "CC-BY-4.0"])
        self.assertTrue(data["authors"])


if __name__ == "__main__":
    unittest.main()
