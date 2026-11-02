"""The maintained docs must describe the code that exists.

PROJECT_REVIEW_2026-07-31.md finding 7: the README pointed users at
docs/API.md, docs/CLI.md and docs/QUERY_GUIDE.md, all of which still taught
`from rex_db import LigandDatabase` — a package removed in the 2026-07-27
cleanup. docs/API.md named `data/sx_ligands.jsonl` as the default database and
docs/QUICKSTART.md claimed it held 99 ligands.

Prose cannot be fully tested, so these check the parts that are mechanically
checkable and were the ones actually wrong: import paths, module paths,
database filenames, and the headline counts.

Historical documents are exempt. `docs/history/` and the dated project
reviews record what *was* true; rewriting them would destroy the record.
"""

import re
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
DOCS = ROOT / "docs"

# These describe the past on purpose.
HISTORICAL = {"PROJECT_REVIEW_2026-07-27.md", "PROJECT_REVIEW_2026-07-31.md",
              "PROJECT_REVIEW_2026-08-21.md", "PROJECT_REVIEW_2026-08-25.md"}

# Linked from README as the way to use the project.
ACTIVE_REFERENCE = ["API.md", "CLI.md", "QUERY_GUIDE.md", "QUICKSTART.md",
                    "DATABASE_STRUCTURE.md", "REST_API.md"]


def maintained_docs():
    return [p for p in sorted(DOCS.glob("*.md")) if p.name not in HISTORICAL]


class TestNoRetiredPackageReferences(unittest.TestCase):

    def test_no_doc_teaches_the_removed_rex_db_package(self):
        offenders = []
        for path in maintained_docs():
            for number, line in enumerate(path.read_text(encoding="utf-8").splitlines(), 1):
                if re.search(r"\brex_db\b", line):
                    offenders.append(f"{path.name}:{number}")
        self.assertEqual(
            offenders, [],
            f"docs still reference the removed rex_db package: {offenders}")

    def test_readme_does_not_reference_rex_db(self):
        self.assertNotRegex((ROOT / "README.md").read_text(encoding="utf-8"),
                            r"\brex_db\b")

    def test_no_doc_points_at_a_superseded_database_generation(self):
        """v2-v13 files live in data/archive/; docs must not name them as current."""
        offenders = []
        for path in maintained_docs() + [ROOT / "README.md"]:
            text = path.read_text(encoding="utf-8")
            for match in re.finditer(
                    r"sx_ligands_merged_v(\d+)[a-z_]*\.jsonl", text):
                if match.group(1) != "14" and "archive" not in text[
                        max(0, match.start() - 80):match.start()]:
                    offenders.append(f"{path.name}: {match.group(0)}")
        self.assertEqual(offenders, [],
                         f"docs name a superseded database: {offenders}")


class TestDocumentedApiExists(unittest.TestCase):
    """Every rel_quest module path a doc names must be a real file."""

    def test_referenced_modules_exist(self):
        missing = []
        for path in maintained_docs():
            text = path.read_text(encoding="utf-8")
            for module in set(re.findall(r"`?(rel_quest/[a-z_]+\.py)`?", text)):
                if not (ROOT / module).exists():
                    missing.append(f"{path.name}: {module}")
        self.assertEqual(missing, [], f"docs reference missing modules: {missing}")

    def test_referenced_scripts_exist(self):
        missing = []
        for path in maintained_docs():
            text = path.read_text(encoding="utf-8")
            # One optional subdirectory: scripts/validation/, scripts/migrations/.
            # Paths under another directory (skills/reel-quest/scripts/) are not root scripts.
            for script in set(re.findall(r"`?(?<![\w/-])(scripts/(?:[a-z_0-9]+/)?[a-z_0-9]+\.py)`?", text)):
                if not (ROOT / script).exists():
                    missing.append(f"{path.name}: {script}")
        self.assertEqual(missing, [], f"docs reference missing scripts: {missing}")

    def test_active_reference_docs_import_rel_quest(self):
        """A usage doc with a Python example must show the current import."""
        for name in ACTIVE_REFERENCE:
            path = DOCS / name
            if not path.exists():  # pragma: no cover
                continue
            text = path.read_text(encoding="utf-8")
            if "LigandDatabase" not in text:
                continue
            with self.subTest(doc=name):
                self.assertIn("rel_quest", text,
                              f"{name} uses LigandDatabase without naming rel_quest")


class TestDocumentedCountsMatchTheDatabase(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        from rel_quest import LigandDatabase

        cls.total = len(LigandDatabase().get_all_ligands())

    def test_no_doc_states_a_stale_ligand_total(self):
        """Any total that is not the live one is stale — derived, not
        enumerated. This used to check a hardcoded list of past totals
        (42/99/176/177), so when the total moved from 166 to 170 the '166'
        left behind in CLI.md sailed through (review 2026-08-25 finding 2).
        A doc quoting an old total as history belongs in HISTORICAL."""
        pattern = re.compile(
            r"(?:database (?:contains|holds)|Total ligands:)\s*(\d+)", re.I)
        offenders = []
        for path in maintained_docs() + [p for p in (ROOT / "README.md", ROOT / "CLAUDE.md") if p.exists()]:
            text = path.read_text(encoding="utf-8")
            for match in pattern.finditer(text):
                if int(match.group(1)) != self.total:
                    offenders.append(f"{path.name}: {match.group(0)!r}")
        self.assertEqual(
            offenders, [],
            f"docs state a ligand total that is not {self.total}: {offenders}")


if __name__ == "__main__":
    unittest.main()
