"""The controlled vocabularies must have exactly one definition.

The seven trend classifications were copied into the validator, the CLI, the
dashboard, the plotting module and the agent tool schema. When trend fitting
moved to a quadratic model and gained Peaked and U_Shaped, five copies were
updated and one was not: `scripts/query_ligands.py --trend` went on offering
five choices, so 18 ligands were queryable through the Python API and the
dashboard but unreachable from the command line
(PROJECT_REVIEW_2026-07-31.md finding 5).

These tests pin the vocabulary to rel_quest.schema and check that no consumer
has quietly grown its own copy again.

They also compile every tracked Python file. pytest only imports what the
tests reach, so a tracked script can be syntactically invalid indefinitely —
which is how scripts/discover_ideal_extractants.py stayed in the repo as an
unterminated docstring (finding 8).
"""

import ast
import re
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent

# Files that legitimately mention trend names for reasons other than defining
# the vocabulary: colour maps keyed by trend, and the classifiers that emit
# one literal 'Insufficient_Data' as an early return.
ALLOWED_LITERAL_FILES = {
    "rel_quest/plotting.py",       # TREND_COLORS colour map
    "app.py",                      # trend colour map
    "rel_quest/ingest_liu.py",
    "rel_quest/ingest_safe.py",
    "rel_quest/ingest_liu_corrected.py",
    "rel_quest/ingest_ideal.py",
    "rel_quest/utils.py",          # merge consensus logic
    "rel_quest/schema.py",         # the definition itself
    "scripts/compute_uq.py",
}


class TestTrendVocabulary(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        from rel_quest.schema import (TREND_CLASSIFICATIONS, TREND_DIRECTIONS,
                                      TREND_DISPLAY_ORDER)
        cls.classifications = TREND_CLASSIFICATIONS
        cls.display_order = TREND_DISPLAY_ORDER
        cls.directions = TREND_DIRECTIONS

    def test_seven_classifications(self):
        self.assertEqual(len(self.classifications), 7)
        self.assertIn("Peaked", self.classifications)
        self.assertIn("U_Shaped", self.classifications)

    def test_display_order_covers_every_classification_plus_unknown(self):
        self.assertEqual(set(self.display_order),
                         set(self.classifications) | {"Unknown"})

    def test_every_classification_has_a_direction(self):
        for trend in self.classifications:
            with self.subTest(trend=trend):
                self.assertIn(trend, self.directions)
                self.assertTrue(self.directions[trend].strip())

    def test_matches_what_the_database_actually_contains(self):
        import json

        db_path = ROOT / "data" / "sx_ligands_merged_v14_flat.jsonl"
        stored = {json.loads(line)["selectivity"]["trend_classification"]
                  for line in db_path.open() if line.strip()}
        self.assertTrue(
            stored.issubset(set(self.classifications)),
            f"database holds trends outside the vocabulary: "
            f"{stored - set(self.classifications)}")

    def test_validator_uses_the_shared_vocabulary(self):
        from rel_quest.validators import validate_ligand_record

        for trend in self.classifications:
            record = {
                "schema_version": "v0.2",
                "ligand": {"Name": "T", "Mechanism": "Solvation",
                           "Extractant_Type": "Diglycolamide"},
                "experimental_data": [{
                    "condition_id": "c1", "reference": {"DOI": "10.1/x"},
                    "distribution_ratios": {"La": {"value": 1.0, "log_value": 0.0}},
                }],
                "selectivity": {"trend_classification": trend,
                                "element_coverage": ["La"]},
                "metadata": {"data_source": "TEST"},
                "performance_summary": {},
                "elements_measured": ["La"],
            }
            with self.subTest(trend=trend):
                is_valid, errors = validate_ligand_record(record)
                self.assertTrue(is_valid, errors)

    def test_agent_tool_schema_offers_every_classification(self):
        from rel_quest import agent_tools

        blob = repr(agent_tools.TOOL_SCHEMAS if hasattr(agent_tools, "TOOL_SCHEMAS")
                    else agent_tools.get_tool_schemas())
        for trend in self.classifications:
            with self.subTest(trend=trend):
                self.assertIn(trend, blob)


class TestCliAcceptsEveryTrend(unittest.TestCase):
    """The regression that started this: --trend rejected two valid values."""

    def _run(self, *args):
        return subprocess.run(
            [sys.executable, str(ROOT / "scripts" / "query_ligands.py"), *args],
            cwd=ROOT, capture_output=True, text=True, timeout=180)

    def test_every_trend_is_an_accepted_choice(self):
        from rel_quest.schema import TREND_CLASSIFICATIONS

        for trend in TREND_CLASSIFICATIONS:
            with self.subTest(trend=trend):
                result = self._run("--trend", trend, "--output", "csv")
                self.assertEqual(
                    result.returncode, 0,
                    f"--trend {trend} rejected: {result.stderr[-300:]}")
                self.assertNotIn("invalid choice", result.stderr)

    def test_mechanism_choices_come_from_the_vocabulary_file(self):
        from rel_quest.schema import canonical_mechanisms

        result = self._run("--help")
        self.assertEqual(result.returncode, 0)
        for mechanism in canonical_mechanisms():
            with self.subTest(mechanism=mechanism):
                self.assertIn(mechanism, result.stdout)

    def test_an_invalid_trend_is_still_rejected(self):
        result = self._run("--trend", "Sideways")
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("invalid choice", result.stderr)


class TestNoDuplicateVocabularyDefinitions(unittest.TestCase):
    """No consumer may re-declare the vocabulary as its own list literal."""

    def test_no_module_redeclares_the_trend_list(self):
        """Walk the AST for list/tuple/set literals of trend-name strings.

        Deliberately AST rather than regex: prose in the agent's system
        prompt names several trends in a sentence, which a text search
        cannot tell from a hardcoded `['Ascending', ...]`.
        """
        from rel_quest.schema import TREND_CLASSIFICATIONS

        names = set(TREND_CLASSIFICATIONS)
        offenders = []
        for path in sorted(ROOT.glob("**/*.py")):
            rel = path.relative_to(ROOT).as_posix()
            if any(part in rel for part in
                   ("data/archive", "docs/history", "venv/", ".claude", "tests/")):
                continue
            if rel in ALLOWED_LITERAL_FILES:
                continue
            try:
                tree = ast.parse(path.read_text(encoding="utf-8", errors="ignore"))
            except SyntaxError:  # pragma: no cover - caught by the compile test
                continue
            for node in ast.walk(tree):
                if not isinstance(node, (ast.List, ast.Tuple, ast.Set)):
                    continue
                literals = {e.value for e in node.elts
                            if isinstance(e, ast.Constant) and isinstance(e.value, str)}
                if len(literals & names) >= 3:
                    offenders.append(
                        f"{rel}:{node.lineno}: {sorted(literals & names)}")
        self.assertEqual(
            offenders, [],
            "these files redeclare the trend vocabulary instead of importing "
            f"it from rel_quest.schema: {offenders}")


class TestEveryTrackedPythonFileCompiles(unittest.TestCase):
    """Finding 8: a tracked script sat broken because nothing imported it."""

    def test_all_tracked_python_parses(self):
        result = subprocess.run(["git", "ls-files", "*.py"],
                                cwd=ROOT, capture_output=True, text=True)
        if result.returncode != 0:  # pragma: no cover
            self.skipTest("not a git checkout")

        broken = []
        for name in result.stdout.split():
            path = ROOT / name
            if not path.exists():
                continue
            try:
                ast.parse(path.read_text(encoding="utf-8", errors="ignore"),
                          filename=name)
            except SyntaxError as exc:
                broken.append(f"{name}: line {exc.lineno}: {exc.msg}")
        self.assertEqual(broken, [],
                         f"{len(broken)} tracked file(s) do not parse: {broken}")

    def test_ci_runs_a_compile_check(self):
        """The test above only helps if CI runs something equivalent."""
        workflow = ROOT / ".github" / "workflows" / "ci.yml"
        if not workflow.exists():  # pragma: no cover
            self.skipTest("no CI workflow")
        self.assertIn("compileall", workflow.read_text(),
                      "CI has no compile step; a broken tracked file could "
                      "merge unnoticed")


if __name__ == "__main__":
    unittest.main()
