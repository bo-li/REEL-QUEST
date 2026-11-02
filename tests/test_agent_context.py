"""The AI assistant must describe the database that is actually loaded.

PROJECT_REVIEW_2026-07-31.md finding 4. Two independent staleness bugs, both
of which fed the model confident falsehoods it had no way to check:

  1. `_build_context()` matched taxonomy questions against a hardcoded map of
     v12 two-layer family labels ("Diglycolamide", "Organophosphorus",
     "Amide_Phosphoryl", ...) while the database had moved to flat v13
     subtypes. Ten of the thirteen entries matched nothing, so "show me
     phosphine oxide ligands" injected "ORGANOPHOSPHORUS LIGANDS (0 total)"
     against 19 real records.
  2. The system prompt hardcoded per-subtype and per-source counts. They had
     drifted — Phosphonate (4) and Phosphoric_acid (2) against a real 5 and
     1, and source counts summing to 217 for a 166-ligand database.

Both are now derived from the live database and vocabulary. These tests fail
if anything reverts to a transcribed number or a hardcoded label list.

No LLM is contacted: `_build_context` and `build_system_prompt` are pure
functions of the database.
"""

import re
import unittest
from pathlib import Path
from rel_quest.schema import PROFILE

ROOT = Path(__file__).resolve().parent.parent


class TestSubtypeAliases(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        from rel_quest import LigandDatabase
        from rel_quest.schema import subtype_aliases

        cls.db = LigandDatabase()
        cls.aliases = subtype_aliases()
        cls.live_subtypes = {l.subtype for l in cls.db.get_all_ligands() if l.subtype}

    def test_every_live_subtype_is_reachable(self):
        """A subtype nobody can ask about is a subtype the assistant hides."""
        unreachable = sorted(self.live_subtypes - set(self.aliases))
        self.assertEqual(unreachable, [],
                         f"subtypes with no alias entry: {unreachable}")

    def test_aliases_resolve_to_real_records(self):
        """Every alias must lead to a subtype the database actually holds.

        This is the check the old hardcoded map would have failed on ten of
        its thirteen entries.
        """
        from rel_quest.agent_tools import search_by_extractant_type

        empty = []
        for subtype in self.aliases:
            if subtype not in self.live_subtypes:
                continue
            found = search_by_extractant_type(self.db, subtype)["total_found"]
            if found == 0:
                empty.append(subtype)
        self.assertEqual(empty, [],
                         f"aliases point at subtypes matching no record: {empty}")

    def test_no_alias_is_ambiguous(self):
        """One phrase must not map to two subtypes."""
        seen = {}
        collisions = []
        for subtype, aliases in self.aliases.items():
            for alias in aliases:
                if alias in seen:
                    collisions.append(f"{alias!r}: {seen[alias]} vs {subtype}")
                seen[alias] = subtype
        self.assertEqual(collisions, [], f"ambiguous aliases: {collisions}")


class TestContextInjection(unittest.TestCase):
    """End-to-end on the real question phrasings the review used."""

    @classmethod
    def setUpClass(cls):
        from rel_quest import LigandDatabase

        cls.db = LigandDatabase()

    def _subtype_block(self, message):
        from rel_quest.agent import _build_context

        match = re.search(r"\n([A-Z0-9_\-]+) LIGANDS \((\d+) total\)",
                          _build_context(self.db, message))
        return (match.group(1), int(match.group(2))) if match else (None, None)

    def test_the_reviews_failing_probe_now_resolves(self):
        """'Show phosphine oxide ligands' returned 0 of 19 before this fix."""
        label, total = self._subtype_block("Show me phosphine oxide ligands")
        self.assertEqual(label, "PHOSPHINE_OXIDE")
        self.assertEqual(total, len(self.db.filter_by_subtype("Phosphine_oxide")))
        self.assertGreater(total, 0)

    def test_representative_questions_reach_the_right_subtype(self):
        expected = {
            "List DGA ligands": "DGA",
            "Which diglycolamide extractants are best?": "DGA",
            "What about CMPO?": "PHOSPHINE_OXIDE",
            "Tell me about BTBP compounds": "MULTIDENT_N_DONOR",
            "crown ether extractants": "CROWN_ETHER",
            "show malonamide ligands": "MALONAMIDE",
            "list TBP data": "PHOSPHATE",
            "phenanthroline diamide ligands": "DAPHEN",
            "show me EDTA chelators": "EDTA",
        }
        for message, label in expected.items():
            with self.subTest(message=message):
                self.assertEqual(self._subtype_block(message)[0], label)

    def test_injected_counts_match_the_database(self):
        from rel_quest.agent import _build_context

        context = _build_context(self.db, "tell me about DGA ligands")
        match = re.search(r"DGA LIGANDS \((\d+) total\)", context)
        self.assertIsNotNone(match)
        self.assertEqual(int(match.group(1)),
                         len(self.db.filter_by_subtype("DGA")))

    def test_no_subtype_block_when_nothing_matches(self):
        """Silence beats a fabricated zero."""
        self.assertEqual(self._subtype_block("what is the weather today")[0], None)


class TestSystemPrompt(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        from rel_quest import LigandDatabase
        from rel_quest.agent import build_system_prompt

        cls.db = LigandDatabase()
        cls.ligands = cls.db.get_all_ligands()
        cls.prompt = build_system_prompt(cls.db)

    def test_ligand_count_is_live(self):
        self.assertIn(f"{len(self.ligands)} unique extractant ligands", self.prompt)

    def test_every_subtype_count_is_correct(self):
        for subtype in {l.subtype for l in self.ligands if l.subtype}:
            with self.subTest(subtype=subtype):
                expected = len(self.db.filter_by_subtype(subtype))
                self.assertIn(f"{subtype} ({expected})", self.prompt)

    @unittest.skipIf(PROFILE == "public", "source labels are in the developer repository only")
    def test_source_counts_sum_to_the_ligand_total(self):
        """The old prompt's counts summed to 217 for 166 ligands."""
        sources = re.search(r"- Data sources: (.+)", self.prompt).group(1)
        counts = [int(n) for n in re.findall(r"\((\d+)\)", sources)]
        self.assertEqual(sum(counts), len(self.ligands))

    def test_lists_all_seven_trends(self):
        from rel_quest.schema import TREND_CLASSIFICATIONS

        for trend in TREND_CLASSIFICATIONS:
            with self.subTest(trend=trend):
                self.assertIn(trend, self.prompt)

    def test_mentions_the_non_lanthanide_elements_present(self):
        """Y and Sc are in the data; a prompt saying 'La-Lu only' misleads."""
        for element in ("Sc", "Y"):
            with self.subTest(element=element):
                self.assertTrue(
                    any(element in l.tested_elements for l in self.ligands),
                    "fixture assumption failed")
        self.assertRegex(self.prompt, r"plus .*(Sc|Y)")

    def test_template_carries_no_hardcoded_counts(self):
        """Numbers in the template are what drifted; they must be placeholders."""
        from rel_quest.agent import SYSTEM_PROMPT_TEMPLATE

        overview = SYSTEM_PROMPT_TEMPLATE.split("## Key Domain Concepts")[0]
        stray = re.findall(r"\(\d+\)", overview)
        self.assertEqual(stray, [],
                         f"template hardcodes counts instead of interpolating: {stray}")


if __name__ == "__main__":
    unittest.main()
