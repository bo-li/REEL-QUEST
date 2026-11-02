"""Accuracy tests for the structure-based Subtype classifier.

The v13 taxonomy was assigned by human review, so the only meaningful check on
an automatic classifier is how well it reproduces those curated labels. This
measures agreement against all 166 records and fails if it regresses.

The classifier is not currently used to build the database — it exists for
ingesting new sources (see scripts/ingest_safe_02102026.py). Its measured
accuracy is what tells a curator how far to trust it.
"""

import json
import unittest
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
DB_PATH = ROOT / "data" / "sx_ligands_merged_v14_flat.jsonl"

# Agreement measured on 2026-07-30: 160/166; 2026-10-07, with the bis-lactam rule
# (BLPhen, BLBPy): 163/167. Raise this floor if the classifier improves; a drop
# means a pattern regressed.
MIN_AGREEMENT = 0.97


class TestSubtypeClassifier(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        try:
            import rdkit  # noqa: F401
        except ImportError:  # pragma: no cover
            raise unittest.SkipTest("RDKit not installed")
        cls.records = [json.loads(line) for line in DB_PATH.open() if line.strip()]

    def test_agreement_with_curated_labels(self):
        from rel_quest.subtype_classifier import classify_subtype

        agree = 0
        mismatches = Counter()
        for record in self.records:
            truth = record["ligand"].get("Subtype")
            predicted, _confidence, _reason = classify_subtype(
                record["ligand"].get("SMILES"))
            if predicted == truth:
                agree += 1
            else:
                mismatches[(truth, predicted)] += 1

        ratio = agree / len(self.records)
        self.assertGreaterEqual(
            ratio, MIN_AGREEMENT,
            f"agreement fell to {ratio:.1%} ({agree}/{len(self.records)}); "
            f"top mismatches: {mismatches.most_common(5)}",
        )
        print(f"✓ Subtype classifier agrees with curated labels on "
              f"{agree}/{len(self.records)} = {ratio:.1%}")

    def test_every_curated_record_is_decidable(self):
        """The classifier must not return None for chemistry already in scope."""
        from rel_quest.subtype_classifier import classify_subtype

        undecided = [
            record["ligand"]["Name"] for record in self.records
            if classify_subtype(record["ligand"].get("SMILES"))[0] is None
        ]
        self.assertEqual(
            undecided, [],
            f"classifier could not decide known-in-scope ligands: {undecided}",
        )

    def test_phosphorus_topology_rule(self):
        """P-bond topology decides for P-containing ligands (the v13 rule)."""
        from rel_quest.subtype_classifier import classify_subtype

        cases = {
            # TBP: three P-OR bonds, no P-C -> Phosphate
            "CCCCOP(=O)(OCCCC)OCCCC": "Phosphate",
            # CMPO-like: P with three P-C bonds -> Phosphine_oxide, despite the
            # amide arm that put it in Amide_Phosphoryl under v12.
            "CCCCCCCCP(=O)(CCCCCCCC)CC(=O)N(CCCC)CCCC": "Phosphine_oxide",
        }
        for smiles, expected in cases.items():
            with self.subTest(smiles=smiles):
                self.assertEqual(classify_subtype(smiles)[0], expected)


if __name__ == "__main__":
    unittest.main()
