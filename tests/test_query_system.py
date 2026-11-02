#!/usr/bin/env python3
"""
Test Query System
Comprehensive tests for REX DB query functionality
"""

import os
import unittest

from rel_quest import LigandDatabase, LigandQuery
from rel_quest.schema import PROFILE


class TestLigandDatabase(unittest.TestCase):
    """Test LigandDatabase class functionality"""

    @classmethod
    def setUpClass(cls):
        """Load database once for all tests"""
        db_path = os.path.join(
            os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
            'data', 'sx_ligands_merged_v14_flat.jsonl'
        )
        cls.db = LigandDatabase(db_path)

    def test_database_loads(self):
        """Test database loads successfully"""
        self.assertIsNotNone(self.db)
        stats = self.db.get_statistics()
        self.assertGreater(stats['total_ligands'], 0)
        print(f"✓ Database loaded: {stats['total_ligands']} ligands")

    def test_query_by_name(self):
        """Test query by ligand name"""
        # Test exact match
        results = self.db.query_by_name('TODGA')
        self.assertGreater(len(results), 0)
        self.assertEqual(results[0].name, 'TODGA')
        print(f"✓ Query by name: Found TODGA")

        # Test fuzzy match
        results = self.db.query_by_name('DGA', fuzzy=True)
        self.assertGreater(len(results), 0)
        print(f"✓ Fuzzy search: Found {len(results)} ligands containing 'DGA'")

    def test_query_by_smiles(self):
        """Test query by SMILES"""
        # Get a SMILES from the database first
        all_ligands = self.db.get_all_ligands()
        test_ligand = all_ligands[0]

        if test_ligand.smiles:
            result = self.db.query_by_smiles(test_ligand.smiles)
            self.assertIsNotNone(result)
            self.assertEqual(result.name, test_ligand.name)
            print(f"✓ Query by SMILES: Found {result.name}")

    def test_filter_by_mechanism(self):
        """Test filter by extraction mechanism"""
        solvation_ligands = self.db.filter_by_mechanism('Solvation')
        self.assertGreater(len(solvation_ligands), 0)

        cation_exchange_ligands = self.db.filter_by_mechanism('Cation_Exchange')
        self.assertGreater(len(cation_exchange_ligands), 0)

        print(f"✓ Filter by mechanism: Solvation={len(solvation_ligands)}, Cation_Exchange={len(cation_exchange_ligands)}")

    def test_filter_by_trend(self):
        """Test filter by selectivity trend"""
        ascending_ligands = self.db.filter_by_trend('Ascending')
        flat_ligands = self.db.filter_by_trend('Flat')

        self.assertGreater(len(ascending_ligands), 0)
        self.assertGreater(len(flat_ligands), 0)

        print(f"✓ Filter by trend: Ascending={len(ascending_ligands)}, Flat={len(flat_ligands)}")

    def test_filter_by_elements(self):
        """Test filter by element coverage"""
        # Find ligands tested on La, Nd, and Lu
        results = self.db.filter_by_elements(['La', 'Nd', 'Lu'])
        self.assertGreater(len(results), 0)

        # Verify all results measured the required elements. The filter reads
        # elements_measured (bounds count as measured); element_coverage lists
        # only the points the trend fit used, so a bound-only Lu is not in it.
        for ligand in results:
            self.assertTrue({'La', 'Nd', 'Lu'} <= set(ligand.tested_elements), ligand.name)

        print(f"✓ Filter by elements: Found {len(results)} ligands with La, Nd, Lu data")

    def test_search(self):
        """Test multi-field search"""
        results = self.db.search('phosph')
        self.assertGreater(len(results), 0)
        print(f"✓ Search 'phosph': Found {len(results)} ligands")


class TestLigandQuery(unittest.TestCase):
    """Test LigandQuery class functionality"""

    @classmethod
    def setUpClass(cls):
        """Load database and get test ligand"""
        db_path = os.path.join(
            os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
            'data', 'sx_ligands_merged_v14_flat.jsonl'
        )
        db = LigandDatabase(db_path)

        # Get TBP as test ligand (has all 15 lanthanides from IDEAL)
        results = db.query_by_name('TBP')
        if results:
            cls.ligand = results[0]
        else:
            # Fallback to any ligand
            cls.ligand = db.get_all_ligands()[0]

    def test_get_d_values(self):
        """Test getting D-values"""
        d_values = self.ligand.get_d_values()
        self.assertGreater(len(d_values), 0)
        self.assertIsInstance(d_values, dict)

        # Verify D-values are positive
        for elem, d_val in d_values.items():
            self.assertGreater(d_val, 0)

        print(f"✓ Get D-values: {len(d_values)} elements for {self.ligand.name}")

    def test_get_log_d_values(self):
        """Test getting log(D) values"""
        log_d_values = self.ligand.get_log_d_values()
        self.assertGreater(len(log_d_values), 0)
        print(f"✓ Get log(D) values: {len(log_d_values)} elements")

    def test_get_trend_summary(self):
        """Test trend classification summary"""
        trend = self.ligand.get_trend_summary()
        self.assertIn('classification', trend)
        self.assertIn('confidence', trend)
        self.assertIn('direction', trend)
        self.assertIn('elements_tested', trend)

        print(f"✓ Trend summary: {self.ligand.name} is {trend['classification']} ({trend['confidence']} confidence)")

    def test_get_statistics(self):
        """Statistics come back in the one canonical shape (2026-08-25)."""
        stats = self.ligand.get_statistics()
        self.assertIn('min_log_d', stats)
        self.assertIn('max_log_d', stats)
        self.assertIn('mean_log_d', stats)
        self.assertIn('n_elements', stats)

        self.assertLessEqual(stats['min_log_d'], stats['mean_log_d'])
        self.assertLessEqual(stats['mean_log_d'], stats['max_log_d'])

        print(f"✓ Statistics: mean_log_d={stats['mean_log_d']:.3f}, std_log_d={stats['std_log_d']:.3f}")

    def test_get_separation_factor(self):
        """Test separation factor calculation"""
        sf = self.ligand.get_separation_factor('La', 'Lu')
        if sf is not None:
            self.assertGreaterEqual(sf, 1.0)
            print(f"✓ Separation factor: SF(La/Lu) = {sf:.2f} for {self.ligand.name}")

    def test_get_condition_list(self):
        """Test getting available conditions"""
        conditions = self.ligand.get_condition_list()
        self.assertGreater(len(conditions), 0)
        print(f"✓ Conditions: {len(conditions)} available for {self.ligand.name}")


class TestDatabaseIntegrity(unittest.TestCase):
    """Test database integrity and data quality"""

    @classmethod
    def setUpClass(cls):
        """Load database"""
        db_path = os.path.join(
            os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
            'data', 'sx_ligands_merged_v14_flat.jsonl'
        )
        cls.db = LigandDatabase(db_path)

    def test_all_ligands_have_required_fields(self):
        """Test all ligands have required fields"""
        all_ligands = self.db.get_all_ligands()

        for ligand in all_ligands:
            self.assertIsNotNone(ligand.name)
            self.assertIsNotNone(ligand.mechanism)
            self.assertIsNotNone(ligand.selectivity_trend)
            self.assertGreater(len(ligand.experimental_conditions), 0)

        print(f"✓ All {len(all_ligands)} ligands have required fields")

    def test_no_duplicate_structures(self):
        """Test no two records describe the same molecule (the dedup invariant)"""
        try:
            from rdkit import Chem, RDLogger
        except ImportError:
            self.skipTest("RDKit not installed")
        RDLogger.DisableLog('rdApp.*')

        canonical = {}
        for lig in self.db.get_all_ligands():
            mol = Chem.MolFromSmiles(lig.smiles) if lig.smiles else None
            if mol is None:
                continue
            smi = Chem.MolToSmiles(mol)
            self.assertNotIn(
                smi, canonical,
                f"Duplicate structure: '{lig.name}' and '{canonical.get(smi)}'",
            )
            canonical[smi] = lig.name

        print(f"✓ No duplicate structures: {len(canonical)} distinct molecules")

    def test_no_duplicate_names(self):
        """Test every ligand Name maps to exactly one record.

        Three names (DMDODGA, DBD1MHDGA, OTPTA) each covered two chemically
        distinct molecules until 2026-07-28; see scripts/fix_name_collisions.py
        and outputs/name_collision_fix_audit.txt. This asserts the invariant
        stays true, since query_by_name() indexes on Name alone.
        """
        from collections import Counter

        names = [lig.name for lig in self.db.get_all_ligands()]
        collisions = sorted(n for n, c in Counter(names).items() if c > 1)

        self.assertEqual(collisions, [], f"Ligand name collision(s): {collisions}")
        print(f"✓ No duplicate names: {len(set(names))} uniquely named ligands")

    def test_alternative_names_are_synonyms(self):
        """Alternative_Names must hold synonyms only, not provenance or classes.

        It previously mixed in source-family labels, variant tags, IUPAC
        duplicates and self-references, all of which app.py renders to users as
        "Also known as" — so DMDOHEMA (a malonamide) advertised the family label
        "Di-n-hexylphosphoric acid (HDHP)" as one of its names. Cleaned up by
        scripts/clean_alternative_names.py; this keeps it clean.
        """
        import re

        variant_tag = re.compile(r"^.+_v\d+$")
        offenders = []

        for record in self.db._records:
            ligand = record['ligand']
            name = ligand.get('Name')
            for entry in (ligand.get('Alternative_Names') or []):
                if entry == name:
                    offenders.append((name, entry, 'self-reference'))
                elif entry == ligand.get('IUPAC_Name'):
                    offenders.append((name, entry, 'IUPAC duplicate'))
                elif ligand.get('Family_Name') and entry == ligand['Family_Name']:
                    offenders.append((name, entry, 'source-family label'))
                elif variant_tag.match(entry):
                    offenders.append((name, entry, 'variant tag'))
                elif entry in (ligand.get('Subtype'), ligand.get('Structural_Class')):
                    offenders.append((name, entry, 'class name'))

        self.assertEqual(offenders, [], f"Non-synonyms in Alternative_Names: {offenders}")

        total = sum(len(r['ligand'].get('Alternative_Names') or []) for r in self.db._records)
        print(f"✓ Alternative_Names holds {total} genuine synonyms only")

    def test_renamed_ligands_resolve(self):
        """Each renamed record is reachable by its new name and keeps the old one."""
        renames = {
            'DMDecODGA': 'DMDODGA',
            'DOTPTA': 'OTPTA',
        }
        # DBD1MHDGA (1-methylhexyl) was Liu's drawing of the paper's 1-methylheptyl compound;
        # v14 verification dropped its two rows as duplicates of SAFE's, kept on DBD1MHpDGA,
        # which still answers to the old name.
        self.assertEqual(self.db.query_by_name('DBD1MHDGA'), [])
        self.assertIn('DBD1MHDGA', self.db.query_by_name('DBD1MHpDGA')[0].alternative_names)
        for new_name, old_name in renames.items():
            hits = self.db.query_by_name(new_name)
            self.assertEqual(len(hits), 1, f"{new_name} should resolve to exactly 1 record")
            self.assertIn(
                old_name, hits[0].alternative_names,
                f"{new_name} should retain '{old_name}' in alternative_names",
            )
            # The old name must still resolve — to the record that kept it.
            kept = self.db.query_by_name(old_name)
            self.assertEqual(len(kept), 1, f"{old_name} should resolve to exactly 1 record")
            self.assertNotEqual(kept[0].smiles, hits[0].smiles)

        print(f"✓ All {len(renames)} renamed ligands resolve unambiguously")

    def test_schema_version_is_current(self):
        """Every record declares the current record-schema version.

        The stamp sat at "v0.1" long after the schema had moved on, because it
        was hardcoded in six writers. rel_quest/schema.py now owns the value.
        """
        from rel_quest.schema import SCHEMA_VERSION

        versions = {r.get('schema_version') for r in self.db._records}
        self.assertEqual(
            versions, {SCHEMA_VERSION},
            f"expected every record at {SCHEMA_VERSION}, found {sorted(versions)}",
        )
        print(f"✓ All records stamped schema_version={SCHEMA_VERSION}")

    def test_schema_version_is_independent_of_taxonomy_version(self):
        """The two version fields track different things and must not be equal.

        schema_version describes record structure; metadata.taxonomy_version
        describes the classification generation. Conflating them is what made
        the old stamp meaningless.
        """
        from rel_quest.schema import SCHEMA_VERSION

        taxonomy = {r.get('metadata', {}).get('taxonomy_version')
                    for r in self.db._records}
        self.assertEqual(taxonomy, {'v13_flat'})
        self.assertNotIn(
            SCHEMA_VERSION, taxonomy,
            "schema_version and taxonomy_version should not share a value",
        )
        print(f"✓ schema_version={SCHEMA_VERSION} distinct from taxonomy_version=v13_flat")

    def test_valid_mechanisms(self):
        """Test all mechanisms are valid"""
        valid_mechanisms = [
            'Solvation', 'Cation_Exchange', 'Anion_Exchange',
            'Ion_Pair', 'Supramolecular', 'Hybrid'
        ]

        stats = self.db.get_statistics()
        for mechanism in stats['mechanisms'].keys():
            self.assertIn(mechanism, valid_mechanisms)

        print(f"✓ All mechanisms valid: {list(stats['mechanisms'].keys())}")

    def test_valid_trends(self):
        """Test all trends are valid"""
        valid_trends = ['Ascending', 'Reverse', 'Flat', 'Mixed', 'Insufficient_Data',
                        'Peaked', 'U_Shaped']

        stats = self.db.get_statistics()
        for trend in stats['trends'].keys():
            self.assertIn(trend, valid_trends)

        print(f"✓ All trends valid: {list(stats['trends'].keys())}")

    @unittest.skipIf(PROFILE == "public", "source labels are in the developer repository only")
    def test_data_sources(self):
        """Test data sources are tracked"""
        all_ligands = self.db.get_all_ligands()

        liu_count = sum(1 for lig in all_ligands if lig.metadata.get('data_source') == 'Liu_2022_JACSAu')
        ideal_count = sum(1 for lig in all_ligands if lig.metadata.get('data_source') == 'OECD_NEA_IDEAL')

        print(f"✓ Data sources: Liu_2022={liu_count}, IDEAL={ideal_count}")
        self.assertGreater(liu_count, 0)
        self.assertGreater(ideal_count, 0)


def main():
    """Run all tests"""
    print("=" * 70)
    print("REX DB Query System Tests")
    print("=" * 70)
    print()

    # Run tests
    unittest.main(verbosity=2)


if __name__ == '__main__':
    main()
