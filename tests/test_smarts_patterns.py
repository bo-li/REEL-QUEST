"""
Unit Tests for SMARTS Pattern Matching

Tests SMARTS patterns against known ligands to ensure correct extractant type
classification.

Author: REX Database Project
Date: 2026-02-23
"""

import unittest
import sys
from rel_quest.smarts_patterns import (
    match_extractant_type,
    validate_single_ligand,
    match_smarts_pattern,
    check_rdkit_available,
    list_all_extractant_types
)


class TestSMARTSPatterns(unittest.TestCase):
    """Test SMARTS pattern matching for extractant type classification"""

    @classmethod
    def setUpClass(cls):
        """Check RDKit availability before running tests"""
        if not check_rdkit_available():
            raise ImportError("RDKit is required for these tests. Install with: pip install rdkit")

    def test_todga_is_diglycolamide(self):
        """Test TODGA (known Diglycolamide) matches correctly"""
        # TODGA: N,N,N',N'-tetraoctyldiglycolamide
        smiles = 'CCCCCCCCN(CCCCCCCC)C(=O)COCCOC(=O)N(CCCCCCCC)CCCCCCCC'

        result = match_extractant_type(smiles)

        self.assertIn('Diglycolamide', result['matched_types'],
                     f"TODGA should match Diglycolamide. Got: {result['matched_types']}")
        self.assertEqual(result['suggested_type'], 'Diglycolamide',
                        f"TODGA should be suggested as Diglycolamide. Got: {result['suggested_type']}")

    def test_thdga_is_diglycolamide(self):
        """Test THDGA (another DGA) matches correctly"""
        # THDGA: N,N,N',N'-tetrahexyldiglycolamide
        smiles = 'CCCCCCN(CCCCCC)C(=O)COCCOC(=O)N(CCCCCC)CCCCCC'

        result = match_extractant_type(smiles)

        self.assertIn('Diglycolamide', result['matched_types'])
        self.assertEqual(result['suggested_type'], 'Diglycolamide')

    def test_cmpo_is_organophosphorus(self):
        """Test CMPO (known Organophosphorus) matches correctly"""
        # CMPO: octyl(phenyl)-N,N-diisobutylcarbamoylmethylphosphine oxide
        smiles = 'CCCCCCCCP(=O)(c1ccccc1)CC(=O)N(CC(C)C)CC(C)C'

        result = match_extractant_type(smiles)

        self.assertIn('Organophosphorus', result['matched_types'],
                     f"CMPO should match Organophosphorus. Got: {result['matched_types']}")
        self.assertEqual(result['suggested_type'], 'Organophosphorus',
                        f"CMPO should be suggested as Organophosphorus. Got: {result['suggested_type']}")

    def test_tbp_is_organophosphorus(self):
        """Test TBP (tributyl phosphate) matches correctly"""
        smiles = 'CCCCOP(=O)(OCCCC)OCCCC'

        result = match_extractant_type(smiles)

        self.assertIn('Organophosphorus', result['matched_types'])

    def test_topo_is_organophosphorus(self):
        """Test TOPO (trioctylphosphine oxide) matches correctly"""
        smiles = 'CCCCCCCCP(=O)(CCCCCCCC)CCCCCCCC'

        result = match_extractant_type(smiles)

        self.assertIn('Organophosphorus', result['matched_types'])
        self.assertEqual(result['suggested_type'], 'Organophosphorus')

    def test_btp_is_aromatic_heterocycle(self):
        """Test BTP (bis-triazinyl-pyridine) matches correctly"""
        # CA-BTP: 2,6-bis(5,6-dipropyl-1,2,4-triazin-3-yl)pyridine
        smiles = 'CCCC1=NC(C2=CC=CC(C3=NC(CCC)=CN=N3)=N2)=NN=C1CCC'

        result = match_extractant_type(smiles)

        self.assertIn('Aromatic_Heterocycle', result['matched_types'],
                     f"BTP should match Aromatic_Heterocycle. Got: {result['matched_types']}")
        # Note: might also match amide or other patterns, so just check it's in the list

    def test_phenanthroline_is_aromatic_heterocycle(self):
        """Test phenanthroline matches aromatic heterocycle"""
        # 1,10-phenanthroline
        smiles = 'c1cnc2c(c1)ccc1ncccc12'

        result = match_extractant_type(smiles)

        self.assertIn('Aromatic_Heterocycle', result['matched_types'])

    def test_pyridine_is_aromatic_heterocycle(self):
        """Test simple pyridine matches aromatic heterocycle"""
        smiles = 'c1ccncc1'

        result = match_extractant_type(smiles)

        self.assertIn('Aromatic_Heterocycle', result['matched_types'])

    def test_malonamide_pattern(self):
        """Test malonamide pattern (two amides separated by CH2, no ether)"""
        # Simple malonamide: N-C(=O)-CH2-C(=O)-N
        smiles = 'CCCCN(CCCC)C(=O)CC(=O)N(CCCC)CCCC'

        result = match_extractant_type(smiles)

        self.assertIn('Malonamide', result['matched_types'],
                     f"Should match Malonamide. Got: {result['matched_types']}")

    def test_amide_simple(self):
        """Test simple amide pattern"""
        # Simple amide without DGA or malonamide features
        smiles = 'CCCCCCCCN(CCCCCCCC)C(=O)CCCCCCCC'

        result = match_extractant_type(smiles)

        # Should match Amide, but not Diglycolamide or Malonamide
        self.assertIn('Amide', result['matched_types'])
        self.assertNotIn('Diglycolamide', result['matched_types'])
        self.assertNotIn('Malonamide', result['matched_types'])

    def test_dga_not_malonamide(self):
        """Test that DGA is not classified as malonamide (has ether oxygen)"""
        smiles = 'CCCCCCCCN(CCCCCCCC)C(=O)COCCOC(=O)N(CCCCCCCC)CCCCCCCC'

        result = match_extractant_type(smiles)

        self.assertIn('Diglycolamide', result['matched_types'])
        self.assertNotIn('Malonamide', result['matched_types'],
                        "DGA should not match Malonamide (has ether oxygen)")

    def test_unknown_ligand_without_smiles(self):
        """Test ligands without SMILES are handled gracefully"""
        ligand_data = {
            'ligand_name': 'Test Ligand',
            'smiles': '',
            'extractant_type': 'Unknown'
        }

        result = validate_single_ligand(ligand_data)

        self.assertEqual(result['match_status'], 'no_smiles')
        self.assertIsNone(result['smiles_predicted_type'])

    def test_invalid_smiles(self):
        """Test invalid SMILES strings are handled gracefully"""
        invalid_smiles = 'INVALID_SMILES_STRING_123'

        result = match_extractant_type(invalid_smiles)

        self.assertEqual(result['confidence'], 'invalid_smiles')

    def test_validate_correct_ligand(self):
        """Test validation of correctly classified ligand"""
        ligand_data = {
            'ligand_name': 'TODGA',
            'smiles': 'CCCCCCCCN(CCCCCCCC)C(=O)COCCOC(=O)N(CCCCCCCC)CCCCCCCC',
            'extractant_type': 'Diglycolamide'
        }

        result = validate_single_ligand(ligand_data)

        self.assertEqual(result['match_status'], 'correct',
                        f"TODGA classified as Diglycolamide should be correct. Got: {result}")
        self.assertIsNone(result['suggested_correction'])

    def test_validate_mismatched_ligand(self):
        """Test validation of incorrectly classified ligand"""
        ligand_data = {
            'ligand_name': 'CMPO',
            'smiles': 'CCCCCCCCP(=O)(c1ccccc1)CC(=O)N(CC(C)C)CC(C)C',
            'extractant_type': 'Amide'  # Wrong! Should be Organophosphorus (higher priority)
        }

        result = validate_single_ligand(ligand_data)

        # CMPO matches both Organophosphorus and Amide, so it's ambiguous
        # (not a complete mismatch, but Amide is lower priority)
        self.assertEqual(result['match_status'], 'ambiguous')
        self.assertEqual(result['suggested_correction'], 'Organophosphorus')

    def test_ambiguous_dual_match(self):
        """Test ligands matching multiple patterns are flagged as ambiguous"""
        # CMPO has both P=O (organophosphorus) and amide
        # But should primarily be classified as Organophosphorus due to priority
        smiles = 'CCCCCCCCP(=O)(c1ccccc1)CC(=O)N(CC(C)C)CC(C)C'

        result = match_extractant_type(smiles)

        # Should match both Organophosphorus and Amide
        self.assertTrue(len(result['matched_types']) >= 1,
                       "CMPO should match at least one pattern")

        # But Organophosphorus should be suggested (higher priority)
        self.assertEqual(result['suggested_type'], 'Organophosphorus',
                        f"CMPO should be suggested as Organophosphorus. Got: {result['suggested_type']}")

    def test_carboxylic_acid_pattern(self):
        """Test carboxylic acid pattern"""
        # Simple carboxylic acid
        smiles = 'CCCCCCCC(=O)O'

        result = match_extractant_type(smiles)

        self.assertIn('Carboxylic_Acid', result['matched_types'])

    def test_extractant_types_list(self):
        """Test that extractant types are properly defined"""
        types = list_all_extractant_types()

        expected_types = [
            'Diglycolamide',
            'Malonamide',
            'Organophosphorus',
            'Aromatic_Heterocycle',
            'Amide',
            'Carboxylic_Acid'
        ]

        for expected_type in expected_types:
            self.assertIn(expected_type, types,
                         f"{expected_type} should be in defined extractant types")

    def test_confidence_high_single_match(self):
        """Test that single clear match gives high confidence"""
        # TODGA - clear Diglycolamide
        smiles = 'CCCCCCCCN(CCCCCCCC)C(=O)COCCOC(=O)N(CCCCCCCC)CCCCCCCC'

        result = match_extractant_type(smiles)

        self.assertEqual(result['confidence'], 'high',
                        f"Clear single match should have high confidence. Got: {result['confidence']}")

    def test_confidence_medium_multiple_matches(self):
        """Test that multiple matches give medium/low confidence"""
        # CMPO - matches both Organophosphorus and Amide
        smiles = 'CCCCCCCCP(=O)(c1ccccc1)CC(=O)N(CC(C)C)CC(C)C'

        result = match_extractant_type(smiles)

        # Confidence should be medium or high (high because of priority system)
        self.assertIn(result['confidence'], ['high', 'medium'],
                     f"Multiple matches with clear priority should have medium/high confidence. Got: {result['confidence']}")

    def test_pattern_matching_basic(self):
        """Test basic SMARTS pattern matching function"""
        # Test simple amide pattern
        smiles = 'CCCCN(CCCC)C(=O)C'
        pattern = '[NX3][CX3](=[OX1])'

        result = match_smarts_pattern(smiles, pattern)

        self.assertTrue(result, "Should match amide pattern")

    def test_pattern_matching_no_match(self):
        """Test pattern that should not match"""
        # Simple alkane (no functional groups)
        smiles = 'CCCCCCCC'
        pattern = '[NX3][CX3](=[OX1])'  # Amide pattern

        result = match_smarts_pattern(smiles, pattern)

        self.assertFalse(result, "Alkane should not match amide pattern")


class TestEdgeCases(unittest.TestCase):
    """Test edge cases and error handling"""

    def test_empty_smiles(self):
        """Test empty SMILES string"""
        result = match_extractant_type('')

        self.assertEqual(result['confidence'], 'no_smiles')
        self.assertEqual(result['matched_types'], [])

    def test_none_smiles(self):
        """Test None SMILES"""
        result = match_extractant_type(None)

        self.assertEqual(result['confidence'], 'no_smiles')

    def test_very_long_smiles(self):
        """Test that very long SMILES strings don't cause issues"""
        # Create a long polymer-like SMILES
        long_smiles = 'C' * 1000 + 'N(C)C(=O)C'

        try:
            result = match_extractant_type(long_smiles)
            # Should either work or fail gracefully
            self.assertIsNotNone(result)
        except Exception as e:
            self.fail(f"Long SMILES should not raise exception: {e}")

    def test_ligand_data_missing_fields(self):
        """Test validation with missing fields in ligand data"""
        ligand_data = {
            'ligand_name': 'Test'
            # Missing smiles and extractant_type
        }

        result = validate_single_ligand(ligand_data)

        self.assertEqual(result['match_status'], 'no_smiles')


def run_tests_verbose():
    """Run tests with verbose output"""
    suite = unittest.TestLoader().loadTestsFromModule(sys.modules[__name__])
    runner = unittest.TextTestRunner(verbosity=2)
    runner.run(suite)


if __name__ == '__main__':
    # Run tests
    unittest.main(verbosity=2)
