"""
SMARTS Pattern Library for Extractant Type Classification

This module provides SMARTS patterns and matching functions for validating
extractant_type classifications in the ligand database based on molecular
structure analysis.

SMARTS (SMiles ARbitrary Target Specification) patterns are used to identify
specific molecular substructures within SMILES strings.

Author: REX Database Project
Date: 2026-02-23
"""

import logging
from typing import Dict, List, Any, Optional, Tuple
from dataclasses import dataclass

try:
    from rdkit import Chem
    from rdkit import RDLogger
    RDKIT_AVAILABLE = True
    # Suppress RDKit warnings
    RDLogger.DisableLog('rdApp.*')
except ImportError:
    RDKIT_AVAILABLE = False
    logging.warning("RDKit not available. SMARTS pattern matching will not work.")


@dataclass
class PatternMatch:
    """Result of a SMARTS pattern match"""
    extractant_type: str
    matched_patterns: List[str]
    confidence: str  # 'high', 'medium', 'low'
    num_matches: int


# SMARTS Pattern Library
# ======================

EXTRACTANT_TYPE_PATTERNS = {
    'Diglycolamide': {
        'description': 'Diglycolamide (DGA) ligands with amide-ether-amide core',
        'patterns': [
            # Core DGA pattern: amide-ether-amide
            # N-C(=O)-C-O-C-C(=O)-N
            '[NX3][CX3](=[OX1])[CX4][OX2][CX4][CX3](=[OX1])[NX3]',

            # Flexible glycolamide moiety (at least one side)
            # Matches single glycolamide unit: N-C(=O)-C-O
            '[NX3][CX3](=[OX1])[CX4][OX2][CX4]',

            # N,N'-disubstituted diglycolamide (both nitrogens have two substituents)
            '[$([NX3]([#6])[#6])][CX3](=[OX1])[CX4][OX2][CX4][CX3](=[OX1])[$([NX3]([#6])[#6])]',
        ],
        'priority': 10,  # High priority for DGA identification
    },

    'Malonamide': {
        'description': 'Malonamide ligands with two amides separated by methylene',
        'patterns': [
            # Malonamide core: N-C(=O)-C-C(=O)-N (no oxygen between carbons)
            '[NX3][CX3](=[OX1])[CX4][CX3](=[OX1])[NX3]',

            # Exclude if it matches DGA pattern (has ether oxygen)
            # This is handled in exclusion logic
        ],
        'priority': 9,  # Slightly lower than DGA
    },

    'Pyridine_Bisamide': {
        'description': 'Pyridine-2,6-dicarboxamide; tridentate NNN coordination; flat/mixed REE selectivity',
        'patterns': [
            # Pyridine-2,6-bis-amide core: amide at both 2 and 6 positions of pyridine
            '[n]1c(C(=O)[NX3])cccc1C(=O)[NX3]',
            # Alternative notation
            'c1cc(C(=O)[NX3])nc(C(=O)[NX3])c1',
        ],
        'priority': 14,  # Highest — very specific pattern
        'exclude_if': [
            'c1cc2ccc3cccnc3c2nc1',   # exclude if phen scaffold present
        ],
    },

    'Bipyridyl_Diamide': {
        'description': '2,2-bipyridyl-6,6-bisamide; tetradentate NNNN coordination; ascending REE selectivity',
        'patterns': [
            # 2,2-bipyridyl with amides at 6,6 positions
            '[n]1c(cccc1-c2cccc(C(=O)[NX3])n2)C(=O)[NX3]',
            # Bipyridyl-diamide via bipy SMARTS + amide check
            'c1ccnc(c1)-c1cccc(C(=O)N)n1',
        ],
        'priority': 14,
        'exclude_if': [
            'c1cc2ccc3cccnc3c2nc1',   # exclude if phen scaffold present
        ],
    },

    'Phenanthroline_Phosphonate': {
        'description': '1,10-phenanthroline + phosphonate/phosphine-oxide pendants; expected REVERSE REE selectivity (POPPhen, MPPhen, BPPhen types)',
        'patterns': [
            'c1cc2ccc3cccnc3c2nc1',   # 1,10-phenanthroline scaffold
        ],
        'priority': 13,  # Between Pyridine_Bisamide (14) and Phenanthroline_Amide (12)
        'require_any': [
            '[PX4](=[OX1])',          # phosphine oxide
            'P(=O)([OX2])[OX2]',     # phosphonate ester
        ],
    },

    'Phenanthroline_Amide': {
        'description': '1,10-phenanthroline + amide pendant groups; ascending REE selectivity (OTPTA, DAPhen, BLPhen types). Note: 3a series (carbazole amide) shows anomalous reverse trend.',
        'patterns': [
            'c1cc2ccc3cccnc3c2nc1',   # 1,10-phenanthroline scaffold
        ],
        'priority': 12,  # Above Organophosphorus (11)
        'require_any': [
            '[CX3](=[OX1])[NX3]',    # amide (non-aromatic N)
            '[CX3](=[OX1])N',         # amide (permissive, catches carbazole amides)
        ],
        'exclude_if': [
            '[PX4](=[OX1])',          # defer to Phenanthroline_Phosphonate if phosphine oxide
            'P(=O)([OX2])[OX2]',     # defer to Phenanthroline_Phosphonate if phosphonate ester
        ],
    },

    'Organophosphorus': {
        'description': 'Organophosphorus ligands including phosphine oxides, phosphates, phosphonates (pure P-donors; phen+P hybrids are Phenanthroline_Phosphonate)',
        'patterns': [
            # Phosphine oxide: P(=O) with 4 bonds
            '[PX4](=[OX1])',

            # Phosphate ester: P(=O)(O)(O)O
            'P(=O)([OX2])([OX2])',

            # CMPO-like (carbamoylmethylphosphine oxide): P(=O)-C-C(=O)-N
            '[PX4](=[OX1])([#6])[CX4][CX3](=[OX1])[NX3]',

            # Phosphonic acid/phosphonate: P(=O)(O)C
            'P(=O)([OH,OX2])[CX4]',
        ],
        'priority': 11,  # Below phenanthroline types (12/13)
        'exclude_if': [
            'c1cc2ccc3cccnc3c2nc1',   # defer to Phenanthroline_* if phen scaffold present
        ],
    },

    'Aromatic_Heterocycle': {
        'description': 'Aromatic heterocycles including pyridine, triazine, phenanthroline, BTP',
        'patterns': [
            # Pyridine ring (6-membered aromatic with N)
            'n1ccccc1',

            # Triazine patterns
            'n1ncncn1',  # 1,3,5-triazine
            'n1nncnc1',  # 1,2,4-triazine
            'c1ncncn1',  # Alternative triazine notation

            # BTP-like (bis-triazinyl-pyridine) - triazine connected to pyridine
            'n1nc(c2ccccn2)ncn1',
            'c1nc(c2ncccn2)ncn1',

            # Phenanthroline (fused ring system with two nitrogens)
            'c1cnc2c(c1)ccc1ncccc12',

            # Bipyridine
            'c1ccncc1-c1ccccn1',

            # Imidazole
            'n1ccnc1',

            # Pyrazole
            'n1nccc1',

            # General aromatic nitrogen-containing heterocycle
            # At least 2 aromatic nitrogens in ring system
            'n1cccnc1',
        ],
        'priority': 8,
    },

    'Amide': {
        'description': 'Simple amide ligands (excluding DGA and malonamide)',
        'patterns': [
            # Generic amide: N-C(=O)-C
            '[NX3][CX3](=[OX1])[CX4]',

            # N,N-disubstituted amide
            '[$([NX3]([#6])[#6])][CX3](=[OX1])',

            # Multiple amides (but not DGA or malonamide pattern)
            '[NX3][CX3](=[OX1])',
        ],
        'priority': 3,  # Low priority (very general)
    },

    'Carboxylic_Acid': {
        'description': 'Carboxylic acids and carboxylates',
        'patterns': [
            # Carboxylic acid: C(=O)OH
            '[CX3](=[OX1])[OX2H1]',

            # Carboxylate: C(=O)[O-]
            '[CX3](=[OX1])[OX1-]',
        ],
        'priority': 7,
    },

    'Crown_Ether': {
        'description': 'Crown ethers and macrocyclic polyethers',
        'patterns': [
            # Multiple ether linkages in ring
            # This is complex - simplified pattern for C-O-C-O-C-O pattern
            '[OX2][CX4][CX4][OX2][CX4][CX4][OX2]',
        ],
        'priority': 9,
    },

    'Calixarene': {
        'description': 'Calixarene macrocycles',
        'patterns': [
            # Phenolic units bridged by methylene
            # Simplified: aromatic-OH with methylene bridges
            'c1ccc(O)cc1[CX4]c1ccc(O)cc1',
        ],
        'priority': 9,
    },
}


# Exclusion patterns - used to differentiate similar types
EXCLUSION_PATTERNS = {
    'Amide': {
        'exclude_if_matches': [
            # Exclude if has glycolamide ether linkage (is DGA)
            '[NX3][CX3](=[OX1])[CX4][OX2]',

            # Exclude if is malonamide
            '[NX3][CX3](=[OX1])[CX4][CX3](=[OX1])[NX3]',
        ]
    },
    'Malonamide': {
        'exclude_if_matches': [
            # Exclude if has ether oxygen (is DGA, not malonamide)
            '[NX3][CX3](=[OX1])[CX4][OX2][CX4][CX3](=[OX1])[NX3]',
        ]
    }
}


def check_rdkit_available() -> bool:
    """Check if RDKit is available"""
    return RDKIT_AVAILABLE


def match_smarts_pattern(smiles: str, pattern: str) -> bool:
    """
    Match SMILES against a single SMARTS pattern

    Args:
        smiles: SMILES string of molecule
        pattern: SMARTS pattern to match

    Returns:
        True if pattern matches, False otherwise
    """
    if not RDKIT_AVAILABLE:
        logging.error("RDKit not available. Cannot perform SMARTS matching.")
        return False

    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            logging.warning(f"Invalid SMILES: {smiles}")
            return False

        smarts = Chem.MolFromSmarts(pattern)
        if smarts is None:
            logging.warning(f"Invalid SMARTS pattern: {pattern}")
            return False

        return mol.HasSubstructMatch(smarts)

    except Exception as e:
        logging.warning(f"SMARTS matching failed for {smiles}: {e}")
        return False


def count_pattern_matches(smiles: str, pattern: str) -> int:
    """
    Count number of times a SMARTS pattern matches in a molecule

    Args:
        smiles: SMILES string
        pattern: SMARTS pattern

    Returns:
        Number of matches (0 if no match)
    """
    if not RDKIT_AVAILABLE:
        return 0

    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return 0

        smarts = Chem.MolFromSmarts(pattern)
        if smarts is None:
            return 0

        matches = mol.GetSubstructMatches(smarts)
        return len(matches)

    except Exception as e:
        logging.warning(f"Pattern counting failed: {e}")
        return 0


def match_extractant_type(smiles: str, verbose: bool = False) -> Dict[str, Any]:
    """
    Match SMILES against all extractant type patterns

    Args:
        smiles: SMILES string to analyze
        verbose: If True, include detailed matching information

    Returns:
        Dictionary with:
            - matched_types: List of extractant types that matched
            - confidence: 'high', 'medium', or 'low'
            - pattern_matches: Dict of type -> list of matched patterns
            - suggested_type: Best guess for extractant type
            - num_matches_per_type: Dict of type -> number of pattern matches
    """
    if not RDKIT_AVAILABLE:
        return {
            'matched_types': [],
            'confidence': 'error',
            'pattern_matches': {},
            'suggested_type': None,
            'error': 'RDKit not available'
        }

    if not smiles:
        return {
            'matched_types': [],
            'confidence': 'no_smiles',
            'pattern_matches': {},
            'suggested_type': None,
        }

    # Validate SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return {
            'matched_types': [],
            'confidence': 'invalid_smiles',
            'pattern_matches': {},
            'suggested_type': None,
            'error': f'Invalid SMILES: {smiles}'
        }

    # Match against all patterns
    results = {}

    for extractant_type, config in EXTRACTANT_TYPE_PATTERNS.items():
        patterns = config['patterns']
        matched_patterns = []
        total_matches = 0

        for pattern in patterns:
            if match_smarts_pattern(smiles, pattern):
                matched_patterns.append(pattern)
                total_matches += count_pattern_matches(smiles, pattern)

        if matched_patterns:
            # Check inline exclude_if list (new-style per-type exclusions)
            should_exclude = False
            for excl_pattern in config.get('exclude_if', []):
                if match_smarts_pattern(smiles, excl_pattern):
                    should_exclude = True
                    if verbose:
                        logging.info(f"Excluding {extractant_type} due to inline exclude_if match: {excl_pattern}")
                    break

            # Check legacy EXCLUSION_PATTERNS dict
            if not should_exclude and extractant_type in EXCLUSION_PATTERNS:
                for excl_pattern in EXCLUSION_PATTERNS[extractant_type]['exclude_if_matches']:
                    if match_smarts_pattern(smiles, excl_pattern):
                        should_exclude = True
                        if verbose:
                            logging.info(f"Excluding {extractant_type} due to exclusion pattern match")
                        break

            if should_exclude:
                continue

            # Check require_any list: at least one must match
            require_any = config.get('require_any', [])
            if require_any:
                has_required = any(match_smarts_pattern(smiles, p) for p in require_any)
                if not has_required:
                    if verbose:
                        logging.info(f"Skipping {extractant_type}: no require_any pattern matched")
                    continue

            results[extractant_type] = {
                'matched_patterns': matched_patterns,
                'num_matches': total_matches,
                'priority': config['priority']
            }

    # Determine matched types
    matched_types = list(results.keys())

    # Determine confidence and suggested type
    if len(matched_types) == 0:
        confidence = 'no_match'
        suggested_type = 'Unknown'
    elif len(matched_types) == 1:
        confidence = 'high'
        suggested_type = matched_types[0]
    else:
        # Multiple matches - use priority to determine best match
        confidence = 'medium'
        sorted_types = sorted(
            matched_types,
            key=lambda t: results[t]['priority'],
            reverse=True
        )
        suggested_type = sorted_types[0]

        # If top two priorities are very close and both have strong matches, lower confidence
        if len(sorted_types) >= 2:
            top_priority = results[sorted_types[0]]['priority']
            second_priority = results[sorted_types[1]]['priority']
            if abs(top_priority - second_priority) <= 2:
                confidence = 'low'

    result = {
        'matched_types': matched_types,
        'confidence': confidence,
        'pattern_matches': {k: v['matched_patterns'] for k, v in results.items()},
        'suggested_type': suggested_type,
        'num_matches_per_type': {k: v['num_matches'] for k, v in results.items()},
    }

    if verbose:
        result['priorities'] = {k: v['priority'] for k, v in results.items()}

    return result


def validate_single_ligand(ligand_data: Dict[str, Any]) -> Dict[str, Any]:
    """
    Validate a single ligand's extractant_type against SMILES structure

    Args:
        ligand_data: Dictionary with keys 'ligand_name', 'smiles', 'extractant_type'

    Returns:
        Validation report with:
            - ligand_name: Name of ligand
            - current_type: Current extractant_type assignment
            - smiles: SMILES string
            - smiles_predicted_type: Type predicted from SMILES
            - matched_types: All types that matched
            - match_status: 'correct' | 'mismatch' | 'ambiguous' | 'no_smiles' | 'invalid_smiles'
            - confidence: Confidence of SMILES prediction
            - suggested_correction: Suggested type if mismatch
            - notes: Additional information
    """
    ligand_name = ligand_data.get('ligand_name', 'Unknown')
    current_type = ligand_data.get('extractant_type', 'Unknown')
    smiles = ligand_data.get('smiles', '')

    # Handle missing SMILES
    if not smiles or smiles == '':
        return {
            'ligand_name': ligand_name,
            'current_type': current_type,
            'smiles': '',
            'smiles_predicted_type': None,
            'matched_types': [],
            'match_status': 'no_smiles',
            'confidence': 'no_smiles',
            'suggested_correction': None,
            'notes': 'No SMILES available for validation'
        }

    # Match against patterns
    match_result = match_extractant_type(smiles)

    # Handle invalid SMILES
    if match_result['confidence'] in ['invalid_smiles', 'error']:
        return {
            'ligand_name': ligand_name,
            'current_type': current_type,
            'smiles': smiles,
            'smiles_predicted_type': None,
            'matched_types': [],
            'match_status': 'invalid_smiles',
            'confidence': match_result['confidence'],
            'suggested_correction': None,
            'notes': match_result.get('error', 'Invalid SMILES')
        }

    predicted_type = match_result['suggested_type']
    matched_types = match_result['matched_types']
    confidence = match_result['confidence']

    # Determine match status
    if confidence == 'no_match':
        match_status = 'no_match'
        suggested_correction = 'Unknown' if current_type != 'Unknown' else None
        notes = 'No SMARTS patterns matched'
    elif current_type == predicted_type:
        # Current type matches the best (highest priority) prediction
        match_status = 'correct'
        suggested_correction = None
        notes = 'Current type matches SMILES structure'
    elif current_type in matched_types:
        # Current type matches a pattern, but not the highest priority one
        match_status = 'ambiguous'
        suggested_correction = predicted_type
        notes = f'Current type matches, but {predicted_type} is higher priority. Matched: {", ".join(matched_types)}'
    elif len(matched_types) > 1:
        match_status = 'ambiguous'
        suggested_correction = predicted_type if predicted_type != current_type else None
        notes = f'Multiple types matched: {", ".join(matched_types)}'
    else:
        match_status = 'mismatch'
        suggested_correction = predicted_type
        notes = f'Current type does not match SMILES structure'

    return {
        'ligand_name': ligand_name,
        'current_type': current_type,
        'smiles': smiles,
        'smiles_predicted_type': predicted_type,
        'matched_types': matched_types,
        'match_status': match_status,
        'confidence': confidence,
        'suggested_correction': suggested_correction,
        'notes': notes
    }


def get_pattern_description(extractant_type: str) -> str:
    """Get description of SMARTS patterns for an extractant type"""
    if extractant_type in EXTRACTANT_TYPE_PATTERNS:
        return EXTRACTANT_TYPE_PATTERNS[extractant_type]['description']
    return f"No pattern definition for {extractant_type}"


def list_all_extractant_types() -> List[str]:
    """Get list of all defined extractant types"""
    return list(EXTRACTANT_TYPE_PATTERNS.keys())


if __name__ == "__main__":
    # Test with example ligands
    print("Testing SMARTS Pattern Library")
    print("=" * 60)

    test_cases = [
        {
            'name': 'TODGA',
            'smiles': 'CCCCCCCCN(CCCCCCCC)C(=O)COCCOC(=O)N(CCCCCCCC)CCCCCCCC',
            'expected': 'Diglycolamide'
        },
        {
            'name': 'CMPO',
            'smiles': 'CCCCCCCCP(=O)(c1ccccc1)CC(=O)N(CC(C)C)CC(C)C',
            'expected': 'Organophosphorus'
        },
        {
            'name': 'BTP',
            'smiles': 'CCCN1N=C(C2=CC=CC(C3=NC(NCCC)=NN3)=N2)C=NC1',
            'expected': 'Aromatic_Heterocycle'
        }
    ]

    for test in test_cases:
        print(f"\nTesting: {test['name']} (Expected: {test['expected']})")
        print(f"SMILES: {test['smiles']}")
        result = match_extractant_type(test['smiles'], verbose=True)
        print(f"Matched types: {result['matched_types']}")
        print(f"Suggested type: {result['suggested_type']}")
        print(f"Confidence: {result['confidence']}")

        if result['suggested_type'] == test['expected']:
            print("✓ PASSED")
        else:
            print("✗ FAILED")
