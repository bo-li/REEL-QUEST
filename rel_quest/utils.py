"""
Utility functions for REX DB
Includes SMILES processing, pattern matching, and helper functions
"""

from .schema import SCHEMA_VERSION
import json
import os
from typing import Dict, Optional, List

# Note: RDKit will be optional dependency
# Code will check for availability and provide graceful degradation

def load_json_file(filepath: str) -> Dict:
    """Load JSON file from relative path"""
    # Get the Database directory (parent of rel_quest)
    db_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    full_path = os.path.join(db_dir, filepath)

    with open(full_path, 'r') as f:
        return json.load(f)


def load_vocabulary(vocab_name: str) -> Dict:
    """
    Load controlled vocabulary file

    Args:
        vocab_name: Name of vocab file (e.g., 'mechanisms', 'elements')

    Returns:
        Dictionary containing vocabulary data
    """
    return load_json_file(f'vocab/{vocab_name}.json')


def canonicalize_smiles(smiles: str) -> Optional[str]:
    """
    Convert SMILES to canonical form using RDKit

    Args:
        smiles: SMILES string

    Returns:
        Canonical SMILES string, or None if RDKit not available or invalid SMILES
    """
    try:
        from rdkit import Chem
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None
        return Chem.MolToSmiles(mol, canonical=True)
    except ImportError:
        # RDKit not available, return original SMILES
        print("Warning: RDKit not available. Using original SMILES without canonicalization.")
        return smiles
    except Exception as e:
        print(f"Error canonicalizing SMILES: {e}")
        return None


def validate_smiles(smiles: str) -> bool:
    """
    Check if SMILES string is valid

    Args:
        smiles: SMILES string to validate

    Returns:
        True if valid, False otherwise
    """
    try:
        from rdkit import Chem
        mol = Chem.MolFromSmiles(smiles)
        return mol is not None
    except ImportError:
        # If RDKit not available, assume valid if non-empty
        return bool(smiles and smiles.strip())
    except:
        return False


def get_atomic_number_from_symbol(symbol: str) -> Optional[int]:
    """
    Get atomic number from element symbol

    Args:
        symbol: Element symbol (e.g., 'La', 'Ce')

    Returns:
        Atomic number, or None if not found
    """
    elements = load_vocabulary('elements')

    # Check lanthanides
    for elem in elements['lanthanides']:
        if elem['symbol'] == symbol:
            return elem['atomic_number']

    # Check other REEs
    for elem in elements.get('other_rees', []):
        if elem['symbol'] == symbol:
            return elem['atomic_number']

    return None


def get_lanthanide_order() -> List[str]:
    """
    Get ordered list of lanthanide symbols (La to Lu)

    Returns:
        List of element symbols in atomic number order
    """
    elements = load_vocabulary('elements')
    return [elem['symbol'] for elem in elements['lanthanides']]


def map_structural_class_to_mechanism(structural_class: str) -> Dict[str, str]:
    """
    Map Liu 2022 Structural_Class to CONTRIBUTING.md mechanism types

    Args:
        structural_class: Structural class from Liu 2022 (e.g., "Diglycolamide (DGA)")

    Returns:
        Dictionary with 'mechanism' and 'extractant_type' keys
    """
    mechanisms = load_vocabulary('mechanisms')
    mapping = mechanisms['mechanism_mapping'].get(structural_class, {})

    if not mapping:
        # Default to Solvation / Other if not found
        return {
            'mechanism': 'Solvation',
            'extractant_type': 'Other'
        }

    return {
        'mechanism': mapping['mechanism'],
        'extractant_type': mapping['extractant_type']
    }


def round_if_not_none(value, decimals: int = 3):
    """Round value if not None, otherwise return None"""
    if value is None or (isinstance(value, float) and (value != value)):  # Check for NaN
        return None
    try:
        return round(float(value), decimals)
    except:
        return None


# Significant figures kept on a stored D-value. Six is well beyond the
# precision of a digitized figure and keeps 10 ** log_D reversible.
D_VALUE_SIG_FIGURES = 6


def round_significant(value, figures: int = D_VALUE_SIG_FIGURES):
    """Round to significant figures, i.e. to a *relative* precision.

    Use this, not round_if_not_none, for distribution ratios. Decimal
    rounding is an absolute precision, so it annihilates small D-values:
    round(0.0001, 3) is 0.0, which then reads as "no measurement" to
    core.get_separation_factor and the plotting layer while log_value still
    says -4.0. That mismatch zeroed 39 measurements before it was corrected
    on 2026-08-13 (scripts/fix_d_value_consistency.py); validators.py now
    rejects it.
    """
    import math

    if value is None or (isinstance(value, float) and value != value):  # NaN
        return None
    try:
        value = float(value)
    except (TypeError, ValueError):
        return None
    if value == 0 or not math.isfinite(value):
        return value
    magnitude = math.floor(math.log10(abs(value)))
    return round(value, -(magnitude - figures + 1))


def collect_elements_measured(record: Dict) -> List[str]:
    """Every element the record actually has a measurement for, Y/Sc included.

    This is what `elements_measured` stores, and it answers a different
    question from `selectivity.element_coverage`: coverage lists the La-Lu
    points the trend fit consumed, so it necessarily excludes Y and Sc (a
    regression of log D on atomic number cannot take Z=39 and Z=21 alongside
    57-71). Filters asking "was this element measured?" must read this one —
    core.filter_by_elements returned nothing for Y and Sc until it did.

    Returns the symbols in periodic order.
    """
    order = {'Sc': 21, 'Y': 39}
    order.update({sym: 57 + i for i, sym in enumerate(get_lanthanide_order())})

    found = set()
    for exp_data in record.get('experimental_data', []):
        for element, data in (exp_data.get('distribution_ratios') or {}).items():
            if isinstance(data, dict) and (
                data.get('value') is not None or data.get('log_value') is not None
            ):
                found.add(element)
    return sorted(found, key=lambda e: (order.get(e, 999), e))


# ============================================================================
# Merge Utilities for Duplicate Ligands
# ============================================================================

def find_duplicate_ligands(db_records: List[Dict]) -> Dict[str, List[Dict]]:
    """
    Identify ligands that appear in multiple sources

    Matching criteria:
    1. Exact name match (case-insensitive)
    2. Canonical SMILES match (if both have SMILES)

    Args:
        db_records: List of ligand records from database

    Returns:
        Dict mapping ligand identifier to list of duplicate records
    """
    from collections import defaultdict

    # Index by name and canonical SMILES
    by_name = defaultdict(list)
    by_smiles = defaultdict(list)

    for record in db_records:
        name = record['ligand']['Name'].upper().strip()
        by_name[name].append(record)

        smiles = record['ligand'].get('SMILES')
        if smiles:
            canonical = canonicalize_smiles(smiles)
            if canonical:
                by_smiles[canonical].append(record)

    # Find duplicates (appears >1 time)
    duplicates = {}

    # Match by name
    for name, records in by_name.items():
        if len(records) > 1:
            duplicates[name] = records

    # Match by SMILES (may find additional duplicates with different names)
    for smiles, records in by_smiles.items():
        if len(records) > 1:
            # Check if not already found by name
            names = [r['ligand']['Name'].upper().strip() for r in records]
            if len(set(names)) > 1:
                # Different names, same SMILES - add with SMILES key
                key = f"SMILES:{smiles[:20]}"
                duplicates[key] = records

    return duplicates


def resolve_smiles_conflict(sources: Dict[str, str]) -> Optional[str]:
    """
    Resolve SMILES conflicts when merging ligands

    Priority order: Liu 2022 > SAFE > IDEAL

    Args:
        sources: Dict mapping data_source to SMILES string

    Returns:
        Best SMILES string, or None if none available

    Raises:
        ValueError: If sources have different canonical SMILES
    """
    # Priority order
    priority = ['Liu_2022_JACSAu', 'SAFE_Database', 'OECD_NEA_IDEAL']

    # Collect SMILES from all sources
    smiles_by_source = {}
    canonical_by_source = {}

    for source, smiles in sources.items():
        if smiles and smiles.strip():
            canonical = canonicalize_smiles(smiles)
            if canonical:
                smiles_by_source[source] = smiles
                canonical_by_source[source] = canonical

    if not smiles_by_source:
        return None

    # Check if all canonical forms match
    unique_canonical = set(canonical_by_source.values())

    if len(unique_canonical) > 1:
        # Conflict - different canonical forms
        conflict_msg = "\n".join([f"  {src}: {can}"
                                  for src, can in canonical_by_source.items()])
        raise ValueError(f"SMILES conflict detected:\n{conflict_msg}")

    # All match or only one source has SMILES - use priority order
    for src in priority:
        if src in smiles_by_source:
            return smiles_by_source[src]

    # Fallback: use any available SMILES
    return list(smiles_by_source.values())[0]


def resolve_mechanism_conflict(sources: Dict[str, tuple]) -> tuple:
    """
    Resolve mechanism conflicts when merging ligands

    Priority order: Liu 2022 > SAFE > IDEAL

    Args:
        sources: Dict mapping data_source to (Mechanism, Extractant_Type) tuple

    Returns:
        (Mechanism, Extractant_Type) tuple
    """
    # Priority order
    priority = ['Liu_2022_JACSAu', 'SAFE_Database', 'OECD_NEA_IDEAL']

    # Check for consensus
    unique_mechs = set(sources.values())

    if len(unique_mechs) == 1:
        # All sources agree
        return list(unique_mechs)[0]

    # Use priority order
    for src in priority:
        if src in sources:
            return sources[src]

    # Fallback
    return list(sources.values())[0]


def resolve_name_conflict(sources: Dict[str, str]) -> str:
    """
    Resolve name conflicts when merging ligands

    Prefer short acronyms (< 15 chars) over long chemical names

    Args:
        sources: Dict mapping data_source to ligand name

    Returns:
        Best ligand name
    """
    names = list(sources.values())

    # Prefer short acronyms (< 15 chars)
    short_names = [n for n in names if len(n) < 15]

    if short_names:
        # Return shortest
        return min(short_names, key=len)

    # Return shortest overall
    return min(names, key=len)


def calculate_merged_selectivity(conditions: List[Dict],
                                 by_source: Dict[str, List[Dict]] = None) -> Dict:
    """
    Calculate selectivity trends across merged conditions from multiple sources

    Args:
        conditions: All experimental conditions (from all sources)
        by_source: Optional dict mapping data_source to its conditions

    Returns:
        Dict with:
        - trend_classification: Overall consensus trend
        - trend_confidence: High/Medium/Low
        - element_coverage: All elements tested
        - trends_by_source: Per-source trend details
    """
    from rel_quest.trends import classify_selectivity_trend

    # Get all elements tested
    all_elements = set()
    for cond in conditions:
        all_elements.update(cond['distribution_ratios'].keys())

    # Only include lanthanides for trend classification
    lanthanides = set(get_lanthanide_order())
    lanthanide_elements = sorted(all_elements & lanthanides)

    # Calculate trends per source if provided
    trends_by_source = {}
    if by_source:
        for source, source_conds in by_source.items():
            # Pass experimental_data list directly to classify_selectivity_trend
            trend_result = classify_selectivity_trend(source_conds)
            trends_by_source[source] = trend_result

    # Determine consensus
    if trends_by_source:
        trend_classifications = [t['trend_classification']
                                for t in trends_by_source.values()]
        unique_trends = set(trend_classifications)

        # Remove Insufficient_Data from consensus calculation
        non_insufficient = [t for t in unique_trends if t != 'Insufficient_Data']

        if len(non_insufficient) == 0:
            # All sources have insufficient data
            consensus = 'Insufficient_Data'
            confidence = 'Low'
        elif len(non_insufficient) == 1:
            # One clear trend (or all agree)
            consensus = non_insufficient[0]
            confidence = 'High' if len(unique_trends) == 1 else 'Medium'
        else:
            # Multiple different trends
            consensus = 'Mixed'
            confidence = 'Low'
    else:
        # No source breakdown, classify overall
        # Pass experimental_data list directly
        trend_result = classify_selectivity_trend(conditions)
        consensus = trend_result['trend_classification']
        confidence = trend_result.get('trend_confidence', 'Medium')
        trends_by_source = {}

    return {
        'trend_classification': consensus,
        'trend_confidence': confidence,
        'element_coverage': lanthanide_elements,
        'trends_by_source': trends_by_source
    }


def merge_ligand_records(records: List[Dict]) -> Dict:
    """
    Merge multiple records for the same ligand from different data sources

    Args:
        records: List of ligand records to merge (must be same ligand)

    Returns:
        Single merged record with combined experimental_data
    """
    from datetime import datetime

    if not records:
        raise ValueError("Cannot merge empty list of records")

    if len(records) == 1:
        return records[0]

    # Use first record as template
    primary = records[0]

    # Resolve conflicts in ligand metadata
    sources_smiles = {r['metadata']['data_source']: r['ligand'].get('SMILES')
                     for r in records}
    sources_mech = {r['metadata']['data_source']:
                   (r['ligand']['Mechanism'], r['ligand']['Extractant_Type'])
                   for r in records}
    sources_name = {r['metadata']['data_source']: r['ligand']['Name']
                   for r in records}

    try:
        resolved_smiles = resolve_smiles_conflict(sources_smiles)
    except ValueError as e:
        print(f"Warning: {e}")
        print("Using SMILES from primary source")
        resolved_smiles = primary['ligand'].get('SMILES')

    resolved_mech = resolve_mechanism_conflict(sources_mech)
    resolved_name = resolve_name_conflict(sources_name)

    # Build merged record
    merged = {
        'schema_version': SCHEMA_VERSION,
        'ligand': {
            'Name': resolved_name,
            'Mechanism': resolved_mech[0],
            'Extractant_Type': resolved_mech[1]
        }
    }

    # Add SMILES if available
    if resolved_smiles:
        merged['ligand']['SMILES'] = resolved_smiles

    # Add Structural_Class if present in primary
    if 'Structural_Class' in primary['ligand']:
        merged['ligand']['Structural_Class'] = primary['ligand']['Structural_Class']

    # Combine experimental_data from all sources
    all_conditions = []
    by_source = {}

    for record in records:
        source = record['metadata']['data_source']
        source_conditions = []

        for cond in record['experimental_data']:
            # Add data_source tag to condition
            cond_copy = cond.copy()
            cond_copy['data_source'] = source
            all_conditions.append(cond_copy)
            source_conditions.append(cond_copy)

        by_source[source] = source_conditions

    merged['experimental_data'] = all_conditions

    # Recalculate selectivity across merged conditions
    merged['selectivity'] = calculate_merged_selectivity(all_conditions, by_source)

    # Recalculate performance summary
    total_measurements = sum(
        len(cond['distribution_ratios'])
        for cond in all_conditions
    )

    # Find best separation pair across all conditions
    best_sf = 0
    best_pair = None

    for cond in all_conditions:
        d_vals = {elem: data['value']
                 for elem, data in cond['distribution_ratios'].items()}

        if len(d_vals) >= 2:
            max_elem = max(d_vals.items(), key=lambda x: x[1])
            min_elem = min(d_vals.items(), key=lambda x: x[1])

            if min_elem[1] > 0:
                sf = max_elem[1] / min_elem[1]
                if sf > best_sf:
                    best_sf = sf
                    best_pair = {
                        'elements': [max_elem[0], min_elem[0]],
                        'SF_max': round(sf, 2)
                    }

    merged['performance_summary'] = {
        'conditions_tested': len(all_conditions),
        'total_measurements': total_measurements
    }

    if best_pair:
        merged['performance_summary']['best_separation_pair'] = best_pair

    # Update metadata
    merged['metadata'] = {
        'data_sources': [r['metadata']['data_source'] for r in records],
        'primary_source': primary['metadata']['data_source'],
        'merge_date': datetime.now().strftime('%Y-%m-%d'),
        'curation_status': 'Merged'
    }

    return merged


def condition_statistics(d_ratios: Dict[str, Dict]) -> Dict:
    """The canonical per-condition `statistics` block.

    One shape for every writer: `mean_log_d`, `std_log_d` (population),
    `min_log_d`, `max_log_d`, `n_elements`, all logs rounded to 3 decimals.
    Until 2026-08-25 the corpus carried two vocabularies — ingest_liu wrote
    D_min/D_max/.../element_count while ingest_safe wrote this shape — and
    178 of 665 conditions had the other one. Both ingest modules now
    delegate here, `scripts/fix_condition_statistics.py` rewrote the stored
    blocks, and tests/test_data_integrity.py pins stored = fresh.

    The block is a cache of the condition's own `distribution_ratios`;
    consumers needing raw-D statistics should read the cells. Returns {}
    when the condition has no finite log D (mirrors what ingest_safe
    always returned).
    """
    import numpy as np

    # A censored cell holds a bound (D > 100), not a value: left out.
    log_values = [v['log_value'] for v in (d_ratios or {}).values()
                  if isinstance(v, dict) and v.get('log_value') is not None
                  and v.get('censor') not in ('>', '<')]
    if not log_values:
        return {}
    # round() is applied to the np.float64 BEFORE float(): numpy's rounding
    # resolves exact .0005 halves to the even digit, and the 487 stored
    # blocks were written that way — float(round(...)) reproduces them,
    # round(float(...)) flips four of them by 0.001.
    return {
        'mean_log_d': float(round(np.mean(log_values), 3)),
        'std_log_d': float(round(np.std(log_values), 3)),
        'min_log_d': round(min(log_values), 3),
        'max_log_d': round(max(log_values), 3),
        'n_elements': len(log_values),
    }
