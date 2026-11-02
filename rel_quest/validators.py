"""
Validation functions for REEL-QUEST
Validates data quality before and after ingestion
"""

import math

import pandas as pd
from typing import Dict, List, Tuple
from .schema import KNOWN_SCHEMA_VERSIONS, SCHEMA_VERSION, TREND_CLASSIFICATIONS, schema_at_least
from .references import check_reference
from .utils import validate_smiles, load_vocabulary

# A stored D and its logarithm may disagree by at most this, in log units.
# The sources round log D to 3-4 decimals, so exact equality is unachievable;
# 0.01 dex (2.3% in D) is loose enough for that yet tight enough to catch the
# decimal-rounding collapse that zeroed 39 measurements before 2026-08-13.
D_LOG_TOLERANCE_DEX = 0.01


def validate_liu_dataframe(df: pd.DataFrame) -> Tuple[bool, List[str]]:
    """
    Validate Liu 2022 Excel data before processing

    Args:
        df: DataFrame loaded from Liu Excel file

    Returns:
        Tuple of (is_valid: bool, error_messages: List[str])
    """
    errors = []

    # Check required columns exist
    required_cols = ['SMILES', 'Ligand_Name', 'Structural_Class', 'log_D',
                     'Atomic Number_metal ', 'acid c.c./M', 'temperature']

    missing_cols = [col for col in required_cols if col not in df.columns]
    if missing_cols:
        errors.append(f"Missing required columns: {missing_cols}")
        return False, errors

    # Check for null SMILES
    null_smiles = df['SMILES'].isna().sum()
    if null_smiles > 0:
        errors.append(f"Found {null_smiles} rows with null SMILES")

    # Validate log_D range
    invalid_logd = df['log_D'].notna() & (~df['log_D'].between(-5, 6))
    if invalid_logd.any():
        count = invalid_logd.sum()
        errors.append(f"Found {count} rows with log_D outside reasonable range [-5, 6]")

    # Validate atomic numbers
    valid_atomic_nums = list(range(57, 72))  # La to Lu
    invalid_atomic = df['Atomic Number_metal '].notna() & (~df['Atomic Number_metal '].isin(valid_atomic_nums))
    if invalid_atomic.any():
        count = invalid_atomic.sum()
        errors.append(f"Found {count} rows with invalid atomic numbers (should be 57-71)")

    # Validate acid concentration
    invalid_acid = df['acid c.c./M'].notna() & (~df['acid c.c./M'].between(0, 15))
    if invalid_acid.any():
        count = invalid_acid.sum()
        errors.append(f"Found {count} rows with acid concentration outside [0, 15] M")

    # Validate temperature
    invalid_temp = df['temperature'].notna() & (~df['temperature'].between(200, 400))
    if invalid_temp.any():
        count = invalid_temp.sum()
        errors.append(f"Found {count} rows with temperature outside [200, 400] K")

    # Check SMILES validity (sample first 10)
    smiles_sample = df['SMILES'].dropna().unique()[:10]
    for smiles in smiles_sample:
        if not validate_smiles(smiles):
            errors.append(f"Invalid SMILES detected: {smiles[:50]}...")
            break

    is_valid = len(errors) == 0
    return is_valid, errors


def check_d_value_consistency(distribution_ratios: Dict) -> List[str]:
    """Check that each stored D agrees with its stored logarithm.

    Both forms are served to consumers, and code disagrees about which to
    read: core.get_separation_factor divides raw D, while the trend fit and
    the plots use log D. When the two disagree the answer depends on which
    one you happen to ask, so the disagreement is itself the defect.

    Returns a list of problem descriptions; empty means consistent.
    """
    problems = []
    for element, data in (distribution_ratios or {}).items():
        if not isinstance(data, dict):
            problems.append(f"{element}: distribution ratio is not an object")
            continue

        value, log_value = data.get('value'), data.get('log_value')
        if value is None and log_value is None:
            continue

        for label, number in (('value', value), ('log_value', log_value)):
            if number is not None and (
                not isinstance(number, (int, float))
                or isinstance(number, bool)
                or not math.isfinite(number)
            ):
                problems.append(f"{element}: {label} is not a finite number ({number!r})")

        if value is None or log_value is None:
            continue
        if not isinstance(value, (int, float)) or not isinstance(log_value, (int, float)):
            continue
        if not math.isfinite(value) or not math.isfinite(log_value):
            continue

        if value <= 0:
            problems.append(
                f"{element}: value is {value} but log_value is {log_value} "
                "(a finite log implies a positive D)")
            continue

        error = abs(math.log10(value) - log_value)
        if error > D_LOG_TOLERANCE_DEX:
            problems.append(
                f"{element}: log10(value)={math.log10(value):.4f} disagrees with "
                f"log_value={log_value} by {error:.4f} dex "
                f"(max {D_LOG_TOLERANCE_DEX})")
    return problems


def check_element_coverage(record: Dict) -> List[str]:
    """Check the two element lists against the measurements themselves.

    `elements_measured` must list exactly what was measured (Y and Sc
    included); `selectivity.element_coverage` lists the La-Lu subset the
    trend fit used and so must be contained in it. Records predating
    elements_measured are skipped rather than failed.
    """
    problems = []
    measured = set()
    for exp_data in record.get('experimental_data', []):
        for element, data in (exp_data.get('distribution_ratios') or {}).items():
            if isinstance(data, dict) and (
                data.get('value') is not None or data.get('log_value') is not None
            ):
                measured.add(element)

    stored = record.get('elements_measured')
    if stored is not None:
        missing = sorted(measured - set(stored))
        extra = sorted(set(stored) - measured)
        if missing:
            problems.append(f"elements_measured omits measured elements: {missing}")
        if extra:
            problems.append(f"elements_measured lists unmeasured elements: {extra}")

    coverage = set(record.get('selectivity', {}).get('element_coverage') or [])
    unmeasured = sorted(coverage - measured)
    if unmeasured:
        problems.append(
            f"element_coverage lists elements with no measurement: {unmeasured}")
    return problems


def validate_ligand_record(record: Dict) -> Tuple[bool, List[str]]:
    """
    Validate a ligand record against schema requirements

    Args:
        record: Ligand dictionary to validate

    Returns:
        Tuple of (is_valid: bool, error_messages: List[str])
    """
    errors = []

    # Check schema version. The value matters, not just its presence — it was
    # stamped "v0.1" long after the record shape had moved on.
    if 'schema_version' not in record:
        errors.append("Missing schema_version")
    elif record['schema_version'] not in KNOWN_SCHEMA_VERSIONS:
        errors.append(
            f"Unknown schema_version {record['schema_version']!r}. "
            f"Must be one of: {list(KNOWN_SCHEMA_VERSIONS)}"
        )

    # Check required ligand fields (from CONTRIBUTING.md)
    ligand = record.get('ligand', {})
    required_ligand_fields = ['Name', 'Mechanism', 'Extractant_Type']

    # SMILES is required for Liu 2022 but optional for IDEAL (can be added later)
    data_source = record.get('metadata', {}).get('data_source', '')
    if data_source == 'Liu_2022_JACSAu' and not ligand.get('SMILES'):
        errors.append(f"Missing required ligand field: SMILES")

    for field in required_ligand_fields:
        if not ligand.get(field):
            errors.append(f"Missing required ligand field: {field}")

    # Validate SMILES
    if ligand.get('SMILES') and not validate_smiles(ligand['SMILES']):
        errors.append(f"Invalid SMILES: {ligand['SMILES'][:50]}...")

    # Validate mechanism is from canonical list
    mechanisms_vocab = load_vocabulary('mechanisms')
    canonical_mechanisms = mechanisms_vocab['canonical_mechanisms']
    if ligand.get('Mechanism') and ligand['Mechanism'] not in canonical_mechanisms:
        errors.append(f"Invalid mechanism '{ligand['Mechanism']}'. Must be one of: {canonical_mechanisms}")

    # Check experimental_data exists and is non-empty
    if 'experimental_data' not in record:
        errors.append("Missing experimental_data field")
    elif not record['experimental_data']:
        errors.append("experimental_data is empty (must have at least 1 condition)")

    # Structured provenance (bare DOI, Locator, Provenance_Level, Aggregator)
    # exists from schema v0.3; archived generations keep their older shape.
    structured_references = schema_at_least(record.get('schema_version'), 'v0.3')

    # Validate each experimental condition
    for i, exp_data in enumerate(record.get('experimental_data', [])):
        # Check reference exists
        if 'reference' not in exp_data:
            errors.append(f"Experimental condition {i} missing reference")
        elif not exp_data['reference'].get('DOI') and data_source == 'Liu_2022_JACSAu':
            # DOI required for Liu 2022, but optional for IDEAL (references in IDEAL database)
            errors.append(f"Experimental condition {i} missing DOI in reference")
        if structured_references and 'reference' in exp_data:
            errors.extend(
                f"Experimental condition {i}: {problem}"
                for problem in check_reference(exp_data['reference'])
            )

        # Check distribution_ratios
        if 'distribution_ratios' not in exp_data:
            errors.append(f"Experimental condition {i} missing distribution_ratios")
        elif len(exp_data['distribution_ratios']) < 1:
            errors.append(f"Experimental condition {i} has no elements (need at least 1)")
        else:
            errors.extend(
                f"Experimental condition {i}: {problem}"
                for problem in check_d_value_consistency(exp_data['distribution_ratios'])
            )

    # Check selectivity classification
    selectivity = record.get('selectivity', {})
    if selectivity.get('trend_classification') not in TREND_CLASSIFICATIONS:
        errors.append(
            f"Invalid trend_classification. Must be one of: "
            f"{list(TREND_CLASSIFICATIONS)}")

    errors.extend(check_element_coverage(record))

    # Check metadata
    if 'metadata' not in record:
        errors.append("Missing metadata field")

    is_valid = len(errors) == 0
    return is_valid, errors


_JSON_SCHEMA_VALIDATOR = None


def validate_against_json_schema(record: Dict) -> List[str]:
    """Problems with `record` against schema/ligand_record.schema.json.

    The JSON Schema describes the current record shape (SCHEMA_VERSION), so
    records from archived generations are not checked against it. It is the
    machine-readable contract published with the data; `validate_ligand_record`
    adds the cross-field checks a schema cannot express (D vs log D,
    coverage vs measurements, a compilation DOI cited as a measurement).
    """
    global _JSON_SCHEMA_VALIDATOR
    if record.get('schema_version') != SCHEMA_VERSION:
        return []
    if _JSON_SCHEMA_VALIDATOR is None:
        import json
        import os
        from jsonschema import Draft202012Validator

        path = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                            'schema', 'ligand_record.schema.json')
        with open(path, encoding='utf-8') as handle:
            _JSON_SCHEMA_VALIDATOR = Draft202012Validator(json.load(handle))
    return [
        f"schema: {'/'.join(str(p) for p in error.absolute_path) or '(record)'}: {error.message[:160]}"
        for error in sorted(_JSON_SCHEMA_VALIDATOR.iter_errors(record), key=lambda e: list(e.absolute_path).__repr__())
    ]


def validate_database_file(jsonl_path: str, max_records: int = None,
                           use_json_schema: bool = True) -> Tuple[int, int, List[str]]:
    """
    Validate entire JSONL database file

    Args:
        jsonl_path: Path to JSONL database file
        max_records: Maximum number of records to validate (None = all)
        use_json_schema: also check current-version records against
            schema/ligand_record.schema.json

    Returns:
        Tuple of (valid_count, invalid_count, error_summary)
    """
    import json

    valid_count = 0
    invalid_count = 0
    error_summary = []

    with open(jsonl_path, 'r') as f:
        for line_num, line in enumerate(f, 1):
            if max_records and line_num > max_records:
                break

            try:
                record = json.loads(line)
                is_valid, errors = validate_ligand_record(record)
                if use_json_schema:
                    schema_errors = validate_against_json_schema(record)
                    errors = errors + schema_errors
                    is_valid = is_valid and not schema_errors

                if is_valid:
                    valid_count += 1
                else:
                    invalid_count += 1
                    ligand_name = record.get('ligand', {}).get('Name', 'Unknown')
                    error_summary.append(f"Line {line_num} ({ligand_name}): {'; '.join(errors)}")

            except json.JSONDecodeError as e:
                invalid_count += 1
                error_summary.append(f"Line {line_num}: Invalid JSON - {e}")

    return valid_count, invalid_count, error_summary
