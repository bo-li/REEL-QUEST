#!/usr/bin/env python3
"""
Solvent property mapping module for Liu 2022 data ingestion

Maps experimental solvent properties (MW, density, BP) to solvent names
using the Liu 2022 solvent reference table.

Functions:
- load_solvent_reference_table(): Load Liu 2022 solvent reference table
- match_solvent_by_properties(): Match solvent properties to name
- format_mixed_solvent(): Format mixed solvent name with ratio
- get_solvent_name(): Extract solvent name(s) from experimental data row
"""

import pandas as pd
import numpy as np
from typing import Optional, Tuple
from pathlib import Path
import math


# Cache for solvent reference table
_SOLVENT_REF_TABLE = None


def load_solvent_reference_table(excel_path: str = None) -> pd.DataFrame:
    """
    Load Liu 2022 solvent reference table from Excel

    Args:
        excel_path: Path to au2c00122_si_002.xlsx (optional, uses default if None)

    Returns:
        DataFrame with columns: name, molar_mass, density, boiling_point
    """
    global _SOLVENT_REF_TABLE

    # Return cached table if available
    if _SOLVENT_REF_TABLE is not None:
        return _SOLVENT_REF_TABLE

    # Default path if not specified
    if excel_path is None:
        excel_path = 'Liu_2022_JACSAu/au2c00122_si_002.xlsx'

    # Load solvent reference sheet
    df = pd.read_excel(excel_path, sheet_name='organic solvents descriptors')

    # Rename columns for easier access
    df_ref = pd.DataFrame({
        'name': df['Unnamed: 0'],
        'molar_mass': df['Molar mass g/mol'],
        'density': df['density g/mL'],
        'boiling_point': df['boiling point K']
    })

    # Cache for future calls
    _SOLVENT_REF_TABLE = df_ref

    return df_ref


def match_solvent_by_properties(mw: float, density: float, bp: float,
                                tolerance_pct: float = 2.0,
                                ref_table: pd.DataFrame = None) -> Optional[str]:
    """
    Match solvent properties to reference table name

    Uses tolerance-based matching (±2% for each property by default)

    Args:
        mw: Molar mass in g/mol
        density: Density in g/mL
        bp: Boiling point in K
        tolerance_pct: Tolerance percentage for matching (default 2%)
        ref_table: Solvent reference table (optional, loads if None)

    Returns:
        Solvent name if match found, None otherwise
    """
    # Handle missing data (all zeros)
    if mw == 0 and density == 0 and bp == 0:
        return None

    # Load reference table if not provided
    if ref_table is None:
        ref_table = load_solvent_reference_table()

    # Match each solvent in reference table
    tolerance = tolerance_pct / 100.0

    for idx, row in ref_table.iterrows():
        ref_mw = row['molar_mass']
        ref_density = row['density']
        ref_bp = row['boiling_point']

        # Check if all properties match within tolerance
        mw_match = abs(mw - ref_mw) / ref_mw <= tolerance if ref_mw > 0 else False
        density_match = abs(density - ref_density) / ref_density <= tolerance if ref_density > 0 else False
        bp_match = abs(bp - ref_bp) / ref_bp <= tolerance if ref_bp > 0 else False

        if mw_match and density_match and bp_match:
            return row['name']

    # No match found - return descriptive unknown
    return f"Unknown_MW{mw:.1f}_D{density:.2f}_BP{bp:.0f}"


def format_mixed_solvent(name_a: str, name_b: str,
                        ratio_a: float, ratio_b: float) -> str:
    """
    Format mixed solvent name with ratio notation

    Args:
        name_a: Primary solvent name
        name_b: Secondary solvent name
        ratio_a: Volume ratio of solvent a
        ratio_b: Volume ratio of solvent b

    Returns:
        Formatted solvent name (e.g., "n-dodecane/toluene (7:3)")
    """
    # Simplify ratio using GCD
    try:
        gcd = math.gcd(int(ratio_a), int(ratio_b))
        simplified_a = int(ratio_a / gcd)
        simplified_b = int(ratio_b / gcd)

        return f"{name_a}/{name_b} ({simplified_a}:{simplified_b})"
    except (ValueError, ZeroDivisionError):
        # Fallback to original ratios if GCD fails
        return f"{name_a}/{name_b} ({ratio_a}:{ratio_b})"


def get_solvent_name(row: pd.Series, ref_table: pd.DataFrame = None) -> Tuple[str, Optional[str]]:
    """
    Extract solvent name(s) from Liu experimental data row

    Handles:
    - Single solvents (ratio_b == 0)
    - Mixed solvents (ratio_a > 0 and ratio_b > 0)
    - Missing data (all zeros)

    Args:
        row: DataFrame row from Liu experimental data
        ref_table: Solvent reference table (optional, loads if None)

    Returns:
        Tuple of (solvent_name, solvent_ratio)
        - solvent_name: Formatted solvent name string
        - solvent_ratio: Ratio string for mixed solvents (e.g., "7:3") or None
    """
    # Load reference table if not provided
    if ref_table is None:
        ref_table = load_solvent_reference_table()

    # Extract volume ratios
    ratio_a = row.get('volume ratio of solvent a', 0)
    ratio_b = row.get('volume ratio of solvent b', 0)

    # Extract properties for solvent a
    mw_a = row.get('Molar mass(a)', 0)
    density_a = row.get('density(a)', 0)
    bp_a = row.get('boiling point(a)', 0)

    # Extract properties for solvent b
    mw_b = row.get('Molar mass(b)', 0)
    density_b = row.get('density(b)', 0)
    bp_b = row.get('boiling point(b)', 0)

    # Match solvent a
    name_a = match_solvent_by_properties(mw_a, density_a, bp_a, ref_table=ref_table)

    # Handle missing data for solvent a
    if name_a is None:
        return 'Unknown', None

    # Single solvent case
    if ratio_b == 0 or (mw_b == 0 and density_b == 0 and bp_b == 0):
        return name_a, None

    # Mixed solvent case
    name_b = match_solvent_by_properties(mw_b, density_b, bp_b, ref_table=ref_table)

    if name_b is None:
        # Fallback if solvent b can't be matched
        return name_a, None

    # Format mixed solvent with ratio
    mixed_name = format_mixed_solvent(name_a, name_b, ratio_a, ratio_b)

    # Calculate simplified ratio for separate field
    try:
        gcd = math.gcd(int(ratio_a), int(ratio_b))
        simplified_a = int(ratio_a / gcd)
        simplified_b = int(ratio_b / gcd)
        ratio_str = f"{simplified_a}:{simplified_b}"
    except (ValueError, ZeroDivisionError):
        ratio_str = f"{ratio_a}:{ratio_b}"

    return mixed_name, ratio_str


# Testing/debugging functions
def test_matcher():
    """Test solvent matching with known values"""
    ref_table = load_solvent_reference_table()

    print("Testing solvent matcher:")
    print("="*70)

    # Test exact matches
    test_cases = [
        (170.34, 0.7495, 489.0, "n-dodecane"),
        (92.141, 0.87, 384.0, "toluene"),
        (119.37, 1.489, 334.3, "chloroform"),
        (266.32, 0.9727, 562.0, "TBP"),
    ]

    for mw, density, bp, expected in test_cases:
        result = match_solvent_by_properties(mw, density, bp, ref_table=ref_table)
        status = "✓" if result == expected else "✗"
        print(f"{status} MW={mw}, D={density}, BP={bp} → {result} (expected: {expected})")

    # Test with tolerance
    print()
    print("Testing with ±2% tolerance:")
    result = match_solvent_by_properties(171.0, 0.75, 490.0, ref_table=ref_table)
    print(f"  MW=171.0, D=0.75, BP=490.0 → {result} (expected: n-dodecane)")

    # Test no match
    print()
    print("Testing unknown solvent:")
    result = match_solvent_by_properties(999.0, 2.5, 1000.0, ref_table=ref_table)
    print(f"  MW=999.0, D=2.5, BP=1000.0 → {result}")

    # Test mixed solvent formatting
    print()
    print("Testing mixed solvent formatting:")
    print(f"  7:3 → {format_mixed_solvent('n-dodecane', 'toluene', 7, 3)}")
    print(f"  5:5 → {format_mixed_solvent('n-dodecane', 'toluene', 5, 5)}")
    print(f"  10:5 → {format_mixed_solvent('kerosene', 'benzene', 10, 5)}")


if __name__ == '__main__':
    test_matcher()
