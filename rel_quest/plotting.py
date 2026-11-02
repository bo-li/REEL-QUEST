"""
REEL-QUEST Plotting Module
Pure-Plotly visualization functions extracted from app.py.
No Streamlit dependencies — returns plotly.graph_objects.Figure objects.
"""

import json
import math
import os
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
import plotly.graph_objects as go
from scipy import stats

from .core import LigandDatabase
from .schema import TREND_DISPLAY_ORDER
from .uncertainty import is_censored, is_mixture

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

LANTHANIDES = [
    'La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd',
    'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu',
]

ATOMIC_NUMBERS = {
    'La': 57, 'Ce': 58, 'Pr': 59, 'Nd': 60, 'Pm': 61, 'Sm': 62, 'Eu': 63,
    'Gd': 64, 'Tb': 65, 'Dy': 66, 'Ho': 67, 'Er': 68, 'Tm': 69, 'Yb': 70,
    'Lu': 71,
}

EXTRACTANT_TYPE_COLORS = {
    # Community-chemistry "Type" labels (the labels surfaced in the UI).
    # Legacy schema-level Extractant_Type keys are kept below for any code path
    # that still passes the raw schema value through this palette.
    'DGA': '#2E86AB',                          # Blue (diglycolamides — largest cluster)
    'TODGA': '#1F6F8B',                        # Deeper blue (canonical DGA reference)
    'TEDGA': '#3A9BC1',                        # Light blue
    'TPDGA': '#5BB3D6',                        # Sky blue
    'TPrDGA': '#75C3DD',                       # Pale blue
    'TDDGA': '#0F4C75',                        # Navy
    'Me-TODGA': '#4A90A4',                     # Slate blue
    'Me-TDDGA': '#3A7080',                     # Dark slate
    'PyranDGA': '#8DB4D2',                     # Periwinkle (oxygen-bridged DGA)
    'Tripodal-DGA': '#A8C8E0',                 # Pale slate
    'DOODA': '#6B8FA8',                        # Steel blue (di-oxa-octanediamide)
    'Bis-DGA': '#476078',                      # Dark steel
    'DAPhen': '#5C4B8A',                       # Deep violet (phenanthroline diamides)
    'PTA': '#7A6BA8',                          # Light violet (PTA/PTdA series)
    'PhenAm': '#9E8FC0',                       # Pale violet
    'BLPhen': '#3F2E6E',                       # Dark violet (bis-lactam phen)
    'POPPhen': '#B07D62',                      # Terracotta/rust (phosphorylated phen)
    'BPPhen': '#C9A07D',                       # Tan (bis-phosphine-oxide phen)
    'MPPhen': '#DDB59A',                       # Pale tan (mono-)
    'BTBP': '#E76F51',                         # Coral (bis-triazinyl-bipyridine)
    'BTPhen': '#EF8068',                       # Salmon
    'BTP': '#F4A582',                          # Light coral (bis-triazinyl-pyridine)
    'BPyDA': '#E9C46A',                        # Amber/gold (bipyridyl diamides)
    'DPA': '#2A9D8F',                          # Teal (dipicolinamides)
    'PyDA': '#3CB5A6',                         # Light teal
    'CMPO': '#F18F01',                         # Orange (carbamoyl-methyl phosphine oxide)
    'NOPOPO': '#D67C00',                       # Dark orange
    'TBP': '#FFB347',                          # Pale orange (tributyl phosphate)
    'HDHP': '#E07B39',                         # Burnt orange (dihexyl phosphoric acid)
    'HEH[EHP]': '#C66A28',                     # Sienna
    'HONTA': '#A85416',                        # Rust
    'Organophosphorus': '#F18F01',             # Orange (fallback)
    'Malonamide': '#C73E1D',                   # Red
    'Amide': '#6A994E',                        # Green
    'DHOA': '#7AAB5E',                         # Light green
    '18-Crown-6': '#BC4749',                   # Crimson (crown ethers)
    'Crown-ether': '#BC4749',                  # Crimson
    'HEDTA': '#386641',                        # Dark green (chelators)
    'NTAamide': '#4F7F4E',                     # Olive green
    'TBSA': '#8AA17C',                         # Sage
    'Aromatic-N-donor': '#A23B72',             # Purple (fallback for aromatic N-donors)
    'Carboxylic-acid': '#386641',              # Dark green
    'Other': '#CCCCCC',                        # Light gray
    'Unknown': '#999999',                      # Dark gray

    # --- Backward-compatible keys for the legacy schema Extractant_Type values ---
    'Diglycolamide': '#2E86AB',
    'Aromatic_Heterocycle': '#A23B72',
    'Phenanthroline_Amide': '#5C4B8A',
    'Phenanthroline_Phosphonate': '#B07D62',
    'Pyridine_Bisamide': '#2A9D8F',
    'Bipyridyl_Diamide': '#E9C46A',
    'Crown_Ether': '#BC4749',
    'Carboxylic_Acid': '#386641',
}

# v13 single-layer palette: 21 Subtype colors, tinted under their donor-chemistry
# group so the pie's slice colors still encode O / N+O / N / P at a glance.
# Records carry the 21-label Subtype directly in both Chemical_Family and
# Subtype fields (donor groups are vocab-side documentation only).
CHEMICAL_FAMILY_COLORS = {
    # --- O-donor subtypes ---
    'DGA':                        '#2E86AB',
    'Tripodal-DGA':               '#4A9BC0',
    'Polyether_Diamide':          '#1F6F8B',
    'Amide':                      '#8DB4D2',
    'Phthalamic_acid':            '#A8C5DC',
    'Malonamide':                 '#5BB3D6',
    'Butanediamide':              '#7AC4E0',
    'Crown_ether':                '#6A994E',
    'Polyphenol':                 '#A7C957',
    # --- N,O-hybrid subtypes ---
    'Nitrilotriacetamide':        '#A23B72',
    'EDTA':                       '#B85A8E',
    'DAPhen':                     '#5C4B8A',
    'MAPhen':                     '#7867A4',
    'BLPhen':                     '#9785BD',
    'Pyridinediamide':            '#9E8FC0',
    'Bipydiamide':                '#B5A8D0',
    # --- N-donor subtypes ---
    'Multident_N_donor':          '#E9C46A',
    # --- P-donor subtypes ---
    'Phosphine_oxide':            '#C73E1D',
    'Phosphonate':                '#E76F51',
    'Phosphate':                  '#F18F01',
    'Phosphoric_acid':            '#9B2D11',
    # --- Fallbacks ---
    'Other':                      '#CCCCCC',
    'Unknown':                    '#999999',
}

SOLVENT_COLORS = {
    'n-dodecane': '#1f77b4',
    'toluene': '#ff7f0e',
    'chloroform': '#2ca02c',
    'kerosene': '#d62728',
    '1-octanol': '#9467bd',
    'Dichloromethane': '#8c564b',
    'meta-nitrobenzotrifluoride': '#e377c2',
    'benzene': '#7f7f7f',
    'nitrobenzene': '#bcbd22',
    'Unknown': '#17becf',
}

SOLVENT_CATEGORY_COLORS = {
    'aliphatic': '#1f77b4',
    'aromatic': '#ff7f0e',
    'chlorinated_hydrocarbon': '#2ca02c',
    'alcohol': '#9467bd',
    'petroleum_mixture': '#d62728',
    'organophosphorus': '#8c564b',
    'ketone': '#e377c2',
    'aromatic_nitro': '#bcbd22',
    'aromatic_fluoro': '#17becf',
    'aromatic_sulfone': '#7f7f7f',
    'Other/Unknown': '#CCCCCC',
    'Unknown': '#999999',
}

LIGAND_ACIDITY_COLORS = {
    'Acidic': '#d62728',
    'Neutral': '#1f77b4',
    'Unknown': '#999999',
}

SELECTIVITY_TREND_COLORS = {
    'Ascending': '#2E86AB',
    'Peaked': '#E76F51',
    'Flat': '#A23B72',
    'U_Shaped': '#457B9D',
    'Reverse': '#F18F01',
    'Mixed': '#C73E1D',
    'Insufficient_Data': '#CCCCCC',
    'Unknown': '#999999',
}

LANTHANIDE_ELEMENT_COLORS = {
    'La': '#E3F2FD', 'Ce': '#BBDEFB', 'Pr': '#90CAF9', 'Nd': '#64B5F6',
    'Pm': '#B9F6CA', 'Sm': '#69F0AE', 'Eu': '#00E676',
    'Gd': '#C8E6C9', 'Tb': '#81C784',
    'Dy': '#FFE082', 'Ho': '#FFD54F', 'Er': '#FFCA28',
    'Tm': '#FFB74D', 'Yb': '#FF9800', 'Lu': '#F57C00',
}

# Default font sizes (non-Streamlit context)
FONT_SIZES = {
    'title': 21,
    'axis_title': 18,
    'tick': 16,
    'hover': 16,
    'legend': 14,
    'annotation': 16,
}

# Aliases used by the dashboard's fullscreen-font toggle
FONT_SIZES_NORMAL = FONT_SIZES

FONT_SIZES_FULLSCREEN = {
    'title': 42,
    'axis_title': 36,
    'tick': 32,
    'hover': 32,
    'legend': 28,
    'annotation': 28,
}


# ---------------------------------------------------------------------------
# Solvent category helpers (standalone, no Streamlit cache)
# ---------------------------------------------------------------------------

_SOLVENT_CATEGORY_CACHE: Optional[Tuple[Dict, Dict]] = None


def load_solvent_categories() -> Tuple[Dict[str, str], Dict]:
    """Load solvent→category mapping from vocab/solvents.json."""
    global _SOLVENT_CATEGORY_CACHE
    if _SOLVENT_CATEGORY_CACHE is not None:
        return _SOLVENT_CATEGORY_CACHE

    vocab_path = os.path.join(
        os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
        'vocab', 'solvents.json',
    )
    with open(vocab_path, 'r') as f:
        solvent_data = json.load(f)

    solvent_to_category: Dict[str, str] = {}
    for info in solvent_data['solvents']:
        solvent_to_category[info['name']] = info['category']
        for alt in info.get('alternative_names', []):
            solvent_to_category[alt] = info['category']

    category_descriptions = solvent_data.get('solvent_categories', {})
    _SOLVENT_CATEGORY_CACHE = (solvent_to_category, category_descriptions)
    return _SOLVENT_CATEGORY_CACHE


def get_solvent_category(solvent_name: str, solvent_to_category: Dict[str, str]) -> str:
    """Map a solvent name to its category."""
    if not solvent_name or solvent_name == 'Unknown':
        return 'Other/Unknown'
    name_lower = solvent_name.lower()
    for known, cat in solvent_to_category.items():
        if known.lower() == name_lower:
            return cat
    if '/' in solvent_name:
        primary = solvent_name.split('/')[0].strip()
        for known, cat in solvent_to_category.items():
            if known.lower() == primary.lower():
                return cat
    if ' with ' in solvent_name.lower():
        primary = solvent_name.lower().split(' with ')[0].strip()
        for known, cat in solvent_to_category.items():
            if known.lower() == primary:
                return cat
    for known, cat in solvent_to_category.items():
        if known.lower() in name_lower:
            return cat
    return 'Other/Unknown'


def classify_ligand_acidity(ligand_name: str, extractant_type: str, smiles: Optional[str] = None) -> str:
    """Classify a ligand as Acidic, Neutral, or Unknown."""
    # Try SMILES-based classification
    if smiles and smiles != '-':
        try:
            from rdkit import Chem
            mol = Chem.MolFromSmiles(smiles)
            if mol:
                acidic_patterns = [
                    '[CX3](=O)[OX2H1]',
                    '[P;X4](=O)[O;H1]',
                    '[S;X4](=O)(=O)[O;H1]',
                ]
                for pat in acidic_patterns:
                    if mol.HasSubstructMatch(Chem.MolFromSmarts(pat)):
                        return 'Acidic'
                has_acidic_oh = any(
                    mol.HasSubstructMatch(Chem.MolFromSmarts(p))
                    for p in acidic_patterns
                )
                if not has_acidic_oh:
                    neutral_patterns = [
                        '[CX3](=O)[NX3]', '[PX4](=O)([CX4,cX3])',
                        '[PX4](=O)(O[CX4])', '[nX2]', '[OX2][CX4][OX2]',
                    ]
                    for pat in neutral_patterns:
                        if mol.HasSubstructMatch(Chem.MolFromSmarts(pat)):
                            return 'Neutral'
        except Exception:
            pass

    if not extractant_type:
        return 'Unknown'
    if extractant_type == 'Carboxylic_Acid':
        return 'Acidic'
    if extractant_type in [
        'Diglycolamide', 'Malonamide', 'Amide', 'Aromatic_Heterocycle',
        'Crown_Ether', 'Phenanthroline_Amide', 'Phenanthroline_Phosphonate',
        'Pyridine_Bisamide', 'Bipyridyl_Diamide',
    ]:
        return 'Neutral'
    if extractant_type == 'Organophosphorus':
        name_lower = (ligand_name or '').lower()
        acidic_kw = ['dehpa', 'd2ehpa', 'hdehp', 'h2mehp', 'cyanex',
                      'phosphoric acid', 'phosphonic acid',
                      'hehehp', 'h2dep', 'h2mep', 'heh[ehp]', 'hdhp']
        if any(k in name_lower for k in acidic_kw):
            return 'Acidic'
        neutral_kw = ['topo', 'cmpo', 'tbpo', 'tbp', 'dmdbtdma',
                       'tppo', 'oxide', 'phen', 'phosphine']
        if any(k in name_lower for k in neutral_kw):
            return 'Neutral'
        return 'Unknown'
    return 'Unknown'


# ---------------------------------------------------------------------------
# Data extraction helpers
# ---------------------------------------------------------------------------

def _safe_float(value, default=0.0):
    """Convert value to float, handling strings with units like '3.0 M'."""
    if value in (None, '', 'N/A'):
        return default
    try:
        return float(value)
    except (ValueError, TypeError):
        try:
            return float(str(value).split()[0])
        except (ValueError, IndexError):
            return default


# --- bounds (censored cells: D > x, D < x) -----------------------------------
#
# A bound is a measurement of its element but not a value. The frames below
# keep it, marked in a 'censor' column; every statistic (box, violin, median,
# fit) is drawn from the values alone, and bounds are drawn as open triangles
# pointing the way the true value lies.

BOUND_SYMBOLS = {'>': 'triangle-up-open', '<': 'triangle-down-open'}
BOUND_LEGEND = 'Bound (▲ greater than, ▼ less than)'
# Near-black, so a bound stands out over a box or violin drawn in its category's colour
BOUND_COLOR = '#222222'


def separation_factor_censor(censor1, censor2) -> Optional[str]:
    """What logSF = log D1 - log D2 is, given each D's censor.

    '=' a value; '>' or '<' a bound (D1 '>' or D2 '<' make it a lower bound,
    D1 '<' or D2 '>' an upper one); None when two bounds point the same way,
    which says nothing about the ratio.
    """
    ways = set()
    if censor1 in ('>', '<'):
        ways.add(censor1)
    if censor2 in ('>', '<'):
        ways.add('<' if censor2 == '>' else '>')
    if not ways:
        return '='
    return ways.pop() if len(ways) == 1 else None


def split_bounds(df: pd.DataFrame) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """(values, bounds) of a frame from extract_element_data / extract_separation_factor_data."""
    if df.empty or 'censor' not in df.columns:
        return df, df.iloc[0:0]
    marked = df['censor'].notna()
    return df[~marked], df[marked]


def category_order(values: pd.DataFrame, bounds: pd.DataFrame, column: str, y: str) -> List:
    """Categories by median of their values, descending; categories holding only bounds last."""
    order = values.groupby(column)[y].median().sort_values(ascending=False).index.tolist()
    return order + sorted(c for c in bounds[column].unique() if c not in order) if not bounds.empty else order


def bound_text(label: str, censor: str, log_value: float) -> str:
    """'D[Eu] > 100 (log 2.00)': how a bound reads in a hover, never as a value."""
    return f"{label} {censor} {10 ** log_value:.3g} (log {log_value:.2f})"


def bound_trace(x, y, censors, hovertexts, color, name=BOUND_LEGEND, showlegend=False) -> go.Scatter:
    """Open triangles at the bounds: up for '>', down for '<'."""
    return go.Scatter(
        x=list(x), y=list(y), mode='markers', name=name,
        marker=dict(symbol=[BOUND_SYMBOLS[c] for c in censors], size=11, color=color,
                    line=dict(width=2, color=color)),
        hovertext=list(hovertexts), hoverinfo='text', showlegend=showlegend,
    )


def extract_element_data(db: LigandDatabase, element: str, all_ligands=None) -> pd.DataFrame:
    """Extract all logD measurements for a specific element.

    Returns a DataFrame with one row per (ligand, condition) measurement,
    including v13 classification columns (community_name, chemical_family,
    subtype) and the ligand-acidity classification. A bound keeps its row,
    with 'censor' '>' or '<' and the bound in 'log_d' (see split_bounds).
    """
    if all_ligands is None:
        all_ligands = db.get_all_ligands()
    data_points = []
    for ligand in all_ligands:
        # tested_elements, not element_coverage: coverage lists only the
        # La-Lu points the trend fit used, so Y and Sc never appear in it.
        if element not in ligand.tested_elements:
            continue
        for cond in ligand.experimental_conditions:
            if is_mixture(cond):   # compares ligands: a mixture's D is not the ligand's alone
                continue
            dist_ratios = cond.get('distribution_ratios', {})
            if element in dist_ratios and dist_ratios[element].get('log_value') is not None:
                acidity = classify_ligand_acidity(ligand.name, ligand.extractant_type, ligand.smiles)
                data_points.append({
                    'ligand_name': ligand.name,
                    'log_d': dist_ratios[element]['log_value'],
                    'd_value': dist_ratios[element]['value'],
                    'acid_conc': _safe_float(cond.get('Acid_Concentration_M'), 0.0),
                    'extractant_conc': _safe_float(cond.get('Extractant_Concentration_M'), 0.0),
                    'temperature': _safe_float(cond.get('Temperature_C'), 25.0),
                    'solvent': cond.get('Solvent_Name', 'Unknown'),
                    'acid_name': cond.get('Acid_Name', 'Unknown'),
                    'mechanism': ligand.mechanism,
                    'structural_class': ligand.structural_class or 'Unknown',
                    'extractant_type': ligand.extractant_type or 'Unknown',
                    'community_name': getattr(ligand, 'community_name', None) or 'Other',
                    'chemical_family': getattr(ligand, 'chemical_family', None) or ligand.extractant_type or 'Unknown',
                    'subtype': getattr(ligand, 'subtype', None) or getattr(ligand, 'community_name', None) or 'Other',
                    'acidity': acidity,
                    'condition_id': cond.get('condition_id', ''),
                    'censor': dist_ratios[element]['censor'] if is_censored(dist_ratios[element]) else None,
                })
    return pd.DataFrame(data_points)


def _sf_safe_float(value):
    """Convert value to float for SF extraction; None when unparseable.

    Unlike _safe_float, this does NOT substitute defaults — a missing acid
    concentration must stay missing (NaN) rather than become a fabricated 0.0.
    """
    if isinstance(value, (int, float)):
        return float(value)
    if isinstance(value, str):
        try:
            numeric_part = ''.join(c for c in value.split()[0] if c.isdigit() or c in '.-')
            return float(numeric_part) if numeric_part else None
        except (ValueError, IndexError):
            return None
    return None


def extract_separation_factor_data(
    db: LigandDatabase, element1: str, element2: str, all_ligands=None,
) -> pd.DataFrame:
    """Extract logSF = logD[elem1] - logD[elem2] for paired measurements.

    A pair with a bound is a bound on logSF ('censor', separation_factor_censor);
    a pair whose two bounds point the same way says nothing and is left out,
    counted in the frame's attrs['no_information'].
    """
    if all_ligands is None:
        all_ligands = db.get_all_ligands()
    data_points, no_information = [], 0
    for ligand in all_ligands:
        # tested_elements, not element_coverage — see the note above.
        tested = ligand.tested_elements
        if element1 not in tested or element2 not in tested:
            continue
        for cond in ligand.experimental_conditions:
            if is_mixture(cond):
                continue
            dr = cond.get('distribution_ratios', {})
            if element1 not in dr or element2 not in dr:
                continue
            log_d1 = dr[element1].get('log_value')
            log_d2 = dr[element2].get('log_value')
            if log_d1 is None or log_d2 is None:
                continue
            c1 = dr[element1]['censor'] if is_censored(dr[element1]) else None
            c2 = dr[element2]['censor'] if is_censored(dr[element2]) else None
            sf_censor = separation_factor_censor(c1, c2)
            if sf_censor is None:
                no_information += 1
                continue
            acidity = classify_ligand_acidity(ligand.name, ligand.extractant_type, ligand.smiles)
            data_points.append({
                'ligand_name': ligand.name,
                'log_sf': log_d1 - log_d2,
                'log_d_elem1': log_d1,
                'log_d_elem2': log_d2,
                'd_elem1': dr[element1].get('value'),
                'd_elem2': dr[element2].get('value'),
                'acid_conc': _sf_safe_float(cond.get('Acid_Concentration_M', cond.get('Acid_Concentration'))),
                'extractant_conc': _sf_safe_float(cond.get('Extractant_Concentration_M', cond.get('Extractant_Concentration'))),
                'temperature': _sf_safe_float(cond.get('Temperature_C', cond.get('Temperature'))),
                'solvent': cond.get('Solvent_Name', 'Unknown'),
                'acid_name': cond.get('Acid_Name', 'Unknown'),
                'mechanism': ligand.mechanism,
                'structural_class': ligand.structural_class or 'Unknown',
                'extractant_type': ligand.extractant_type or 'Unknown',
                'community_name': getattr(ligand, 'community_name', None) or 'Other',
                'chemical_family': getattr(ligand, 'chemical_family', None) or ligand.extractant_type or 'Unknown',
                'subtype': getattr(ligand, 'subtype', None) or getattr(ligand, 'community_name', None) or 'Other',
                'acidity': acidity,
                'condition_id': cond.get('condition_id', ''),
                'censor': None if sf_censor == '=' else sf_censor,
                'censor_elem1': c1,
                'censor_elem2': c2,
            })
    df = pd.DataFrame(data_points)
    df.attrs['no_information'] = no_information
    return df


# ---------------------------------------------------------------------------
# Plot functions — each returns a go.Figure
# ---------------------------------------------------------------------------

def plot_selectivity_trend(
    db: LigandDatabase,
    ligand_name: str,
    condition_id: Optional[str] = None,
) -> Optional[go.Figure]:
    """
    Plot logD vs atomic number for a ligand with linear fit + CI bands.

    Returns (figure, summary_text) or (None, error_text).
    """
    # Find ligand
    matches = [l for l in db.get_all_ligands() if l.name.lower() == ligand_name.lower()]
    if not matches:
        matches = [l for l in db.get_all_ligands() if ligand_name.lower() in l.name.lower()]
    if not matches:
        return None

    ligand = matches[0]

    if condition_id:
        conditions = [c for c in ligand.experimental_conditions if c['condition_id'] == condition_id]
        if not conditions:
            return None
        d_ratios = conditions[0].get('distribution_ratios', {})
        title_suffix = f"<br><sub>Condition: {condition_id[:60]}...</sub>"
    else:
        d_ratios = ligand.get_aggregated_d_values(aggregation='median', lanthanides_only=True)
        title_suffix = "<br><sub>Aggregated (median across all conditions)</sub>"

    elements, z_vals, d_vals, log_d_vals, sigma_log_d_vals = [], [], [], [], []
    bounds = []   # (Z, log bound, censor, element): drawn, never fitted
    for elem in LANTHANIDES:
        if elem in d_ratios:
            entry = d_ratios[elem]
            if is_censored(entry):
                bounds.append((ATOMIC_NUMBERS[elem], entry['log_value'], entry['censor'], elem))
                continue
            d_val = entry.get('value', 0) if isinstance(entry, dict) else entry
            sigma = entry.get('sigma_log_D', 0.033) if isinstance(entry, dict) else 0.033
            if isinstance(d_val, (int, float)) and d_val > 0:
                elements.append(elem)
                z_vals.append(ATOMIC_NUMBERS[elem])
                d_vals.append(d_val)
                log_d_vals.append(np.log10(d_val))
                sigma_log_d_vals.append(sigma if sigma else 0.033)

    if len(z_vals) < 2:
        return None

    slope, intercept, r_value, p_value, std_err = stats.linregress(z_vals, log_d_vals)
    r_squared = r_value ** 2

    sel = getattr(ligand, 'selectivity', None) or {}
    stored_trend = sel.get('trend_classification')
    if stored_trend in ('Peaked', 'U_Shaped', 'Ascending', 'Reverse', 'Flat', 'Mixed'):
        trend = stored_trend
    elif r_squared < 0.6:
        trend = 'Flat/Mixed'
    elif slope > 0.02:
        trend = 'Ascending'
    elif slope < -0.02:
        trend = 'Reverse'
    else:
        trend = 'Flat'

    trend_colors = {
        'Ascending': 'green', 'Reverse': 'red',
        'Peaked': '#E76F51', 'U_Shaped': '#457B9D',
    }
    marker_color = trend_colors.get(trend, 'gray')

    fig = go.Figure()
    z_range = np.linspace(min(z_vals) - 0.5, max(z_vals) + 0.5, 60)

    # Trend line and CI band both pivot on the centroid of the plotted
    # points. The stored slope is the mean of per-condition fits and the
    # stored CI is its bootstrap interval over conditions
    # (rel_quest.uncertainty.bootstrap_trend_ci), so the aggregated view
    # draws that slope and that band; a single-condition view shows only
    # its own OLS fit, because the stored CI does not describe it. (The
    # band was previously drawn as ci * Z + intercept, pivoting at Z = 0,
    # so its height at Z ~ 64 was 64 x the CI width.)
    z_bar = float(np.mean(z_vals))
    y_bar = float(np.mean(log_d_vals))
    stored_slope = sel.get('slope')
    ci_lower = sel.get('slope_CI_lower')
    ci_upper = sel.get('slope_CI_upper')
    aggregated = condition_id is None

    if aggregated and ci_lower is not None and ci_upper is not None:
        y_lower = y_bar + ci_lower * (z_range - z_bar)
        y_upper = y_bar + ci_upper * (z_range - z_bar)
        fig.add_trace(go.Scatter(
            x=np.concatenate([z_range, z_range[::-1]]).tolist(),
            y=np.concatenate([y_upper, y_lower[::-1]]).tolist(),
            fill='toself', fillcolor='rgba(31,119,180,0.15)',
            line=dict(color='rgba(0,0,0,0)'),
            name='95% CI of trend slope (bootstrap over conditions)',
            showlegend=True, hoverinfo='skip',
        ))

    if aggregated and stored_slope is not None:
        line_slope, line_name = stored_slope, 'Trend slope (mean of per-condition fits)'
    else:
        line_slope, line_name = slope, 'Linear fit'  # OLS through the centroid
    fig.add_trace(go.Scatter(
        x=z_range.tolist(),
        y=(y_bar + line_slope * (z_range - z_bar)).tolist(),
        mode='lines', line=dict(color=marker_color, width=2, dash='dash'),
        name=line_name, hoverinfo='skip',
    ))

    peak_elem = sel.get('peak_element')
    if trend in ('Peaked', 'U_Shaped') and peak_elem and peak_elem in elements:
        peak_idx = elements.index(peak_elem)
        label = 'Peak' if trend == 'Peaked' else 'Trough'
        fig.add_trace(go.Scatter(
            x=[z_vals[peak_idx]], y=[log_d_vals[peak_idx]],
            mode='markers+text',
            marker=dict(size=16, color=marker_color, symbol='diamond',
                        line=dict(width=2, color='white')),
            text=[peak_elem], textposition='top center',
            name=f'{label} ({peak_elem})',
            hovertemplate=f'<b>{label}: {peak_elem}</b><br>log(D) = {log_d_vals[peak_idx]:.3f}<extra></extra>',
        ))

    fig.add_trace(go.Scatter(
        x=z_vals, y=log_d_vals, mode='markers',
        marker=dict(size=12, color=marker_color, line=dict(width=1, color='white')),
        error_y=dict(type='data', array=sigma_log_d_vals, visible=True,
                     color=marker_color, thickness=1.5, width=4),
        text=elements, customdata=d_vals,
        hovertemplate='<b>%{text}</b><br>Atomic #: %{x}<br>D = %{customdata:.3f}<br>log(D) = %{y:.3f}<extra></extra>',
        name='Data',
    ))
    if bounds:
        fig.add_trace(bound_trace([b[0] for b in bounds], [b[1] for b in bounds], [b[2] for b in bounds],
                                  [f"<b>{b[3]}</b><br>{bound_text('D', b[2], b[1])}<br>not in the fit"
                                   for b in bounds], marker_color, showlegend=True))

    fig.update_layout(
        title=f"Selectivity Trend: {ligand.name} ({trend}){title_suffix}",
        xaxis_title="Atomic Number", yaxis_title="log(D)",
        xaxis=dict(dtick=1, range=[56.5, 71.5]),
        hovermode='closest', showlegend=True, height=450,
        template="plotly_white",
    )

    return fig


def get_selectivity_trend_summary(
    db: LigandDatabase, ligand_name: str,
) -> str:
    """Return a text summary of the selectivity trend for a ligand."""
    matches = [l for l in db.get_all_ligands() if l.name.lower() == ligand_name.lower()]
    if not matches:
        matches = [l for l in db.get_all_ligands() if ligand_name.lower() in l.name.lower()]
    if not matches:
        return f"Ligand '{ligand_name}' not found in the database."

    ligand = matches[0]
    agg = ligand.get_aggregated_d_values(aggregation='median', lanthanides_only=True)

    z_vals, log_d_vals = [], []
    for elem in LANTHANIDES:
        if elem in agg:
            d_val = agg[elem].get('value', 0) if isinstance(agg[elem], dict) else agg[elem]
            if isinstance(d_val, (int, float)) and d_val > 0:
                z_vals.append(ATOMIC_NUMBERS[elem])
                log_d_vals.append(np.log10(d_val))

    sel = getattr(ligand, 'selectivity', None) or {}
    trend = sel.get('trend_classification', 'Unknown')

    if len(z_vals) >= 2:
        slope, _, r_value, _, _ = stats.linregress(z_vals, log_d_vals)
        r_squared = r_value ** 2
        logd_min = min(log_d_vals)
        logd_max = max(log_d_vals)
        return (
            f"{ligand.name} shows a {trend} trend (slope={slope:.4f}, "
            f"R²={r_squared:.3f}). logD ranges from {logd_min:.2f} to "
            f"{logd_max:.2f} across {len(z_vals)} measured lanthanides."
        )
    else:
        return f"{ligand.name}: insufficient data for trend analysis ({len(z_vals)} elements measured)."


def plot_element_distribution(
    db: LigandDatabase,
    element: str,
    group_by: str = "ligand_type",
) -> Optional[go.Figure]:
    """
    Box+violin plot of logD for one element, grouped by category.

    group_by: "ligand_type" | "solvent" | "solvent_category" | "ligand_acidity"
    """
    df = extract_element_data(db, element)
    if df.empty:
        return None

    # Determine grouping column and color map
    group_map = {
        "ligand_type": ("extractant_type", "Ligand Type", EXTRACTANT_TYPE_COLORS),
        "solvent": ("solvent", "Solvent", SOLVENT_COLORS),
        "solvent_category": (None, "Solvent Category", SOLVENT_CATEGORY_COLORS),
        "ligand_acidity": ("acidity", "Ligand Acidity", LIGAND_ACIDITY_COLORS),
    }

    if group_by not in group_map:
        group_by = "ligand_type"

    x_column, x_label, color_map = group_map[group_by]

    plot_df = df.copy()
    if group_by == "solvent_category":
        solvent_to_category, _ = load_solvent_categories()
        plot_df['category'] = plot_df['solvent'].apply(
            lambda s: get_solvent_category(s, solvent_to_category)
        )
        x_column = 'category'

    # Statistics from values only; bounds drawn on top as triangles
    plot_df, bound_df = split_bounds(plot_df)
    categories_sorted = category_order(plot_df, bound_df, x_column, 'log_d')

    fig = go.Figure()
    for category in categories_sorted:
        cat_bounds = bound_df[bound_df[x_column] == category]
        cat_color = color_map.get(category, '#999999')
        if not cat_bounds.empty:
            fig.add_trace(bound_trace(
                [category] * len(cat_bounds), cat_bounds['log_d'], cat_bounds['censor'],
                [f"<b>{r.ligand_name}</b><br>{bound_text(f'D[{element}]', r.censor, r.log_d)}"
                 for r in cat_bounds.itertuples()], BOUND_COLOR))
        cat_data = plot_df[plot_df[x_column] == category]
        if cat_data.empty:
            continue
        n = len(cat_data)
        median_val = cat_data['log_d'].median()

        # Violin (underneath)
        fig.add_trace(go.Violin(
            y=cat_data['log_d'], name=category,
            line_color=cat_color, fillcolor=cat_color, opacity=0.3,
            box_visible=False, meanline_visible=False, points=False,
            scalemode='width', width=0.8, showlegend=False, hoverinfo='skip',
        ))
        # Box (on top)
        fig.add_trace(go.Box(
            y=cat_data['log_d'], name=category,
            marker=dict(color=cat_color, line=dict(width=1.5, color='white')),
            line=dict(color=cat_color, width=2), boxmean=True, boxpoints='outliers',
            hovertemplate=(
                f"<b>{category}</b><br>"
                f"Median: {median_val:.2f}<br>"
                f"Count: {n}<extra></extra>"
            ),
            showlegend=False, fillcolor='rgba(255,255,255,0.8)',
        ))

    fig.update_layout(
        title=dict(text=f"Distribution Coefficient for {element}",
                   font=dict(size=FONT_SIZES['title'])),
        xaxis=dict(title=dict(text=x_label, font=dict(size=FONT_SIZES['axis_title'])),
                   tickfont=dict(size=FONT_SIZES['tick']), tickangle=-45,
                   categoryorder='array', categoryarray=categories_sorted),
        yaxis=dict(title=dict(text=f"log(D[{element}])", font=dict(size=FONT_SIZES['axis_title'])),
                   tickfont=dict(size=FONT_SIZES['tick'])),
        template="plotly_white", height=500, hovermode='closest',
        plot_bgcolor='rgba(250,250,250,0.5)',
    )
    return fig


def plot_separation_factor(
    db: LigandDatabase,
    element1: str,
    element2: str,
    group_by: str = "ligand_type",
) -> Optional[go.Figure]:
    """
    Box+violin plot of logSF[elem1/elem2], grouped by category.

    group_by: "ligand_type" | "solvent" | "solvent_category" | "ligand_acidity"
    """
    df = extract_separation_factor_data(db, element1, element2)
    if df.empty:
        return None

    group_map = {
        "ligand_type": ("extractant_type", "Ligand Type", EXTRACTANT_TYPE_COLORS),
        "solvent": ("solvent", "Solvent", SOLVENT_COLORS),
        "solvent_category": (None, "Solvent Category", SOLVENT_CATEGORY_COLORS),
        "ligand_acidity": ("acidity", "Ligand Acidity", LIGAND_ACIDITY_COLORS),
    }
    if group_by not in group_map:
        group_by = "ligand_type"

    x_column, x_label, color_map = group_map[group_by]

    plot_df = df.copy()
    if group_by == "solvent_category":
        solvent_to_category, _ = load_solvent_categories()
        plot_df['category'] = plot_df['solvent'].apply(
            lambda s: get_solvent_category(s, solvent_to_category)
        )
        x_column = 'category'

    # Statistics from values only; bounds drawn on top as triangles
    plot_df, bound_df = split_bounds(plot_df)
    categories_sorted = category_order(plot_df, bound_df, x_column, 'log_sf')

    fig = go.Figure()
    for category in categories_sorted:
        cat_bounds = bound_df[bound_df[x_column] == category]
        cat_color = color_map.get(category, '#999999')
        if not cat_bounds.empty:
            fig.add_trace(bound_trace(
                [category] * len(cat_bounds), cat_bounds['log_sf'], cat_bounds['censor'],
                [f"<b>{r.ligand_name}</b><br>{bound_text(f'SF[{element1}/{element2}]', r.censor, r.log_sf)}"
                 for r in cat_bounds.itertuples()], BOUND_COLOR))
        cat_data = plot_df[plot_df[x_column] == category]
        if cat_data.empty:
            continue
        n = len(cat_data)
        median_val = cat_data['log_sf'].median()

        fig.add_trace(go.Violin(
            y=cat_data['log_sf'], name=category,
            line_color=cat_color, fillcolor=cat_color, opacity=0.3,
            box_visible=False, meanline_visible=False, points=False,
            scalemode='width', width=0.8, showlegend=False, hoverinfo='skip',
        ))
        fig.add_trace(go.Box(
            y=cat_data['log_sf'], name=category,
            marker=dict(color=cat_color, line=dict(width=1.5, color='white')),
            line=dict(color=cat_color, width=2), boxmean=True, boxpoints='outliers',
            hovertemplate=(
                f"<b>{category}</b><br>"
                f"Median logSF: {median_val:.2f}<br>"
                f"Count: {n}<extra></extra>"
            ),
            showlegend=False, fillcolor='rgba(255,255,255,0.8)',
        ))

    # Reference line at y=0 (no selectivity)
    fig.add_hline(y=0, line_dash="dash", line_color="gray", opacity=0.5,
                  annotation_text="No selectivity", annotation_position="bottom right")

    fig.update_layout(
        title=dict(text=f"Separation Factor: logSF[{element1}/{element2}]",
                   font=dict(size=FONT_SIZES['title'])),
        xaxis=dict(title=dict(text=x_label, font=dict(size=FONT_SIZES['axis_title'])),
                   tickfont=dict(size=FONT_SIZES['tick']), tickangle=-45,
                   categoryorder='array', categoryarray=categories_sorted),
        yaxis=dict(title=dict(text=f"logSF[{element1}/{element2}]",
                              font=dict(size=FONT_SIZES['axis_title'])),
                   tickfont=dict(size=FONT_SIZES['tick'])),
        template="plotly_white", height=500, hovermode='closest',
        plot_bgcolor='rgba(250,250,250,0.5)',
    )
    return fig


def plot_database_overview(
    db: LigandDatabase,
    chart_type: str = "extractant_types",
) -> Optional[go.Figure]:
    """
    Pie/donut chart for database overview.

    chart_type: "extractant_types" | "selectivity_trends" | "element_coverage"
    """
    ligands = db.get_all_ligands()
    if not ligands:
        return None

    if chart_type == "extractant_types":
        counts: Dict[str, int] = {}
        for lig in ligands:
            t = lig.extractant_type or 'Unknown'
            counts[t] = counts.get(t, 0) + 1
        sorted_items = sorted(counts.items(), key=lambda x: x[1], reverse=True)
        labels = [i[0] for i in sorted_items]
        values = [i[1] for i in sorted_items]
        colors = [EXTRACTANT_TYPE_COLORS.get(l, '#999999') for l in labels]
        title = f"Extractant Type Distribution ({len(ligands)} total)"

    elif chart_type == "selectivity_trends":
        counts = {}
        for lig in ligands:
            t = lig.selectivity_trend or 'Unknown'
            counts[t] = counts.get(t, 0) + 1
        trend_order = list(TREND_DISPLAY_ORDER)
        labels, values = [], []
        for t in trend_order:
            if t in counts:
                labels.append(t)
                values.append(counts[t])
        for t, c in counts.items():
            if t not in trend_order:
                labels.append(t)
                values.append(c)
        colors = [SELECTIVITY_TREND_COLORS.get(l, '#999999') for l in labels]
        title = f"Selectivity Trend Distribution ({len(ligands)} total)"

    elif chart_type == "element_coverage":
        from collections import Counter
        elem_counts = Counter()
        for lig in ligands:
            seen = set()
            for cond in lig.experimental_conditions:
                for elem in LANTHANIDES:
                    if elem in cond.get('distribution_ratios', {}):
                        seen.add(elem)
            for elem in seen:
                elem_counts[elem] += 1
        labels = [e for e in LANTHANIDES if e in elem_counts]
        values = [elem_counts[e] for e in labels]
        colors = [LANTHANIDE_ELEMENT_COLORS.get(e, '#999999') for e in labels]
        title = "Measured Lanthanides Distribution"
    else:
        return None

    fig = go.Figure(data=[go.Pie(
        labels=labels, values=values, hole=0.3,
        hovertemplate='<b>%{label}</b><br>Count: %{value}<br>Percentage: %{percent}<extra></extra>',
        textposition='auto', textinfo='label+percent',
        marker=dict(colors=colors, line=dict(color='white', width=2)),
    )])
    fig.update_layout(
        title=dict(text=title, x=0.5, xanchor='center'),
        height=450, template="plotly_white", showlegend=True,
        legend=dict(orientation="v", yanchor="middle", y=0.5, xanchor="left", x=1.05),
    )
    return fig


def plot_ligand_comparison(
    db: LigandDatabase,
    ligand_names: List[str],
) -> Optional[go.Figure]:
    """
    Overlay selectivity trends for multiple ligands on one plot.
    """
    all_ligands = db.get_all_ligands()

    # Color palette for comparison
    palette = ['#2E86AB', '#E76F51', '#2A9D8F', '#F4A261', '#264653',
               '#E9C46A', '#A23B72', '#6A994E']

    fig = go.Figure()
    found_any = False

    for idx, name in enumerate(ligand_names[:8]):
        name_lower = name.lower()
        matches = [l for l in all_ligands if l.name.lower() == name_lower]
        if not matches:
            matches = [l for l in all_ligands if name_lower in l.name.lower()]
        if not matches:
            continue

        ligand = matches[0]
        agg = ligand.get_aggregated_d_values(aggregation='median', lanthanides_only=True)

        elems, z, logd = [], [], []
        for elem in LANTHANIDES:
            if elem in agg:
                d_val = agg[elem].get('value', 0) if isinstance(agg[elem], dict) else agg[elem]
                if isinstance(d_val, (int, float)) and d_val > 0:
                    elems.append(elem)
                    z.append(ATOMIC_NUMBERS[elem])
                    logd.append(np.log10(d_val))

        if len(z) < 2:
            continue
        found_any = True

        color = palette[idx % len(palette)]
        sel = getattr(ligand, 'selectivity', None) or {}
        trend = sel.get('trend_classification', '?')

        fig.add_trace(go.Scatter(
            x=z, y=logd, mode='lines+markers',
            marker=dict(size=9, color=color, line=dict(width=1, color='white')),
            line=dict(color=color, width=2),
            name=f"{ligand.name} ({trend})",
            text=elems,
            hovertemplate='<b>%{text}</b><br>log(D) = %{y:.3f}<extra>' + ligand.name + '</extra>',
        ))

    if not found_any:
        return None

    fig.update_layout(
        title="Ligand Comparison: logD vs Atomic Number",
        xaxis_title="Atomic Number", yaxis_title="log(D)",
        xaxis=dict(dtick=1, range=[56.5, 71.5]),
        hovermode='closest', showlegend=True, height=450,
        template="plotly_white",
    )
    return fig
