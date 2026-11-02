"""Selectivity-trend classification across the lanthanide series.

`classify_selectivity_trend` labels a record Ascending / Reverse / Flat /
Mixed / Peaked / U_Shaped / Insufficient_Data from its conditions: a linear fit
of log D on atomic number per condition, overridden by
`peak_position_heuristic` when a clear interior extremum exists, then votes
aggregated across conditions. It is the function behind every stored
`selectivity` block (rel_quest.derived.recompute_record).

Moved here from rel_quest/ingest_liu.py (2026-10-21), which re-exports both
names; it was never specific to one source.
"""

from collections import Counter
from typing import Dict, List, Tuple

import numpy as np

from .uncertainty import is_mixture
from .utils import get_lanthanide_order, round_if_not_none


def peak_position_heuristic(points: List[Tuple[int, float]],
                             min_elements: int = 4,
                             drop_threshold: float = 0.3) -> str:
    """
    Detect Peaked or U_Shaped profiles using the position of the max/min log_D.

    A profile is Peaked when:
      - max log_D is at an *interior* element (not the first or last measured)
      - log_D drops by >= drop_threshold on BOTH sides of that maximum

    A profile is U_Shaped when:
      - min log_D is at an interior element
      - log_D rises by >= drop_threshold on BOTH sides of that minimum

    Args:
        points: List of (atomic_number, log_D) pairs, any order.
        min_elements: Minimum distinct elements required (default 4, the same
            minimum the linear fit uses, so a 4-element condition with a clear
            interior maximum is not forced into a linear class).
        drop_threshold: Minimum log_D change required on each side (default 0.3).

    Returns:
        'Peaked', 'U_Shaped', or 'Linear' (= no clear arch/trough detected).
    """
    if len(points) < min_elements:
        return 'Linear'

    # Sort by atomic number
    sorted_pts = sorted(points, key=lambda p: p[0])
    log_d = [p[1] for p in sorted_pts]

    # --- Peaked check ---
    peak_idx = int(np.argmax(log_d))
    if 0 < peak_idx < len(log_d) - 1:
        drop_left  = log_d[peak_idx] - log_d[0]
        drop_right = log_d[peak_idx] - log_d[-1]
        if drop_left >= drop_threshold and drop_right >= drop_threshold:
            return 'Peaked'

    # --- U_Shaped check ---
    trough_idx = int(np.argmin(log_d))
    if 0 < trough_idx < len(log_d) - 1:
        rise_left  = log_d[0]   - log_d[trough_idx]
        rise_right = log_d[-1]  - log_d[trough_idx]
        if rise_left >= drop_threshold and rise_right >= drop_threshold:
            return 'U_Shaped'

    return 'Linear'


def classify_selectivity_trend(experimental_data: List[Dict], threshold_r2: float = 0.6,
                               threshold_slope: float = 0.02,
                               peak_min_elements: int = 4,
                               peak_drop_threshold: float = 0.3) -> Dict:
    """
    Classify selectivity trend as Ascending/Reverse/Flat/Mixed/Peaked/U_Shaped.

    Algorithm (per condition):
    1. Extract (atomic_number, log_D) pairs
    2. Linear regression → candidate label (Ascending/Reverse/Flat)
    3. Peak-position heuristic → if clear arch/trough, override with Peaked/U_Shaped
    4. Aggregate votes across conditions

    The peak-position heuristic fires when:
    - >= peak_min_elements distinct lanthanides measured
    - The max (or min) log_D sits at an interior element (not first/last)
    - Both sides drop (or rise) by >= peak_drop_threshold log units from that extremum

    Args:
        experimental_data: List of experimental condition dictionaries
        threshold_r2: Minimum R² for a monotone linear trend (default 0.6)
        threshold_slope: Minimum |slope| for Ascending/Reverse (default 0.02)
        peak_min_elements: Minimum elements needed for peak heuristic (default 4)
        peak_drop_threshold: Min log_D drop on each side of peak (default 0.3)

    Returns:
        Dict with trend_classification, confidence, element_coverage, slope, r_squared.
        When Peaked/U_Shaped: also peak_element.
    """
    from scipy.stats import linregress

    lanthanides_ordered = get_lanthanide_order()
    _Z_TO_ELEM = {57 + i: sym for i, sym in enumerate(lanthanides_ordered)}

    trend_votes = []
    all_elements = set()
    slope_values = []
    r_squared_values = []
    peak_element_votes = []

    for exp_data in experimental_data:
        if is_mixture(exp_data):   # a mixed-extractant D is not this ligand's alone
            continue
        d_ratios = exp_data.get('distribution_ratios', {})

        points = []
        for elem in lanthanides_ordered:
            if elem in d_ratios and d_ratios[elem].get('log_value') is not None \
                    and d_ratios[elem].get('censor') not in ('>', '<'):   # a bound is not a point
                atomic_num = 57 + lanthanides_ordered.index(elem)
                points.append((atomic_num, d_ratios[elem]['log_value']))
                all_elements.add(elem)

        if len(points) < 4:
            continue

        # Linear regression
        try:
            x_arr = [p[0] for p in points]
            y_arr = [p[1] for p in points]
            slope, intercept, r_value, p_value, std_err = linregress(x_arr, y_arr)
            r_squared = r_value ** 2
            slope_values.append(slope)
            r_squared_values.append(r_squared)
        except Exception as e:
            print(f"Warning: Linear regression failed - {e}")
            continue

        # Linear label
        if r_squared < threshold_r2:
            label = 'Flat'
        elif slope > threshold_slope:
            label = 'Ascending'
        elif slope < -threshold_slope:
            label = 'Reverse'
        else:
            label = 'Flat'

        # Peak-position heuristic (may override linear label)
        peak_shape = peak_position_heuristic(
            points,
            min_elements=peak_min_elements,
            drop_threshold=peak_drop_threshold,
        )
        if peak_shape in ('Peaked', 'U_Shaped'):
            label = peak_shape
            # Record which element is the peak/trough
            sorted_pts = sorted(points, key=lambda p: p[0])
            log_d = [p[1] for p in sorted_pts]
            if peak_shape == 'Peaked':
                extremum_idx = int(np.argmax(log_d))
            else:
                extremum_idx = int(np.argmin(log_d))
            peak_z = sorted_pts[extremum_idx][0]
            peak_elem = _Z_TO_ELEM.get(peak_z)
            if peak_elem:
                peak_element_votes.append(peak_elem)

        trend_votes.append(label)

    # Conditions measured with a co-extractant are not in the fit; say so, or a
    # mixture-only record would read as one with nothing measured.
    n_mixed = sum(1 for c in experimental_data if is_mixture(c))
    mixed_note = (f'; {n_mixed} condition(s) measured with a co-extractant (mixed solvent) are not in the fit'
                  if n_mixed else '')

    # Insufficient data
    if not trend_votes:
        if n_mixed and n_mixed == len(experimental_data):
            insufficient = (f'All {n_mixed} condition(s) were measured with a co-extractant (mixed solvent); '
                            'none enters the fit')
        else:
            insufficient = 'Less than 4 elements tested per condition' + mixed_note
        return {
            'trend_classification': 'Insufficient_Data',
            'trend_confidence': 'None',
            'element_coverage': sorted(list(all_elements), key=lambda x: lanthanides_ordered.index(x) if x in lanthanides_ordered else 99),
            'slope': None,
            'r_squared': None,
            'notes': insufficient
        }

    # Aggregate votes
    vote_counts = Counter(trend_votes)
    dominant_trend = vote_counts.most_common(1)[0][0]
    dominant_count = vote_counts.most_common(1)[0][1]

    if len(vote_counts) == 1:
        confidence = 'High'
    elif dominant_count / len(trend_votes) > 0.7:
        confidence = 'Medium'
    else:
        confidence = 'Low'
        dominant_trend = 'Mixed'

    avg_slope = np.mean(slope_values) if slope_values else None
    avg_r_squared = np.mean(r_squared_values) if r_squared_values else None

    result = {
        'trend_classification': dominant_trend,
        'trend_confidence': confidence,
        'element_coverage': sorted(list(all_elements), key=lambda x: lanthanides_ordered.index(x) if x in lanthanides_ordered else 99),
        'slope': round_if_not_none(avg_slope, 4),
        'r_squared': round_if_not_none(avg_r_squared, 3),
        'notes': f'Based on {len(trend_votes)} experimental condition(s) with linear regression + peak heuristic'
                 + mixed_note
    }

    if dominant_trend in ('Peaked', 'U_Shaped') and peak_element_votes:
        peak_counter = Counter(peak_element_votes)
        result['peak_element'] = peak_counter.most_common(1)[0][0]

    return result
