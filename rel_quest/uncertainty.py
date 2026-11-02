"""
Uncertainty Quantification (UQ) for REEL-QUEST ligand database.

Three sources of uncertainty:
1. Bootstrap CI on trend slopes
2. Measurement uncertainty (sigma_log_D per D-value)
3. Outlier detection via Z-score on per-condition regression residuals
"""

import math
from typing import Dict, List, Optional, Tuple

import numpy as np
from scipy.stats import linregress, zscore

LANTHANIDES = ['La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu']
ATOMIC_NUMBER = {sym: 57 + i for i, sym in enumerate(LANTHANIDES)}


# A condition contributes to the trend fit only if it holds at least this
# many lanthanide log D values. Mirrors the `len(points) < 4: continue` in
# rel_quest.trends.classify_selectivity_trend; tests pin the two together.
MIN_POINTS_PER_CONDITION = 4


def is_censored(cell) -> bool:
    """True when a cell holds a bound (D > x or D < x), not a measured value.

    `value`/`log_value` of a censored cell carry the bound, so every consumer
    that treats a cell as a measurement (trend fits, outliers, statistics,
    separation factors) must skip it.
    """
    return isinstance(cell, dict) and cell.get('censor') in ('>', '<')


def is_mixture(condition) -> bool:
    """True when a condition was measured with other extractants in the organic phase.

    Its D belongs to the mixture, not to the record's ligand alone: it stays on
    the record (a real measurement, shown in the ligand's condition table) but
    is left out of everything pooled across a ligand's conditions - the trend
    fit and its bootstrap, outliers, aggregates and the dashboard's viewers.
    """
    return isinstance(condition, dict) and bool(condition.get('Co_Extractants'))


def condition_trend_points(condition: Dict) -> List[Tuple[int, float]]:
    """(atomic_number, log_value) pairs for one condition, lanthanides only,
    in periodic order. Y and Sc are excluded exactly as the classifier
    excludes them: the fit regresses on Z across 57-71. A mixed-extractant
    condition (is_mixture) gives no points: no fit consumes it."""
    if is_mixture(condition):
        return []
    d_ratios = condition.get('distribution_ratios', {}) or {}
    points = []
    for elem in LANTHANIDES:
        data = d_ratios.get(elem)
        if isinstance(data, dict) and data.get('log_value') is not None and not is_censored(data):
            points.append((ATOMIC_NUMBER[elem], data['log_value']))
    return points


def per_condition_slopes(experimental_data: List[Dict],
                         min_points: int = MIN_POINTS_PER_CONDITION
                         ) -> List[Tuple[float, List[Tuple[int, float]]]]:
    """OLS slope of log D on Z for every eligible condition.

    These are the units `selectivity.slope` averages over
    (`classify_selectivity_trend` stores the mean of exactly these slopes),
    so anything that quantifies the uncertainty of `slope` must resample
    them, not the pooled points.

    Returns a list of (slope, points) pairs, one per eligible condition.
    """
    fits = []
    for cond in experimental_data:
        points = condition_trend_points(cond)
        if len(points) < min_points:
            continue
        x = [p[0] for p in points]
        y = [p[1] for p in points]
        if len(set(x)) < 2:
            continue
        slope, *_ = linregress(x, y)
        fits.append((float(slope), points))
    return fits


def _within_condition_slopes(points: List[Tuple[int, float]], n_draws: int,
                             rng: np.random.Generator) -> np.ndarray:
    """`n_draws` OLS slopes from case-resampling one condition's points.

    Vectorised: slope = Sxy / Sxx per resample. A resample in which every
    drawn point has the same Z has no slope; those rows are redrawn (with
    four or more points this is < 2% of draws).
    """
    x = np.array([p[0] for p in points], dtype=float)
    y = np.array([p[1] for p in points], dtype=float)
    n = len(points)
    out = np.empty(n_draws)
    pending = np.arange(n_draws)
    while pending.size:
        idx = rng.integers(0, n, size=(pending.size, n))
        xs, ys = x[idx], y[idx]
        xm = xs.mean(axis=1, keepdims=True)
        ym = ys.mean(axis=1, keepdims=True)
        sxx = ((xs - xm) ** 2).sum(axis=1)
        ok = sxx > 0
        out[pending[ok]] = ((xs - xm) * (ys - ym)).sum(axis=1)[ok] / sxx[ok]
        pending = pending[~ok]
    return out


def bootstrap_trend_ci(experimental_data: List[Dict], n_bootstrap: int = 1000,
                       seed: int = 42) -> Optional[Dict]:
    """
    95% bootstrap confidence interval for `selectivity.slope`.

    The stored slope is the mean of per-condition regression slopes, so the
    sampling unit here is the *condition*, in a two-stage (cluster, then
    within-cluster) bootstrap:

    1. Draw `n_cond` conditions with replacement from the eligible ones
       (those with at least MIN_POINTS_PER_CONDITION lanthanides — the
       same ones the classifier averaged).
    2. For each drawn condition, case-resample its (Z, log D) points and
       refit its slope.
    3. The replicate is the mean of those slopes. Repeat `n_bootstrap`
       times; the 2.5th / 97.5th percentiles are the interval and the
       standard deviation is `slope_std`.

    Stage 1 carries the between-condition spread — different acidities,
    diluents, sources — which is the dominant real uncertainty whenever
    there are several conditions. Stage 2 keeps the within-condition fit
    uncertainty as a floor: with only two or three conditions a
    conditions-only bootstrap collapses to the range of their slopes, and
    to zero width when two conditions happen to agree, which would be
    read as certainty on two data points. With a single eligible condition
    stage 1 is the identity and the result is that condition's own case
    bootstrap — the same code path, no special case. Resampling at both
    levels slightly overstates variance (Davison & Hinkley 1997, §3.8);
    with this few clusters that is the conservative side to err on.

    No eligible condition (every one has fewer than
    MIN_POINTS_PER_CONDITION lanthanides) returns None — the record is
    `Insufficient_Data` and carries no slope.

    Until 2026-08-21 this function pooled every point from every condition
    into one regression and bootstrapped that. Pooling answers a different
    question from the per-condition mean whenever conditions differ in
    intercept (acid concentration shifts log D bodily) or direction
    (`Mixed`), and on eight ligands the stored slope lay outside its own
    stored interval. See docs/PROJECT_REVIEW_2026-08-21.md finding 2.

    Args:
        experimental_data: List of experimental condition dicts
        n_bootstrap: Number of bootstrap replicates (default 1000)
        seed: Seed for numpy's Generator, so stored values are reproducible

    Returns:
        Dict with slope_CI_lower, slope_CI_upper, slope_std, or None
    """
    fits = per_condition_slopes(experimental_data)
    if not fits:
        return None

    rng = np.random.default_rng(seed)
    n_cond = len(fits)

    # Stage 2 pools: n_bootstrap within-condition resampled slopes for each
    # condition, drawn up front so stage 1 can index them.
    pools = np.vstack([_within_condition_slopes(points, n_bootstrap, rng)
                       for _, points in fits])            # (n_cond, n_bootstrap)

    # Stage 1: each replicate draws n_cond conditions with replacement and,
    # for every draw, an independent within-condition resample.
    cond_draws = rng.integers(0, n_cond, size=(n_bootstrap, n_cond))
    pool_draws = rng.integers(0, n_bootstrap, size=(n_bootstrap, n_cond))
    boot = pools[cond_draws, pool_draws].mean(axis=1)

    return {
        'slope_CI_lower': round(float(np.percentile(boot, 2.5)), 6),
        'slope_CI_upper': round(float(np.percentile(boot, 97.5)), 6),
        'slope_std': round(float(np.std(boot)), 6),
    }


# sigma(log D) given to a value without a measured sigma, by how it was read.
#   figure  measured: the median over calibrated figures of the reading uncertainty
#           found on them (figure_prior_evidence; tests/test_sigma_priors.py keeps the
#           two within 0.005 dex). 7.5 % on D (0.033 dex) until 2026-10-23.
#   table   assumed: 3.5 % on D for rounding and repeatability. Too few papers report
#           an experimental sigma to measure it.
TABLE_PRIOR_ON_D = 0.035
# The floor vector-figure readings were clamped to (calibrations' sigma_floor): the
# table prior when they were read, kept fixed here so the evidence below does not
# move if the table prior ever does.
VECTOR_FIGURE_FLOOR = 0.0152
SIGMA_PRIORS = {'table': round(TABLE_PRIOR_ON_D / math.log(10), 4), 'figure': 0.0715}
# The same priors as relative uncertainty on D, the unit compute_sigma_log_d takes.
QUALITY_PRIORS = {quality: sigma * math.log(10) for quality, sigma in SIGMA_PRIORS.items()}


def figure_prior_evidence(records: List[Dict]) -> Dict:
    """The reading uncertainty measured on calibrated figures, one median per figure.

    Counts measured sigmas of digitised, uncensored cells. A reading at exactly
    VECTOR_FIGURE_FLOOR is left out: vector figures clamp their readings to that
    floor (digitisation.digitisation_sigma), so it says "at least the floor", not
    what reading the figure cost. A figure is a (DOI, Locator) pair; one whose
    readings are all at the floor counts in n_figures_calibrated only.
    """
    figures: Dict[Tuple, List[float]] = {}
    calibrated = set()
    clamped = 0
    for record in records:
        for condition in record.get('experimental_data', []):
            reference = condition.get('reference') or {}
            key = (reference.get('DOI'), reference.get('Locator'))
            for cell in (condition.get('distribution_ratios') or {}).values():
                if not isinstance(cell, dict) or cell.get('sigma_source') != 'measured' \
                        or cell.get('value_origin') != 'digitised' or is_censored(cell):
                    continue
                calibrated.add(key)
                if cell.get('sigma_log_D') == VECTOR_FIGURE_FLOOR:
                    clamped += 1
                    continue
                figures.setdefault(key, []).append(cell['sigma_log_D'])
    per_figure = [float(np.median(sigmas)) for sigmas in figures.values()]
    q1, median, q3 = (np.percentile(per_figure, [25, 50, 75]) if per_figure else (None,) * 3)
    return {
        'n_cells': sum(len(sigmas) for sigmas in figures.values()),
        'n_figures': len(per_figure),
        'n_figures_calibrated': len(calibrated),
        'n_floor_excluded': clamped,
        'median_dex': None if median is None else round(float(median), 4),
        'q1_dex': None if q1 is None else round(float(q1), 4),
        'q3_dex': None if q3 is None else round(float(q3), 4),
    }


def classify_measurement_quality(condition: Dict) -> Tuple[str, float]:
    """
    Heuristic classification of measurement quality based on condition metadata.

    Returns:
        (source_quality, sigma_rel) where source_quality is 'table'|'figure'|'unknown'
        and sigma_rel is the relative uncertainty fraction on D.
    """
    condition_id = condition.get('condition_id', '')
    reference = condition.get('reference', {}) or {}
    notes = str(reference.get('Notes', '') or '')

    if condition_id.startswith(('SAFE_', 'Lit_')):
        # Tabulated only when every location named is a table (SAFE, and the
        # conditions this project read directly from a paper: 'Lit_'). Reads the
        # structured reference.Locator (schema v0.3); older records fall back
        # to parsing the notes. The previous rule — 'Table' anywhere in the
        # notes — missed SAFE's lowercase "table1" (6 conditions) and called
        # "Figure 10 and Table 1" tabulated (2), where the tighter prior is
        # not justified.
        locator = reference.get('Locator')
        if locator is None:
            from .references import parse_locator
            locator = parse_locator(notes) or ''
        kinds = {part.split()[0] for part in locator.split('; ') if part}
        if kinds == {'Table'}:
            return ('table', QUALITY_PRIORS['table'])
        return ('figure', QUALITY_PRIORS['figure'])

    # Pramanik 2024: D-values transcribed from SI Tables S1/S2 (triplicate
    # ICP measurements with reported standard deviations), not digitized.
    if condition_id.startswith('Pramanik2024_'):
        return ('table', QUALITY_PRIORS['table'])

    # Liu 2022 and IDEAL: digitized from paper figures
    if condition_id.startswith(('Liu2022_', 'IDEAL_')):
        return ('figure', QUALITY_PRIORS['figure'])

    # Any other id: the public release renames conditions to neutral ids
    # (scripts/export_public_repo.py), so the class the cells were stamped
    # with at ingest is the input. One class per condition, as written above.
    stamped = {cell.get('source_quality')
               for cell in (condition.get('distribution_ratios') or {}).values()}
    if len(stamped) == 1 and stamped <= set(QUALITY_PRIORS):
        quality = stamped.pop()
        return (quality, QUALITY_PRIORS[quality])
    return ('figure', QUALITY_PRIORS['figure'])


def compute_sigma_log_d(sigma_rel: float) -> float:
    """
    Propagate relative D uncertainty to log10(D) uncertainty.

    sigma_log_D = sigma_rel / ln(10)

    Args:
        sigma_rel: Relative uncertainty in D (e.g. 0.035 for 3.5%)

    Returns:
        Uncertainty in log10(D)
    """
    return round(sigma_rel / math.log(10), 4)


def detect_outliers(experimental_data: List[Dict], threshold: float = 3.0) -> List[Dict]:
    """
    Detect outlier D-value measurements: Z-score on *per-condition* regression
    residuals.

    Every eligible condition (>= MIN_POINTS_PER_CONDITION lanthanides, the
    same threshold the trend classifier uses) gets its own OLS fit of log D
    on atomic number; a point's residual is its distance from *its own
    condition's* line. Residuals from all eligible conditions are then
    pooled — each condition's fit has already removed its intercept and
    slope, so they are comparable — and Z-scored together, which gives a
    stable noise scale even when single conditions are small.

    An outlier is therefore "a point that deviates from its condition's own
    trend". Conditions too small to define a trend contribute no residuals
    and cannot be flagged; with a single 4-point condition no |Z| can reach
    3, which is honest — a 2-parameter fit through 4 points has no room to
    call any of them wrong.

    The per-condition model is linear, like the classifier's, so on a
    strongly Peaked or U_Shaped condition the extremes carry systematic
    residuals and can be flagged: such a flag means "off the linear trend",
    which for curved conditions includes curvature, not only error.

    Until 2026-08-25 this function fitted one pooled regression across every
    condition, so on multi-condition records the residuals were dominated by
    between-condition offsets (acid concentration shifts log D bodily) and
    by slope disagreement (`Mixed` records), not by measurement error — the
    same estimator mismatch fixed in `bootstrap_trend_ci` on 2026-08-21.
    See docs/PROJECT_REVIEW_2026-08-21.md finding 2.

    Args:
        experimental_data: List of experimental condition dicts
        threshold: Z-score threshold for outlier detection (default 3.0)

    Returns:
        List of outlier flag dicts, or empty list if insufficient data
    """
    records = []
    residuals = []
    for cond in experimental_data:
        cond_id = cond.get('condition_id', '')
        points = condition_trend_points(cond)
        if len(points) < MIN_POINTS_PER_CONDITION:
            continue
        x = np.array([p[0] for p in points], dtype=float)
        y = np.array([p[1] for p in points], dtype=float)
        slope, intercept, *_ = linregress(x, y)
        y_pred = slope * x + intercept
        for i, (z, log_value) in enumerate(points):
            records.append({
                'element': LANTHANIDES[z - 57],
                'condition_id': cond_id,
                'log_value': log_value,
                'predicted': float(y_pred[i]),
            })
            residuals.append(log_value - float(y_pred[i]))

    if len(records) < 4:
        return []

    residuals = np.array(residuals)
    if np.std(residuals) == 0:
        return []

    z_scores = zscore(residuals)

    flags = []
    for i, rec in enumerate(records):
        z_s = float(z_scores[i])
        if abs(z_s) > threshold:
            severity = 'High' if abs(z_s) > 4 else 'Medium'
            flags.append({
                'element': rec['element'],
                'condition_id': rec['condition_id'],
                'log_D_measured': round(rec['log_value'], 4),
                'log_D_predicted': round(rec['predicted'], 4),
                'residual': round(float(residuals[i]), 4),
                'z_score': round(z_s, 3),
                'severity': severity,
            })

    return flags


def compute_quality_score(record: Dict, ci_result: Optional[Dict],
                          outlier_flags: List[Dict]) -> Tuple[int, str]:
    """
    Compute overall data quality score (0-100) and grade.

    Factors:
      1. CI width (30 pts)
      2. Element coverage (30 pts)
      3. Measurement quality (20 pts)
      4. Outlier penalty (20 pts)

    Args:
        record: Full ligand record dict
        ci_result: Output of bootstrap_trend_ci (or None)
        outlier_flags: Output of detect_outliers

    Returns:
        (score: int, grade: str) where grade is 'High'|'Medium'|'Low'
    """
    score = 0

    # --- Factor 1: CI width (30 pts) ---
    if ci_result is not None:
        ci_width = ci_result['slope_CI_upper'] - ci_result['slope_CI_lower']
        if ci_width < 0.02:
            score += 30
        elif ci_width < 0.05:
            score += 20
        elif ci_width < 0.1:
            score += 10
        # else 0

    # --- Factor 2: Element coverage (30 pts) ---
    n_elements = len(record.get('selectivity', {}).get('element_coverage', []))
    if n_elements >= 10:
        score += 30
    elif n_elements >= 7:
        score += 20
    elif n_elements >= 4:
        score += 10
    else:
        score += 5

    # --- Factor 3: Measurement quality (20 pts) ---
    total_d = 0
    table_d = 0
    for cond in record.get('experimental_data', []):
        sq, _ = classify_measurement_quality(cond)
        n = len(cond.get('distribution_ratios', {}))
        total_d += n
        if sq == 'table':
            table_d += n

    if total_d > 0:
        table_fraction = table_d / total_d
        if table_fraction > 0.75:
            score += 20
        elif table_fraction > 0.50:
            score += 15
        elif table_fraction > 0.25:
            score += 10
        else:
            score += 5
    else:
        score += 5

    # --- Factor 4: Outlier penalty (20 pts) ---
    n_outliers = len(outlier_flags)
    has_high = any(f['severity'] == 'High' for f in outlier_flags)
    if has_high:
        score += 0
    elif n_outliers == 0:
        score += 20
    elif n_outliers == 1:
        score += 15
    elif n_outliers == 2:
        score += 10
    else:
        score += 5

    if score > 75:
        grade = 'High'
    elif score >= 50:
        grade = 'Medium'
    else:
        grade = 'Low'

    return (score, grade)
