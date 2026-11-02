# Three Sources of Uncertainty in RE-QUEST

**Date**: March 4, 2026

---

## Overview

You identified the three main sources of uncertainty correctly:

1. **Trend Uncertainty** - Statistical uncertainty in slope/classification due to limited data coverage
2. **Measurement Uncertainty** - Precision of individual D-value measurements (figure vs. table extraction)
3. **Cross-Source Discrepancy** - Inconsistency when the same ligand appears in multiple databases

These are **distinct but related** - they combine to give the overall uncertainty in our predictions.

---

## 1. Trend Uncertainty (Statistical Coverage)

### What It Is
Uncertainty in the **slope** and **trend classification** due to:
- Limited number of lanthanides measured (e.g., 5 out of 15)
- Limited number of experimental conditions
- Natural variability in trends across different conditions
- Random sampling variation

### Even If Individual Measurements Were Perfect...
We'd still have trend uncertainty because we don't measure all 15 lanthanides at all possible conditions.

### How We Quantify It
**Two-stage bootstrap** (1000 replicates, `rel_quest.uncertainty.bootstrap_trend_ci`):
- The stored slope is the *mean of per-condition regression slopes*, so the
  sampling unit is the condition
- Each replicate: draw the conditions with replacement, then within each
  drawn condition resample its (atomic_number, log_D) pairs and refit;
  average the slopes
- 2.5th / 97.5th percentiles → 95% confidence interval; standard deviation → `slope_std`
- Between-condition spread (acidity, diluent, source) dominates when there
  are several conditions; the within-condition stage keeps the fit
  uncertainty as a floor so two agreeing conditions do not read as certainty.
  With one condition the result is that condition's own case bootstrap.
- (Before 2026-08-21 the pairs from all conditions were pooled into one
  regression and resampled — a different estimator from `slope`.)

### Example: TODGA
```
Data available:
  15 elements × 12 conditions = 180 (Z, log_D) pairs

Bootstrap result:
  Slope: 0.0847 ± 0.0062
  95% CI: [0.0785, 0.0909]
  R²: 0.94

Interpretation:
  - Trend: Ascending (slope > 0)
  - Confidence: HIGH (CI doesn't overlap 0, narrow interval)
  - If we had only 5 elements instead of 15:
    → Wider CI, maybe MEDIUM confidence
```

### Factors That Reduce This Uncertainty
✅ Measure more lanthanides (especially missing ones)
✅ Measure more experimental conditions
✅ Ensure lanthanides span the series (La, Eu, Lu better than La, Ce, Pr)

---

## 2. Measurement Uncertainty (Data Source Quality)

### What It Is
Precision/accuracy of individual D-value measurements due to:
- **Figure extraction**: Digitizing graphs → ±5-10% error
- **Table data**: Rounding, significant figures → ±2-5% error
- **Experimental variation**: If replicates reported → use actual std dev

### This Affects Individual Data Points
Unlike trend uncertainty (which is about the *pattern*), this is about the *quality* of each measurement.

### How We Quantify It
**Error propagation**:
```python
# For D-value from table, without a reported error
D = 52.0
σ_D = 52.0 × 0.035 = 1.8  # ±3.5% assumed (rounding, repeatability)

# Propagate to log D
log_D = log₁₀(52.0) = 1.716
σ_log_D = σ_D / (D × ln(10)) = 1.8 / (52.0 × 2.303) = 0.015

# Report
log_D = 1.716 ± 0.015
```

A value read from a figure gets σ_log_D = 0.0715 directly (about ±17% on D).
That prior is measured, not assumed: it is the median, over the 64 figures
calibrated during the v14 verification, of the reading uncertainty found on
each (one marker radius in dex; `uncertainty.figure_prior_evidence`, pinned by
`tests/test_sigma_priors.py`). Until 2026-10-23 it was an assumed 7.5% on D
(0.033 dex), half what reading the calibrated figures actually costs.

### Example: Same Ligand, Different Source Quality
```
Ligand X extracting Eu at 1M HNO₃:

Source A (Table in paper):
  D = 15.2 ± 0.5  (3.3% error from reported uncertainty)
  log_D = 1.182 ± 0.014

Source B (Digitized from Figure 3):
  D = 14.8 ± 2.4  (figure prior, 0.0715 dex)
  log_D = 1.170 ± 0.072

Both measure ~15, but Source A is more precise (narrower error bars)
```

### How This Combines with Trend Uncertainty
When doing bootstrap resampling, we can **weight** data points by their measurement quality:
- High-quality (table) measurements get more weight
- Low-quality (figure) measurements get less weight

```python
# Weighted bootstrap
def weighted_bootstrap(Z, log_D, weights, n_bootstrap=1000):
    """
    Resample with probability proportional to measurement quality

    weights = 1 / σ²_log_D  (inverse variance weighting)
    """
    slopes = []
    for i in range(n_bootstrap):
        # Resample with weights
        indices = np.random.choice(len(Z), size=len(Z),
                                  replace=True, p=weights/weights.sum())
        # ... compute slope
    return slope_mean, slope_CI
```

### Factors That Reduce This Uncertainty
✅ Prefer table data over figures when available
✅ Use papers that report error bars (rare in old literature)
✅ Conduct new experiments with replicate measurements

---

## 3. Cross-Source Discrepancy (Consistency)

### What It Is
When the **same ligand** appears in multiple databases (Liu 2022, SAFE, IDEAL), do they report consistent D-values for the same conditions?

### This Is a Quality Check, Not Independent Data
If Liu reports log_D = 1.5 and SAFE reports log_D = 1.8 for the **exact same** ligand/element/condition:
- They can't both be right
- One (or both) might have transcription errors
- We need to decide which to trust or average them

### How We Quantify It
**Consistency score**:
```python
# For each (element, acid, temp, concentration) combination:
#   - Collect log_D from all sources
#   - Calculate standard deviation
#   - Flag if max_difference > 0.5 log units

Example: TDDGA (merged from Liu + SAFE)

Condition: Eu, 1M HNO₃, 25°C, 0.1M extractant
  Liu 2022:  log_D = 1.45
  SAFE:      log_D = 1.52
  Difference: 0.07 (within 0.5 threshold ✓)

Condition: Lu, 3M HNO₃, 25°C, 0.05M extractant
  Liu 2022:  log_D = 2.31
  SAFE:      log_D = 2.85
  Difference: 0.54 (exceeds threshold! ⚠️)

Consistency score: 1 agreement / 2 comparisons = 50%
→ Flag for manual review
```

### How This Relates to Other Uncertainties
- If sources disagree, we have **systematic uncertainty** (which is correct?)
- If sources agree, we have **confirmation** (reduces uncertainty)
- We can use agreement level as a quality weight

### Resolution Strategies
1. **Weighted average** (if difference < threshold):
   ```
   Liu: log_D = 1.45 ± 0.03 (table, high quality, weight=10)
   SAFE: log_D = 1.52 ± 0.05 (figure, medium quality, weight=4)

   Weighted mean = (10×1.45 + 4×1.52) / (10+4) = 1.47
   Combined uncertainty = sqrt(1/(10+4)) = 0.027
   ```

2. **Flag and exclude** (if difference > threshold):
   ```
   Lu discrepancy of 0.54 log units → exclude from analysis
   Mark in database: "Conflicting measurements - needs verification"
   ```

3. **Manual verification**:
   ```
   Check original papers:
   - Liu 2022 SI Table S2: log_D = 2.31 ✓
   - SAFE record links to: Wrong paper? Transcription error?
   → Keep Liu value, flag SAFE as potential error
   ```

### Factors That Improve Consistency
✅ Careful data curation (double-check transcription)
✅ Standardized extraction protocols (reduces experimental variation)
✅ Direct data sharing between databases (avoid re-digitization)

---

## How The Three Combine: Complete Uncertainty Budget

### Example: TODGA extracting Eu

#### Individual Measurement (Source 2)
```
Source: Liu 2022 Figure 3a
D = 54.0 (digitized from graph)
Measurement uncertainty: σ_D = 54.0 × 0.075 = 4.0
→ log_D = 1.732 ± 0.032
```

#### Trend Across Lanthanide Series
```
All measurements (15 elements, 12 conditions):
Two-stage bootstrap (resample conditions, then points within each):
  Slope = 0.0847 ± 0.0062

At Z=63 (Eu):
  Predicted: log_D = 0.0847 × 63 + intercept = 1.69
  Trend uncertainty: ±0.15 (from slope uncertainty propagated)
```

#### Cross-Source Comparison
```
Liu 2022:  Eu at 1M HNO₃ → log_D = 1.732
SAFE:      Eu at 1M HNO₃ → log_D = 1.745
Difference: 0.013 (excellent agreement ✓)

Consistency check: PASS
Combined value: (1.732 + 1.745)/2 = 1.739
```

#### **Total Uncertainty**
Combining all three sources:
```
1. Measurement: ±0.032 (from digitization)
2. Trend:       ±0.15  (from limited coverage)
3. Cross-source: ±0.007 (std of two measurements)

Combined (root-sum-square):
σ_total = sqrt(0.032² + 0.15² + 0.007²) = 0.153

Final result:
log_D[Eu] = 1.74 ± 0.15
```

**Which dominates?** Trend uncertainty (±0.15) >> Measurement (±0.03)
→ Measuring more lanthanides helps more than improving individual measurement precision!

---

## Practical Implications

### What This Means for Users

#### High-Quality Ligand (Low Uncertainty)
```
Ligand: CyMe4BTBP
✓ 15/15 lanthanides measured
✓ 20 experimental conditions
✓ All from tables (no figure digitization)
✓ Liu + SAFE data agree (consistency = 95%)

Result:
  Slope: 0.112 ± 0.008
  Trend: Ascending
  Confidence: HIGH

→ Use immediately for process design
→ High confidence in performance predictions
```

#### Low-Quality Ligand (High Uncertainty)
```
Ligand: Obscure_DGA_variant
✗ Only 4/15 lanthanides measured (La, Eu, Gd, Lu)
✗ Only 2 experimental conditions
✗ All from figures (±0.07 dex, about ±17%)
✗ Only in IDEAL database (no cross-validation)

Result:
  Slope: 0.045 ± 0.089  (huge uncertainty!)
  Trend: Insufficient_Data
  Confidence: N/A

→ DO NOT use for process design
→ Priority candidate for active learning (measure more elements)
→ Recommend measuring: Pr, Nd, Dy, Er (fill gaps)
```

---

## How They Interact in Active Learning

### Scenario: Recommend Next Experiment

```python
# Ligand X current state:
Current_data = {
    "elements_measured": ["La", "Ce", "Eu", "Lu"],  # 4/15
    "conditions": 3,
    "sources": ["Liu2023"],  # Only one source
    "measurement_quality": "Medium" (mix of table/figure)
}

# Uncertainties:
Trend_uncertainty = 0.052  # Large (only 4 elements)
Measurement_uncertainty = 0.028  # Medium (mixed quality)
Cross_source_uncertainty = None  # Only one source (can't check)

# Which uncertainty should we reduce first?
Priority 1: Trend uncertainty (largest contributor)
  → Measure more lanthanides: Recommend Pr, Dy, Er
  → Expected reduction: 0.052 → 0.025 (50% reduction!)

Priority 2: Cross-source validation
  → Find this ligand in SAFE or conduct new experiment
  → Check if trend agrees

Priority 3: Measurement quality
  → Re-measure from tables instead of figures
  → Expected reduction: 0.028 → 0.015 (smaller impact)

Recommendation:
  "Measure Pr, Dy, Er for ligand X at 1M HNO₃, 25°C
   Expected uncertainty reduction: 0.052 → 0.025
   Priority: HIGH"
```

---

## Summary Table

| Source of Uncertainty | What It Measures | How to Quantify | How to Reduce |
|----------------------|------------------|-----------------|---------------|
| **1. Trend Uncertainty** | Statistical uncertainty in slope due to limited data coverage | Bootstrap resampling (95% CI for slope) | Measure more lanthanides, more conditions |
| **2. Measurement Uncertainty** | Precision of individual D-values | Error propagation (σ_D → σ_log_D) | Use table data, get replicates, new experiments |
| **3. Cross-Source Discrepancy** | Consistency when same ligand in multiple databases | Standard deviation across sources, flag if >0.5 log units | Manual verification, weighted averaging |

### Which Dominates?
- **For well-studied ligands** (15 elements, many conditions): Measurement uncertainty dominates
- **For poorly-studied ligands** (few elements): Trend uncertainty dominates ← Most common case!
- **For merged ligands**: Cross-source discrepancy can dominate if sources disagree

### Combined Uncertainty
```
σ_total = sqrt(σ_trend² + σ_measurement² + σ_cross_source²)
```

---

## Database Schema Enhancement

To track all three, add to each ligand record:

```json
{
  "ligand": {
    "Name": "TODGA",
    "SMILES": "..."
  },
  "selectivity": {
    "trend_classification": "Ascending",
    "slope": 0.0847,
    "slope_uncertainty": {
      "trend_CI_lower": 0.0785,
      "trend_CI_upper": 0.0909,
      "trend_std": 0.0062,
      "source": "bootstrap_1000_iterations"
    },
    "confidence": "High",
    "R2": 0.94
  },
  "experimental_data": [
    {
      "condition_id": "...",
      "distribution_ratios": {
        "Eu": {
          "value": 54.0,
          "log_value": 1.732,
          "measurement_uncertainty": {
            "sigma_log_D": 0.032,
            "source_quality": "figure",
            "error_percent": 7.5
          },
          "n_replicates": 1
        }
      }
    }
  ],
  "data_quality": {
    "cross_source_consistency": {
      "score": 0.95,  # 95% agreement
      "sources_compared": ["source A", "source B"],
      "flagged_conditions": [
        {
          "element": "Lu",
          "condition": "3M_HNO3_25C",
          "discrepancy": 0.54,
          "action": "manual_review_needed"
        }
      ]
    },
    "overall_quality": "High",  # Based on all three factors
    "recommended_measurements": ["Tm", "Dy"]  # Active learning
  }
}
```

---

## Next Steps

1. **Validate the assumptions**:
   - ~~Check if ±7.5% for figures is reasonable~~ — it was not; the figure prior is now measured (0.0715 dex, 2026-10-23)
   - Compare estimated uncertainties to papers that report error bars
   - Test bootstrap CI calibration (do 95% of new measurements fall in 95% CI?)

2. **Implement in priority order**:
   - Month 1: Bootstrap trend uncertainty (biggest impact for most ligands)
   - Month 2: Measurement uncertainty classification
   - Month 3: Cross-source consistency checks

3. **Visualize in dashboard**:
   - Show all three uncertainties in separate columns
   - Color-code quality: Green (high), Yellow (medium), Red (low)
   - Recommend measurements based on which uncertainty dominates

Would you like me to create a validation notebook that tests these methods on 5 representative ligands (high/medium/low quality)?
