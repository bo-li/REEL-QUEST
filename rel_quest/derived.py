"""Every stored field that is a function of the measurements, in one place.

A record stores two kinds of data:

* **source-transcribed** — identifiers, names, conditions, each cell's `value`
  and `log_value`, and `reference`; these can only be checked against the
  source papers;
* **derived** — everything below, which must reproduce exactly from the
  record's own measurements.

Each derived field already had a canonical function, and `tests/test_data_integrity.py`
pinned most of them one by one. `recompute_record` applies all of them in the
order the ingest writers use, so

    python scripts/recompute_derived_fields.py --check

is a single, one-command reproducibility check for the whole database, and an
edit to the measurements (a correction, or dropping one source's conditions
for a release bundle) can be followed by one call instead of seven.

Derived fields, and the function that owns each:

    cell.sigma_log_D, cell.source_quality   uncertainty.classify_measurement_quality
                                            + uncertainty.compute_sigma_log_d
                                            (sigma_log_D only unless sigma_source is
                                            "measured", when it is source data)
    condition.statistics                    utils.condition_statistics
    elements_measured                       utils.collect_elements_measured
    selectivity (trend, slope, notes, ...)  trends.classify_selectivity_trend
    selectivity.slope_CI_*                  uncertainty.bootstrap_trend_ci (seeded)
    data_quality.outlier_flags              uncertainty.detect_outliers
    data_quality.overall_quality_score,
    data_quality.quality_grade              uncertainty.compute_quality_score
    performance_summary                     performance_summary() below
"""

from __future__ import annotations

import copy
from typing import Dict, List

from .trends import classify_selectivity_trend
from .uncertainty import (
    bootstrap_trend_ci,
    classify_measurement_quality,
    compute_quality_score,
    compute_sigma_log_d,
    detect_outliers,
)
from .utils import collect_elements_measured, condition_statistics

CI_KEYS = ("slope_CI_lower", "slope_CI_upper", "slope_std")


def performance_summary(record: Dict) -> Dict:
    """Condition and cell counts plus the measured elements.

    Until v1.0.0 this block was written only by some ingest paths and never
    refreshed: 62 records had it empty and 34 carried counts that no longer
    matched their conditions after the v11 merge — and the REST API returns
    it verbatim.
    """
    conditions = record.get("experimental_data", [])
    return {
        "conditions_tested": len(conditions),
        "total_measurements": sum(len(c.get("distribution_ratios") or {}) for c in conditions),
        "elements_tested": list(record.get("elements_measured") or []),
    }


def recompute_record(record: Dict) -> Dict:
    """Return a deep copy of `record` with every derived field recomputed.

    Source-transcribed fields are copied unchanged.
    """
    fresh = copy.deepcopy(record)
    conditions = fresh.get("experimental_data", [])

    for condition in conditions:
        quality, sigma_rel = classify_measurement_quality(condition)
        sigma = compute_sigma_log_d(sigma_rel)
        for cell in (condition.get("distribution_ratios") or {}).values():
            if cell.get("sigma_source") != "measured":   # a measured sigma is source data
                cell["sigma_log_D"] = sigma
            cell["source_quality"] = quality
        condition["statistics"] = condition_statistics(condition.get("distribution_ratios"))

    fresh["elements_measured"] = collect_elements_measured(fresh)

    selectivity = classify_selectivity_trend(conditions)
    ci = bootstrap_trend_ci(conditions) if selectivity.get("slope") is not None else None
    for key in CI_KEYS:
        selectivity[key] = ci[key] if ci else None
    fresh["selectivity"] = selectivity

    flags = detect_outliers(conditions)
    score, grade = compute_quality_score(fresh, ci, flags)
    fresh["data_quality"] = {
        "overall_quality_score": score,
        "quality_grade": grade,
        "outlier_flags": flags,
    }

    fresh["performance_summary"] = performance_summary(fresh)
    return fresh


def derived_differences(stored: Dict, fresh: Dict) -> List[str]:
    """Paths of derived fields whose stored value differs from `fresh`."""
    diffs = []
    for key in ("elements_measured", "performance_summary"):
        if stored.get(key) != fresh.get(key):
            diffs.append(key)
    for key in sorted(set(stored.get("selectivity") or {}) | set(fresh["selectivity"])):
        if (stored.get("selectivity") or {}).get(key) != fresh["selectivity"].get(key):
            diffs.append(f"selectivity.{key}")
    for key in ("overall_quality_score", "quality_grade", "outlier_flags"):
        if (stored.get("data_quality") or {}).get(key) != fresh["data_quality"][key]:
            diffs.append(f"data_quality.{key}")
    for old, new in zip(stored.get("experimental_data", []), fresh["experimental_data"]):
        if old.get("statistics") != new["statistics"]:
            diffs.append(f"{old.get('condition_id')}.statistics")
        for element, cell in (old.get("distribution_ratios") or {}).items():
            new_cell = new["distribution_ratios"][element]
            for key in ("sigma_log_D", "source_quality"):
                if cell.get(key) != new_cell[key]:
                    diffs.append(f"{old.get('condition_id')}.{element}.{key}")
    return diffs
