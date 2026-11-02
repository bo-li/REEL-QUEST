"""Cross-field integrity invariants for the canonical database.

The schema validator used to check containers, vocabularies and SMILES, but
nothing that related one stored number to another. Two defects lived in that
gap until 2026-08-13 (PROJECT_REVIEW_2026-07-31.md findings 1 and 2), and
both were invisible to a passing validation run:

  1. 407 measurements stored a D-value that contradicted its own logarithm,
     39 of them as a flat 0.0 next to a finite log — so an answer depended
     on whether the caller read `value` or `log_value`.
  2. Twelve selectivity blocks were left over from before the v11 merge
     appended more conditions, carrying stale trends, slopes and coverage.

These tests pin the repaired invariants against the live database. They are
deliberately recomputed from the record's own measurements rather than
compared to a stored snapshot, so they keep meaning as data is added.
"""

import json
import math
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
DB_PATH = ROOT / "data" / "sx_ligands_merged_v14_flat.jsonl"

# Matches validators.D_LOG_TOLERANCE_DEX; sources round log D to 3-4
# decimals, so exact equality is not achievable.
TOLERANCE_DEX = 0.01

LANTHANIDES = ["La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd",
               "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu"]


def load_records():
    return [json.loads(line) for line in DB_PATH.open() if line.strip()]


def measured_elements(record):
    found = set()
    for condition in record.get("experimental_data", []):
        for element, data in (condition.get("distribution_ratios") or {}).items():
            if isinstance(data, dict) and (
                data.get("value") is not None or data.get("log_value") is not None
            ):
                found.add(element)
    return found


class TestDValueConsistency(unittest.TestCase):
    """Finding 1: the two stored forms of a measurement must agree."""

    @classmethod
    def setUpClass(cls):
        cls.records = load_records()

    def test_no_zero_d_alongside_a_finite_log(self):
        """A finite log D implies a strictly positive D."""
        offenders = []
        for record in self.records:
            for condition in record.get("experimental_data", []):
                for element, data in (condition.get("distribution_ratios") or {}).items():
                    if not isinstance(data, dict):
                        continue
                    value, log_value = data.get("value"), data.get("log_value")
                    if log_value is None or value is None:
                        continue
                    if math.isfinite(log_value) and value <= 0:
                        offenders.append(
                            f"{record['ligand']['Name']}/{element}: "
                            f"value={value} log_value={log_value}")
        self.assertEqual(
            offenders, [],
            f"{len(offenders)} measurement(s) store a non-positive D next to a "
            f"finite log D: {offenders[:5]}")

    def test_d_and_log_d_agree_within_tolerance(self):
        worst = (0.0, None)
        offenders = []
        for record in self.records:
            for condition in record.get("experimental_data", []):
                for element, data in (condition.get("distribution_ratios") or {}).items():
                    if not isinstance(data, dict):
                        continue
                    value, log_value = data.get("value"), data.get("log_value")
                    if value is None or log_value is None or value <= 0:
                        continue
                    error = abs(math.log10(value) - log_value)
                    if error > worst[0]:
                        worst = (error, f"{record['ligand']['Name']}/{element}")
                    if error > TOLERANCE_DEX:
                        offenders.append(
                            f"{record['ligand']['Name']}/{element}: {error:.4f} dex")
        self.assertEqual(
            offenders, [],
            f"{len(offenders)} measurement(s) exceed {TOLERANCE_DEX} dex; "
            f"worst overall {worst[0]:.5f} dex at {worst[1]}")

    def test_every_d_value_is_a_finite_number(self):
        for record in self.records:
            for condition in record.get("experimental_data", []):
                for element, data in (condition.get("distribution_ratios") or {}).items():
                    with self.subTest(ligand=record["ligand"]["Name"], element=element):
                        self.assertIsInstance(data, dict)
                        for key in ("value", "log_value"):
                            number = data.get(key)
                            if number is not None:
                                self.assertIsInstance(number, (int, float))
                                self.assertFalse(isinstance(number, bool))
                                self.assertTrue(math.isfinite(number))


class TestElementCoverage(unittest.TestCase):
    """Finding 2: the element lists must match the measurements."""

    @classmethod
    def setUpClass(cls):
        cls.records = load_records()

    def test_every_record_declares_elements_measured(self):
        missing = [r["ligand"]["Name"] for r in self.records
                   if r.get("elements_measured") is None]
        self.assertEqual(missing, [], f"records without elements_measured: {missing}")

    def test_elements_measured_matches_the_measurements(self):
        for record in self.records:
            with self.subTest(ligand=record["ligand"]["Name"]):
                self.assertEqual(set(record["elements_measured"]),
                                 measured_elements(record))

    def test_element_coverage_is_a_subset_of_what_was_measured(self):
        """Coverage lists fitted points, so it cannot name unmeasured elements."""
        for record in self.records:
            coverage = set(record.get("selectivity", {}).get("element_coverage") or [])
            with self.subTest(ligand=record["ligand"]["Name"]):
                self.assertTrue(
                    coverage.issubset(measured_elements(record)),
                    f"coverage names unmeasured elements: "
                    f"{sorted(coverage - measured_elements(record))}")

    def test_coverage_holds_every_measured_lanthanide(self):
        """The fit must use all La-Lu points available, not a stale subset.

        This is what the v11 merge broke: conditions were appended without
        re-running classification, leaving coverage frozen at the pre-merge
        element set. Records with too few points to fit are exempt.
        """
        from rel_quest.uncertainty import is_censored, is_mixture
        for record in self.records:
            # what the fit can use: values (not bounds) from conditions without a co-extractant
            fittable = {"experimental_data": [c for c in record.get("experimental_data", []) if not is_mixture(c)]}
            measured_ln = {el for c in fittable["experimental_data"]
                           for el, d in (c.get("distribution_ratios") or {}).items()
                           if not is_censored(d) and d.get("log_value") is not None} & set(LANTHANIDES)
            if len(measured_ln) < 4:
                continue
            coverage = set(record.get("selectivity", {}).get("element_coverage") or [])
            with self.subTest(ligand=record["ligand"]["Name"]):
                self.assertEqual(
                    measured_ln - coverage, set(),
                    "coverage omits measured lanthanides (stale selectivity block?)")

    def test_y_and_sc_are_findable(self):
        """filter_by_elements must read elements_measured, not coverage.

        Y and Sc can never appear in element_coverage — the trend fit
        regresses on atomic number across La-Lu — so a filter built on that
        field silently returned nothing for them.
        """
        from rel_quest import LigandDatabase

        db = LigandDatabase()
        for element in ("Y", "Sc"):
            with self.subTest(element=element):
                expected = sum(1 for r in self.records
                               if element in measured_elements(r))
                self.assertGreater(expected, 0, "fixture assumption failed")
                self.assertEqual(len(db.filter_by_elements([element])), expected)


class TestSelectivityIsReproducible(unittest.TestCase):
    """Every stored selectivity block must be refittable from its own data."""

    @classmethod
    def setUpClass(cls):
        cls.records = load_records()

    def test_stored_selectivity_matches_a_fresh_classification(self):
        from rel_quest.trends import classify_selectivity_trend

        drifted = []
        for record in self.records:
            stored = record.get("selectivity") or {}
            fresh = classify_selectivity_trend(record.get("experimental_data", []))
            # `notes` and `peak_element` are included on purpose: the
            # 2026-08-13 refit preserved them, and twelve records then
            # described a fit that no longer existed (review 2026-08-21).
            for key in ("trend_classification", "trend_confidence",
                        "element_coverage", "slope", "r_squared",
                        "notes", "peak_element"):
                if stored.get(key) != fresh.get(key):
                    drifted.append(
                        f"{record['ligand']['Name']}.{key}: "
                        f"{stored.get(key)!r} != {fresh.get(key)!r}")
        self.assertEqual(
            drifted, [],
            f"{len(drifted)} selectivity field(s) are not reproducible from the "
            f"record's own measurements: {drifted[:5]}")

    def test_stored_ci_matches_a_fresh_bootstrap(self):
        """The interval must come from the current estimator, seeded."""
        from rel_quest.uncertainty import bootstrap_trend_ci

        drifted = []
        for record in self.records:
            sel = record.get("selectivity") or {}
            fresh = bootstrap_trend_ci(record.get("experimental_data", []))
            if sel.get("slope") is None:
                # Insufficient_Data carries no slope and must carry no interval.
                if sel.get("slope_CI_lower") is not None or fresh is not None:
                    drifted.append(f"{record['ligand']['Name']}: interval without a slope")
                continue
            for key in ("slope_CI_lower", "slope_CI_upper", "slope_std"):
                if sel.get(key) != (fresh or {}).get(key):
                    drifted.append(
                        f"{record['ligand']['Name']}.{key}: "
                        f"{sel.get(key)!r} != {(fresh or {}).get(key)!r}")
        self.assertEqual(drifted, [], drifted[:5])

    def test_slope_lies_inside_its_own_interval(self):
        """Review 2026-08-21 finding 2: the CI used to describe a pooled
        regression while `slope` is the mean of per-condition fits, and eight
        ligands' slopes sat outside their own interval."""
        outside = []
        for record in self.records:
            sel = record.get("selectivity") or {}
            slope, lo, hi = sel.get("slope"), sel.get("slope_CI_lower"), sel.get("slope_CI_upper")
            if slope is None or lo is None or hi is None:
                continue
            if not (lo - 1e-9 <= slope <= hi + 1e-9):
                outside.append(f"{record['ligand']['Name']}: {slope} not in [{lo}, {hi}]")
        self.assertEqual(outside, [], outside)

    def test_stored_condition_statistics_match_fresh_computation(self):
        """The `statistics` block is a cache of the condition's own cells and
        must reproduce from them in the one canonical vocabulary
        (`rel_quest.utils.condition_statistics`). Until 2026-08-25 the corpus
        carried two shapes — ingest_liu's D_min/.../element_count on 178
        conditions, ingest_safe's mean_log_d/.../n_elements on 487."""
        from rel_quest.utils import condition_statistics

        drifted = []
        for record in self.records:
            for cond in record.get("experimental_data", []):
                fresh = condition_statistics(cond.get("distribution_ratios"))
                if cond.get("statistics") != fresh:
                    drifted.append(f"{record['ligand']['Name']}:"
                                   f"{cond.get('condition_id')}")
        self.assertEqual(drifted, [], drifted[:5])

    def test_stored_outlier_flags_match_fresh_detection(self):
        """Flags must come from the per-condition detector: a flag means
        'off its own condition's trend'. Until 2026-08-25 they came from a
        pooled regression, whose residuals on multi-condition records were
        between-condition offsets, not measurement error."""
        from rel_quest.uncertainty import detect_outliers

        drifted = []
        for record in self.records:
            stored = (record.get("data_quality") or {}).get("outlier_flags") or []
            fresh = detect_outliers(record.get("experimental_data", []))
            if stored != fresh:
                drifted.append(f"{record['ligand']['Name']}: "
                               f"{len(stored)} stored vs {len(fresh)} fresh")
        self.assertEqual(drifted, [], drifted[:5])

    def test_quality_score_reproduces_from_the_stored_interval_and_flags(self):
        """`compute_quality_score` spends 30 points on CI width, so a CI
        change that leaves the score behind is a stale record."""
        from rel_quest.uncertainty import compute_quality_score

        stale = []
        for record in self.records:
            sel = record.get("selectivity") or {}
            dq = record.get("data_quality") or {}
            ci = (None if sel.get("slope_CI_lower") is None else
                  {k: sel[k] for k in ("slope_CI_lower", "slope_CI_upper", "slope_std")})
            score, grade = compute_quality_score(record, ci, dq.get("outlier_flags") or [])
            if (dq.get("overall_quality_score"), dq.get("quality_grade")) != (score, grade):
                stale.append(f"{record['ligand']['Name']}: stored "
                             f"{dq.get('overall_quality_score')}/{dq.get('quality_grade')} "
                             f"vs {score}/{grade}")
        self.assertEqual(stale, [], stale[:5])


class TestEveryDerivedFieldReproduces(unittest.TestCase):
    """One check over every derived field (rel_quest/derived.py).

    The per-field tests above say *which* invariant broke; this one says
    nothing was missed. It is what `scripts/recompute_derived_fields.py
    --check` runs in CI, and it caught `performance_summary`: empty on 62
    records, stale on 34, and served verbatim by the REST API.
    """

    def test_recompute_is_a_fixed_point(self):
        from rel_quest.derived import derived_differences, recompute_record

        drifted = {}
        for record in load_records():
            diffs = derived_differences(record, recompute_record(record))
            if diffs:
                drifted[record["ligand"]["Name"]] = diffs[:3]
        self.assertEqual(drifted, {}, list(drifted.items())[:5])

    def test_performance_summary_is_uniform(self):
        from rel_quest.derived import performance_summary

        for record in load_records():
            with self.subTest(ligand=record["ligand"]["Name"]):
                self.assertEqual(record["performance_summary"], performance_summary(record))


class TestValidatorCatchesTheseDefects(unittest.TestCase):
    """The validator must fail these records, not just the tests above."""

    def _record(self, ratios, coverage=("La", "Ce", "Pr", "Nd")):
        from rel_quest.schema import SCHEMA_VERSION

        return {
            "schema_version": SCHEMA_VERSION,
            "ligand": {"Name": "T", "Mechanism": "Solvation",
                       "Extractant_Type": "Diglycolamide", "SMILES": "CCO"},
            "experimental_data": [{
                "condition_id": "c1",
                "reference": {"DOI": "10.1000/x", "Locator": None,
                              "Provenance_Level": "primary", "Aggregator": None},
                "distribution_ratios": ratios,
            }],
            "selectivity": {"trend_classification": "Flat",
                            "element_coverage": list(coverage)},
            "metadata": {"data_source": "TEST"},
            "performance_summary": {},
        }

    def test_zero_d_with_finite_log_is_rejected(self):
        from rel_quest.validators import validate_ligand_record

        record = self._record({"La": {"value": 0.0, "log_value": -4.0}})
        is_valid, errors = validate_ligand_record(record)
        self.assertFalse(is_valid)
        self.assertTrue(any("positive D" in e for e in errors), errors)

    def test_contradictory_d_and_log_is_rejected(self):
        from rel_quest.validators import validate_ligand_record

        record = self._record({"La": {"value": 10.0, "log_value": -4.0}})
        is_valid, errors = validate_ligand_record(record)
        self.assertFalse(is_valid)
        self.assertTrue(any("disagrees" in e for e in errors), errors)

    def test_coverage_naming_unmeasured_elements_is_rejected(self):
        from rel_quest.validators import validate_ligand_record

        record = self._record({"La": {"value": 1.0, "log_value": 0.0}},
                              coverage=("La", "Lu"))
        is_valid, errors = validate_ligand_record(record)
        self.assertFalse(is_valid)
        self.assertTrue(any("no measurement" in e for e in errors), errors)

    def test_wrong_elements_measured_is_rejected(self):
        from rel_quest.validators import validate_ligand_record

        record = self._record({"La": {"value": 1.0, "log_value": 0.0}},
                              coverage=("La",))
        record["elements_measured"] = ["La", "Eu"]
        is_valid, errors = validate_ligand_record(record)
        self.assertFalse(is_valid)
        self.assertTrue(any("unmeasured elements" in e for e in errors), errors)

    def test_a_consistent_record_passes(self):
        from rel_quest.validators import validate_ligand_record

        record = self._record({"La": {"value": 0.0001, "log_value": -4.0}},
                              coverage=("La",))
        record["elements_measured"] = ["La"]
        is_valid, errors = validate_ligand_record(record)
        self.assertTrue(is_valid, errors)


class TestRoundSignificant(unittest.TestCase):
    """The helper the ingest pipelines now use instead of decimal rounding."""

    def test_small_values_survive(self):
        from rel_quest.utils import round_significant

        self.assertEqual(round_significant(10 ** -4.0), 0.0001)
        self.assertGreater(round_significant(10 ** -5.5), 0)

    def test_round_trips_through_log_within_tolerance(self):
        from rel_quest.utils import round_significant

        for log_d in (-6.0, -4.4881, -2.757, -0.5, 0.0, 1.5, 3.9):
            with self.subTest(log_d=log_d):
                value = round_significant(10 ** log_d)
                self.assertAlmostEqual(math.log10(value), log_d, places=5)

    def test_handles_none_and_nan(self):
        from rel_quest.utils import round_significant

        self.assertIsNone(round_significant(None))
        self.assertIsNone(round_significant(float("nan")))


if __name__ == "__main__":
    unittest.main()
