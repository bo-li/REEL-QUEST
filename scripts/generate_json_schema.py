"""Generate the machine-readable schema and data dictionary for a record.

Writes, deterministically:

  schema/ligand_record.schema.json   JSON Schema (draft 2020-12) for one line of
                                     data/sx_ligands_merged_v14_flat.jsonl
  schema/data_dictionary.csv         one row per field: path, type, required,
                                     units, provenance, description
  schema/contribution.schema.json    JSON Schema for one contribution file
                                     (contributions/<paper>.json; rendered from
                                     rel_quest.contributions.contribution_schema)

Every enumeration is injected from the code or the vocabularies it already
reads (rel_quest.schema, rel_quest.references, vocab/*.json), so the schema
cannot drift from the validator. `tests/test_json_schema.py` checks that the
files on disk equal this script's output, that every stored record validates,
and that every field in the data appears in the dictionary.

Each leaf carries two annotation keywords (ignored by validators):
  x-units       physical unit, or "dimensionless" / "log10" / "dex"
  x-provenance  source-transcribed | curated | derived: <function> | computed-at-ingest

Run:
  python scripts/generate_json_schema.py [--check]

In the public repository (rel_quest/_release.py PROFILE = 'public') the schema
and dictionary describe only rel_quest.schema.PUBLIC_FIELDS (`public_view`).
"""

from __future__ import annotations

import argparse
import copy
import csv
import io
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from rel_quest.references import DOI_PATTERN  # noqa: E402
from rel_quest.schema import (  # noqa: E402
    PROFILE, PROVENANCE_LEVELS, PUBLIC_FIELDS, SCHEMA_VERSION, TREND_CLASSIFICATIONS,
    canonical_mechanisms)

SCHEMA_PATH = ROOT / "schema" / "ligand_record.schema.json"
DICTIONARY_PATH = ROOT / "schema" / "data_dictionary.csv"
CONTRIBUTION_PATH = ROOT / "schema" / "contribution.schema.json"
SCHEMA_ID = "https://github.com/bo-li/REEL-QUEST/blob/main/schema/ligand_record.schema.json"
PUBLIC_SCHEMA_ID = "https://github.com/bo-li/REEL-QUEST/blob/main/schema/ligand_record.schema.json"
# Descriptions that name a compilation, reworded for the public schema.
PUBLIC_DESCRIPTIONS = {
    ("experimental_data[]", "condition_id"): "Unique, stable condition identifier",
    ("experimental_data[]", "Acid_Name"): "Aqueous acid (mostly HNO3; 'Hcl' is hydrochloric acid "
                                          "as spelled in the source)",
}

ST, CUR = "source-transcribed", "curated"


def vocab(name):
    return json.loads((ROOT / "vocab" / f"{name}.json").read_text(encoding="utf-8"))


SUBTYPES = sorted(vocab("chemical_families")["subtypes"])
LANTHANIDES = ["La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy",
               "Ho", "Er", "Tm", "Yb", "Lu"]
ELEMENTS = LANTHANIDES + ["Y", "Sc"]
DATA_SOURCES = ["Liu_2022_JACSAu", "SAFE_01262026", "OECD_NEA_IDEAL", "Pramanik_2024_JACS", "Literature_curated"]
AGGREGATOR_NAME_PATTERN = "^(Liu_2022_JACSAu|OECD_NEA_IDEAL|SAFE_\\d{8})$"
DATE = "^\\d{4}-\\d{2}-\\d{2}$"


def leaf(types, description, units=None, provenance=ST, **extra):
    """A documented leaf. `types` is a JSON type or list of types."""
    node = {"type": types, "description": description}
    node.update(extra)
    if units:
        node["x-units"] = units
    node["x-provenance"] = provenance
    return node


def obj(properties, required=(), additional=False, description=None, nullable=False):
    node = {"type": ["object", "null"] if nullable else "object"}
    if description:
        node["description"] = description
    node["properties"] = properties
    if required:
        node["required"] = list(required)
    node["additionalProperties"] = additional
    return node


def build_schema():
    element = {"type": "string", "enum": ELEMENTS}
    lanthanide = {"type": "string", "enum": LANTHANIDES}

    cell = obj({
        "value": leaf(["number", "null"], "Distribution ratio D = [M]org/[M]aq, 6 significant figures "
                      "(null only where the source gave no usable value)", "dimensionless",
                      exclusiveMinimum=0),
        "log_value": leaf(["number", "null"], "log10(D) as reported or digitised; agrees with "
                          "`value` within 0.01 dex (validators.check_d_value_consistency)", "log10"),
        "n_replicates": leaf("integer", "Number of source values behind this cell", minimum=1),
        "sigma_log_D": leaf("number", "Uncertainty of log D (1 sigma). When `sigma_source` is "
                            "'measured', computed for this value (e.g. from the figure's marker size "
                            "and axis scale); otherwise a prior by source type (uncertainty.SIGMA_PRIORS: "
                            "0.0715 dex for figure-digitised, the median reading uncertainty of the "
                            "calibrated figures; 3.5 % relative, 0.0152 dex, for tabulated)",
                            "dex", "derived: uncertainty.classify_measurement_quality + "
                            "compute_sigma_log_d, unless sigma_source is 'measured' (source data)",
                            minimum=0),
        "sigma_source": leaf("string", "Where sigma_log_D comes from: 'prior' by source type, or "
                             "'measured' for this value; absent means 'prior'", enum=["prior", "measured"]),
        "value_origin": leaf("string", "How the value was obtained: 'reported' (a table or the "
                             "text), 'digitised' (read off a figure marker), 'interpolated' (read off "
                             "a fitted line or between points, or at an inferred condition), 'derived' "
                             "(computed from other quantities), 'censored' (only a bound is known; see "
                             "`censor`). Absent means reported or digitised, per source_quality",
                             enum=["reported", "digitised", "interpolated", "derived", "censored"]),
        "censor": leaf(["string", "null"], "For a censored value, the direction of the bound: '>' "
                       "(D exceeds `value`) or '<' (D is below it); `value`/`log_value` hold the "
                       "bound. Censored cells are left out of trend fits, statistics and "
                       "separation factors", enum=[">", "<", None]),
        "source_quality": leaf("string", "Whether the value was read from a figure or a table",
                               provenance="derived: uncertainty.classify_measurement_quality",
                               enum=["figure", "table"]),
    }, required=("value", "log_value", "n_replicates", "sigma_log_D", "source_quality"))
    # A censored origin needs its bound direction; a bound needs the censored origin.
    cell["allOf"] = [
        {"if": {"properties": {"value_origin": {"const": "censored"}}, "required": ["value_origin"]},
         "then": {"properties": {"censor": {"enum": [">", "<"]}}, "required": ["censor"]}},
        {"if": {"properties": {"censor": {"enum": [">", "<"]}}, "required": ["censor"]},
         "then": {"properties": {"value_origin": {"const": "censored"}}}},
    ]

    statistics = obj({
        "mean_log_d": leaf("number", "Mean log D over the condition's elements", "log10",
                           "derived: utils.condition_statistics"),
        "std_log_d": leaf("number", "Population standard deviation of log D", "dex",
                          "derived: utils.condition_statistics"),
        "min_log_d": leaf("number", "Smallest log D", "log10", "derived: utils.condition_statistics"),
        "max_log_d": leaf("number", "Largest log D", "log10", "derived: utils.condition_statistics"),
        "n_elements": leaf("integer", "Number of elements with a log D", provenance=
                           "derived: utils.condition_statistics", minimum=0),
    }, description="Cache of the condition's own distribution_ratios, measured values only (censored "
                   "bounds are left out); empty when every value of the condition is a bound")
    # Either every statistic, or none: a condition of bounds only has no measured value to summarise.
    statistics["anyOf"] = [{"maxProperties": 0},
                           {"required": ["mean_log_d", "std_log_d", "min_log_d", "max_log_d", "n_elements"]}]

    aggregator = obj({
        "name": leaf("string", "Compilation the value was taken from (matches "
                     "metadata.data_source vocabulary)", pattern=AGGREGATOR_NAME_PATTERN),
        "DOI": leaf(["string", "null"], "The compilation's own DOI", pattern=DOI_PATTERN.pattern),
        "URL": leaf(["string", "null"], "The compilation's web address", format="uri"),
        "entry": leaf(["string", "null"], "The row's identifier inside the compilation: Liu SI "
                      "reference number, SAFE exp_id, or IDEaL '<name> ref [n]'"),
        "curation_DOI": leaf(["string", "null"], "Publication the compilation curated this row "
                             "for (SAFE: 10.1021/jacs.5c19738, SAFE-MolGen). Not a measurement "
                             "source", pattern=DOI_PATTERN.pattern),
    }, required=("name", "DOI", "URL", "entry", "curation_DOI"), nullable=True,
        description="Set for aggregator rows, null for primary rows")

    reference = obj({
        "DOI": leaf(["string", "null"], "Bare, lower-case DOI of the paper that measured the "
                    "values; never a compilation DOI", pattern=DOI_PATTERN.pattern),
        "URL": leaf("string", "Source link when the paper has no DOI", format="uri"),
        "Locator": leaf(["string", "null"], "Figure or table in that paper, e.g. 'Figure 3', "
                        "'Table S2'", provenance="curated (parsed by references.parse_locator)"),
        "Provenance_Level": leaf("string", "primary: transcribed by this project from the paper; "
                                 "aggregator: taken from a compilation citing it",
                                 provenance=CUR, enum=list(PROVENANCE_LEVELS)),
        "Aggregator": aggregator,
        "Source_Type": leaf("string", "Legacy source label (superseded by Provenance_Level)",
                            enum=["Experimental", "Database"]),
        "Data_Used": leaf("string", "Legacy free-text description of what was taken"),
        "Reference_Number": leaf(["integer", "null"], "Liu 2022 SI reference number (Liu rows)"),
        "Notes": leaf("string", "Source free text; for SAFE rows the full comments_description "
                      "(data location, title, authors, year)"),
    }, required=("DOI", "Locator", "Provenance_Level", "Aggregator"),
        description="Where the condition's numbers came from (schema v0.3, "
                    "rel_quest.references.upgrade_reference)")
    reference["allOf"] = [
        {"if": {"properties": {"DOI": {"type": "null"}}},
         "then": {"required": ["URL"]}},
        {"if": {"properties": {"Provenance_Level": {"const": "aggregator"}}},
         "then": {"properties": {"Aggregator": {"type": "object"}}}},
        {"if": {"properties": {"Provenance_Level": {"const": "primary"}}},
         "then": {"properties": {"Aggregator": {"type": "null"}}}},
    ]

    condition = obj({
        "condition_id": leaf("string", "Unique condition identifier; its prefix names the source "
                             "(Liu2022_, SAFE_, IDEAL_, Pramanik2024_)", provenance=CUR),
        "Acid_Name": leaf(["string", "null"], "Aqueous acid (mostly HNO3; SAFE spells hydrochloric "
                          "acid 'Hcl')"),
        "Acid_Concentration_M": leaf(["number", "null"], "Initial aqueous acid concentration",
                                     "mol L-1", minimum=0),
        "Temperature_C": leaf("number", "Temperature", "degC"),
        "Solvent_Name": leaf("string", "Organic diluent"),
        "Solvent_Ratio": leaf(["string", "null"], "Volume ratio for mixed diluents, e.g. '7:3'"),
        "Extractant_Concentration_M": leaf(["number", "null"], "Extractant concentration in the "
                                           "organic phase", "mol L-1", minimum=0),
        "Phase_Modifier_Name": leaf("string", "Organic phase modifier (empty string if none)"),
        "Phase_Modifier_Concentration_M": leaf(["number", "null"], "Phase modifier concentration",
                                               "mol L-1"),
        "Contact_Time": leaf("string", "Contact/equilibration time with unit, e.g. '60 min'"),
        "Additives": leaf("string", "Other aqueous components, free text"),
        "Co_Extractants": {
            "type": ["array", "null"],
            "description": "Other extractants in the organic phase (a mixed solvent). The D belongs to the "
                           "mixture: the condition stays on the main ligand's record but is left out of "
                           "everything pooled across conditions (uncertainty.is_mixture)",
            "items": obj({
                "Name": leaf("string", "Co-extractant name"),
                "SMILES": leaf(["string", "null"], "Co-extractant structure"),
                "Concentration_M": leaf(["number", "null"], "Co-extractant concentration in the organic "
                                        "phase", "mol L-1", minimum=0),
            }, required=["Name"]),
        },
        "data_source": leaf("string", "Per-condition source label where it differs from the "
                            "record's metadata.data_source"),
        "distribution_ratios": {
            "type": "object",
            "description": "D-values keyed by element symbol",
            "propertyNames": element,
            "additionalProperties": cell,
            "minProperties": 1,
        },
        "statistics": statistics,
        "reference": reference,
    }, required=("condition_id", "distribution_ratios", "statistics", "reference"))

    selectivity = obj({
        "trend_classification": leaf("string", "Selectivity trend across La-Lu (see "
                                     "docs/DATABASE_STRUCTURE.md)", provenance=
                                     "derived: trends.classify_selectivity_trend",
                                     enum=list(TREND_CLASSIFICATIONS)),
        "trend_confidence": leaf("string", "Agreement of per-condition votes: High (unanimous), "
                                 "Medium (>70 %), Low; 'None' for Insufficient_Data",
                                 provenance="derived: trends.classify_selectivity_trend",
                                 enum=["High", "Medium", "Low", "None"]),
        "element_coverage": {"type": "array", "items": lanthanide, "uniqueItems": True,
                             "description": "La-Lu elements the trend fit used (never Y/Sc)",
                             "x-provenance": "derived: trends.classify_selectivity_trend"},
        "slope": leaf(["number", "null"], "Mean of per-condition slopes of log D vs atomic number",
                      "dex per unit Z", "derived: trends.classify_selectivity_trend"),
        "r_squared": leaf(["number", "null"], "Mean per-condition R^2 of that fit",
                          "dimensionless", "derived: trends.classify_selectivity_trend",
                          minimum=0, maximum=1),
        "notes": leaf("string", "Plain-language basis of the classification",
                      provenance="derived: trends.classify_selectivity_trend"),
        "peak_element": {**lanthanide, "description": "Position of the extremum for Peaked / "
                         "U_Shaped", "x-provenance": "derived: trends.classify_selectivity_trend"},
        "slope_CI_lower": leaf(["number", "null"], "Lower bound of the 95 % two-stage bootstrap "
                               "interval of `slope`", "dex per unit Z",
                               "derived: uncertainty.bootstrap_trend_ci (seeded)"),
        "slope_CI_upper": leaf(["number", "null"], "Upper bound of that interval",
                               "dex per unit Z", "derived: uncertainty.bootstrap_trend_ci (seeded)"),
        "slope_std": leaf(["number", "null"], "Bootstrap standard deviation of `slope`",
                          "dex per unit Z", "derived: uncertainty.bootstrap_trend_ci (seeded)",
                          minimum=0),
    }, required=("trend_classification", "trend_confidence", "element_coverage", "slope",
                 "r_squared", "notes", "slope_CI_lower", "slope_CI_upper", "slope_std"),
        description="Selectivity trend across the lanthanide series, fitted from the record's "
                    "own conditions")

    outlier = obj({
        "condition_id": leaf("string", "Condition holding the flagged value",
                             provenance="derived: uncertainty.detect_outliers"),
        "element": {**element, "description": "Element of the flagged value",
                    "x-provenance": "derived: uncertainty.detect_outliers"},
        "log_D_measured": leaf("number", "Stored log D", "log10", "derived: uncertainty.detect_outliers"),
        "log_D_predicted": leaf("number", "Value on the condition's own linear fit", "log10",
                                "derived: uncertainty.detect_outliers"),
        "residual": leaf("number", "Measured minus predicted", "dex",
                         "derived: uncertainty.detect_outliers"),
        "z_score": leaf("number", "Residual z-score across the record's eligible conditions",
                        "dimensionless", "derived: uncertainty.detect_outliers"),
        "severity": leaf("string", "|z| > 4 High, > 3 Medium",
                         provenance="derived: uncertainty.detect_outliers", enum=["High", "Medium"]),
    }, required=("condition_id", "element", "log_D_measured", "log_D_predicted", "residual",
                 "z_score", "severity"))

    data_quality = obj({
        "overall_quality_score": leaf("integer", "Uncalibrated ordinal index 0-100 (CI width 30, "
                                      "coverage 30, measurement quality 20, outliers 20)",
                                      provenance="derived: uncertainty.compute_quality_score",
                                      minimum=0, maximum=100),
        "quality_grade": leaf("string", "Grade from the score", provenance=
                              "derived: uncertainty.compute_quality_score",
                              enum=["High", "Medium", "Low"]),
        "outlier_flags": {"type": "array", "items": outlier,
                          "description": "Measurements off their own condition's linear trend",
                          "x-provenance": "derived: uncertainty.detect_outliers"},
    }, required=("overall_quality_score", "quality_grade", "outlier_flags"),
        description="Outlier flags and an uncalibrated quality index")

    performance = obj({
        "conditions_tested": leaf("integer", "Number of conditions", provenance=
                                  "derived: derived.performance_summary", minimum=0),
        "total_measurements": leaf("integer", "Number of D-values", provenance=
                                   "derived: derived.performance_summary", minimum=0),
        "elements_tested": {"type": "array", "items": element,
                            "description": "Same list as elements_measured",
                            "x-provenance": "derived: derived.performance_summary"},
    }, required=("conditions_tested", "total_measurements", "elements_tested"),
        description="Cache of counts (rel_quest.derived.performance_summary)")

    ligand = obj({
        "Name": leaf("string", "Curated short name; unique across the database", provenance=CUR,
                     minLength=1),
        "SMILES": leaf("string", "RDKit-canonical SMILES; the primary key (no duplicates)",
                       provenance=CUR, minLength=1),
        "InChI": leaf("string", "Standard InChI", provenance="derived at ingest (RDKit)"),
        "InChI_Key": leaf("string", "Standard InChIKey", provenance="derived at ingest (RDKit)",
                          pattern="^[A-Z]{14}-[A-Z]{10}-[A-Z]$"),
        "IUPAC_Name": leaf(["string", "null"], "IUPAC name", provenance=CUR),
        "CAS_Number": leaf(["string", "null"], "CAS Registry Number", provenance=CUR,
                           pattern="^\\d{2,7}-\\d{2}-\\d$"),
        "PubChem_CID": leaf(["integer", "null"], "PubChem compound ID", provenance=CUR, minimum=1),
        "Molecular_Formula": leaf("string", "Hill formula", provenance="derived at ingest (RDKit)"),
        "Molecular_Weight": leaf(["number", "null"], "Average molecular weight", "g mol-1",
                                 "derived at ingest (RDKit)", exclusiveMinimum=0),
        "Mechanism": leaf("string", "Extraction mechanism (vocab/mechanisms.json)",
                          provenance=CUR, enum=list(canonical_mechanisms())),
        "Subtype": leaf("string", "Chemical subtype, 21-value flat taxonomy "
                        "(vocab/chemical_families.json)", provenance=CUR, enum=SUBTYPES),
        "Chemical_Family": leaf("string", "Same label as Subtype (kept for compatibility)",
                                provenance=CUR, enum=SUBTYPES),
        "Extractant_Type": leaf("string", "Legacy SMARTS class; prefer Subtype", provenance=CUR),
        "Structural_Class": leaf(["string", "null"], "Deprecated legacy free text from the "
                                 "sources; uncontrolled, use Subtype", deprecated=True),
        "Family_Name": leaf("string", "Source grouping key from the Liu SI; not a chemical claim"),
        "Alternative_Names": {"type": "array", "items": {"type": "string"}, "uniqueItems": True,
                              "description": "Genuine synonyms only", "x-provenance": CUR},
        "Community_Name": leaf("string", "Informal family nickname; uncurated, not for analysis",
                               provenance=CUR),
        "Community_Name_Source": leaf("string", "How Community_Name was assigned", provenance=CUR),
        "Name_Source": leaf("string", "How IUPAC_Name was obtained", provenance=CUR),
        "Name_Confidence": leaf("string", "Confidence in IUPAC_Name", provenance=CUR),
        "Name_Provenance": obj({
            "doi": leaf("string", "Paper the name was taken from"),
            "paper_compound_code": leaf("string", "Compound code in that paper"),
            "pdf": leaf(["string", "null"], "Local PDF path used"),
            "section": leaf("string", "Section of the paper"),
        }, description="Where a PDF-extracted name came from"),
        "Classification_Source": leaf("string", "How Subtype was assigned", provenance=CUR),
        "Classification_Confidence": leaf("string", "Confidence of an automatic Subtype",
                                          provenance=CUR),
        "Classification_Reason": leaf("string", "Rule that assigned an automatic Subtype",
                                      provenance=CUR),
        "Extraction_Method": leaf("string", "How the structure was extracted from the source",
                                  provenance=CUR),
        "Origin": leaf("string", "'LLM_generated' marks ligands proposed by an LLM design study "
                       "(no experimental origin)", provenance=CUR),
        "Notes": leaf("string", "Curator notes", provenance=CUR),
    }, required=("Name", "SMILES", "Mechanism", "Extractant_Type", "Subtype", "Chemical_Family",
                 "Alternative_Names"),
        description="Identity, structure, names and classification of the extractant")

    small = lambda props, req, text: obj(props, required=req, description=text)  # noqa: E731
    metadata = obj({
        "data_source": leaf("string", "Source that first contributed the ligand",
                            provenance=CUR, enum=DATA_SOURCES),
        "added_date": leaf("string", "Date the record was added", provenance=CUR, pattern=DATE),
        "last_updated": leaf("string", "Date the record's content last changed", provenance=CUR,
                             pattern=DATE),
        "curation_status": leaf("string", "Curation state", provenance=CUR,
                                enum=["Validated", "Imported", "Needs_Review", "Incomplete",
                                      "Deprecated"]),
        "curator": leaf("string", "Who or what curated the record", provenance=CUR),
        "taxonomy_version": leaf("string", "Classification generation (independent of "
                                 "schema_version)", provenance=CUR),
        "reclassified_at": leaf("string", "When the v13 taxonomy was applied", provenance=CUR),
        "prior_classification": small({
            "Subtype": leaf("string", "Subtype before the v13 reclassification"),
            "Chemical_Family": leaf("string", "Chemical_Family before the v13 reclassification"),
        }, ("Subtype", "Chemical_Family"), "Classification before the v13 taxonomy"),
        "merged_from": {"type": "array", "description": "Lineage of the v11 SMILES dedup",
                        "items": obj({
                            "record": leaf("string", "Merged record identifier"),
                            "name": leaf("string", "Its name before the merge"),
                            "n_experimental_data": leaf("integer", "Conditions it contributed"),
                            "merged_on": leaf("string", "Merge date", pattern=DATE),
                            "reason": leaf("string", "Why it was merged"),
                        }, required=("record", "name", "n_experimental_data", "merged_on"))},
        "name_corrections": {"type": "array", "description": "Renames resolving name collisions",
                             "items": obj({
                                 "previous_name": leaf("string", "Old name"),
                                 "new_name": leaf("string", "New name"),
                                 "collided_with": leaf("string", "Record it collided with"),
                                 "date": leaf("string", "Date", pattern=DATE),
                                 "reason": leaf("string", "Why"),
                             }, required=("previous_name", "new_name", "date"))},
        "subtype_corrections": {"type": "array", "description": "Manual Subtype corrections",
                                "items": obj({
                                    "from": leaf("string", "Old Subtype"),
                                    "to": leaf("string", "New Subtype"),
                                    "date": leaf("string", "Date", pattern=DATE),
                                    "reason": leaf("string", "Why"),
                                }, required=("from", "to", "date"))},
        "source_variant": leaf("string", "Liu SI variant tag (X_vN), provenance only"),
        "data_sources": {"type": "array", "items": {"type": "string"},
                         "description": "Legacy list of contributing sources (v6 enrichment)"},
        "enrichment_status": leaf("string", "Legacy v6 enrichment marker"),
        "ideal_data_summary": obj({
            "has_references": leaf("boolean", "IDEaL entry lists references"),
            "total_d_values": leaf("integer", "D-values listed in IDEaL"),
            "total_separation_factors": leaf("integer", "Separation factors listed in IDEaL"),
            "total_decontamination_factors": leaf("integer", "Decontamination factors in IDEaL"),
        }, description="Counts reported by the IDEaL entry at scrape time"),
        "pubchem_cid": leaf("null", "Legacy placeholder (use ligand.PubChem_CID)"),
        "rdkit_descriptors_available": leaf("boolean", "Legacy ingest flag"),
        "validation_notes": leaf("string", "Curator notes on validation"),
    }, required=("data_source", "added_date", "last_updated", "curation_status", "curator",
                 "taxonomy_version"),
        description="Curation state, dates and lineage")

    return {
        "$schema": "https://json-schema.org/draft/2020-12/schema",
        "$id": SCHEMA_ID,
        "title": "REEL-QUEST ligand record",
        "description": (f"One line of data/sx_ligands_merged_v14_flat.jsonl, schema {SCHEMA_VERSION}. "
                        "Generated by scripts/generate_json_schema.py; see "
                        "docs/DATABASE_STRUCTURE.md for prose and schema/data_dictionary.csv "
                        "for units and provenance."),
        "type": "object",
        "properties": {
            "schema_version": {"const": SCHEMA_VERSION, "description": "Record-shape version"},
            "ligand": ligand,
            "selectivity": selectivity,
            "experimental_data": {"type": "array", "items": condition, "minItems": 1,
                                  "description": "Experimental conditions"},
            "performance_summary": performance,
            "metadata": metadata,
            "data_quality": data_quality,
            "elements_measured": {"type": "array", "items": element, "uniqueItems": True,
                                  "description": "Every element with a measurement, incl. Y/Sc",
                                  "x-provenance": "derived: utils.collect_elements_measured"},
        },
        "required": ["schema_version", "ligand", "selectivity", "experimental_data",
                     "performance_summary", "metadata", "data_quality", "elements_measured"],
        "additionalProperties": False,
    }


def _rule_keys(rule):
    """Property names an if/then rule mentions."""
    keys = set()
    for part in (rule.get("if") or {}, rule.get("then") or {}):
        keys |= set(part.get("properties") or {}) | set(part.get("required") or [])
    return keys


def public_view(schema):
    """The schema of the public release: only rel_quest.schema.PUBLIC_FIELDS.

    Levels the whitelist does not name (selectivity, data_quality, cells, ...)
    ship whole; every listed key must exist in the full schema.
    """
    view = copy.deepcopy(schema)
    condition = view["properties"]["experimental_data"]["items"]
    nodes = {"": view, "ligand": view["properties"]["ligand"],
             "metadata": view["properties"]["metadata"], "experimental_data[]": condition,
             "experimental_data[].reference": condition["properties"]["reference"]}
    if set(nodes) != set(PUBLIC_FIELDS):
        raise ValueError(f"PUBLIC_FIELDS levels {sorted(PUBLIC_FIELDS)} != {sorted(nodes)}")
    for path, keep in PUBLIC_FIELDS.items():
        node = nodes[path]
        unknown = set(keep) - set(node["properties"])
        if unknown:
            raise ValueError(f"PUBLIC_FIELDS[{path!r}] names fields the schema lacks: {sorted(unknown)}")
        node["properties"] = {k: v for k, v in node["properties"].items() if k in keep}
        if "required" in node:
            node["required"] = [k for k in node["required"] if k in keep]
        if "allOf" in node:
            node["allOf"] = [rule for rule in node["allOf"] if _rule_keys(rule) <= set(keep)]
    for (path, key), text in PUBLIC_DESCRIPTIONS.items():
        nodes[path]["properties"][key]["description"] = text
    view["$id"] = PUBLIC_SCHEMA_ID
    return view


def dictionary_rows(node, path="", required=True):
    """Flatten the schema into (path, type, required, units, provenance, description)."""
    rows = []
    kind = node.get("type")
    kind_text = "|".join(kind) if isinstance(kind, list) else (kind or "")
    if path:
        rows.append((path, kind_text, "yes" if required else "no", node.get("x-units", ""),
                     node.get("x-provenance", ""), node.get("description", "")))
    props = node.get("properties") or {}
    for name, child in props.items():
        rows += dictionary_rows(child, f"{path}.{name}" if path else name,
                                name in (node.get("required") or []))
    items = node.get("items")
    if isinstance(items, dict) and items.get("properties"):
        for name, child in items["properties"].items():
            rows += dictionary_rows(child, f"{path}[].{name}", name in (items.get("required") or []))
    extra = node.get("additionalProperties")
    if isinstance(extra, dict) and extra.get("properties"):
        for name, child in extra["properties"].items():
            rows += dictionary_rows(child, f"{path}.<element>.{name}",
                                    name in (extra.get("required") or []))
    return rows


def render(profile=PROFILE):
    schema = build_schema()
    if profile == "public":
        schema = public_view(schema)
    schema_text = json.dumps(schema, indent=2, ensure_ascii=False) + "\n"
    buffer = io.StringIO()
    writer = csv.writer(buffer, lineterminator="\n")
    writer.writerow(["path", "type", "required", "units", "provenance", "description"])
    writer.writerows(dictionary_rows(schema))
    return schema_text, buffer.getvalue()


def render_contribution() -> str:
    """schema/contribution.schema.json as text (the module owns the schema; this writes it)."""
    from rel_quest.contributions import contribution_schema
    return json.dumps(contribution_schema(), indent=2, ensure_ascii=False) + "\n"


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--check", action="store_true",
                        help="exit 1 if the files on disk differ from a fresh render")
    args = parser.parse_args()
    schema_text, dictionary_text = render()
    contribution_text = render_contribution()
    if args.check:
        stale = [p.relative_to(ROOT) for p, text in ((SCHEMA_PATH, schema_text),
                                                      (DICTIONARY_PATH, dictionary_text),
                                                      (CONTRIBUTION_PATH, contribution_text))
                 if not p.exists() or p.read_text(encoding="utf-8") != text]
        if stale:
            print(f"stale: {stale} — run scripts/generate_json_schema.py")
            return 1
        print("schema, data dictionary and contribution schema are current")
        return 0
    SCHEMA_PATH.parent.mkdir(parents=True, exist_ok=True)
    SCHEMA_PATH.write_text(schema_text, encoding="utf-8")
    DICTIONARY_PATH.write_text(dictionary_text, encoding="utf-8")
    CONTRIBUTION_PATH.write_text(contribution_text, encoding="utf-8")
    print(f"wrote {SCHEMA_PATH.relative_to(ROOT)}, {DICTIONARY_PATH.relative_to(ROOT)} "
          f"({dictionary_text.count(chr(10)) - 1} fields) and {CONTRIBUTION_PATH.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
