"""Record-schema version for the ligand database.

`schema_version` describes the **shape of a record** — which fields exist and
what they hold. It is not the taxonomy generation: that is
`metadata.taxonomy_version` (e.g. "v13_flat"), which tracks how ligands are
*classified* and changes independently.

    schema_version          record structure      this module
    metadata.taxonomy_version  classification     vocab/chemical_families.json

Kept in one place because it was previously hardcoded as "v0.1" in six separate
writers, which is precisely why it never got bumped as the schema grew.

Versioning scheme (see docs/DATABASE_STRUCTURE.md):
    v0.x  developing; additive changes expected
    v1.x  stable; backwards-compatible additions only
    v2.x  breaking changes

History:
    v0.1  original Liu 2022 ingest: ligand / selectivity / experimental_data /
          performance_summary / metadata
    v0.2  current. Adds the top-level `data_quality` block; naming fields on
          `ligand` (IUPAC_Name, InChI, InChI_Key, CAS_Number, PubChem_CID,
          Community_Name, Name_Source, Name_Confidence, Name_Provenance);
          classification fields (Subtype, Chemical_Family, Family_Name,
          Classification_Source); selectivity confidence intervals
          (slope_CI_lower, slope_CI_upper, slope_std, peak_element); and
          metadata provenance (taxonomy_version, prior_classification,
          merged_from, source_variant, name_corrections). Distribution ratios
          gained sigma_log_D and source_quality.
    v0.3  current. Structured reference provenance on every condition
          (rel_quest/references.py): `reference.DOI` is a bare lower-case DOI
          of the measuring paper (it held comma-joined URL pairs on 214
          conditions); new `URL`, `Locator`, `Provenance_Level` and
          `Aggregator` {name, DOI, URL, entry, curation_DOI}. Legacy keys
          Source_Type, Data_Used, Reference_Number and Notes are unchanged.
"""

SCHEMA_VERSION = "v0.3"

# Versions a record may legitimately declare. Older versions are retained so
# archived generations under data/archive/ still validate against their own
# schema; checks introduced later apply only from the version that added them.
KNOWN_SCHEMA_VERSIONS = ("v0.1", "v0.2", "v0.3")


def schema_at_least(version, minimum):
    """True if schema `version` (e.g. "v0.3") is `minimum` or newer."""
    def parts(v):
        try:
            return tuple(int(x) for x in str(v).lstrip("v").split("."))
        except ValueError:
            return (0,)
    return parts(version) >= parts(minimum)


# Where a condition's numbers came from (reference.Provenance_Level, v0.3+).
#   primary     transcribed by this project from the measuring paper itself
#   aggregator  taken from a compilation (Liu 2022 SI, SAFE, IDEAL) that
#               cites the measuring paper; see reference.Aggregator
PROVENANCE_LEVELS = ("primary", "aggregator")


# ---------------------------------------------------------------------------
# Public release profile
# ---------------------------------------------------------------------------
# The public repository (scripts/export_public_repo.py) ships every record
# projected onto these fields. Each condition keeps the paper that measured it
# (reference.DOI + Locator); nothing names the compilation it came through
# (metadata.data_source, reference.Aggregator, the condition-id prefix, curator
# notes), which stays in this developer repository. A whitelist, so a field
# added later stays private until someone decides it is public. Keys are paths
# into the record; a level not listed (selectivity, data_quality, cells, ...)
# is derived and ships whole.
from ._release import PROFILE  # noqa: E402  ('dev' here, 'public' in the export)

PUBLIC_FIELDS = {
    "": ("schema_version", "ligand", "selectivity", "experimental_data",
         "performance_summary", "metadata", "data_quality", "elements_measured"),
    "ligand": ("Name", "SMILES", "InChI", "InChI_Key", "IUPAC_Name", "CAS_Number",
               "PubChem_CID", "Molecular_Formula", "Molecular_Weight", "Mechanism",
               "Subtype", "Chemical_Family", "Extractant_Type", "Alternative_Names",
               "Community_Name"),
    "metadata": ("added_date", "last_updated", "taxonomy_version"),
    "experimental_data[]": ("condition_id", "Acid_Name", "Acid_Concentration_M",
                            "Temperature_C", "Solvent_Name", "Solvent_Ratio",
                            "Extractant_Concentration_M", "Phase_Modifier_Name",
                            "Phase_Modifier_Concentration_M", "Contact_Time", "Additives",
                            "Co_Extractants", "distribution_ratios", "statistics",
                            "reference"),
    "experimental_data[].reference": ("DOI", "URL", "Locator"),
}


def project_public(record, rename=None):
    """`record` reduced to PUBLIC_FIELDS, key order kept; `rename` maps condition ids.

    scripts/export_public_repo.py passes the developer->public id map; the
    contribution check in the public repository projects with ids unchanged.
    """
    import copy

    rename = rename or (lambda condition_id: condition_id)

    def keep(block, allowed):
        return {k: copy.deepcopy(v) for k, v in block.items() if k in allowed}

    out = keep(record, PUBLIC_FIELDS[""])
    out["ligand"] = keep(record["ligand"], PUBLIC_FIELDS["ligand"])
    out["metadata"] = keep(record["metadata"], PUBLIC_FIELDS["metadata"])
    conditions = []
    for condition in record["experimental_data"]:
        public = keep(condition, PUBLIC_FIELDS["experimental_data[]"])
        public["condition_id"] = rename(condition["condition_id"])
        public["reference"] = keep(condition["reference"], PUBLIC_FIELDS["experimental_data[].reference"])
        conditions.append(public)
    out["experimental_data"] = conditions
    for flag in (out.get("data_quality") or {}).get("outlier_flags") or []:
        flag["condition_id"] = rename(flag["condition_id"])
    return out


# ---------------------------------------------------------------------------
# Controlled vocabularies
# ---------------------------------------------------------------------------
# Kept here for the same reason as SCHEMA_VERSION above: these lists were
# copied into the validator, the CLI, the dashboard, the plotting module and
# the agent tool schema, and they drifted. When the peak-position heuristic
# was added to the linear trend fit and the vocabulary gained Peaked and
# U_Shaped, five of six copies were
# updated. The sixth was `scripts/query_ligands.py --trend`, which went on
# rejecting two valid classifications covering 18 ligands — queryable through
# the Python API but unreachable from the command line.

TREND_CLASSIFICATIONS = (
    "Ascending",
    "Reverse",
    "Flat",
    "Mixed",
    "Peaked",
    "U_Shaped",
    "Insufficient_Data",
)

# Display order for charts and pickers: monotone increase, through the two
# curved shapes, to monotone decrease, then the non-committal labels. This is
# a presentation concern, so it also carries "Unknown" for records with no
# stored classification. Do not use it for validation.
TREND_DISPLAY_ORDER = (
    "Ascending",
    "Peaked",
    "Flat",
    "U_Shaped",
    "Reverse",
    "Mixed",
    "Insufficient_Data",
    "Unknown",
)

# Plain-language reading of each classification, for API responses and UI.
TREND_DIRECTIONS = {
    "Ascending": "La→Lu increasing (light to heavy REE preference)",
    "Reverse": "Lu→La increasing (heavy to light REE preference)",
    "Flat": "No clear selectivity trend across lanthanide series",
    "Mixed": "Multiple conflicting trends across different conditions",
    "Peaked": "Maximum at a mid-series element (see selectivity.peak_element)",
    "U_Shaped": "Minimum at a mid-series element (see selectivity.peak_element)",
    "Insufficient_Data": "Not enough data points to determine trend",
}


# Phrasings a researcher is likely to type for each v13 Subtype, beyond the
# label itself and the subtype's `representative` compound (both of which
# subtype_aliases() derives automatically from vocab/chemical_families.json).
#
# These are trade names, abbreviations and older family words — the things a
# question actually contains. "Show me phosphine oxide ligands" has to reach
# Phosphine_oxide, and "what about CMPO" has to reach it too, because CMPO was
# a v12 family of its own before the P-topology rule folded it in.
_EXTRA_SUBTYPE_ALIASES = {
    "DGA": ["diglycolamide", "tehdga", "tddga", "tdddga", "dodga", "dmdodga"],
    "Tripodal-DGA": ["tripodal", "tripodal dga"],
    "Phosphine_oxide": ["phosphine oxide", "cmpo", "topo", "amide phosphoryl",
                        "carbamoylmethylphosphine", "organophosphorus"],
    "Phosphonate": ["phosphonate", "amidophosphonate", "dhdecmp", "dhdecep"],
    "Phosphate": ["tbp", "tributyl phosphate"],
    "Phosphoric_acid": ["phosphoric acid", "hdehp", "d2ehpa", "dhpa"],
    "Malonamide": ["malonamide", "dmdohema", "diamide"],
    "Butanediamide": ["butanediamide", "succinamide", "tbbda"],
    "Amide": ["monoamide", "dhoa"],
    "DAPhen": ["daphen", "phenanthroline diamide", "phenanthroline amide"],
    "MAPhen": ["maphen", "phenanthroline monoamide"],
    "BLPhen": ["blphen", "bi-loop phenanthroline"],
    "Multident_N_donor": ["btp", "btbp", "btphen", "bpyrphen", "multident",
                          "n-donor", "n donor", "polyaza"],
    "Bipydiamide": ["bipydiamide", "bipyridine diamide", "bipyridyl"],
    "Pyridinediamide": ["pyridinediamide", "pyridine diamide", "dipicolinamide"],
    "Nitrilotriacetamide": ["nitrilotriacetamide", "honta", "n-tripod", "tripod chelator"],
    "EDTA": ["edta", "hedta", "cydta", "aminopolycarboxylate", "dta"],
    "Crown_ether": ["crown ether", "crown", "18-crown-6"],
    "Polyether_Diamide": ["polyether diamide", "dooda", "polyether"],
    "Polyphenol": ["polyphenol", "quercetin", "catechol", "calixarene"],
    "Phthalamic_acid": ["phthalamic", "phthalamic acid"],
}


def subtype_aliases():
    """Map each v13 Subtype to the phrases that should resolve to it.

    Built from vocab/chemical_families.json — the subtype label itself (with
    underscores read as spaces) and its `representative` compound — plus the
    curated phrasings above. Deriving the base from the vocabulary means a
    new subtype is reachable the moment it is added, rather than whenever
    someone remembers to update a keyword table.

    That is the failure this replaces: the AI assistant matched on a
    hardcoded map of v12 two-layer family names ("Diglycolamide",
    "Organophosphorus", "Amide_Phosphoryl", ...) while the database had moved
    to flat v13 subtypes. Ten of its thirteen entries matched nothing, so
    "show me phosphine oxide ligands" injected "ORGANOPHOSPHORUS LIGANDS
    (0 total)" against 19 real records.

    Returns:
        dict of {Subtype: [alias, ...]}, aliases lowercased and deduplicated,
        longest first so a specific phrase wins over a substring of it.
    """
    import json
    import os

    path = os.path.join(
        os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
        "vocab", "chemical_families.json",
    )
    with open(path) as handle:
        subtypes = json.load(handle)["subtypes"]

    index = {}
    for subtype, entry in subtypes.items():
        aliases = {subtype.lower(), subtype.lower().replace("_", " ")}
        representative = (entry or {}).get("representative")
        if representative:
            aliases.add(str(representative).lower())
        aliases.update(a.lower() for a in _EXTRA_SUBTYPE_ALIASES.get(subtype, []))
        index[subtype] = sorted(aliases, key=lambda a: (-len(a), a))
    return index


def canonical_mechanisms():
    """The mechanism vocabulary, read from vocab/mechanisms.json.

    A function rather than a constant so the JSON file stays the single
    source of truth; callers that need a static list for an argparse
    `choices=` can call this at parser-construction time.
    """
    import json
    import os

    path = os.path.join(
        os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
        "vocab", "mechanisms.json",
    )
    with open(path) as handle:
        return tuple(json.load(handle)["canonical_mechanisms"])


def write_skip_log(skipped_records, output_path, source):
    """Persist skipped records so they survive past the console.

    Ingest pipelines used to only print their skips, which meant a dropped
    record vanished into scrollback. Writing a sibling log next to the output
    keeps the omission auditable — and, since the file is written only when
    something was actually skipped, its *absence* is meaningful evidence that
    nothing was.

    Args:
        skipped_records: list of dicts describing each skip (must say why)
        output_path: the ingest output path; the log is written alongside it
        source: label for the pipeline, recorded on every entry

    Returns:
        Path to the log, or None if there was nothing to write.
    """
    import json
    from datetime import date
    from pathlib import Path

    if not skipped_records:
        return None

    log_path = Path(output_path).with_name(
        Path(output_path).stem + "_skipped.jsonl"
    )
    log_path.parent.mkdir(parents=True, exist_ok=True)
    with log_path.open("w") as handle:
        for record in skipped_records:
            handle.write(json.dumps(
                {"source": source, "logged": date.today().isoformat(), **record}
            ) + "\n")
    print(f"  {len(skipped_records)} skipped record(s) logged to: {log_path}")
    return log_path
