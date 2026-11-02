"""Build the data-only release bundle for Zenodo (and the paper's Data Records).

Writes dist/reel-quest-data-v<version>/ and a zip of it:

  sx_ligands_merged_v14_flat.jsonl   the database (one record per line)
  reel_quest_measurements.csv        one row per D-value, with its condition,
                                     ligand identifiers and full provenance
  reel_quest_ligands.csv             one row per ligand: identity, class, trend
  schema/                            JSON Schema + data dictionary (units, provenance)
  vocab/                             controlled vocabularies
  splits/                            benchmark split definitions, if present
  DATABASE_STRUCTURE.md, DATA_PROVENANCE_AND_TERMS.md, LICENSE-DATA,
  CITATION.cff, references.bib, README.md (generated), zenodo_dataset.json (metadata for the
  manual dataset upload), SHA256SUMS

The zip is byte-reproducible: entries are sorted and carry a fixed timestamp.

--exclude-aggregator NAME drops every condition taken from that compilation
(`reference.Aggregator.name`, e.g. SAFE_01262026), removes records left with
no conditions, and recomputes every derived field of the rest with
rel_quest.derived.recompute_record — the fallback if a source's terms do not
allow redistribution (docs/DATA_PROVENANCE_AND_TERMS.md). The canonical
database file is never modified.

Run:
  python scripts/export_release_bundle.py [--exclude-aggregator SAFE_01262026] [--out dist]
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import shutil
import sys
import zipfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

import rel_quest  # noqa: E402
from rel_quest.derived import recompute_record  # noqa: E402
from rel_quest.schema import PROFILE  # noqa: E402
from rel_quest.validators import validate_against_json_schema, validate_ligand_record  # noqa: E402

DB_PATH = ROOT / "data" / "sx_ligands_merged_v14_flat.jsonl"
ATOMIC_NUMBER = {"Sc": 21, "Y": 39, "La": 57, "Ce": 58, "Pr": 59, "Nd": 60, "Pm": 61,
                 "Sm": 62, "Eu": 63, "Gd": 64, "Tb": 65, "Dy": 66, "Ho": 67, "Er": 68,
                 "Tm": 69, "Yb": 70, "Lu": 71}
FIXED_ZIP_TIME = (2026, 1, 1, 0, 0, 0)

MEASUREMENT_COLUMNS = [
    "ligand_name", "smiles", "inchi_key", "subtype", "mechanism", "condition_id",
    "acid_name", "acid_concentration_M", "extractant_concentration_M", "solvent",
    "solvent_ratio", "phase_modifier", "temperature_C", "contact_time", "element",
    "atomic_number", "D", "log_D", "sigma_log_D", "source_quality", "doi", "url",
    "locator", "provenance_level", "aggregator", "aggregator_entry", "curation_doi",
]
LIGAND_COLUMNS = [
    "ligand_name", "smiles", "inchi_key", "iupac_name", "subtype", "mechanism",
    "trend_classification", "trend_confidence", "slope", "slope_CI_lower",
    "slope_CI_upper", "r_squared", "peak_element", "n_conditions", "n_measurements",
    "elements_measured", "data_source", "quality_grade",
]


# Columns that name the compilation a value came through. The public release
# carries no source labels (rel_quest.schema.PUBLIC_FIELDS), so its CSVs omit them.
SOURCE_COLUMNS = {"provenance_level", "aggregator", "aggregator_entry", "curation_doi", "data_source"}
if PROFILE == "public":
    MEASUREMENT_COLUMNS = [c for c in MEASUREMENT_COLUMNS if c not in SOURCE_COLUMNS]
    LIGAND_COLUMNS = [c for c in LIGAND_COLUMNS if c not in SOURCE_COLUMNS]


def load_records():
    return [json.loads(line) for line in DB_PATH.open(encoding="utf-8") if line.strip()]


def exclude(records, names):
    """Drop conditions from the named compilations; recompute what remains."""
    if not names:
        return records, []
    kept, dropped = [], []
    for record in records:
        conditions = [c for c in record["experimental_data"]
                      if ((c["reference"].get("Aggregator") or {}).get("name")) not in names]
        if not conditions:
            dropped.append(record["ligand"]["Name"])
            continue
        if len(conditions) != len(record["experimental_data"]):
            record = dict(record, experimental_data=conditions)
            record = recompute_record(record)
        kept.append(record)
    return kept, dropped


def measurement_rows(records):
    for record in records:
        ligand = record["ligand"]
        for condition in record["experimental_data"]:
            reference = condition["reference"]
            aggregator = reference.get("Aggregator") or {}
            for element, cell in condition["distribution_ratios"].items():
                yield {
                    "ligand_name": ligand["Name"],
                    "smiles": ligand["SMILES"],
                    "inchi_key": ligand.get("InChI_Key"),
                    "subtype": ligand["Subtype"],
                    "mechanism": ligand["Mechanism"],
                    "condition_id": condition["condition_id"],
                    "acid_name": condition.get("Acid_Name"),
                    "acid_concentration_M": condition.get("Acid_Concentration_M"),
                    "extractant_concentration_M": condition.get("Extractant_Concentration_M"),
                    "solvent": condition.get("Solvent_Name"),
                    "solvent_ratio": condition.get("Solvent_Ratio"),
                    "phase_modifier": condition.get("Phase_Modifier_Name") or None,
                    "temperature_C": condition.get("Temperature_C"),
                    "contact_time": condition.get("Contact_Time"),
                    "element": element,
                    "atomic_number": ATOMIC_NUMBER[element],
                    "D": cell.get("value"),
                    "log_D": cell.get("log_value"),
                    "sigma_log_D": cell.get("sigma_log_D"),
                    "source_quality": cell.get("source_quality"),
                    "doi": reference.get("DOI"),
                    "url": reference.get("URL"),
                    "locator": reference.get("Locator"),
                    "provenance_level": reference.get("Provenance_Level"),
                    "aggregator": aggregator.get("name"),
                    "aggregator_entry": aggregator.get("entry"),
                    "curation_doi": aggregator.get("curation_DOI"),
                }


def ligand_rows(records):
    for record in records:
        ligand, sel = record["ligand"], record["selectivity"]
        yield {
            "ligand_name": ligand["Name"],
            "smiles": ligand["SMILES"],
            "inchi_key": ligand.get("InChI_Key"),
            "iupac_name": ligand.get("IUPAC_Name"),
            "subtype": ligand["Subtype"],
            "mechanism": ligand["Mechanism"],
            "trend_classification": sel["trend_classification"],
            "trend_confidence": sel["trend_confidence"],
            "slope": sel.get("slope"),
            "slope_CI_lower": sel.get("slope_CI_lower"),
            "slope_CI_upper": sel.get("slope_CI_upper"),
            "r_squared": sel.get("r_squared"),
            "peak_element": sel.get("peak_element"),
            "n_conditions": record["performance_summary"]["conditions_tested"],
            "n_measurements": record["performance_summary"]["total_measurements"],
            "elements_measured": ";".join(record["elements_measured"]),
            "data_source": record["metadata"].get("data_source"),
            "quality_grade": record["data_quality"]["quality_grade"],
        }


def write_csv(path, columns, rows):
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, lineterminator="\n", extrasaction="ignore")
        writer.writeheader()
        count = 0
        for row in rows:
            writer.writerow({k: ("" if v is None else v) for k, v in row.items()})
            count += 1
    return count


def readme(version, records, n_cells, excluded, dropped):
    conditions = sum(len(r["experimental_data"]) for r in records)
    lines = [
        f"# REEL-QUEST data release v{version}",
        "",
        "Curated liquid-liquid extraction distribution ratios for the lanthanides (La-Lu),",
        "Y and Sc, with per-measurement provenance and per-ligand selectivity trends.",
        "",
        f"- Ligands: {len(records)}",
        f"- Experimental conditions: {conditions}",
        f"- D-values: {n_cells}",
        "",
        "## Files",
        "",
        "| File | Content |",
        "|---|---|",
        "| `sx_ligands_merged_v14_flat.jsonl` | The database, one ligand record per line |",
        "| `reel_quest_measurements.csv` | One row per D-value with its condition and provenance |",
        "| `reel_quest_ligands.csv` | One row per ligand: identifiers, subtype, trend, counts |",
        "| `schema/ligand_record.schema.json` | JSON Schema (draft 2020-12) for a record |",
        "| `schema/data_dictionary.csv` | Every field with type, units and provenance |",
        "| `vocab/` | Controlled vocabularies (subtypes, mechanisms, solvents, ...) |",
        "| `DATABASE_STRUCTURE.md` | Schema documentation |",
        "| `DATA_PROVENANCE_AND_TERMS.md` | Sources, counts and third-party terms |",
        "| `references.bib` | Every paper the values come from, with the manuscript's other references (BibTeX) |",
        "| `SHA256SUMS` | Checksums of every file |",
        "",
        "## Provenance",
        "",
        "Every condition cites the paper that measured it (`reference.DOI`, or `doi` in the",
        "CSV), with the figure or table where recorded (`Locator`).",
    ]
    if PROFILE == "dev":
        lines += [
            "Values taken through a compilation name it (`Aggregator`): the Liu et al. 2022",
            "JACS Au Supporting Information, the SAFE database, or the OECD-NEA IDEaL database.",
        ]
    else:
        lines += [
            "The values were compiled from the Liu et al. 2022 JACS Au Supporting Information,",
            "the SAFE database, the OECD-NEA IDEaL database and papers read directly, and every",
            "condition was checked against the paper that measured it.",
        ]
    lines += [""]
    if excluded:
        lines += [
            "## Excluded sources",
            "",
            f"This bundle excludes every condition taken from: {', '.join(excluded)}.",
            f"{len(dropped)} ligands had no other conditions and are omitted; derived fields",
            "of the remaining records were recomputed from what remains.",
            "",
        ]
    lines += [
        "## Licence and citation",
        "",
        "CC BY 4.0 (see `LICENSE-DATA`). Cite the release (see `CITATION.cff`) and the papers",
        "behind the measurements you use (all in `references.bib`).",
        "",
        "Code: https://github.com/bo-li/REEL-QUEST (MIT).",
        "",
    ]
    return "\n".join(lines)


def zenodo_metadata(version, excluded):
    zenodo = json.loads((ROOT / ".zenodo.json").read_text(encoding="utf-8"))
    zenodo.update({
        "upload_type": "dataset",
        "license": "cc-by-4.0",
        "title": f"REEL-QUEST data release v{version}: lanthanide solvent-extraction "
                 f"distribution ratios with provenance and selectivity trends",
        "notes": "Data under CC BY 4.0. Code under MIT at https://github.com/bo-li/REEL-QUEST."
                 + (f" Excludes conditions from: {', '.join(excluded)}." if excluded else ""),
    })
    return zenodo


def build(out_dir: Path, excluded=()):
    version = rel_quest.__version__
    records, dropped = exclude(load_records(), set(excluded))
    problems = []
    for record in records:
        ok, errors = validate_ligand_record(record)
        errors = errors + validate_against_json_schema(record)
        if errors:
            problems.append((record["ligand"]["Name"], errors[:2]))
    if problems:
        raise SystemExit(f"refusing: {len(problems)} records fail validation: {problems[:3]}")

    suffix = "" if not excluded else "-without-" + "-".join(sorted(excluded))
    name = f"reel-quest-data-v{version}{suffix}"
    bundle = out_dir / name
    if bundle.exists():
        shutil.rmtree(bundle)
    (bundle / "schema").mkdir(parents=True)

    with (bundle / "sx_ligands_merged_v14_flat.jsonl").open("w", encoding="utf-8", newline="\n") as handle:
        for record in records:
            handle.write(json.dumps(record) + "\n")
    n_cells = write_csv(bundle / "reel_quest_measurements.csv", MEASUREMENT_COLUMNS,
                        measurement_rows(records))
    write_csv(bundle / "reel_quest_ligands.csv", LIGAND_COLUMNS, ligand_rows(records))
    for source in ("ligand_record.schema.json", "data_dictionary.csv"):
        shutil.copy2(ROOT / "schema" / source, bundle / "schema" / source)
    shutil.copytree(ROOT / "vocab", bundle / "vocab")
    if (ROOT / "data" / "splits").is_dir():
        shutil.copytree(ROOT / "data" / "splits", bundle / "splits")
    for source in ("docs/DATABASE_STRUCTURE.md", "docs/DATA_PROVENANCE_AND_TERMS.md",
                   "LICENSE-DATA", "CITATION.cff", "paper/references.bib"):
        shutil.copy2(ROOT / source, bundle / Path(source).name)
    (bundle / "README.md").write_text(readme(version, records, n_cells, sorted(excluded), dropped),
                                      encoding="utf-8")
    (bundle / "zenodo_dataset.json").write_text(
        json.dumps(zenodo_metadata(version, sorted(excluded)), indent=2) + "\n", encoding="utf-8")

    files = sorted(p for p in bundle.rglob("*") if p.is_file())
    sums = [f"{hashlib.sha256(p.read_bytes()).hexdigest()}  {p.relative_to(bundle).as_posix()}"
            for p in files]
    (bundle / "SHA256SUMS").write_text("\n".join(sums) + "\n", encoding="utf-8")

    archive = out_dir / f"{name}.zip"
    with zipfile.ZipFile(archive, "w", compression=zipfile.ZIP_DEFLATED) as zf:
        for path in sorted(p for p in bundle.rglob("*") if p.is_file()):
            info = zipfile.ZipInfo(f"{name}/{path.relative_to(bundle).as_posix()}", FIXED_ZIP_TIME)
            info.compress_type = zipfile.ZIP_DEFLATED
            info.external_attr = 0o644 << 16
            zf.writestr(info, path.read_bytes())
    return bundle, archive, records, n_cells, dropped


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--exclude-aggregator", action="append", default=[],
                        help="drop conditions from this compilation (repeatable)")
    parser.add_argument("--out", type=Path, default=ROOT / "dist")
    args = parser.parse_args()
    bundle, archive, records, n_cells, dropped = build(args.out, args.exclude_aggregator)
    conditions = sum(len(r["experimental_data"]) for r in records)
    print(f"Bundle:  {bundle}")
    print(f"Archive: {archive}  ({archive.stat().st_size / 1e6:.2f} MB, "
          f"sha256 {hashlib.sha256(archive.read_bytes()).hexdigest()[:16]}...)")
    print(f"Ligands {len(records)}, conditions {conditions}, D-values {n_cells}")
    if dropped:
        print(f"Omitted {len(dropped)} ligands with no remaining conditions")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
