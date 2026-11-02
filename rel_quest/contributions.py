"""One contribution = one paper = one file, applied to the released database.

contributions/<doi-slug>.json (schema/contribution.schema.json, rendered from
`contribution_schema()` by scripts/generate_json_schema.py):

  DOI           the paper that measured the values, bare ("10.1021/jacs.4c07332"),
                never a compilation's (SAFE, SAFE-MolGen: rel_quest.references)
  citation      one line, as the paper cites itself
  proposed_by   who read the paper, when and how ("A. Chemist, 2026-11-02, Table 2")
  reviewed_by   who read the paper against the file before it was applied (the review
                rule, GOVERNANCE.md §4); scripts/apply_contributions.py refuses a file
                without it unless told --unreviewed
  new_ligands   ligands no record names yet: Name, SMILES, Mechanism, Extractant_Type,
                Subtype (vocab/chemical_families.json), reason; optionally
                Chemical_Family (= Subtype), IUPAC_Name, Alternative_Names,
                Structural_Class, Family_Name, Classification_Reason
  conditions    one per experiment: id (A0001 ...), ligand (an existing Name or one of
                new_ligands), conditions (schema fields: Acid_Name, Acid_Concentration_M,
                Extractant_Concentration_M, Solvent_Name, Temperature_C, Contact_Time,
                Phase_Modifier_Name, Phase_Modifier_Concentration_M, Solvent_Ratio,
                Additives), additives (free text, appended), co_extractants
                ([{Name, SMILES, Concentration_M}] for a mixed organic phase),
                reference {Locator: "Table 2" | "Figure 3" | "Text p. 5"}, cells
                {El: {log_D | D, value_origin, sigma_log_D?, n_replicates?} or
                {censor: ">" | "<", log_bound}}, why, note
  status        written by scripts/apply_contributions.py: "applied", with applied_on,
                applied_to and condition_ids

A cell's `sigma_log_D` is given only when the source states one (a printed
error, a calibrated figure); otherwise rel_quest.derived.recompute_record
assigns the prior by source type. A bound (censor) stays out of every fit,
statistic and separation factor. Each condition becomes a `Lit_<ligand>_<doi>_<id>`
condition with primary provenance (rel_quest.references.upgrade_reference) on
its ligand's record; a new ligand's record is created the first time a
condition names it, with its identity fields computed from the SMILES.

This is the grammar scripts/apply_verification.py used for the literature
additions of v14 (outputs/verification/additions/), so the v14 rebuild and the
live route share this module: the cell and condition editors, the resolution
grammar, the structure index and the dating rule live here and nowhere else.
"""

from __future__ import annotations

import copy
import json
import math
import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Callable, Dict, List, Optional, Tuple

from .derived import recompute_record
from .references import DOI_PATTERN, NON_MEASUREMENT_DOIS, normalize_doi, upgrade_reference
from .schema import PROFILE, SCHEMA_VERSION, canonical_mechanisms, project_public
from .uncertainty import SIGMA_PRIORS
from .utils import load_vocabulary, round_significant
from .validators import validate_against_json_schema, validate_ligand_record

ELEMENTS = ("La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Y", "Sc")
TABLE_PRIOR = SIGMA_PRIORS["table"]
FIGURE_PRIOR = SIGMA_PRIORS["figure"]
# Condition fields a resolution or a contribution may set (schema/ligand_record.schema.json)
RESOLVABLE_FIELDS = {"Acid_Name", "Acid_Concentration_M", "Temperature_C", "Solvent_Name", "Solvent_Ratio",
                     "Extractant_Concentration_M", "Phase_Modifier_Name", "Phase_Modifier_Concentration_M",
                     "Contact_Time", "Additives"}
NUMERIC_FIELDS = {"Acid_Concentration_M", "Temperature_C", "Extractant_Concentration_M",
                  "Phase_Modifier_Concentration_M"}
NOT_AMENDED = ("move_to", "merge_into", "recover", "parts")
ADDITION_ID = re.compile(r"^A\d{4}$")
LITERATURE = "Literature_curated"
VALUE_ORIGINS = ("reported", "digitised", "interpolated", "derived")
STATUS_APPLIED = "applied"
KNOWN_TOP_KEYS = {"DOI", "citation", "proposed_by", "reviewed_by", "new_ligands", "conditions",
                  "status", "applied_on", "applied_to", "condition_ids"}
NEW_LIGAND_REQUIRED = ("Name", "SMILES", "Mechanism", "Extractant_Type", "Subtype", "reason")
NEW_LIGAND_OPTIONAL = ("Chemical_Family", "IUPAC_Name", "Alternative_Names", "Structural_Class", "Family_Name",
                       "Classification_Reason", "CAS_Number", "PubChem_CID")
# What goes from a declaration into ligand.* (the rest is curation metadata)
NEW_LIGAND_PRIVATE = ("reason", "proposed_by")
DATE_STAMPS = ("last_updated", "merged_on", "added_date")
RECOVERED_ON = re.compile(r"\bRecovered \d{4}-\d{2}-\d{2}:")


# --- cell and condition editors ---------------------------------------------------


def set_cell(cells: dict, element: str, log_d: float, **extra) -> None:
    cell = dict(cells.get(element) or {"n_replicates": 1, "sigma_log_D": TABLE_PRIOR, "source_quality": "figure"})
    cell.update(log_value=round(float(log_d), 4), value=round_significant(10 ** float(log_d)))
    for k in ("censor",):
        cell.pop(k, None)
    cell.update(extra)
    cells[element] = cell


def set_bound(cells: dict, element: str, direction: str, log_bound: float) -> None:
    cell = dict(cells.get(element) or {"n_replicates": 1, "sigma_log_D": TABLE_PRIOR, "source_quality": "figure"})
    cell.update(log_value=round(float(log_bound), 4), value=round_significant(10 ** float(log_bound)),
                censor=direction, value_origin="censored", sigma_source="prior")
    cells[element] = cell


def add_additive(condition: dict, text: str) -> None:
    have = condition.get("Additives") or ""
    if text not in have:
        condition["Additives"] = f"{have}; {text}" if have else text


def set_reference(condition: dict, doi=None, locator=None) -> None:
    """Change a condition's source through rel_quest.references (DOI normalised, aggregator kept)."""
    ref = dict(condition["reference"])
    if doi:
        ref["DOI"] = doi
        ref.pop("URL", None)
    upgraded = upgrade_reference(ref, condition["condition_id"])
    if locator:
        upgraded["Locator"] = locator    # curator-verified; upgrade_reference would re-derive it from the notes
    elif ref.get("Locator"):
        upgraded["Locator"] = ref["Locator"]
    condition["reference"] = upgraded


# --- the resolution grammar (a cell, a condition field, an additive, a reference) ----


def resolution_problems(spec: dict) -> list[str]:
    """What makes a resolution unusable as written (checked before anything is applied)."""
    bad = []
    if spec.get("action") not in ("replace", "split", "drop", "resolved", "amend"):
        bad.append(f"unknown action {spec.get('action')!r}")
    if spec.get("action") == "amend" and any(spec.get(k) for k in NOT_AMENDED):
        bad.append("an amend edits the condition in place: no " + ", ".join(NOT_AMENDED))
    if not spec.get("why"):
        bad.append("no 'why'")
    parts = spec.get("parts") or []
    if spec.get("action") == "split":
        ids = [p.get("id") for p in parts]
        if len(parts) < 2 or len(set(ids)) != len(ids) or not all(ids):
            bad.append("a split needs two or more parts with distinct ids")
    for part in [spec] + parts:
        for field_ in part.get("conditions") or {}:
            if field_ not in RESOLVABLE_FIELDS:
                bad.append(f"condition field {field_!r} cannot be set")
        for el, c in (part.get("cells") or {}).items():
            if el not in ELEMENTS:
                bad.append(f"unknown element {el!r}")
            elif c is not None and not (c.get("censor") in (">", "<") and isinstance(c.get("log_bound"), (int, float))
                                        or isinstance(c.get("log_D"), (int, float))
                                        and isinstance(c.get("sigma_log_D", 0.0), (int, float)) and c.get("value_origin")):
                bad.append(f"{el}: a cell needs log_D and value_origin (and sigma_log_D when the source states one), "
                           f"or censor and log_bound")
    return bad


def apply_resolution(condition: dict, spec: dict, own_cells_only: bool = False) -> list[str]:
    """Write one resolution (or split part) into a condition; returns what changed."""
    changes = []
    cells = {} if own_cells_only else condition["distribution_ratios"]
    for el, c in (spec.get("cells") or {}).items():
        if c is None:
            if cells.pop(el, None) is not None:
                changes.append(f"{el} removed")
            continue
        if c.get("censor") in (">", "<"):
            set_bound(cells, el, c["censor"], c["log_bound"])
            changes.append(f"{el}{c['censor']}{c['log_bound']}")
        elif c.get("sigma_log_D") is not None:
            set_cell(cells, el, c["log_D"], sigma_log_D=round(float(c["sigma_log_D"]), 4), sigma_source="measured",
                     value_origin=c["value_origin"])
            changes.append(f"{el}={c['log_D']}±{c['sigma_log_D']}")
        else:
            # the source states no uncertainty: the condition's prior applies (recompute_record sets it)
            set_cell(cells, el, c["log_D"], value_origin=c["value_origin"])
            cells[el].pop("sigma_source", None)
            changes.append(f"{el}={c['log_D']}")
        if c.get("n_replicates"):
            cells[el]["n_replicates"] = int(c["n_replicates"])
    condition["distribution_ratios"] = cells
    for field_, value in (spec.get("conditions") or {}).items():
        condition[field_] = value
        changes.append(f"{field_}={value}")
    if spec.get("additives"):
        add_additive(condition, spec["additives"])
        changes.append(f"Additives+={spec['additives']}")
    if spec.get("co_extractants"):
        condition["Co_Extractants"] = spec["co_extractants"]
        changes.append("Co_Extractants=" + ", ".join(x["Name"] for x in spec["co_extractants"]))
    ref = spec.get("reference") or {}
    if ref:
        set_reference(condition, ref.get("DOI"), ref.get("Locator"))
        changes.append(f"reference={condition['reference'].get('DOI')} {condition['reference'].get('Locator')}")
    if spec.get("note"):
        notes = condition["reference"].get("Notes") or ""
        condition["reference"]["Notes"] = f"{notes}; {spec['note']}" if notes else spec["note"]
    return changes


def addition_problems(paper: dict, a: dict) -> list[str]:
    """What makes one condition of a paper file unusable as written (no database needed)."""
    bad = [f"the file has no {k}" for k in ("DOI", "proposed_by") if not paper.get(k)]
    if not ADDITION_ID.match(str(a.get("id") or "")):
        bad.append(f"id {a.get('id')!r} is not A followed by four digits")
    if not (a.get("reference") or {}).get("Locator"):
        bad.append("no Locator")
    if not a.get("cells"):
        bad.append("no cells")
    cells, cell_bad = normalise_cells(a.get("cells") or {})
    bad += cell_bad
    spec = {"action": "replace", "why": a.get("why"), "cells": cells, "conditions": a.get("conditions")}
    bad += resolution_problems(spec)
    for field_, value in (a.get("conditions") or {}).items():
        if field_ in NUMERIC_FIELDS and not isinstance(value, (int, float)):
            bad.append(f"{field_} must be a number, not {value!r}")
    return bad


def normalise_cells(cells: dict) -> tuple[dict, list[str]]:
    """Cells with `D` turned into `log_D` (a contributor may give either); problems found doing so."""
    out, bad = {}, []
    for el, c in cells.items():
        if not isinstance(c, dict):
            out[el] = c
            continue
        c = dict(c)
        if "D" in c:
            d = c.pop("D")
            if not isinstance(d, (int, float)) or d <= 0:
                bad.append(f"{el}: D must be a positive number, not {d!r}")
            else:
                log_d = round(math.log10(d), 4)
                if isinstance(c.get("log_D"), (int, float)) and abs(c["log_D"] - log_d) > 0.01:
                    bad.append(f"{el}: D {d} and log_D {c['log_D']} disagree (log10 D = {log_d})")
                c.setdefault("log_D", log_d)
        out[el] = c
    return out, bad


def cell_signature(condition: dict) -> tuple:
    # A cell the paper marks as not determined stores no value (Pramanik 2024, Table S2).
    return tuple(sorted((el, None if c.get("log_value") is None else round(c["log_value"], 3), c.get("censor") or "")
                        for el, c in condition["distribution_ratios"].items()))


# --- structures ---------------------------------------------------------------------


def inchikey(smiles: str):
    from rdkit import Chem
    mol = Chem.MolFromSmiles(smiles) if smiles else None
    return Chem.MolToInchiKey(mol) if mol else None


def structure_index(records) -> tuple[dict, dict]:
    """({InChIKey: Name}, {connectivity layer: [Name, ...]}); the first record with a key wins."""
    by_key, by_conn = {}, {}
    for r in records:
        key = inchikey(r["ligand"]["SMILES"])
        if key and key not in by_key:
            by_key[key] = r["ligand"]["Name"]
            by_conn.setdefault(key.split("-")[0], []).append(r["ligand"]["Name"])
    return by_key, by_conn


def record_for(smiles: str, by_key: dict, by_conn: dict):
    """The record holding this structure: the full InChIKey (stereo included) first, so that
    stereoisomers (RS/SS-DMeTODGA share connectivity) stay apart; else the connectivity, if unique."""
    key = inchikey(smiles)
    if not key:
        return None
    if key in by_key:
        return by_key[key]
    names = by_conn.get(key.split("-")[0], [])
    return names[0] if len(names) == 1 else None


def new_ligand_record(lig: dict, curator: str, today: str, data_source: Optional[str] = None) -> dict:
    """An empty record for a declared ligand; identity fields come from the SMILES, derived fields
    from recompute_record once conditions have moved in. Raises ValueError on an unparseable SMILES."""
    from rdkit import Chem
    from rdkit.Chem.Descriptors import MolWt
    from rdkit.Chem.rdMolDescriptors import CalcMolFormula
    mol = Chem.MolFromSmiles(lig["SMILES"])
    if mol is None:
        raise ValueError(f"{lig.get('Name')}: unparseable SMILES {lig['SMILES']!r}")
    ligand = {k: v for k, v in lig.items() if k not in NEW_LIGAND_PRIVATE}
    ligand.setdefault("Chemical_Family", ligand.get("Subtype"))
    ligand.update(InChI=Chem.MolToInchi(mol), InChI_Key=Chem.MolToInchiKey(mol),
                  Molecular_Formula=CalcMolFormula(mol), Molecular_Weight=round(MolWt(mol), 1),
                  Alternative_Names=lig.get("Alternative_Names", []), Name_Source="curator",
                  Name_Confidence="high", Classification_Source="curator")
    return {
        "schema_version": SCHEMA_VERSION, "ligand": ligand,
        "selectivity": {}, "performance_summary": {}, "data_quality": {}, "elements_measured": [],
        "metadata": {"data_source": data_source, "added_date": today, "last_updated": today,
                     "curation_status": "Validated", "curator": curator,
                     "validation_notes": lig.get("reason", ""), "taxonomy_version": "v13_flat"},
        "experimental_data": []}


# --- dating -------------------------------------------------------------------------


def undated(obj):
    """obj without the build's own date stamps (what decides whether content changed)."""
    if isinstance(obj, dict):
        return {k: undated(v) for k, v in obj.items() if k not in DATE_STAMPS}
    if isinstance(obj, list):
        return [undated(x) for x in obj]
    if isinstance(obj, str):
        return RECOVERED_ON.sub("Recovered:", obj)
    return obj


def keep_released_dates(fresh: dict, released: Optional[dict], today: str) -> dict:
    """The released record when the content is unchanged; else fresh, re-dated today, with the
    released dates of its unchanged conditions and merge entries."""
    if released is None:
        return fresh
    if undated(fresh) == undated(released):
        return released
    old = {c["condition_id"]: c for c in released["experimental_data"]}
    fresh["experimental_data"] = [old[c["condition_id"]] if c["condition_id"] in old
                                  and undated(c) == undated(old[c["condition_id"]]) else c
                                  for c in fresh["experimental_data"]]
    metadata = dict(fresh["metadata"], last_updated=today)
    if released["metadata"].get("added_date"):
        metadata["added_date"] = released["metadata"]["added_date"]
    old_merges = released["metadata"].get("merged_from") or []
    if metadata.get("merged_from"):
        metadata["merged_from"] = [next((o for o in old_merges if undated(o) == undated(m)), m)
                                   for m in metadata["merged_from"]]
    fresh["metadata"] = metadata
    return fresh


# --- the contribution file ----------------------------------------------------------


def doi_slug(doi: str) -> str:
    return re.sub(r"[^A-Za-z0-9.]+", "_", doi)


def condition_id(ligand: str, doi: str, aid: str) -> str:
    return f"Lit_{ligand}_{doi_slug(doi)}_{aid}"


def load_contribution(path) -> dict:
    return json.loads(Path(path).read_text(encoding="utf-8"))


def save_contribution(path, paper: dict) -> None:
    Path(path).write_text(json.dumps(paper, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")


def pending_files(directory) -> list[Path]:
    """The files of a contributions directory not yet applied, in name order."""
    out = []
    for path in sorted(Path(directory).glob("*.json")):
        try:
            paper = load_contribution(path)
        except (OSError, ValueError):
            out.append(path)          # unreadable: the check will say so
            continue
        if isinstance(paper, dict) and paper.get("status") != STATUS_APPLIED:
            out.append(path)
    return out


def load_records(path) -> list[dict]:
    with Path(path).open(encoding="utf-8") as h:
        return [json.loads(line) for line in h if line.strip()]


def write_records(path, records) -> None:
    with Path(path).open("w", encoding="utf-8", newline="\n") as h:
        for r in records:
            h.write(json.dumps(r) + "\n")


def contribution_schema() -> dict:
    """JSON Schema (draft 2020-12) for one contribution file."""
    number = {"type": "number"}
    string = {"type": "string"}
    nullable_string = {"type": ["string", "null"]}
    condition_fields = {f: (number if f in NUMERIC_FIELDS else string) for f in sorted(RESOLVABLE_FIELDS)}
    cell = {
        "type": "object",
        "description": "A measured value (log_D or D, with value_origin, and sigma_log_D only when the source "
                       "states one) or a bound (censor and log_bound)",
        "properties": {
            "log_D": {"type": "number", "description": "log10 of the distribution ratio"},
            "D": {"type": "number", "exclusiveMinimum": 0, "description": "The distribution ratio itself; "
                  "either log_D or D"},
            "value_origin": {"type": "string", "enum": list(VALUE_ORIGINS),
                             "description": "reported (a table or the text), digitised (read off a figure), "
                                            "interpolated, derived"},
            "sigma_log_D": {"type": "number", "minimum": 0, "description": "1 sigma of log D, only when the "
                            "source states it (a printed error, a calibrated figure)"},
            "n_replicates": {"type": "integer", "minimum": 1},
            "censor": {"type": "string", "enum": [">", "<"], "description": "The value is a bound: D exceeds "
                       "(>) or is below (<) 10**log_bound"},
            "log_bound": {"type": "number"},
        },
        "additionalProperties": False,
    }
    co_extractant = {
        "type": "object",
        "properties": {"Name": string, "SMILES": nullable_string, "Concentration_M": {"type": ["number", "null"]}},
        "required": ["Name"], "additionalProperties": False,
    }
    condition = {
        "type": "object",
        "properties": {
            "id": {"type": "string", "pattern": ADDITION_ID.pattern, "description": "A0001, A0002, ... unique "
                   "in the file; part of the condition_id"},
            "ligand": {"type": "string", "description": "An existing record's Name, or one of new_ligands"},
            "conditions": {"type": "object", "properties": condition_fields, "additionalProperties": False,
                           "description": "Schema condition fields (units: M, degC)"},
            "additives": nullable_string,
            "co_extractants": {"type": ["array", "null"], "items": co_extractant},
            "reference": {"type": "object", "properties": {"Locator": string}, "required": ["Locator"],
                          "additionalProperties": False,
                          "description": "Where in the paper the values are: 'Table 2', 'Figure 3', 'Text p. 5'"},
            "cells": {"type": "object", "propertyNames": {"enum": list(ELEMENTS)},
                      "additionalProperties": cell, "minProperties": 1},
            "why": {"type": "string", "description": "Why these values and conditions: what the paper says"},
            "note": nullable_string,
        },
        "required": ["id", "ligand", "reference", "cells", "why"],
        "additionalProperties": False,
    }
    new_ligand = {
        "type": "object",
        "properties": {
            "Name": string, "SMILES": string,
            "Mechanism": {"type": "string", "enum": list(canonical_mechanisms())},
            "Extractant_Type": {"type": "string", "description": "Legacy SMARTS class label, e.g. Diglycolamide"},
            "Subtype": {"type": "string", "enum": sorted(load_vocabulary("chemical_families")["subtypes"])},
            "Chemical_Family": {"type": "string", "description": "Same value as Subtype"},
            "IUPAC_Name": nullable_string, "Alternative_Names": {"type": "array", "items": string},
            "Structural_Class": nullable_string, "Family_Name": nullable_string,
            "Classification_Reason": nullable_string, "CAS_Number": nullable_string,
            "PubChem_CID": {"type": ["integer", "string", "null"]},
            "reason": {"type": "string", "description": "Where the paper identifies the molecule (formula, "
                       "name, table)"},
        },
        "required": list(NEW_LIGAND_REQUIRED),
        "additionalProperties": False,
    }
    return {
        "$schema": "https://json-schema.org/draft/2020-12/schema",
        "$id": "https://github.com/bo-li/REEL-QUEST/blob/main/schema/contribution.schema.json",
        "title": "REEL-QUEST contribution: one paper's measurements",
        "description": "One file per paper under contributions/, checked by scripts/check_contribution.py and "
                       "applied by scripts/apply_contributions.py (rel_quest/contributions.py)",
        "type": "object",
        "properties": {
            "DOI": {"type": "string", "pattern": DOI_PATTERN.pattern,
                    "description": "The paper that measured the values; never a compilation's DOI"},
            "citation": string,
            "proposed_by": {"type": "string", "description": "Who read the paper, when, and from what"},
            "reviewed_by": {"type": "string", "description": "Who read the paper against this file before it "
                            "was applied (GOVERNANCE.md §4); required by apply_contributions.py"},
            "new_ligands": {"type": "array", "items": new_ligand},
            "conditions": {"type": "array", "items": condition, "minItems": 1},
            "status": {"type": "string", "enum": [STATUS_APPLIED]},
            "applied_on": {"type": "string", "pattern": "^\\d{4}-\\d{2}-\\d{2}$"},
            "applied_to": string,
            "condition_ids": {"type": "array", "items": string},
        },
        "required": ["DOI", "citation", "proposed_by", "conditions"],
        "additionalProperties": True,
    }


def contribution_problems(paper: dict) -> list[str]:
    """What makes a file unusable before the database is even looked at: the schema, the DOI,
    the ids, the new-ligand declarations, each condition's grammar."""
    from jsonschema import Draft202012Validator
    bad = []
    if not isinstance(paper, dict):
        return ["the file is not a JSON object"]
    for err in sorted(Draft202012Validator(contribution_schema()).iter_errors(paper), key=lambda e: list(e.path)):
        where = "/".join(str(p) for p in err.path) or "file"
        bad.append(f"schema: {where}: {err.message[:160]}")
    doi = paper.get("DOI")
    if isinstance(doi, str):
        if normalize_doi(doi) is None:
            bad.append(f"DOI {doi!r} is not a DOI (10.xxxx/...)")
        elif normalize_doi(doi) != doi:
            bad.append(f"DOI must be bare and lower-case: {normalize_doi(doi)!r}")
        elif doi in NON_MEASUREMENT_DOIS:
            bad.append(f"DOI {doi} is a compilation's ({NON_MEASUREMENT_DOIS[doi]}); cite the paper that "
                       f"measured the values")
    ids = [str(a.get("id")) for a in paper.get("conditions") or [] if isinstance(a, dict)]
    for aid in sorted({i for i in ids if ids.count(i) > 1}):
        bad.append(f"id {aid} is used by more than one condition")
    declared = [lig for lig in paper.get("new_ligands") or [] if isinstance(lig, dict)]
    names = [lig.get("Name") for lig in declared]
    for name in sorted({n for n in names if names.count(n) > 1}, key=str):
        bad.append(f"new ligand {name!r} is declared twice")
    for lig in declared:
        if lig.get("SMILES") and inchikey(lig["SMILES"]) is None:
            bad.append(f"new ligand {lig.get('Name')}: SMILES {lig['SMILES']!r} does not parse")
        if lig.get("Chemical_Family") not in (None, lig.get("Subtype")):
            bad.append(f"new ligand {lig.get('Name')}: leave Chemical_Family out; it is set to Subtype "
                       f"({lig.get('Subtype')!r}), not {lig.get('Chemical_Family')!r}")
        used = any(isinstance(a, dict) and a.get("ligand") == lig.get("Name") for a in paper.get("conditions") or [])
        if not used:
            bad.append(f"new ligand {lig.get('Name')} has no condition in this file")
    for a in paper.get("conditions") or []:
        if isinstance(a, dict):
            bad += [f"{a.get('id')}: {p}" for p in addition_problems(paper, a)]
    return sorted(set(bad), key=bad.index)


@dataclass
class CheckResult:
    path: str
    doi: str = ""
    problems: List[str] = field(default_factory=list)
    warnings: List[str] = field(default_factory=list)
    rows: List[dict] = field(default_factory=list)              # one per condition
    new_ligands: List[Tuple[str, Optional[str]]] = field(default_factory=list)   # (name, class after)
    label_changes: List[Tuple[str, str, str]] = field(default_factory=list)      # (name, before, after)
    touched: Dict[str, dict] = field(default_factory=dict)      # name -> record after a dry-run apply

    @property
    def ok(self) -> bool:
        return not self.problems


def people(paper: dict) -> str:
    """'proposed by X; reviewed by Y' (what is known), for the record's curator and the condition's notes."""
    parts = [f"proposed by {paper['proposed_by']}" if paper.get("proposed_by") else "",
             f"reviewed by {paper['reviewed_by']}" if paper.get("reviewed_by") else ""]
    return "; ".join(x for x in parts if x)


def build_condition(paper: dict, a: dict) -> dict:
    """The `Lit_` condition one entry of a paper file describes (cells normalised, reference upgraded)."""
    locator = a["reference"]["Locator"]
    cells, _ = normalise_cells(a["cells"])
    c = {"condition_id": condition_id(a["ligand"], paper["DOI"], a["id"])}
    c.update(a.get("conditions") or {})
    notes = "; ".join(x for x in (paper.get("citation") or "", people(paper)) if x)
    c.update(distribution_ratios={}, reference={"DOI": paper["DOI"], "Source_Type": "Experimental",
                                                "Data_Used": locator, "Notes": notes})
    apply_resolution(c, {"cells": cells, "additives": a.get("additives"), "co_extractants": a.get("co_extractants"),
                         "reference": {"Locator": locator}, "note": a.get("note")})
    return c


def check_contribution(paper: dict, records: list[dict], today: str, path: str = "",
                       doi_exists: Optional[Callable[[str], Optional[bool]]] = None) -> CheckResult:
    """Everything that must hold before a file is applied, and what applying it would do.

    `doi_exists(doi)` -> True | False | None (could not tell); False is a problem, None a warning.
    The database and the file are not changed; `touched` holds the records a dry-run apply
    produced, which must pass the validators and the JSON Schema.
    """
    result = CheckResult(path=str(path), doi=str(paper.get("DOI") or "") if isinstance(paper, dict) else "")
    result.problems += contribution_problems(paper)
    if not isinstance(paper, dict):
        return result
    for key in sorted(set(paper) - KNOWN_TOP_KEYS):
        result.warnings.append(f"top-level key {key!r} is not part of the format (kept, ignored)")
    if paper.get("status") == STATUS_APPLIED:
        result.problems.append(f"the file is already applied ({paper.get('applied_on')}, "
                               f"{len(paper.get('condition_ids') or [])} conditions)")
    doi = paper.get("DOI")
    doi_ok = isinstance(doi, str) and normalize_doi(doi) == doi and doi not in NON_MEASUREMENT_DOIS
    if doi_exists is not None and doi_ok:
        known = doi_exists(doi)
        if known is False:
            result.problems.append(f"DOI {doi} does not resolve at doi.org")
        elif known is None:
            result.warnings.append(f"could not confirm DOI {doi} at doi.org (offline?)")

    names = {r["ligand"]["Name"]: r for r in records}
    stored_ids = {c["condition_id"] for r in records for c in r["experimental_data"]}
    aliases = {alt: r["ligand"]["Name"] for r in records for alt in (r["ligand"].get("Alternative_Names") or [])}
    by_key, by_conn = structure_index(records)
    declared = {}
    for lig in paper.get("new_ligands") or []:
        if not isinstance(lig, dict) or not lig.get("Name"):
            continue
        name = lig["Name"]
        if name in names:
            result.problems.append(f"new ligand {name}: a record of that name exists; add conditions to it instead")
        elif name in aliases:
            result.warnings.append(f"new ligand {name}: {aliases[name]} lists that name as a synonym")
        key = inchikey(lig.get("SMILES") or "")
        if key and key in by_key:
            result.problems.append(f"new ligand {name}: same structure as {by_key[key]} (InChIKey {key}); "
                                   f"add conditions to that record instead")
        elif key and by_conn.get(key.split("-")[0]):
            result.warnings.append(f"new ligand {name}: same connectivity as "
                                   f"{', '.join(by_conn[key.split('-')[0]])} (a stereoisomer? say so in reason)")
        declared[name] = lig

    signatures = {}
    for r in records:
        for c in r["experimental_data"]:
            signatures.setdefault(cell_signature(c), []).append((r["ligand"]["Name"], c["condition_id"]))
    seen_in_file = {}
    for a in paper.get("conditions") or []:
        if not isinstance(a, dict):
            continue
        aid, name = str(a.get("id")), a.get("ligand")
        row = {"id": aid, "ligand": name, "condition_id": "", "cells": 0, "bounds": 0, "measured_sigma": 0,
               "locator": (a.get("reference") or {}).get("Locator") or "", "outcome": ""}
        result.rows.append(row)
        if addition_problems(paper, a) or not doi_ok:
            row["outcome"] = "refused (see problems)"
            continue
        if name not in names and name not in declared:
            result.problems.append(f"{aid}: no record and no new_ligands entry named {name!r}")
            row["outcome"] = "refused: unknown ligand"
            continue
        cid = condition_id(name, doi, aid)
        row["condition_id"] = cid
        if cid in stored_ids:
            result.problems.append(f"{aid}: {cid} is already in the database (applied before, or the id reused)")
            row["outcome"] = "refused: already in the database"
            continue
        try:
            c = build_condition(paper, a)
        except (ValueError, KeyError, TypeError) as e:
            result.problems.append(f"{aid}: cannot be built as written: {e}")
            row["outcome"] = "refused: cannot be built"
            continue
        cells = c["distribution_ratios"]
        row["cells"] = len(cells)
        row["bounds"] = sum(1 for v in cells.values() if v.get("censor"))
        row["measured_sigma"] = sum(1 for v in cells.values() if v.get("sigma_source") == "measured")
        sig = cell_signature(c)
        same = signatures.get(sig, [])
        own = [other for owner, other in same if owner == name]
        others = [f"{owner}'s {other}" for owner, other in same if owner != name]
        if own:
            result.problems.append(f"{aid}: duplicate of {own[0]} (same elements and log D on {name})")
            row["outcome"] = f"refused: duplicate of {own[0]}"
            continue
        if others:
            result.warnings.append(f"{aid}: same cells as {others[0]} (another ligand; a transcription slip?)")
        if sig in seen_in_file:
            result.problems.append(f"{aid}: same cells as {seen_in_file[sig]} in this file")
            row["outcome"] = f"refused: repeats {seen_in_file[sig]}"
            continue
        seen_in_file[sig] = aid
        row["outcome"] = f"add to {name}" + (" (new record)" if name in declared else "")

    # GOVERNANCE.md §1: single-element series (Eu in an Am/Eu study) are in scope for ligands
    # already held; a new record needs at least one condition with two or more elements.
    for name in declared:
        sizes = [len(a.get("cells") or {}) for a in paper.get("conditions") or []
                 if isinstance(a, dict) and a.get("ligand") == name]
        if sizes and max(sizes) < 2:
            result.warnings.append(f"new ligand {name}: every condition holds a single element; GOVERNANCE.md "
                                   f"§1 admits such series only for ligands already in the database")

    if result.problems:
        return result
    applied = apply_contribution(paper, records, today)
    result.problems += applied.problems
    result.touched = applied.touched
    result.label_changes = applied.label_changes
    result.new_ligands = applied.new_ligands
    return result


@dataclass
class ApplyResult:
    records: List[dict]
    touched: Dict[str, dict]
    condition_ids: List[str]
    new_ligands: List[Tuple[str, Optional[str]]]
    label_changes: List[Tuple[str, str, str]]
    problems: List[str] = field(default_factory=list)


def apply_contribution(paper: dict, records: list[dict], today: str) -> ApplyResult:
    """The records with the paper's conditions added (untouched records are the same objects;
    touched ones are recomputed copies re-dated today; new records are created once a condition
    names them). Validates every touched record; `problems` lists what fails. The caller runs
    check_contribution first: this applies what it is given."""
    names = {r["ligand"]["Name"]: r for r in records}
    curator = people(paper)
    new = {lig["Name"]: new_ligand_record(lig, curator, today, LITERATURE)
           for lig in paper.get("new_ligands") or []}
    touched: Dict[str, dict] = {}
    cids = []
    for a in paper["conditions"]:
        name = a["ligand"]
        target = touched.get(name)
        if target is None:
            target = copy.deepcopy(names[name]) if name in names else new[name]
            touched[name] = target
        c = build_condition(paper, a)
        if target["metadata"].get("data_source") != LITERATURE:
            c["data_source"] = LITERATURE
        target["experimental_data"].append(c)
        cids.append(c["condition_id"])
    out, problems, changes, created = [], [], [], []
    for r in records:
        name = r["ligand"]["Name"]
        if name not in touched:
            out.append(r)
            continue
        fresh = recompute_record(touched[name])
        fresh["metadata"] = dict(fresh["metadata"], last_updated=today)
        touched[name] = fresh
        out.append(fresh)
        before, after = r["selectivity"].get("trend_classification"), fresh["selectivity"].get("trend_classification")
        changes.append((name, before, after))
    for name, r in new.items():
        if name in touched:
            fresh = recompute_record(touched[name])
            touched[name] = fresh
            out.append(fresh)
            created.append((name, fresh["selectivity"].get("trend_classification")))
    if PROFILE == "public":   # the public schema allows only PUBLIC_FIELDS (no source labels)
        touched = {name: project_public(r) for name, r in touched.items()}
        out = [touched.get(r["ligand"]["Name"], r) for r in out]
    for name, r in touched.items():
        ok, errors = validate_ligand_record(r)
        errors += validate_against_json_schema(r)
        problems += [f"{name}: {e}" for e in errors[:5]]
    return ApplyResult(out, touched, cids, created, changes, problems)


def report_lines(result: CheckResult, paper: dict) -> list[str]:
    """The check report (Markdown): verdict, problems, warnings, what applying would do."""
    stem = Path(result.path).stem if result.path else "contribution"
    lines = [f"# Contribution check: {stem}", "",
             f"- DOI: {result.doi or '(none)'}",
             f"- Citation: {paper.get('citation') if isinstance(paper, dict) else ''}",
             f"- Proposed by: {paper.get('proposed_by') if isinstance(paper, dict) else ''}",
             f"- Reviewed by: {(paper.get('reviewed_by') if isinstance(paper, dict) else None) or '(not yet)'}",
             f"- Verdict: **{'No problems' if result.ok else f'{len(result.problems)} problem(s)'}**"
             + (f", {len(result.warnings)} warning(s)" if result.warnings else ""), ""]
    if result.problems:
        lines += ["## Problems", ""] + [f"- {p}" for p in result.problems] + [""]
    if result.warnings:
        lines += ["## Warnings", ""] + [f"- {w}" for w in result.warnings] + [""]
    if result.rows:
        lines += ["## Conditions", "", "| id | ligand | condition_id | locator | cells | bounds | measured σ | outcome |",
                  "|---|---|---|---|---:|---:|---:|---|"]
        lines += [f"| {r['id']} | {r['ligand']} | {r['condition_id'] or '—'} | {r['locator']} | {r['cells']} | "
                  f"{r['bounds']} | {r['measured_sigma']} | {r['outcome']} |" for r in result.rows]
        lines.append("")
    if result.ok:
        lines += ["## What applying would do", ""]
        for name, cls in result.new_ligands:
            lines.append(f"- {name}: new record ({cls})")
        for name, before, after in result.label_changes:
            flag = "" if before == after else " **(label changes)**"
            lines.append(f"- {name}: {before} -> {after}{flag}")
        lines.append("")
    return lines
