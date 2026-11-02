# Contributing Guide  
## Rare Earth Solvent Extraction Database (SX)

Thank you for contributing to the **Rare Earth Separation Database**.  
This document defines **how data should be added, validated, and reviewed**.
Who decides, the review rule, versioning and retirement are in `GOVERNANCE.md`.

This database is designed to be:
- machine-readable,
- chemistry-faithful,
- citation-driven,
- and extensible without breaking backward compatibility.

Both **human contributors and AI agents** are expected to follow these rules.
The short version: one paper, one file under `contributions/`, checked by
`scripts/check_contribution.py`, reviewed against the paper, applied by a
maintainer (`contributions/README.md`).

---

## 1. Data Sources

The database was compiled from the Liu et al. 2022 JACS Au Supporting
Information (`10.1021/jacsau.2c00122`), the SAFE database
(https://safe.lanl.gov/; cite `10.1080/07366299.2025.2564381`), the OECD-NEA
IDEaL database (https://www.oecd-nea.org/ideal/) and papers read directly.
Every condition cites the paper that measured the values: a bare DOI in
`reference.DOI` (never a compilation's DOI) and the figure or table in
`reference.Locator`. New data enter from the measuring paper itself, one file
per paper (§4); `docs/DATA_PROVENANCE_AND_TERMS.md` has the sources' terms.

---

## 2. What Can Be Contributed

Valid contributions include:

- New **ligand (or ligand system)** entries
- Additional **performance data** for existing ligands
- Improved **experimental condition metadata**
- Additional **references** clarifying data provenance
- Minor corrections (typos, units, formatting)

Invalid contributions include:
- Uncited inferred data
- Hypothetical or simulated results without labeling
- Aggregated trends without source attribution
- Process flowsheets (not supported)

---

## 2. Fundamental Rule

> **One database record = one ligand (or ligand system).**

Do NOT:
- merge multiple ligands into one entry
- split one ligand into multiple entries unless chemically distinct

---

## 3. Required Files to Read Before Contributing

All contributors (human or AI) MUST read:

1. `DATABASE_STRUCTURE.md`  
2. Controlled vocabulary files in `/vocab/`

Contributions that violate the schema will be rejected.

---

## 4. How to Add a Paper's Data

One paper, one file: `contributions/<doi-slug>.json` (format and a worked
example in `contributions/README.md`; schema in
`schema/contribution.schema.json`). The file names the paper (its DOI, never a
compilation's), who read it, the ligands it adds (only those no record holds),
and one condition per experiment with its locator (table, figure, text page)
and its cells.

### Step 1 — Find the record

```bash
python scripts/query_ligands.py --name TODGA
```

Search by name and by structure (the check refuses a SMILES whose InChIKey is
already in the database). A ligand with a record takes its conditions under
`"ligand": "<Name>"`; one without is declared once in `new_ligands` with
`Name`, `SMILES`, `Mechanism`, `Extractant_Type`, `Subtype`
(`vocab/chemical_families.json`) and a `reason` saying where the paper
identifies the molecule.

### Step 2 — One condition per experiment

Each condition carries the schema fields the paper states
(`Acid_Name`, `Acid_Concentration_M`, `Extractant_Concentration_M`,
`Solvent_Name`, `Temperature_C`, `Contact_Time`, …), the aqueous additives as
text, `reference.Locator`, and its cells: `log_D` or `D` with `value_origin`
(`reported` from a table or the text, `digitised` from a figure),
`sigma_log_D` only when the source states an uncertainty, a bound as
`censor` + `log_bound`. Do not pool experiments into one condition, do not
average, do not fill a field the paper leaves blank.

### Step 3 — Say why

`why` states what in the paper supports the condition and the values (the
row, the series, the caption). The reviewer reads it first.

### Step 4 — Check

```bash
python scripts/check_contribution.py contributions/<doi-slug>.json
```

The report (`outputs/contributions/<doi-slug>_check.md`, created on first run) lists problems (the
file is refused), warnings, and what applying it would do, including any
trend label that changes.

---

## 5. Partial Data Is Allowed (If Explicit)

The database allows incomplete records IF:

- missing data is clearly stated
- no values are fabricated
- references specify what was and was not extracted

Example:
- Allowed: selectivity trend without D values
- Not allowed: estimated D values without citation

---

## 6. Units and Conventions

- Distribution ratios (D): dimensionless
- Separation factors (SF): dimensionless
- Acid concentration: mol/L (`Acid_Concentration_M`)
- Temperature: Celsius (`Temperature_C`)
- REEs: use chemical symbols (`La`, `Nd`, `Dy`, etc.)

Do NOT mix units.

---

## 7. Validation Before Submission

`scripts/check_contribution.py` runs the same checks the CI job runs on a
pull request: the contribution schema; a bare DOI that resolves and is not a
compilation's; names and structures against the database; the cell grammar;
duplicates of stored conditions; and a dry-run apply whose records must pass
the validator and the JSON Schema. A file with a problem is refused; fix it
and run the check again.

For code changes, the usual gates apply: `python scripts/validate_database.py`,
`python scripts/recompute_derived_fields.py --check` and `pytest`
(`.github/workflows/ci.yml`).

---

## 8. Contribution Workflow

1. Fork or branch; add `contributions/<doi-slug>.json`.
2. Run `python scripts/check_contribution.py contributions/<doi-slug>.json`
   until it reports no problems.
3. Open a pull request (the template asks for the DOI, the locators and the
   licence line). The `Contributions` workflow posts the check report as the
   job summary.
4. A reviewer reads the paper against the file (`GOVERNANCE.md` §4), asks
   for changes in the pull request, and adds their name as `reviewed_by`.
5. A maintainer applies it in the maintainers' repository, which recomputes
   every derived field of the records it touches and records what was done,
   and the next release publishes it here: the file is marked
   `"status": "applied"` with the condition ids it created, and the change is
   listed in `CHANGELOG.md`.

No JSON? Open an issue with the DOI and the table or figure ("Add a paper's
data"), or send a filled verification sheet; a maintainer drafts the file and
names you in `proposed_by`. A wrong stored value: the "A stored value is
wrong" issue template.


---

## 9. Rules for AI Agents

AI agents contributing to this database MUST:

- Never hallucinate data
- Never infer missing values
- Always include explicit references
- Prefer descriptive completeness over numerical completeness
- Respect controlled vocabularies exactly

If uncertain, leave fields empty and document uncertainty.

---

## 10. Review Criteria

Contributions are evaluated based on:

- Schema compliance
- Data provenance clarity
- Chemical correctness
- Consistency with existing entries
- Minimal redundancy

---

## 11. Versioning Policy

- Current schema version: **v0.3** (`rel_quest/schema.py`; history in its
  docstring and in `docs/DATABASE_STRUCTURE.md`)
- Contributions must not break backward compatibility
- Structural changes require maintainer approval

---

## 12. Licence of contributions

Contributed code is released under MIT (`LICENSE`); contributed data under
CC BY 4.0 (`LICENSE-DATA`). Contribute only values you may redistribute, and
cite the paper that reported them.

## 13. Questions and Extensions

If you are unsure:
- open an issue
- leave a comment in `Notes`
- flag the entry as incomplete

Future schema extensions (v1+) will be documented separately.

---

Thank you for contributing to a **clean, reliable, and AI-ready** rare earth separation database.
