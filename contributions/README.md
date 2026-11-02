# contributions/ — one paper, one file

This directory is how measurements enter REEL-QUEST after v14. Each file holds
the distribution ratios one paper reports, read from its tables or figures, and
is applied to the released database by a maintainer once it has been checked
and reviewed. The format is `schema/contribution.schema.json`; the rules are in
`rel_quest/contributions.py` and the process in `GOVERNANCE.md`.

## The file

`contributions/<doi-slug>.json`, where the slug is the DOI with every character
other than letters, digits and dots replaced by `_` (`10.1021/jacs.4c07332` →
`10.1021_jacs.4c07332.json`).

```json
{
  "DOI": "10.1021/jacs.4c07332",
  "citation": "Pramanik et al., J. Am. Chem. Soc. 2024, 146, 25669",
  "proposed_by": "A. Chemist, 2026-11-02, SI Tables S1 and S2",
  "reviewed_by": "B. Reviewer, 2026-11-09, against the SI",
  "new_ligands": [
    {
      "Name": "BLBPy",
      "SMILES": "...",
      "Mechanism": "Solvation",
      "Extractant_Type": "Bipyridine_Amide",
      "Subtype": "Bipydiamide",
      "reason": "SI p. S4: formula C.., name ..; Table S1 row 2"
    }
  ],
  "conditions": [
    {
      "id": "A0001",
      "ligand": "BLBPy",
      "conditions": {
        "Acid_Name": "HNO3",
        "Acid_Concentration_M": 1.0,
        "Extractant_Concentration_M": 0.01,
        "Solvent_Name": "1,2-dichloroethane",
        "Temperature_C": 25.0,
        "Contact_Time": "60 min"
      },
      "additives": "1 M NaNO3",
      "reference": {"Locator": "Table S1"},
      "cells": {
        "La": {"D": 0.12, "value_origin": "reported"},
        "Eu": {"log_D": 0.45, "value_origin": "reported"},
        "Nd": {"log_D": 0.10, "sigma_log_D": 0.03, "value_origin": "digitised"},
        "Lu": {"censor": ">", "log_bound": 2.0}
      },
      "why": "Table S1, the 1 M HNO3 row; the SI gives D directly",
      "note": null
    }
  ]
}
```

| Key | Meaning |
|---|---|
| `DOI` | The paper that **measured** the values, bare and lower-case. Never a compilation's DOI (SAFE, Liu 2022, IDEaL): those are where a value was found, not where it was measured. |
| `proposed_by` | Who read the paper, when, from what. It is kept with the curation record of the conditions it adds. |
| `reviewed_by` | Who read the paper against the file before it was applied (the review rule, `GOVERNANCE.md` §4). The reviewer adds it in the pull request; the apply step refuses a file without it unless run with `--unreviewed`, which its audit records. |
| `new_ligands` | Only for molecules no record holds, and only with at least one condition of two or more elements (a single-element series, such as Eu alone, can only add to a stored ligand; `GOVERNANCE.md` §1). Check first with `python scripts/query_ligands.py --name <name>` and by structure: the check refuses a SMILES whose InChIKey is already in the database and warns on a stereoisomer. `Subtype` comes from `vocab/chemical_families.json`; `Mechanism` from `vocab/mechanisms.json`. Leave `Chemical_Family` out: it is set to `Subtype`. |
| `conditions[].id` | `A0001`, `A0002`, … unique within the file; the applied file lists the condition id each one became. |
| `conditions[].conditions` | Schema condition fields with units in the names (`_M` molar, `Temperature_C`). Leave out what the paper does not state; do not guess. |
| `reference.Locator` | `Table 2`, `Figure 3`, `Text p. 5`: where the values are. A table gives `value_origin` `reported`; a figure `digitised`. |
| `cells` | One entry per element: `log_D` or `D` with `value_origin`; `sigma_log_D` only when the source states an uncertainty (a printed error, a calibrated figure); a bound as `censor` (`>` or `<`) with `log_bound`. Bounds stay out of every fit, statistic and separation factor. |
| `why` | What in the paper supports these conditions and values. The reviewer reads this first. |

A mixed organic phase (two extractants) goes in `co_extractants`
(`[{"Name": ..., "SMILES": ..., "Concentration_M": ...}]`); the condition then
stays on the main ligand's record but out of pooled statistics.

## Values read off a figure

Keep the calibration you read them with:
`contributions/figures/<doi-slug>/<locator>.json` (format in
`.agents/skills/reel-quest-contribute/references/figures.md`). It records
the page, the axis ticks and each marker in pixels, so a reviewer with the
paper can redraw your readings on the figure. The paper's PDF is not
committed; the calibration names it as `papers/<doi-slug>.pdf`.

## With an AI agent

The repository ships a skill, `.agents/skills/reel-quest-contribute/`, that
drafts this file from a paper's PDF or does the second read of one against
its paper. Codex and Gemini CLI find it there; Claude Code finds it through
`.claude/skills/reel-quest-contribute/`. Open the agent in your clone and ask,
for example, "draft a contribution from papers/<file>.pdf". The file it
writes goes through the same check and review as any other.

## Before opening a pull request

```bash
python scripts/check_contribution.py contributions/<doi-slug>.json
```

The report lands in `outputs/contributions/<doi-slug>_check.md`. It refuses a
file for a DOI that is not a DOI or does not resolve, a structure or name the
database already has, a cell without `value_origin`, a condition that repeats
one already stored, a field the schema has not got, or applied records that do
not validate; it warns on a stereoisomer of a stored record, cells identical to
another ligand's, and keys it does not know. The same check runs in CI on every
pull request that touches this directory (`.github/workflows/contributions.yml`).

## What happens next

A maintainer reads the paper against the file (the review rule in
`GOVERNANCE.md`: a second read, as every v14 condition had) and asks for
changes in the pull request. Once it is reviewed, the maintainer applies it in
the maintainers' repository, which recomputes every derived field of the
records it touches and writes an audit, and the next release publishes it
here, with the file marked `"status": "applied"` and the condition ids it
created. The file stays here as the record of where those conditions came from.

## No JSON? Send the paper

If writing the file is the obstacle, open an issue with the DOI and the table
or figure the values are in. A maintainer, or an AI agent, drafts the file;
you are named in `proposed_by` and the review is the same.
