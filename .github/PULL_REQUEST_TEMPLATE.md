<!-- For a data contribution (a file under contributions/), keep the checklist; for code, delete it. -->

## What this adds

- Paper (DOI):
- Where the values come from (table, figure, text page):
- Ligands: existing records / new (declared in `new_ligands`):

## Checklist for a data contribution

- [ ] One file per paper under `contributions/`, named by the DOI slug (`contributions/README.md`)
- [ ] `DOI` is the paper that measured the values, not a compilation
- [ ] Every cell says whether it was `reported` or `digitised`; `sigma_log_D` only where the source states one; bounds as `censor`
- [ ] `why` says what in the paper supports each condition
- [ ] `python scripts/check_contribution.py contributions/<file>.json` reports no problems (the CI job repeats it)
- [ ] I may redistribute these values under CC BY 4.0 (`LICENSE-DATA`): they are the paper's reported numbers with attribution, nothing copied from a database with its own terms
