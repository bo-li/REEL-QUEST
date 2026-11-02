# Governance

How REEL-QUEST is maintained: who decides what enters the database, how an
entry is checked, how releases are versioned and cited, what happens to data
that turns out to be wrong, and on what terms the data are redistributed.
The contribution mechanics are in `contributions/README.md` and
`CONTRIBUTING.md`; the data terms in `docs/DATA_PROVENANCE_AND_TERMS.md`.

## 1. Scope

The database holds liquid–liquid extraction distribution ratios of the
lanthanides (La–Lu) plus Y and Sc, one record per extractant molecule, one
condition per experiment, every value traced to the paper that measured it and
the table or figure it was read from. Out of scope: actinide-only data,
process flowsheets, simulated or predicted values, and values whose source
cannot be named. The 21-subtype taxonomy (`vocab/chemical_families.json`) is
the chemistry vocabulary; a molecule outside it (an S-donor, say) needs a
taxonomy decision before a record.

A series that reports one element only (Eu in an Am/Eu study, say) is in
scope for a ligand the database already holds: it adds conditions, often the
acid dependence, to a record characterised elsewhere. A new ligand needs at
least one condition with two or more elements; one known only from
single-element series would be a record no trend or separation factor can
use, the reason the sparse Am/Eu structures of SAFE's 2026-02-10 snapshot
were not ingested. `scripts/check_contribution.py` warns on such a new
ligand; the maintainer decides.

## 2. Roles

| Role | Who | Does |
|---|---|---|
| Maintainer | Bo Li (repository owner); further maintainers are named in `CITATION.cff` | Reviews and applies contributions, decides taxonomy and schema changes, cuts releases, answers issues |
| Reviewer | A maintainer or a named chemist | Reads the paper against a contribution file (§4) |
| Contributor | Anyone, human or AI agent | Proposes a paper's data (§3), reports a wrong value |

Decisions that change the shape of the data (a schema field, a subtype, a
mechanism label, a source policy) are made by a maintainer and written into
`CHANGELOG.md`.

## 3. How data enter

Two intakes, one format (`schema/contribution.schema.json`), one apply path:

1. **A contribution file**: one paper, one file under `contributions/`, by a
   pull request. The CI job (`.github/workflows/contributions.yml`) runs
   `scripts/check_contribution.py` on it: schema, DOI, names and structures
   against the database, the cell grammar, duplicates, and a dry-run apply
   whose records must validate. For developers and AI agents.
2. **A paper, or a filled sheet**: an issue ("Add a paper's data") with the
   DOI and the tables or figures. A maintainer, or an AI agent, drafts the
   file; the proposer is named in `proposed_by`. For chemists.

Either way the file is reviewed (§4) and applied by a maintainer in the
maintainers' repository, which recomputes every derived field of the records
it touches, re-dates them, marks the file applied with the condition ids it
created and writes an audit; the next release publishes the result here. The
file stays in the repository as the provenance of those conditions.

What cannot be contributed: values without a source, values copied from a
database with terms that forbid redistribution, pooled or averaged conditions
presented as one experiment, and edits to derived fields (they are recomputed).

## 4. The review rule

Every condition is read against its paper by someone other than the person who
transcribed it before it is applied: a second reader with the PDF, or an AI
check with page-level evidence followed by a human spot-check, as every v14
condition had. The reviewer confirms
the molecule, the acid and its concentration, the extractant concentration,
the diluent, the locator, and the values; a figure-read value is checked on
the rendered figure. Disagreements go in the pull request; what the reviewer
decides, and why, goes in the file's `why` and `note` fields.

A contribution from an author of the paper still gets the second read.

The reviewer is named in the file (`reviewed_by`), and the apply step refuses
a file that names none unless run with `--unreviewed`, which its audit
records.

## 5. Integrity

The checks that run on every push and pull request (`.github/workflows/ci.yml`):

- every record passes the validator and the JSON Schema
  (`scripts/validate_database.py`);
- every derived field reproduces from the stored measurements
  (`scripts/recompute_derived_fields.py --check`): trends, bootstrap intervals,
  statistics, per-cell σ, outlier flags, quality grades, element lists;
- a stored D agrees with its own log D; `element_coverage` and
  `elements_measured` mean what `docs/DATABASE_STRUCTURE.md` says they mean;
- every reference carries a bare DOI of the measuring paper and a locator;
- the numbers in the README and the manuscript come from `paper/stats.py`.

Nothing edits the database by hand. Every change is made by a script in the
maintainers' repository that writes an audit, and each release here is
exported from it (`RELEASE_PROVENANCE.json` names the build). Access to write
is by pull request; the maintainer merges.

## 6. Versions and releases

- The data file is versioned by build (`v14` is the verified build of
  2026-10-06). Contributions are applied to it in
  place and listed in `CHANGELOG.md`; a release tags the state.
- Releases follow semantic versioning on the data: a patch release adds or
  corrects conditions, a minor release adds ligands or a field, a major release
  changes the meaning of a field or the taxonomy. Each release is a GitHub
  Release and a Zenodo deposit with its own DOI; `CITATION.cff` carries the
  concept DOI.
- The schema has its own version (`schema_version`, `rel_quest/schema.py`);
  a schema change is additive within a major version and documented in
  `docs/DATABASE_STRUCTURE.md`.

## 7. Retirement

A value, condition or record is never deleted silently.

- A wrong value is corrected through a contribution or a repair script with an
  audit; the previous value stays in the previous release and the audit.
- A condition that turns out not to be in its paper, or to duplicate another,
  is dropped by a script that lists it, and the release notes say so.
- A record that loses every condition is retired to a listing with its
  structure and the reason.
- A source whose terms turn out to forbid redistribution is removed, every
  derived field is recomputed from what remains, and the release notes say so.

## 8. Access and terms

The code is MIT (`LICENSE`); the data, vocabularies and schema are CC BY 4.0
(`LICENSE-DATA`). Everything is public: the repository, the Zenodo bundle, the
dashboard, the REST API and the MCP server. Citing: `CITATION.cff`, plus the
papers that reported the values used (`reference.DOI`) and the SAFE paper,
from which part of the data were compiled. The provenance of every source and the basis on which
its values are redistributed is in `docs/DATA_PROVENANCE_AND_TERMS.md`; a
contributor may add only values they may redistribute.

## 9. Sustainment

The database is maintained by the group that built it, as part of its
research on rare-earth separations; the cost of keeping it is the CI and the
review time of §4. The code depends on pinned packages
(`requirements.lock`) and on nothing hosted but GitHub, Zenodo and the
dashboard host. If maintenance stops, the last release remains citable on
Zenodo and the repository can be forked under its licences.

## 10. Contact

Open an issue (templates: a wrong value; a paper to add), or write to the
maintainer named in `CITATION.cff`.
