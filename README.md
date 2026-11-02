# REEL-QUEST — Rare-Earth Extraction Ligand Database

[![CI](https://github.com/bo-li/REEL-QUEST/actions/workflows/ci.yml/badge.svg)](https://github.com/bo-li/REEL-QUEST/actions/workflows/ci.yml)

**REEL-QUEST** (Rare-Earth Extraction Ligands — Quantitative Understanding of Elemental Selectivity Trends) is a curated, queryable database of liquid–liquid extraction data for rare-earth-element (REE) separation, with mechanism classification, selectivity-trend analysis with uncertainty, and access for people and AI agents: an interactive dashboard with an AI assistant, an MCP server, an agent skill, a REST API, a Python package and a CLI.

## Quick start

**Interactive dashboard:**

```bash
pip install -e .            # installs deps and puts rel_quest on the path
streamlit run app.py        # then open http://localhost:8501
```

Install **editable** (`-e`): the database and vocabularies live at the repository
root and are resolved relative to the package, so a non-editable install imports
but cannot find `data/`. `pip install -r requirements.txt` is equivalent for the
dependencies alone (it is what the Dockerfile uses).

**CLI queries:**

```bash
python scripts/query_ligands.py --stats
python scripts/query_ligands.py --name TODGA --format full
python scripts/query_ligands.py --trend Ascending --output csv
```

**Python API:**

```python
from rel_quest import LigandDatabase

db = LigandDatabase()                 # finds data/sx_ligands_merged_v14_flat.jsonl
ligand = db.query_by_name('TODGA')[0]
print(ligand.get_d_values())          # {'La': 0.5, ..., 'Lu': 15.7}
print(ligand.get_trend_summary()['classification'])
```

**REST API:**

```bash
uvicorn rel_quest.api:app             # docs at http://127.0.0.1:8000/api/docs
curl 'localhost:8000/api/v1/ligands?subtype=DGA&elements=Eu,La'
```

See `docs/REST_API.md` for endpoints and guarantees. The Hugging Face Space
serves the same API at `/api` (once the Space is public).

**MCP server** (query the database from Claude Code, Claude Desktop or any MCP client):

```bash
pip install -e ".[mcp]"
python -m rel_quest.mcp_server        # stdio
```

See `docs/MCP_SERVER.md` for the tools and client configuration, and for the
Hugging Face Space (`space/`): a Gradio demo that also serves the MCP tools at `/mcp`.
For agents that run code there is also a Claude Code skill (`skills/reel-quest/`,
`docs/CLAUDE_CODE_SKILL.md`): `python scripts/build_skill.py --out ~/.claude/skills/reel-quest`.

**AI assistant:** the dashboard's assistant runs on OpenRouter (or any
OpenAI-compatible provider). Put a key in `.streamlit/secrets.toml` — copy
`.streamlit/secrets.toml.example` — and visitors can chat on free models
without supplying a key of their own; they can also enter their own key in
⚙️ Settings. A local Ollama server, when running, is offered as a third option.

## Deploying

`docs/DEPLOYMENT_GUIDE.md` covers Streamlit Community Cloud (main file `app.py`,
Python 3.11, secrets pasted into the app's settings), Docker, EC2 and Cloud Run.

## Contributing data

One paper, one file: `contributions/<doi-slug>.json` with the paper's DOI, the
ligands it adds and one condition per experiment, each cell marked reported or
digitised, bounds kept as bounds. `python scripts/check_contribution.py` checks
it against the database, and the CI job repeats the check on every pull
request. A reviewer reads the paper against it, and a maintainer applies it to
the next release. Format: `contributions/README.md`; rules and roles:
`CONTRIBUTING.md`, `GOVERNANCE.md`. No JSON? Open an issue with the DOI and the
table or figure.

## Tests and reproducibility

```bash
pip install -e ".[dev]"                          # adds pytest
pytest                                           # tests, incl. a dashboard smoke test
python scripts/validate_database.py              # validator + JSON Schema, every record
python scripts/recompute_derived_fields.py --check   # every derived field reproduces
python paper/stats.py --check                    # the manuscript's numbers
python paper/figures.py                          # the manuscript's figures (pip install -e ".[paper]")
```

CI runs the suite on Python 3.10 and 3.11, and again against the fully pinned
`requirements.lock` (Python 3.11), which is also what the Dockerfile installs.

What reproduces from this repository, and what does not: every derived field of
every record, every figure, and every number in `paper/stats.json` except those
that describe how the database was assembled (counts per source compilation,
the cross-source comparison of Fig. S2, and the verification tallies). Those
need the curation records, which stay with the maintainers; `paper/README.md`
lists them, and this repository carries their values as computed there.
`RELEASE_PROVENANCE.json` names the build this release was exported from.

## Data availability

- **Database:** `data/sx_ligands_merged_v14_flat.jsonl` (JSON Lines, one ligand per line), with flat CSVs of all D-values and ligands in the release bundle (`python scripts/export_release_bundle.py` → `dist/`).
- **Schema and data dictionary:** `schema/`.
- **Benchmark splits:** `data/splits/benchmark_splits.json` — leakage-controlled 5-fold and leave-subtype-out splits for machine-learning use.
- **Archive:** Zenodo, https://doi.org/10.5281/zenodo.23819219 (v1.0.0).

## Licence

Code (`rel_quest/`, `app.py`, `scripts/`, `tests/`, `paper/`) is MIT (`LICENSE`).
The database, vocabularies and schema are CC BY 4.0 (`LICENSE-DATA`).

## How to cite

Use GitHub's "Cite this repository" (from `CITATION.cff`). When you rely on
specific measurements, also cite the papers in their `reference.DOI`;
`paper/references.bib` has every one of them in BibTeX.

## Repository layout

```
REEL-QUEST/
├── app.py                  # Streamlit dashboard
├── rel_quest/              # Python package: queries, trends, uncertainty, validation, APIs, agents
├── scripts/                # CLI tools, checks, release bundle, skill and Space builders
├── data/                   # The database and benchmark splits
├── schema/                 # JSON Schema + data dictionary (generated)
├── vocab/                  # Controlled vocabularies (taxonomy, mechanisms, solvents…)
├── paper/                  # The manuscript's statistics and figures (generated)
├── contributions/          # One file per contributed paper
├── skills/, space/         # Claude Code skill; Hugging Face Space
├── docs/                   # Documentation (schema, API, CLI, MCP, deployment)
└── tests/                  # pytest suite, run in CI
```

## Documentation

- `docs/DATABASE_STRUCTURE.md` — record schema (prose); `schema/` — machine-readable
- `docs/API.md`, `docs/CLI.md`, `docs/QUICKSTART.md`, `docs/USAGE_EXAMPLES.md` — usage
- `docs/REST_API.md`, `docs/MCP_SERVER.md`, `docs/CLAUDE_CODE_SKILL.md` — programmatic and agent access
- `docs/LIGAND_FAMILY_CLASSIFICATION_MEMO.md` — taxonomy rationale (21 subtypes)
- `docs/DATA_PROVENANCE_AND_TERMS.md` — sources and their terms
- `CONTRIBUTING.md`, `GOVERNANCE.md`, `CHANGELOG.md`
