# scripts/

| Script | What it does |
|---|---|
| `query_ligands.py` | Query the database from the command line (`--stats`, `--name`, `--trend`, `--mechanism`, ...); see `docs/CLI.md` |
| `validate_database.py` | Every record against the validator and the JSON Schema (CI) |
| `recompute_derived_fields.py` | `--check`: every derived field reproduces from the stored measurements (CI) |
| `generate_json_schema.py` | Writes `schema/` (JSON Schema, data dictionary, contribution schema); `--check` in CI |
| `make_benchmark_splits.py` | Writes `data/splits/benchmark_splits.json` |
| `export_release_bundle.py` | The data release bundle (JSONL, CSVs, schema, vocabularies, licence, checksums) → `dist/` |
| `check_contribution.py` | Checks `contributions/<paper>.json` against the database (CI on data pull requests) |
| `build_skill.py` | Assembles the Claude Code skill (`skills/reel-quest/`) |
| `build_space.py` | Assembles the Hugging Face Space bundle (`space/`) |

The manuscript's numbers and figures: `paper/stats.py`, `paper/figures.py`
(`paper/README.md`).
