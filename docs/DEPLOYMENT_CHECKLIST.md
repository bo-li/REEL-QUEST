# Quick Deployment Checklist

Getting the REEL-QUEST dashboard (and optionally the REST API) running on a
new host. For production hardening — systemd, nginx, HTTPS, cloud — see
`DEPLOYMENT_GUIDE.md`.

## Option A: Docker (fewest moving parts)

```bash
docker compose up -d          # http://localhost:8501
```

`docker-compose.yml` also wires the dashboard to an Ollama container on the
`llm-server_default` network for the AI assistant. Drop the `networks:`
block if you are not running one — the rest of the dashboard works without.

## Option B: The helper script

```bash
git clone <repo> && cd Database
./deploy.sh                              # http://localhost:8501
./deploy.sh --port 8080 --host 0.0.0.0   # external access
./deploy.sh --skip-checks                # skip validation + tests
```

It checks Python (3.10+ required), verifies the required files are present,
creates and activates `venv/`, installs `requirements.txt` plus the package
itself editable, confirms the database loads and validates, runs the test
suite, then launches Streamlit.

## Option C: By hand

```bash
python3 -m venv venv
source venv/bin/activate                 # Windows: venv\Scripts\activate
pip install -e .                         # editable — see below
streamlit run app.py
```

**Install editable (`-e`).** The database and vocabularies live at the repo
root and are resolved relative to the package, so a non-editable install
imports fine but raises `FileNotFoundError` from `LigandDatabase()`.

## What has to be present

The dashboard needs the whole checkout, not a file subset — `app.py` imports
shared constants and plot helpers from `rel_quest/`. If you are copying by
hand rather than cloning, these are what `deploy.sh` refuses to start
without:

```
app.py
requirements.txt
pyproject.toml
data/sx_ligands_merged_v14_flat.jsonl     (~3 MB, 169 ligands)
rel_quest/                                (the package)
vocab/                                    (elements, mechanisms, families, acids, solvents)
```


## Verification

After launch, open the dashboard and confirm:

- [ ] Page loads without a traceback
- [ ] Sidebar filters render
- [ ] Sidebar reports **170 ligands**
- [ ] Database Overview plots draw
- [ ] Ligand Browser table loads and a row opens the detail view

From the command line:

```bash
python scripts/validate_database.py     # expect 170 valid, 0 invalid
pytest -q                               # full suite
python scripts/query_ligands.py --stats
```

## REST API (optional)

```bash
uvicorn rel_quest.api:app --host 0.0.0.0 --port 8000
```

Docs at `/api/docs`; endpoints and guarantees in `docs/REST_API.md`. It is
read-only and unauthenticated — put it behind a reverse proxy before
exposing it publicly.

## Troubleshooting

**`ModuleNotFoundError: rel_quest`** — run `pip install -e .`, or launch from
the repo root so the CWD is on `sys.path`.

**`FileNotFoundError` from `LigandDatabase()`** — a non-editable install. See
the note above; reinstall with `-e`, or pass an explicit `db_path`.

**`./deploy.sh: /bin/bash^M: bad interpreter`** — the checkout rewrote the
script to CRLF. `.gitattributes` pins `*.sh` to LF; verify with
`file deploy.sh`, and re-clone with `core.autocrlf=input` if it persists.

**`Permission denied` running `./deploy.sh`** — `chmod +x deploy.sh`, or run
`bash deploy.sh`.

**Port already in use** — `./deploy.sh --port 8080`.

**RDKit missing (non-critical)** — molecular structure images disappear;
everything else works. `conda install -c conda-forge rdkit` if pip fails.

## System requirements

|  | Minimum | Recommended |
|---|---|---|
| Python | 3.10 | 3.11 |
| RAM | 2 GB | 4 GB |
| Disk | 1 GB | 2 GB |

Dependencies total roughly 500 MB installed, dominated by RDKit and
Streamlit's transitive tree.
