# REST API

Read-only HTTP access to the REEL-QUEST database (`rel_quest/api.py`), built
with FastAPI on top of the same `LigandDatabase` the dashboard and CLI use.
It exists so the database can be queried without cloning the repo — the first
step toward a public resource in the style of unicornedb.org.

## Hosted

The Hugging Face Space (`space/app.py`, `docs/MCP_SERVER.md`) mounts this app
at `/api` of its Gradio server, on the database it already loads: the same
paths, docs and CORS policy as below, at `https://<the Space's host>/api/v1/`
and `https://<the Space's host>/api/docs`. Not deployed yet; public once the
manuscript is published. A free Space sleeps when idle, so the first request
after a pause waits for it to start. `tests/test_mcp_server.py`
(`TestRESTAPIInSpaceApp`) checks the mounted API against this one.

`create_app(db=..., mount="")` is how another app mounts it: `mount` is the
path the API lives under in the app that serves it (`/api` standalone, empty
when mounted at `/api`).

## Run it

```bash
pip install -e .                       # fastapi + uvicorn are core deps
uvicorn rel_quest.api:app              # http://127.0.0.1:8000
```

Interactive documentation is generated from the code:

| URL | What |
|---|---|
| `/api/docs` | Swagger UI (try requests in the browser) |
| `/api/redoc` | ReDoc reference |
| `/api/v1/openapi.json` | OpenAPI 3 spec |

`/` redirects to `/api/docs`.

## Endpoints

All endpoints are GET and return JSON. There is no authentication and no
write path — the API serves the same canonical
`data/sx_ligands_merged_v14_flat.jsonl` that ships with the checkout.

### `GET /api/v1/health`

Liveness plus live record counts:

```json
{"status": "ok", "version": "0.1.0", "ligands": 170, "measurements": 665}
```

### `GET /api/v1/stats`

`LigandDatabase.get_statistics()` verbatim — totals and distributions by
mechanism, trend, and structural class.

### `GET /api/v1/vocab`

Every value accepted by the list filters, with ligand counts. Use this for
discoverability instead of hardcoding taxonomy labels:

```json
{
  "subtypes":   {"DGA": 75, "Phosphine_oxide": 19, ...},
  "mechanisms": {"Solvation": 120, "Cation_Exchange": 50},
  "trends":     {"Insufficient_Data": 54, ..., "Peaked": 10, "U_Shaped": 8},
  "sources":    {...},
  "elements":   ["Sc", "Y", "La", "Ce", ..., "Lu"]
}
```

### `GET /api/v1/ligands`

Filtered, paginated summaries. All filters combine with AND; omit a parameter
to leave it unfiltered.

| Parameter | Meaning |
|---|---|
| `q` | case-insensitive name substring |
| `subtype` | v13 Subtype (validated against `/vocab`) |
| `mechanism` | e.g. `Solvation` |
| `trend` | full 7-value vocabulary — `Peaked` and `U_Shaped` included |
| `source` | source compilation, e.g. `Liu_2022_JACSAu`, where records carry one (the public release does not: `sources` in `/vocab` is then empty) |
| `elements` | comma-separated symbols the ligand must have measurements for, e.g. `Eu,La`. **Y and Sc work** (the filter reads `elements_measured`, not the fit-only `element_coverage` field) |
| `logd_min` / `logd_max` | keep ligands with at least one log D in range |
| `limit` / `offset` | pagination (default 50, max 200) |

```bash
curl 'localhost:8000/api/v1/ligands?subtype=DGA&trend=Ascending&elements=Eu,La'
```

Response envelope: `{"total", "offset", "limit", "count", "results": [...]}`.
An invalid filter value returns 422 with the allowed values listed.

### `GET /api/v1/ligands/{name}`

The full record: identity, selectivity block, performance summary,
data-quality block, metadata, and all experimental conditions.

Name resolution, in order (`matched_by` in the response says which applied):

1. exact `Name` (unique across the database),
2. exact `Alternative_Names` entry — pre-2026-07-28 names like `C5BTBP`
   still resolve after renames/merges,
3. unambiguous substring match.

Unknown names return 404 with up to 10 suggestions.

## Guarantees pinned by tests (`tests/test_api.py`)

- Counts always agree with the Python API (asserted against a live
  `LigandDatabase`, not hardcoded numbers).
- All 7 trend values are queryable — unlike `scripts/query_ligands.py
  --trend`, which still rejects `Peaked`/`U_Shaped`.
- Y/Sc element filtering works (regression-tested against the stale
  `element_coverage` bug, PROJECT_REVIEW_2026-07-31.md finding 2).
- Old ligand names resolve via `Alternative_Names`.

## Caveats

- The database loads once at startup; restart the server after rebuilding
  the JSONL.
- CORS is open (`*`) for GET — the data is public and read-only.
