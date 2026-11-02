# MCP server

`rel_quest/mcp_server.py` serves the database over the
[Model Context Protocol](https://modelcontextprotocol.io), so Claude Code,
Claude Desktop or any MCP client can query it directly. The tools are the
read-only query functions of `rel_quest/agent_tools.py`, the same ones the
dashboard's AI assistant calls.

## Install and run (from a checkout)

```bash
git clone https://github.com/bo-li/REEL-QUEST.git && cd REEL-QUEST
pip install -e ".[mcp]"          # the MCP SDK (mcp>=2.0) is an optional extra
python -m rel_quest.mcp_server   # or: reel-quest-mcp
```

The server speaks stdio: a client launches it as a subprocess and talks to it
over stdin/stdout. Running it by hand just waits for input, so use Ctrl-C to stop it.
`RELQUEST_DB=<file.jsonl>` serves a different database file, e.g. a candidate
build, as the dashboard does.

## Connect a client

**Claude Code.** The repository's `.mcp.json` registers the server for this
project. On the next `claude` start in the checkout, approve `reel-quest` when
asked. The command is `python`, so the interpreter on `PATH` must have
the `mcp` extra installed. Elsewhere:

```bash
claude mcp add reel-quest -- /path/to/python -m rel_quest.mcp_server
```

`/path/to/python` is the interpreter that ran `pip install -e ".[mcp]"`. The
editable install is what lets it import `rel_quest` and find `data/` from any
working directory.

**Claude Desktop** (`claude_desktop_config.json`, Settings → Developer):

```json
{"mcpServers": {"reel-quest": {
  "command": "/path/to/python",
  "args": ["-m", "rel_quest.mcp_server"]
}}}
```

On Windows with the checkout and Python inside WSL, launch through `wsl.exe`:

```json
{"mcpServers": {"reel-quest": {
  "command": "wsl.exe",
  "args": ["bash", "-lc", "cd /mnt/c/.../Database && python -m rel_quest.mcp_server"]
}}}
```

## Tools

| Tool | Returns |
|---|---|
| `get_database_statistics` | counts by family, mechanism, trend, source, element |
| `search_ligands_by_name` | ligands whose name contains a string |
| `get_ligand_details` | identity, subtype, trend, median log D per lanthanide, 3 sample conditions |
| `get_ligand_conditions` | every condition: acid, concentrations, diluent, temperature, D per element (bounds marked), DOI and table/figure |
| `search_by_trend` | ligands with one of the 7 trend classifications |
| `search_by_mechanism` | Solvation / Cation_Exchange |
| `search_by_subtype` | one of the populated subtypes (DGA, Phosphine_oxide, ...) |
| `get_best_separators` | ligands ranked by median log SF over conditions that measured both elements |
| `compare_ligands` | side-by-side summary of 2-5 ligands |
| `get_ligands_for_element` | ligands with data for one element, Y and Sc included |

The resource `reel-quest://vocab` lists the valid enum values with their counts.
The trend, mechanism, subtype and element parameters are enums in the tool
schemas, so a client sees the vocabulary up front. Mechanism and subtype values
are read from the loaded file, so they always match the data.
All tools are marked read-only. The `plot_*` tools are left out because they
return Plotly figures, which MCP clients cannot render.

The server's `instructions` give the client model the conventions it needs:
- bounds (`censor` `>`/`<`) and mixed-extractant conditions stay out of medians;
- separation factors are paired within one condition;
- DGA names follow SAFE (TODGA C8, TDDGA C10, TDdDGA C12);
- cite `reference.DOI`, since the data are CC BY 4.0.

## Run over HTTP (hosting)

```bash
python -m rel_quest.mcp_server --transport http                 # http://127.0.0.1:8000/mcp
python -m rel_quest.mcp_server --transport http --host 0.0.0.0 --port 7860   # in a container
```

`--port` defaults to `$PORT`, then 8000. The connector URL is `<base>/mcp`, and
`GET /` serves a short page that names it. To try the HTTP mode from Claude Code:
`claude mcp add --transport http reel-quest-http http://127.0.0.1:8000/mcp`.

- **Stateless, plain JSON.** Every request stands alone and gets a JSON reply,
  not a stream. A restarted server (for example a free Space waking from sleep)
  answers a client's old session id normally, and no proxy can hold back a stream.
- **Host header.** Bound to `127.0.0.1` (the default), the SDK refuses any `Host`
  but localhost (421), which guards a laptop against DNS rebinding. A container
  binds `0.0.0.0` so that its public host name is accepted.
- **Behind a TLS proxy**, the landing page reads `X-Forwarded-Proto`, so it shows
  the `https://` URL users should paste.
- There is no authentication; the server is read-only over public data.

### Hugging Face Space

`space/app.py` is a Gradio app for a Space. It serves a demo of the database:
separation factors for an element pair, every condition of a ligand with its
source papers, and the ligands measured for one element. It also mounts this
server's `/mcp` route in the same app, with the same ten tools, and the REST
API (`docs/REST_API.md`) at `/api`.
`scripts/build_space.py` copies `space/`, `rel_quest/`, `vocab/`, the database
and the licence files into a bundle, in the checkout's own layout, and can
upload it. The Space repo is generated: change `space/`, not the Space.

**Why Gradio on ZeroGPU, not Docker.** As of October 2026, Hugging Face
charges for Docker and Gradio Spaces on `cpu-basic`: they need PRO (about
$9 a month). A free personal account can host up to two ZeroGPU Spaces, and
those are Gradio only
([ZeroGPU docs](https://huggingface.co/docs/hub/spaces-zerogpu)). ZeroGPU will
not start a Space without a `@spaces.GPU` function. The app has one that is
never called, so it spends no GPU quota. On PRO, the same app runs on
`cpu-basic`, where the `spaces` package does nothing.

**How the server is mounted.**
- Gradio's built-in MCP mode is not used. It needs `mcp<2`, and this server
  is written for `mcp` 2.x. Against 2.x, Gradio 6.30 cannot build its own
  server and quietly turns the mode off; the app also passes
  `mcp_server=False`.
- The SDK's `/mcp` route goes into Gradio's app through `launch(app_kwargs=...)`,
  along with the SDK's lifespan, which starts the session manager.
- `ssr_mode=False`, and `GRADIO_SSR_MODE=false` set when the Space is created.
  With SSR on, Node fronts the app and `/mcp` would get its page instead.
- The DNS-rebinding guard follows `GRADIO_SERVER_NAME`, as `--host` does for
  the standalone server.

**Status (2026-10-13): tested locally, not deployed.** Local tests used
Python 3.12, Gradio 6.30.0 and mcp 2.3.0 (see Testing). Creating
`bolivu/reel-quest` as a **private** ZeroGPU Space on the free account was
refused: *402 "You must be subscribed to PRO to host Spaces with ZeroGPU. If
you recently created your account, please wait 30 days or request a community
grant."* The account meets the documented free conditions: verified email,
created 2026-05-01. HF's docs don't say whether private ZeroGPU Spaces need
PRO. So either subscribe to PRO, or create the Space public at release.

**Deploy:**

```bash
hf auth login                                   # pip install huggingface_hub
hf repos create <owner>/reel-quest --type space --space-sdk gradio \
    --flavor zero-a10g --env GRADIO_SSR_MODE=false [--private]   # PRO: cpu-basic works too
python scripts/build_space.py --out build/space --upload <owner>/reel-quest
hf spaces logs <owner>/reel-quest --build --follow
hf spaces info <owner>/reel-quest --expand runtime              # wait for RUNNING
```

`--upload` refuses a Space that does not exist yet. `hf upload` would
otherwise create it, public and without the settings above.

Then run the test battery against the live Space (Testing below). The
connector URL is `https://<owner>-reel-quest.hf.space/mcp`; the demo's
"Use with an AI assistant" tab shows it. Paste it into Claude → Settings →
Connectors → Add custom connector, or run
`claude mcp add --transport http reel-quest <url>`.

Before relying on it:
- **First deploy: check which `mcp` pip installed.** Gradio's `mcp` extra pins
  `mcp<2`, and the docs don't say whether HF's Gradio image installs it. If it
  does, the build log shows a conflict with `mcp>=2.0`, or the run log shows
  `from mcp.server.mcpserver import MCPServer` failing. The fix would be an
  mcp 1.x fallback import in `rel_quest/mcp_server.py` (`FastMCP`,
  `mcp.server.fastmcp.exceptions.ToolError`, `mcp.types.ToolAnnotations`).
  None of that is written yet.
- **Private Space.** Every request needs `Authorization: Bearer <HF token>`.
  Claude Code can send it (`claude mcp add ... --header "Authorization: Bearer hf_..."`);
  claude.ai custom connectors need a public Space.
- **Sleep.** A free Space sleeps after 48 h without traffic. The first
  request after that waits for the app to start, which can exceed a client's
  timeout. Opening the Space page wakes it.
- **Data updates.** Run `scripts/build_space.py --upload` again. Every upload
  replaces the database and code, and the commit message names the source
  revision.
- **Elsewhere.** `--transport http --host 0.0.0.0` (the section above) serves
  the same tools from any container host. The Space is the free option.

## Installing without the checkout (not yet supported)

A non-editable install (`pip install "reel-quest[mcp] @ git+https://github.com/bo-li/REEL-QUEST"`,
or `uvx`) imports fine, but the data are not package data (see `pyproject.toml`).
`RELQUEST_DB` can name the JSONL from the Zenodo bundle. However, `LigandDatabase` and
the plotting helpers still read `vocab/` (elements, solvents) from beside the
package, so most tools fail. `tests/test_mcp_server.py` pins this as an
expected failure. Closing the gap needs one of:

1. **A data-directory setting.** Point the package at a directory holding the
   JSONL and `vocab/`, which is exactly what the unzipped Zenodo bundle contains.
   This keeps code (MIT) and data (CC BY 4.0) in separate artifacts.
2. **Ship `vocab/` and the JSONL as package data**, so that
   `uvx --from git+https://github.com/bo-li/REEL-QUEST reel-quest-mcp` works in one line.
   This reverses the packaging decision recorded in `pyproject.toml`.

The base dependencies (streamlit, playwright, rdkit, ...) make either install
heavy; a slim MCP-only install would need them moved into extras.

## Testing

`tests/test_mcp_server.py` launches the server as a subprocess and drives it
with the MCP client SDK in the classic `initialize` handshake. It checks
every answer against the raw JSONL, and that stdout carries nothing but JSON-RPC.
The same checks then rerun over HTTP. Further HTTP checks cover:
- the 2026-07-28 discovery handshake;
- a stale session id;
- JSON (not SSE) replies;
- the Host-header rules;
- the landing page's https URL, through the machine's non-loopback address
  as a proxy would reach it.

The Space checks (failure modes S1-S10 in the test file) need a Python like
the Space's, named by `RELQUEST_SPACE_PYTHON`; without it they skip:

```bash
uv venv --python 3.12 build/space-venv       # build/ is gitignored
uv pip install --python build/space-venv/bin/python gradio==6.30.0 spaces -r space/requirements.txt
RELQUEST_SPACE_PYTHON=build/space-venv/bin/python python -m pytest tests/test_mcp_server.py
```

Keep the Gradio version equal to the card's `sdk_version`; the tests check
that. They:
- build the bundle with `scripts/build_space.py` and check it holds exactly
  what the app reads;
- check the card and requirements;
- launch `app.py` from the bundle as a Space would (bound to 0.0.0.0, with
  `SPACE_HOST` and a hostile `GRADIO_MCP_SERVER=true`);
- rerun the whole battery against its `/mcp`;
- check that the Space host name is accepted, that `/mcp` belongs to the MCP
  transport, and that the page shows the https connector URL;
- check that each demo tab agrees with the MCP tool it is built on.

After a deploy, the same battery runs against the live Space, through
Hugging Face's proxy:

```bash
RELQUEST_SPACE_URL=https://<owner>-reel-quest.hf.space HF_TOKEN=$(hf auth token) \
    python -m pytest tests/test_mcp_server.py -k LiveSpace
```

That is the only check that ZeroGPU accepted the Space (S10) and that SSR is
really off on Hugging Face (S4).

Each run writes the request/response transcripts to
`outputs/mcp_e2e_transcript*.json` (stdio, http, space, live_space).
