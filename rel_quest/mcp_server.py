"""MCP server for the REEL-QUEST ligand database.

Exposes the read-only query tools of rel_quest.agent_tools (the same ones the
dashboard's AI assistant calls) over the Model Context Protocol, so Claude
Desktop, Claude Code or any MCP client can query the database directly.

Run (stdio, the transport local clients launch):
    python -m rel_quest.mcp_server          # or the `reel-quest-mcp` script

Run over HTTP (hosted, e.g. a Hugging Face Space; connector URL <base>/mcp):
    python -m rel_quest.mcp_server --transport http --host 0.0.0.0 --port 7860

HTTP mode is stateless and answers in plain JSON, so a restart (a free Space
waking from sleep) never strands a client's session and no proxy buffers a
stream. Bound to 127.0.0.1 (the default) the SDK refuses any Host header but
localhost; a container binds 0.0.0.0 to accept its public host name.

RELQUEST_DB names a database file other than the canonical one in data/ — for
example the JSONL from the Zenodo bundle when the package is installed without
the repository checkout.

The plot_* tools are not exposed: they return Plotly figures, which MCP
clients cannot render. See docs/MCP_SERVER.md.
"""

import argparse
import html
import json
import os
from typing import Any, Dict, List, Literal, Optional

from mcp.server.mcpserver import MCPServer
from mcp.server.mcpserver.exceptions import ToolError
from mcp_types import ToolAnnotations
from starlette.requests import Request
from starlette.responses import HTMLResponse

from . import __version__, agent_tools
from .core import LigandDatabase
from .schema import TREND_CLASSIFICATIONS, TREND_DIRECTIONS

INSTRUCTIONS = """\
REEL-QUEST: curated liquid-liquid extraction data for the rare-earth elements
(La-Lu, Y, Sc) from the published literature. Each ligand record holds
experimental conditions; each condition holds distribution ratios
D = [M]org / [M]aq for the elements measured under it.

Reading the data:
- A cell with censor '>' or '<' is a bound (D > x, D < x), not a value. Bounds
  are left out of every median and separation factor.
- A condition with co_extractants (is_mixture) was measured with other
  extractants present; it is left out of everything pooled across conditions.
- Medians pool conditions that differ in acid, concentration and diluent. For a
  real comparison use get_ligand_conditions and compare like with like.
- Separation factors pair both elements within one condition.
- Selectivity trends: Ascending = D rises La->Lu (heavy-REE preference),
  Reverse = falls La->Lu (light-REE preference), Peaked / U_Shaped = a mid-series
  maximum / minimum, Mixed = conditions disagree.
- DGA naming follows SAFE: TODGA = tetra-n-octyl (C8), TDDGA = tetra-n-decyl
  (C10), TDdDGA = tetra-n-dodecyl (C12). Many papers call the C12 compound
  "TDDGA"; that one is TDdDGA here.

Citing: data are CC BY 4.0. Cite the DOI of the paper that measured a value
(reference.DOI from get_ligand_conditions), and the REEL-QUEST database.
"""

# Periodic (Z) order: Sc 21, Y 39, then La-Lu.
ELEMENTS = ("Sc", "Y") + tuple(agent_tools.LANTHANIDES)
Element = Literal[ELEMENTS]  # type: ignore[valid-type]
Trend = Literal[TREND_CLASSIFICATIONS]  # type: ignore[valid-type]

READ_ONLY = ToolAnnotations(read_only_hint=True, destructive_hint=False,
                            idempotent_hint=True, open_world_hint=False)


def _checked(result):
    """Turn agent_tools' {"error": ...} into an MCP tool error."""
    if isinstance(result, dict) and "error" in result:
        raise ToolError(result["error"])
    return result


def _counted(values) -> Dict[str, int]:
    counts: Dict[str, int] = {}
    for value in values:
        if value:
            counts[value] = counts.get(value, 0) + 1
    return dict(sorted(counts.items(), key=lambda kv: (-kv[1], kv[0])))


def create_server(db_path: Optional[str] = None,
                  db: Optional[LigandDatabase] = None) -> MCPServer:
    """Build the server over one database file (default: RELQUEST_DB, then data/).

    `db` shares a database already loaded (the Space's demo, space/app.py).
    """
    if db is None:
        db = LigandDatabase(db_path or os.environ.get("RELQUEST_DB") or None)
    ligands = db.get_all_ligands()

    # Enums come from the loaded file, as in rel_quest.api: they can never
    # name a value the data lacks, nor miss one it holds.
    vocab = {
        "trends": {t: TREND_DIRECTIONS[t] for t in TREND_CLASSIFICATIONS},
        "subtypes": _counted(l.subtype for l in ligands),
        "mechanisms": _counted(l.mechanism for l in ligands),
        "data_sources": _counted(s for l in ligands for s in l.data_sources),
        "elements": [e for e in ELEMENTS if any(e in l.tested_elements for l in ligands)],
    }
    Subtype = Literal[tuple(vocab["subtypes"])]  # type: ignore[valid-type]
    Mechanism = Literal[tuple(vocab["mechanisms"])]  # type: ignore[valid-type]

    server = MCPServer(
        name="reel-quest",
        title="REEL-QUEST rare-earth ligand database",
        version=__version__,
        instructions=INSTRUCTIONS,
    )
    tool = server.tool

    @tool(annotations=READ_ONLY)
    def get_database_statistics() -> Dict[str, Any]:
        """Counts of ligands by family, mechanism, selectivity trend, data source and element."""
        return agent_tools.get_database_statistics(db)

    @tool(annotations=READ_ONLY)
    def search_ligands_by_name(name: str, limit: int = 15) -> List[Dict[str, Any]]:
        """Find ligands whose name contains `name` (case-insensitive), e.g. 'TODGA', 'BTBP', 'CMPO'."""
        return agent_tools.search_ligands_by_name(db, name, limit)

    @tool(annotations=READ_ONLY)
    def get_ligand_details(name: str) -> Dict[str, Any]:
        """One ligand: SMILES, subtype, mechanism, selectivity trend and median log D per lanthanide.

        Medians pool all of the ligand's conditions; use get_ligand_conditions
        for the individual measurements and their sources.
        """
        return _checked(agent_tools.get_ligand_details(db, name))

    @tool(annotations=READ_ONLY)
    def get_ligand_conditions(name: str, element: Optional[Element] = None,
                              limit: int = 50) -> Dict[str, Any]:
        """Every experimental condition of one ligand: acid, concentrations, diluent,
        temperature, D values (with bounds marked) and the DOI and table/figure of the source paper.

        `element` keeps only conditions that measured that element.
        """
        return _checked(agent_tools.get_ligand_conditions(db, name, element, limit))

    @tool(annotations=READ_ONLY)
    def search_by_trend(trend: Trend, limit: int = 10) -> Dict[str, Any]:
        """Ligands with a given selectivity trend across the lanthanide series
        (Ascending = heavy-REE preference, Reverse = light-REE preference)."""
        return agent_tools.search_by_trend(db, trend, limit)

    @tool(annotations=READ_ONLY)
    def search_by_mechanism(mechanism: Mechanism, limit: int = 10) -> Dict[str, Any]:
        """Ligands by extraction mechanism (Solvation = neutral donors, Cation_Exchange = acidic)."""
        return agent_tools.search_by_mechanism(db, mechanism, limit)

    @tool(annotations=READ_ONLY)
    def search_by_subtype(subtype: Subtype, limit: int = 10) -> Dict[str, Any]:
        """Ligands of one chemical subtype (DGA, Phosphine_oxide, DAPhen, Malonamide, ...)."""
        return agent_tools.search_by_extractant_type(db, subtype, limit)

    @tool(annotations=READ_ONLY)
    def get_best_separators(element1: Element, element2: Element,
                            limit: int = 10) -> Dict[str, Any]:
        """Ligands ranked by separation factor SF = D(element1)/D(element2).

        log SF is taken within each condition that measured both elements; the
        median over a ligand's paired conditions is reported with their count.
        """
        return agent_tools.get_best_separators(db, element1, element2, limit)

    @tool(annotations=READ_ONLY)
    def compare_ligands(names: List[str]) -> Dict[str, Any]:
        """Side-by-side summary of 2-5 ligands: family, mechanism, trend, median log D at key elements."""
        return agent_tools.compare_ligands(db, names)

    @tool(annotations=READ_ONLY)
    def get_ligands_for_element(element: Element, limit: int = 10) -> Dict[str, Any]:
        """Ligands with measurements for one element (Y and Sc included), highest median log D first."""
        return agent_tools.get_ligands_for_element(db, element, limit)

    @server.resource("reel-quest://vocab", name="vocabulary", mime_type="application/json",
                     description="Valid trend, subtype, mechanism, source and element values, with counts")
    def vocabulary() -> str:
        return json.dumps(vocab, ensure_ascii=False)

    @server.custom_route("/", methods=["GET"], include_in_schema=False)
    async def landing(request: Request) -> HTMLResponse:
        # Behind a TLS-terminating proxy (Hugging Face Spaces) the request
        # arrives as http; the URL to paste into a client is the public one.
        scheme = request.headers.get("x-forwarded-proto", request.url.scheme)
        host = request.headers.get("host", request.url.netloc)
        endpoint = html.escape(f"{scheme}://{host}/mcp")
        n_conditions = sum(len(l.experimental_conditions) for l in ligands)
        return HTMLResponse(LANDING.format(
            version=html.escape(__version__), n_ligands=len(ligands),
            n_conditions=n_conditions, endpoint=endpoint))

    return server


LANDING = """<!doctype html>
<html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>REEL-QUEST MCP server</title>
<style>body{{font:16px/1.5 system-ui,sans-serif;max-width:42rem;margin:2rem auto;padding:0 1rem}}
code{{background:#f2f2f2;padding:.1rem .3rem;border-radius:3px;word-break:break-all}}</style>
</head><body>
<h1>REEL-QUEST MCP server</h1>
<p>Rare-earth-element solvent-extraction data ({n_ligands} ligands, {n_conditions}
experimental conditions) for AI assistants, over the Model Context Protocol.
Version {version}.</p>
<p>Connector URL: <code>{endpoint}</code></p>
<p>In Claude: Settings &rarr; Connectors &rarr; Add custom connector, and paste the URL.
Any MCP client that speaks streamable HTTP works the same way.</p>
<p>Data are CC BY 4.0. Cite the paper behind each value (its DOI comes with
every condition) and the REEL-QUEST database.</p>
</body></html>
"""


def main(argv: Optional[List[str]] = None) -> None:
    parser = argparse.ArgumentParser(
        prog="reel-quest-mcp", description="REEL-QUEST MCP server (see docs/MCP_SERVER.md)")
    parser.add_argument("--transport", choices=("stdio", "http"), default="stdio",
                        help="stdio for a local client (default); http to host it")
    parser.add_argument("--host", default="127.0.0.1",
                        help="HTTP bind address (default 127.0.0.1; 0.0.0.0 in a container)")
    parser.add_argument("--port", type=int, default=int(os.environ.get("PORT", 8000)),
                        help="HTTP port (default $PORT, else 8000)")
    args = parser.parse_args(argv)

    server = create_server()
    if args.transport == "stdio":
        server.run()
    else:
        server.run("streamable-http", host=args.host, port=args.port,
                   stateless_http=True, json_response=True)


if __name__ == "__main__":
    main()
