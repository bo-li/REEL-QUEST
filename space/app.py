"""REEL-QUEST on a Hugging Face Space.

A Gradio demo over the ligand database, with the MCP server
(rel_quest/mcp_server.py, the same ten tools the local stdio server offers)
mounted at /mcp of the same app, and the read-only REST API
(rel_quest/api.py) at /api.

scripts/build_space.py copies this file beside rel_quest/, vocab/ and the
database, in the repository's own layout, and uploads the result. Testing and
deployment: docs/MCP_SERVER.md.

Why it is put together this way:
- A free account can host a Gradio Space only on ZeroGPU, and ZeroGPU will not
  start one without a @spaces.GPU function. Nothing here needs a GPU, so the
  one decorated function is never called and no quota is spent.
- ZeroGPU hooks gr.Blocks.launch, so the app starts through demo.launch(),
  not through a uvicorn of its own.
- Gradio's built-in MCP mode needs mcp<2, and the server here is written for
  mcp 2.x. Gradio's mode stays off; the SDK's /mcp route goes into Gradio's
  app instead, with the SDK's lifespan (it starts the session manager).
- The REST API is mounted whole at /api, not route by route, so it keeps its
  own CORS middleware and docs; Gradio serves nothing under /api/.
"""

import spaces  # first: on ZeroGPU it patches CUDA and Gradio at import

import os

import gradio as gr
import pandas as pd

from rel_quest import __version__, agent_tools, plotting
from rel_quest.api import create_app as create_rest_api
from rel_quest.core import LigandDatabase
from rel_quest.mcp_server import ELEMENTS, create_server
from starlette.routing import Mount


@spaces.GPU(duration=1)
def _zerogpu_requires_one_gpu_function() -> None:
    """Never called; see the module docstring."""


db = LigandDatabase()
LIGANDS = db.get_all_ligands()
NAMES = sorted((l.name for l in LIGANDS), key=str.lower)
MEASURED = [e for e in ELEMENTS if any(e in l.tested_elements for l in LIGANDS)]
N_CONDITIONS = sum(len(l.experimental_conditions) for l in LIGANDS)
N_VALUES = sum(len(c.get("distribution_ratios", {}))
               for l in LIGANDS for c in l.experimental_conditions)

# --- the MCP server ------------------------------------------------------------
# Stateless, with plain JSON replies, as `--transport http` serves it: a Space
# restarted from sleep strands no client session, and no proxy buffers a
# stream. The SDK's DNS-rebinding guard follows Gradio's bind address: on for
# loopback (a laptop), off for 0.0.0.0 (the Space, whose proxy sends the
# public host name).
mcp_app = create_server(db=db).streamable_http_app(
    host=os.environ.get("GRADIO_SERVER_NAME", "127.0.0.1"),
    stateless_http=True, json_response=True)
MCP_ROUTES = [route for route in mcp_app.routes if getattr(route, "path", None) == "/mcp"]
if len(MCP_ROUTES) != 1:
    raise RuntimeError(f"expected one /mcp route from the MCP SDK, got {MCP_ROUTES}")

# --- the REST API ----------------------------------------------------------------
# The same database object, mounted at /api: its paths are /api/v1/..., its
# docs /api/docs, exactly as `uvicorn rel_quest.api:app` serves them.
REST_API = Mount("/api", app=create_rest_api(db=db, mount=""))


def base_url() -> str:
    host = os.environ.get("SPACE_HOST")  # set by Hugging Face on a Space
    if host:
        return f"https://{host}"
    return f"http://127.0.0.1:{os.environ.get('GRADIO_SERVER_PORT', '7860')}"


def connector_url() -> str:
    return f"{base_url()}/mcp"


# --- the demo --------------------------------------------------------------------

def _checked(result):
    """Turn agent_tools' {"error": ...} into a Gradio error."""
    if isinstance(result, dict) and "error" in result:
        raise gr.Error(result["error"])
    return result


SEPARATOR_COLUMNS = ["Ligand", "Subtype", "median log SF", "SF", "min log SF", "max log SF",
                     "paired conditions", "bound pairs", "prefers"]


def best_separators(element1: str, element2: str, limit: int = 25):
    """Ligands ranked by separation factor SF = D(element1) / D(element2).

    log SF is taken within each condition that measured both elements, and the
    median over a ligand's paired conditions is shown with their count. Bounds
    (D > x, D < x) and mixed-extractant conditions stay out of the median.

    Args:
        element1: first element, e.g. "Eu".
        element2: second element, e.g. "La".
        limit: how many ligands to list.
    """
    if element1 == element2:
        raise gr.Error("Pick two different elements.")
    result = agent_tools.get_best_separators(db, element1, element2, int(limit))
    table = pd.DataFrame(
        [[r["name"], r["ligand_family"], r["median_log_SF"], r["SF"], r["min_log_SF"],
          r["max_log_SF"], r["n_paired_conditions"], r["n_bound_pairs"], r["prefers"]]
         for r in result["top_separators"]],
        columns=SEPARATOR_COLUMNS)
    summary = (f"**{result['total_with_both_elements']} ligands** were measured for both "
               f"{element1} and {element2} under at least one condition. {result['note']}")
    figure = plotting.plot_separation_factor(db, element1, element2) if len(table) else None
    return summary, table, figure


CONDITION_COLUMNS = ["Acid", "[acid] / M", "[extractant] / M", "Diluent", "T / °C", "Mixture"]


def _log_d(cell) -> str:
    if cell["log_D"] is None:
        return ""
    value = f"{cell['log_D']:.2f}"
    return f"{cell['censor']} {value}" if cell["censor"] else value


def ligand(name: str, element: str = "All"):
    """Every experimental condition of one ligand: acid, concentrations, diluent,
    temperature, log D per element (bounds marked > or <) and the source paper.

    Args:
        name: ligand name, e.g. "TODGA"; a partial name picks the closest match.
        element: "All", or an element symbol to keep only the conditions that measured it.
    """
    details = _checked(agent_tools.get_ligand_details(db, name))
    result = _checked(agent_tools.get_ligand_conditions(
        db, details["name"], None if element == "All" else element, limit=1_000_000))
    conditions = result["conditions"]
    elements = [e for e in ELEMENTS if any(e in c["D"] for c in conditions)]

    rows, papers = [], {}
    for cond in conditions:
        doi = cond["reference"]["DOI"]
        rows.append([
            cond["acid"], cond["acid_conc_M"], cond["extractant_conc_M"], cond["solvent"],
            cond["temperature_C"], "yes" if cond["is_mixture"] else "",
            *[_log_d(cond["D"][e]) if e in cond["D"] else "" for e in elements],
            doi or "", cond["reference"]["Locator"] or "", cond["condition_id"],
        ])
        if doi:
            papers[doi] = papers.get(doi, 0) + 1
    columns = CONDITION_COLUMNS + elements + ["DOI", "Where in the paper", "Condition"]
    table = pd.DataFrame(rows, columns=columns)

    asked = name.strip()
    matched = "" if details["name"].lower() == asked.lower() else f" (closest match to “{asked}”)"
    trend = f"selectivity trend **{details['selectivity_trend']}**"
    if details.get("trend_confidence"):
        trend += f" ({details['trend_confidence']} confidence)"
    shown = (f"{len(conditions)} conditions" if element == "All" else
             f"{len(conditions)} of {details['total_conditions']} conditions measured {element}")
    others = (f" Other names that match: {', '.join(details['other_matches'])}."
              if details["other_matches"] else "")
    summary = (
        f"### {details['name']}{matched}\n"
        f"{details['subtype']} · {details['mechanism']} · {trend}"
        + (f" · record from {', '.join(details['data_sources'])}" if details.get("data_sources") else "")
        + "\n\n"
        f"SMILES `{details['smiles']}`\n\n"
        f"{shown}, log D per element. `>` and `<` mark bounds (D > x, D < x), which no "
        f"median uses; *Mixture* marks a condition measured with other extractants "
        f"present.{others}\n\n"
        f"**Sources** (cite these; conditions from each): " + "; ".join(
            f"[{doi}](https://doi.org/{doi}) ({n})" for doi, n in papers.items())
    )
    figure = plotting.plot_selectivity_trend(db, details["name"])
    return summary, table, figure


def ligands_for_element(element: str, limit: int = 50):
    """Ligands with measurements for one element (Y and Sc included), highest median log D first.

    Medians pool conditions that differ in acid, concentration and diluent;
    compare like with like in the Ligand tab.

    Args:
        element: element symbol, e.g. "Y".
        limit: how many ligands to list.
    """
    result = agent_tools.get_ligands_for_element(db, element, int(limit))
    key = f"median_logD_{result['element']}"
    table = pd.DataFrame(
        [[r["name"], r["ligand_family"], r["selectivity_trend"], r[key]] for r in result["results"]],
        columns=["Ligand", "Subtype", "Trend", f"median log D({element})"])
    summary = f"**{result['total_ligands_with_data']} ligands** have measurements for {element}."
    return summary, table


HEADER = f"""# REEL-QUEST
Curated liquid-liquid extraction data for the rare-earth elements from the
published literature: **{len(LIGANDS)} ligands**, {N_CONDITIONS:,} experimental
conditions and {N_VALUES:,} distribution ratios D = [M]org / [M]aq. Version {__version__}.

Data are CC BY 4.0. Cite the paper behind each value (the Ligand tab gives
every condition's DOI) and the REEL-QUEST database.
"""

CONNECT = f"""
### Connector URL

`{connector_url()}`

This Space is also an MCP server, with the same ten read-only tools as the
local one: search by name, trend, mechanism or subtype; a ligand's details
and conditions; best separators for an element pair; compare ligands; the
ligands measured for an element.

- **Claude**: Settings → Connectors → Add custom connector, and paste the URL.
- **Claude Code**: `claude mcp add --transport http reel-quest {connector_url()}`
- Any client that speaks streamable HTTP works the same way.

If this Space is private, clients must also send
`Authorization: Bearer <your Hugging Face token>` (Claude Code:
`--header "Authorization: Bearer hf_..."`).

### REST API

Plain HTTP and JSON, for scripts in any language: `{base_url()}/api/v1/`.
Interactive documentation: `{base_url()}/api/docs`. For example,
`{base_url()}/api/v1/ligands?subtype=DGA&elements=Eu,La` lists the
diglycolamides with both Eu and La measured, and
`{base_url()}/api/v1/ligands/TODGA` gives one ligand's full record.

### Reading the data

- A bound (D > x, D < x) is not a value. It is left out of every median and
  separation factor.
- A condition measured with other extractants present (a mixture) is left
  out of everything pooled across conditions.
- Separation factors pair both elements within one condition.
- DGA names follow SAFE: TODGA = tetra-n-octyl (C8), TDDGA = tetra-n-decyl
  (C10), TDdDGA = tetra-n-dodecyl (C12), the compound many papers call "TDDGA".
"""

sf_initial = best_separators("Eu", "La")
ligand_initial = ligand("TODGA")
element_initial = ligands_for_element("Y")

with gr.Blocks(title="REEL-QUEST") as demo:
    gr.Markdown(HEADER)

    with gr.Tab("Separation factors"):
        with gr.Row():
            sf_element1 = gr.Dropdown(MEASURED, value="Eu", label="Element 1")
            sf_element2 = gr.Dropdown(MEASURED, value="La", label="Element 2")
            sf_limit = gr.Slider(5, len(LIGANDS), value=25, step=1, label="Ligands to list")
        sf_button = gr.Button("Rank ligands", variant="primary")
        sf_summary = gr.Markdown(sf_initial[0])
        sf_table = gr.DataFrame(sf_initial[1], interactive=False, max_height=420)
        sf_plot = gr.Plot(sf_initial[2])
        gr.Examples([["Eu", "La"], ["Nd", "Pr"], ["Dy", "Tb"], ["Lu", "Yb"]],
                    inputs=[sf_element1, sf_element2])
        sf_button.click(best_separators, [sf_element1, sf_element2, sf_limit],
                        [sf_summary, sf_table, sf_plot], api_name="best_separators")

    with gr.Tab("Ligand"):
        with gr.Row():
            lig_name = gr.Dropdown(NAMES, value="TODGA", label="Ligand", filterable=True,
                                   allow_custom_value=True)
            lig_element = gr.Dropdown(["All"] + MEASURED, value="All",
                                      label="Only conditions that measured")
        lig_button = gr.Button("Show conditions", variant="primary")
        lig_summary = gr.Markdown(ligand_initial[0])
        lig_plot = gr.Plot(ligand_initial[2])
        lig_table = gr.DataFrame(ligand_initial[1], interactive=False, max_height=520)
        gr.Examples([["TODGA", "All"], ["CyMe4-BTBP", "All"], ["TBP", "Eu"], ["DBPhthA", "All"]],
                    inputs=[lig_name, lig_element])
        lig_button.click(ligand, [lig_name, lig_element], [lig_summary, lig_table, lig_plot],
                         api_name="ligand")

    with gr.Tab("By element"):
        with gr.Row():
            el_element = gr.Dropdown(MEASURED, value="Y", label="Element")
            el_limit = gr.Slider(5, len(LIGANDS), value=50, step=1, label="Ligands to list")
        el_button = gr.Button("List ligands", variant="primary")
        el_summary = gr.Markdown(element_initial[0])
        el_table = gr.DataFrame(element_initial[1], interactive=False, max_height=520)
        el_button.click(ligands_for_element, [el_element, el_limit], [el_summary, el_table],
                        api_name="ligands_for_element")

    with gr.Tab("Use with an AI assistant"):
        gr.Markdown(CONNECT)


if __name__ == "__main__":
    demo.launch(
        mcp_server=False,  # Gradio's MCP mode needs mcp<2, even if GRADIO_MCP_SERVER is set
        ssr_mode=False,    # Python serves every route; with SSR, Node fronts it
        app_kwargs={"routes": [*MCP_ROUTES, REST_API], "lifespan": mcp_app.router.lifespan_context},
    )
