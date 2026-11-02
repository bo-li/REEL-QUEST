---
title: REEL-QUEST
emoji: 🧪
colorFrom: indigo
colorTo: blue
sdk: gradio
sdk_version: 6.30.0
python_version: "3.12"
app_file: app.py
short_description: Rare-earth solvent-extraction data, with an MCP server
license: cc-by-4.0
pinned: false
---

# REEL-QUEST

Curated liquid-liquid extraction data for the rare-earth elements (La-Lu, Y,
Sc) from the published literature: ligands, the experimental conditions they
were measured under, and distribution ratios D = [M]org / [M]aq.

This Space has three faces:

- **The demo** (this page): separation factors for any element pair, every
  condition of a ligand with its source paper, and the ligands measured for
  one element.
- **An MCP server** at `https://<this Space's host>/mcp`, with the same ten
  read-only query tools the local server offers. In Claude: Settings →
  Connectors → Add custom connector, and paste that URL. The demo's
  "Use with an AI assistant" tab shows the exact URL.
- **A REST API** at `https://<this Space's host>/api/v1/`, plain JSON for
  scripts in any language; interactive documentation at `/api/docs`.

Data are CC BY 4.0: cite the paper behind each value (its DOI comes with every
condition) and the REEL-QUEST database (CITATION.cff). Code is MIT (LICENSE).

The files here are generated: `scripts/build_space.py` in the source
repository copies `space/` and the `rel_quest` package, `vocab/` and the
database beside it. Change them there, not here.
