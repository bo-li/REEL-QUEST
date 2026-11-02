# paper/

Everything the manuscript quotes or shows, generated from the database. Do not
edit the outputs by hand.

| File | Made by | Checked by |
|---|---|---|
| `stats.json` | `python paper/stats.py` | `tests/test_paper_stats.py` (must equal a fresh computation) |
| `figures/*.pdf`, `figures/*.png`, `figures/CAPTIONS.md` | `python paper/figures.py` (needs the `paper` extra: `pip install -e ".[paper]"`) | `tests/test_paper_stats.py` builds every figure |

`references.bib` is every reference the manuscript and its ESI cite, in BibTeX,
including each paper the database's values come from.
`tests/test_citation_metadata.py` checks that every `reference.DOI` in the
database is in it.


`stats.json` also carries the cross-source summary
(`outputs/validation/cross_source_agreement.json`).

**What this repository cannot recompute.** The database here carries no source
labels (which compilation each value came through) and no curation records, so
these `stats.json` keys are kept as computed in the maintainers' repository and
carried over by `paper/stats.py` (`DEV_ONLY_KEYS`): `ligands_by_data_source`,
`conditions_by_origin`, `d_values_by_origin`, `provenance_level_conditions`,
`n_merged_records`, `safe_conditions_without_concentrations` and `verification`.
Fig. S2 is drawn from the paired cells in
`outputs/validation/cross_source_agreement.csv`, which was written there.
Everything else is recomputed and compared by `tests/test_paper_stats.py`.

Main-text figures:

| Figure | Content |
|---|---|
| `fig1_pipeline` | What the database holds: data inputs (generic ligand cores, counts) → REEL-QUEST (one JSONL record per ligand, its contents, six access routes: dashboard, AI assistant, MCP server, agent skill, REST API, Python/CLI) → analytical capabilities |
| `fig2_platform` | The platform in use: (a) the dashboard's ligand browser, (b) its AI assistant answering a question, (c) one MCP tool call. Composed from `figures/platform/`, which `python paper/capture_platform.py` captures from the running dashboard and MCP server (needs Playwright and the `mcp` SDK; the assistant uses the free OpenRouter key in `.streamlit/secrets.toml`). Recapture after a data or software change; the captions read the model and date from `capture.json` |
| `fig3_overview` | Overview: ligands by subtype (coloured by donor chemistry, each with its generic core structure), trend-class donut, lanthanide coverage |
| `fig4_trend_gallery` | One representative ligand per trend class (Ascending, Descending, Flat, Peaked, U-shaped, Mixed); conditions aligned by their fitted offsets, median line with IQR and full-range bands; Mixed shows one band per vote with its acid range; picked by rule |
| `fig5_eu_la_showcase` | Violin-box plots of (a) log D(Eu), (b) log D(La) and (c) SF(Eu/La) (log axis) per condition, by subtype, each with its core structure; one paper's values stored on several variants count once |
| `fig6_uncertainty` | A ligand whose trend class and adjacent-pair SF signs its own σ can overturn (Monte Carlo, `stats.py` `perturbed_draws`); picked by rule |
| `fig7_trend_by_subtype` | Share of each subtype's ligands in each trend class (100 % bars, Fig. 3a row order, core structures) |
| `fig8_data_gaps` | Per subtype, where more experiments would change the answer: ligands without a class, lanthanides and conditions per ligand, class stability under σ, adjacent SF within 2σ of 1 (`stats.json` `by_subtype`) |

Supporting information: `figS1_ci_width` (slope CI width vs conditions),
`figS2_cross_source` (the same condition reported by two compilations),
`figS3_sf_viewer` (the dashboard's separation-factor viewer set to Nd/Pr, captured like Fig. 2:
`python paper/capture_platform.py --only sf_viewer`).

All plots have a black frame and no grid. "Descending" is the display name
of the stored class `Reverse`.

`stats.py` runs a seeded Monte Carlo over every fitted ligand
(`trend_label_stability`), so `tests/test_paper_stats.py` takes a couple of
minutes.

`figures/CAPTIONS.md` holds draft captions, written by the same run; their
numbers are computed from the records, never typed.


The trend palette is a colour-vision-deficiency-validated categorical set (see
the docstring of `figures.py`); Mixed and Insufficient data are neutral greys.

Regenerate both after any change to the database:

```bash
python paper/stats.py && python paper/figures.py
```
