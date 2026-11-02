---
name: reel-quest
description: Query and analyse REEL-QUEST, a curated database of rare-earth liquid-liquid (solvent) extraction data from the literature. It holds extractants (ligands) with distribution ratios for the lanthanides La-Lu plus Y and Sc under each experimental condition, with the source paper's DOI. Use it for questions about rare-earth or lanthanide separation, extractant performance or selectivity (TODGA and other diglycolamides, BTBP, CMPO, phosphine oxides, phenanthroline diamides, ...), a distribution ratio, a separation factor between two elements such as Eu/La or Nd/Pr, selectivity trends across the lanthanide series, or any custom analysis of this data.
---

# REEL-QUEST

Distribution ratios D = [M]org / [M]aq from published papers: one record per
extractant (ligand), each with its experimental conditions. Everything here
reads `SKILL_DIR/data/sx_ligands_merged_v14_flat.jsonl` through the bundled
`rel_quest` package, the code that produced the database's published
statistics. Use that code rather than parsing the JSONL yourself: the rules
at the end are built into it.

**`SKILL_DIR`** below means this skill's base directory, shown when the skill
loads. Replace it in every command and script.

## Setup

Needs Python 3.10+ with pandas, numpy, scipy and plotly. Check with
`python -c "import pandas, numpy, scipy, plotly"`. If that fails, ask the user
before running `pip install -r SKILL_DIR/requirements.txt`.

## Quick answers: query.py

```bash
python SKILL_DIR/scripts/query.py get_best_separators '{"element1": "Eu", "element2": "La"}'
```

It prints JSON. The tools are the same as the REEL-QUEST MCP server's:

| Tool | Arguments (defaults) |
|---|---|
| `get_database_statistics` | none |
| `search_ligands_by_name` | `name`, `limit` (15) |
| `get_ligand_details` | `name`: SMILES, subtype, trend, median log D per element |
| `get_ligand_conditions` | `name`, `element` (all), `limit` (50): every condition with D values and its DOI |
| `search_by_trend` | `trend`, `limit` (10) |
| `search_by_mechanism` | `mechanism`, `limit` (10) |
| `search_by_subtype` | `subtype` (DGA, Phosphine_oxide, DAPhen, ...), `limit` (10) |
| `get_best_separators` | `element1`, `element2`, `limit` (10): median log SF paired within conditions |
| `compare_ligands` | `names`, a list of 2-5 |
| `get_ligands_for_element` | `element` (Y and Sc too), `limit` (10) |

`python SKILL_DIR/scripts/query.py vocab` lists the valid trends, subtypes,
mechanisms and elements. A ligand that is not found exits 1; a bad tool or
value exits 2. In both cases the reason goes to stderr.

## Custom analysis

When a question needs a filter or a statistic those tools don't offer, write
Python on rel_quest's helpers. Start with the `sys.path` line. Short scripts
can run as `python - <<'EOF' ... EOF` from any directory. For longer ones, or
if your agent holds a heredoc for approval (some shell checks flag braces
next to quotes), save the script as a `.py` file and run `python file.py`. Each
recipe below runs on its own. `LigandDatabase()` logs "Loaded N ligands" to
stderr.

```python
# Recipe: Look up a ligand
import sys; sys.path.insert(0, "SKILL_DIR")
from rel_quest import LigandDatabase

db = LigandDatabase()
lig = db.query_by_name("TODGA")[0]          # exact name; query_by_name(x, fuzzy=True) for partial
print(lig.name, lig.subtype, lig.mechanism, lig.selectivity_trend)
print(f"{len(lig.experimental_conditions)} conditions; elements measured: {lig.tested_elements}")
cond = lig.experimental_conditions[0]       # one condition = one experiment in one paper
print(cond["Acid_Name"], cond["Acid_Concentration_M"], cond["Extractant_Concentration_M"],
      cond["Solvent_Name"], cond["Temperature_C"], cond["reference"]["DOI"])
element, cell = next(iter(cond["distribution_ratios"].items()))
print(element, cell["value"], cell["log_value"], cell.get("censor"))   # censor '>'/'<': a bound
```

```python
# Recipe: Paired separation factors under chosen conditions
import sys; sys.path.insert(0, "SKILL_DIR")
from rel_quest import LigandDatabase
from rel_quest.plotting import extract_separation_factor_data, split_bounds

db = LigandDatabase()
sf = extract_separation_factor_data(db, "Eu", "La")  # one row per condition that measured both
values, bounds = split_bounds(sf)                    # bounds: pairs where log SF is only > or <
pick = lambda df: df[(df.ligand_name == "TODGA") & (df.acid_name == "HNO3") & (df.acid_conc >= 3)]
rows = pick(values)
print(f"median log SF(Eu/La) = {rows.log_sf.median():.3f} over {len(rows)} conditions")
print(f"left out: {len(pick(bounds))} pairs that give only a bound")
# Pairs whose two bounds point the same way say nothing about log SF and are
# dropped before the frame is built; sf.attrs["no_information"] counts them (all ligands).
```

```python
# Recipe: Cite the papers behind a result
import sys; sys.path.insert(0, "SKILL_DIR")
from rel_quest import LigandDatabase
from rel_quest.plotting import extract_separation_factor_data, split_bounds

db = LigandDatabase()
values, _ = split_bounds(extract_separation_factor_data(db, "Eu", "La"))
rows = values[(values.ligand_name == "TODGA") & (values.acid_name == "HNO3") & (values.acid_conc >= 3)]
# condition_id is unique across the database; the frames carry it, not the DOI
reference = {c["condition_id"]: c["reference"]
             for lig in db.get_all_ligands() for c in lig.experimental_conditions}
for doi, n in rows.condition_id.map(lambda cid: reference[cid]["DOI"]).value_counts().items():
    print(f"https://doi.org/{doi}  ({n} conditions)")
```

```python
# Recipe: One element across ligands
import sys; sys.path.insert(0, "SKILL_DIR")
from rel_quest import LigandDatabase
from rel_quest.plotting import extract_element_data, split_bounds

db = LigandDatabase()
values, bounds = split_bounds(extract_element_data(db, "Y"))  # Y and Sc work; mixtures excluded
summary = values.groupby("ligand_name").log_d.agg(["median", "count"]).sort_values("median", ascending=False)
for name, row in summary.iterrows():
    print(f"{name:24s} {row['median']:.2f}  (median log D over {int(row['count'])} conditions)")
print("only bounds for Y:", sorted(set(bounds.ligand_name) - set(values.ligand_name)))
```

```python
# Recipe: Recompute a selectivity trend on a subset
import sys; sys.path.insert(0, "SKILL_DIR")
from rel_quest import LigandDatabase
from rel_quest.trends import classify_selectivity_trend

db = LigandDatabase()
conditions = db.query_by_name("TODGA")[0].experimental_conditions
strong = [c for c in conditions
          if c.get("Acid_Name") == "HNO3" and (c.get("Acid_Concentration_M") or 0) >= 3]
# the stored label comes from this call on all conditions, at these defaults
print("all conditions:", classify_selectivity_trend(conditions)["trend_classification"])
print("HNO3 >= 3 M:", classify_selectivity_trend(strong)["trend_classification"])
```

## Rules

- **Bounds.** A cell with `censor` `>` or `<` means D > x or D < x; it is
  not a value. Keep bounds out of medians, means, fits and separation factors
  (`split_bounds`), and say how many you left out.
- **Mixtures.** A condition with `Co_Extractants`
  (`rel_quest.uncertainty.is_mixture`) was measured with other extractants
  present. The `extract_*` helpers drop it. Keep it out of anything pooled.
- **Pair within a condition.** A separation factor is log D1 - log D2 inside
  one condition. Never divide medians taken from different conditions or
  papers.
- **Pooled medians mix conditions** that differ in acid, concentration and
  diluent. Say so, or filter to comparable conditions.
- **Missing parameters.** `extract_element_data` puts 0.0 for a missing acid
  concentration and 25 °C for a missing temperature. To filter on those, use
  `extract_separation_factor_data` (NaN when missing) or the condition dict.
- **Which elements were measured.** Use `lig.tested_elements`, which includes
  Y and Sc. `selectivity.element_coverage` lists only the La-Lu points the
  trend fit used.
- **Chemistry.** Classify by `lig.subtype` (21 subtypes by donor chemistry).
- **DGA names follow SAFE:** TODGA = tetra-n-octyl (C8), TDDGA = tetra-n-decyl
  (C10), TDdDGA = tetra-n-dodecyl (C12). Papers often call the C12 compound
  "TDDGA"; here it is TDdDGA.
- **Trends.** Ascending = D rises La to Lu (heavy-REE preference); Reverse =
  falls; Peaked / U_Shaped = a mid-series maximum / minimum; Mixed =
  conditions disagree.

## Citing

Data are CC BY 4.0. Cite the DOI of the paper that measured each value
(`condition["reference"]["DOI"]`; `Locator` gives its table or figure), and
the REEL-QUEST database (`SKILL_DIR/CITATION.cff`).

## Reference

`SKILL_DIR/reference/data_dictionary.csv` lists every field of a record: its
path, type, units, provenance and meaning. Read it when you need a field not
shown here.
