# Drafting a contribution from a paper

The result is `contributions/<doi-slug>.json` that passes
`python scripts/check_contribution.py contributions/<doi-slug>.json`, holds
every lanthanide (and Y, Sc) distribution ratio the paper itself measured,
and says in its `why` and `note` fields where each came from and what is
uncertain. The format, key by key, is `contributions/README.md`; the schema
is `schema/contribution.schema.json`. Read the README before you start.

## 1. The paper

- Find the DOI on the first page. It must be this paper's, bare and lower
  case (`10.1021/jacs.4c07332`), never a compilation's or a review's.
- Get the text: `python .agents/skills/reel-quest-contribute/scripts/figure.py text <pdf>`.
  Read the experimental section: reagents, how D was measured, temperature,
  contact time, diluent, acid, extractant concentration.
- If the paper's tables or figures are in supporting information you do not
  have, say so to the person and ask for the SI file before going further.

## 2. An inventory, before any value

List every place the paper reports D for La-Lu, Y or Sc: each table, each
figure panel, each sentence of text. For each: the locator ("Table 2",
"Figure 2b", "Text p. 6"), the extractant(s), the conditions that vary, and
whether this paper measured it (rule 3 of SKILL.md).

A series with a single element (Eu in an Am/Eu study, say) is in scope only
for a ligand the database already holds (step 3 tells you). A new ligand
needs at least one condition with two or more elements; if the paper gives a
new ligand only single-element series, leave that ligand out and list it in
the file's first condition `note` (GOVERNANCE.md §1). Leave out:

- D of other metals (Am, Cm, U, Fe ...), percentages extracted with no D
  (unless D can be derived exactly, below), and separation factors alone;
- isotherms, slope analyses and loading curves whose x axis is not a
  condition field (log D against log [L] gives D at stated [L]: those points
  are conditions; say so in `why`);
- data from cited papers.

One measurement shown twice (a figure point at the same condition as a table
value, a slope plot repeating a curve's point) is one condition: take the
table, or the clearer figure, and say in `note` where else it appears.

## 3. Each molecule

For every extractant in the inventory:

1. Write its SMILES from the paper's structure and name. Check it:
   `python .agents/skills/reel-quest-contribute/scripts/ligand.py formula '<SMILES>'`
   and compare the formula and mass with what the paper prints (elemental
   analysis, HRMS m/z for [M+H]+ is mass + 1.007). Record the comparison in
   the ligand's `reason`.
2. Is it already in the database?
   `ligand.py lookup --smiles '<SMILES>'` and `ligand.py lookup --name '<the paper's name>'`.
   - Found under `same_structure`: use that record's `Name` as the
     condition's `ligand`; do not declare it.
   - Found only under `stereoisomers`: a different molecule only if the paper
     states the other configuration; otherwise ask the person.
   - Found by name with a different structure: names collide in this field
     (TDDGA is the decyl compound here, TDdDGA the dodecyl one). Trust the
     structure; mention the collision in `note`.
   - Not found: declare it in `new_ligands`.
3. A new ligand needs `Name` (a short code in the database's style, e.g.
   DMDOHEMA, not the paper's "L3"; `ligand.py lookup --name` shows
   siblings), `SMILES`, `Subtype`, `Mechanism` and `Extractant_Type`
   (`ligand.py subtypes` gives the vocabulary and what records of each
   Subtype carry), and `reason`: the paper's compound number and name, and
   the formula check. Put the paper's own code in `Alternative_Names`.

## 4. Each condition

One entry per experiment (rule 4), ids `A0001`, `A0002`, ... in the order of
the inventory.

- `conditions`: only fields the paper states (rule 5), with the schema's
  names and units: `Acid_Name` (`HNO3`, `HCl`, ...), `Acid_Concentration_M`,
  `Extractant_Concentration_M`, `Solvent_Name` (the diluent, as the paper
  names it, with its abbreviation explained in `why`), `Solvent_Ratio`
  (only for a mixed diluent: its volume ratio, "7:3"; the organic:aqueous
  phase ratio has no field and goes in `why`), `Temperature_C`, `Contact_Time` ("30 min"), `Phase_Modifier_Name`,
  `Phase_Modifier_Concentration_M`. Salts or complexants in the aqueous
  phase go in `additives` as text ("1 M NaNO3").
- `reference.Locator`: the table or figure ("Table 2", "Figure 2", "Figure S3").
- `why`: which column, row, curve or symbol, on which page, and anything you
  converted. The reviewer reads this first.
- `note`: what you were unsure of, and anything odd in the paper (a value
  that disagrees with its own separation factor, a caption that contradicts
  the text). Null if nothing.

## 5. The values

| Source | Cell | `value_origin` |
|---|---|---|
| A table or the text printing D | `{"D": 1.1, "value_origin": "reported"}` (as printed; do not round) | `reported` |
| A table printing log D | `{"log_D": 0.04, "value_origin": "reported"}` | `reported` |
| A figure | from `figure.py read` (`references/figures.md`) | `digitised` |
| Computed exactly from what is printed (D from %E at a stated phase ratio) | `log_D`, the formula in `why` | `derived` |
| "< 0.01", "n.d.", a point on the axis floor, an arrow | `{"censor": "<", "log_bound": -2.0}` | none |

`sigma_log_D` only when the paper prints an uncertainty for that value, or
it comes from a figure calibration. A blanket "errors are about 10 %" is
not a per-value uncertainty: mention it in `note` instead.

Check every table value twice: read the rendered table, not the text layer,
and compare each element with its neighbours. A value that breaks the
smooth trend across the series is the most common transcription error;
confirm it on the page, and if the paper really prints it, say so in `note`.

## 6. Write and check

Build the file with a short Python script (`json.dump(..., indent=1,
ensure_ascii=False)`), with `proposed_by` naming the person, the date and
what was read ("A. Chemist, 2026-11-02, Table 2 and Figure 4"). Leave out
`reviewed_by` and `status`.

```bash
python scripts/check_contribution.py contributions/<doi-slug>.json
```

Add `--offline` when there is no network (the DOI is then not looked up).
Fix every problem it reports and run it again until it says OK. Read its
warnings: each one is either fixed or explained in a `note`. The report is
in `outputs/contributions/<doi-slug>_check.md`.

Then tell the person what the file holds (conditions, new ligands, which
tables and figures), what you were unsure of, and that the next step is a
pull request adding the file (and its `contributions/figures/` calibrations),
where it will be reviewed.
