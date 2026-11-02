# Reviewing a contribution against its paper

GOVERNANCE.md's review rule: every condition is read against its paper by
someone other than its transcriber before it is applied, and an AI read is
followed by a person's spot-check. You are that AI read. You give each
condition and each new ligand a verdict with page-level evidence; the person
checks your verdicts and signs `reviewed_by`. **A wrong "Correct" is far
worse than an honest "Unsure"**: be conservative.

Never edit the contribution file. Your output is
`outputs/contributions/<doi-slug>_review.json`.

## Before you start

- Run `python scripts/check_contribution.py contributions/<doi-slug>.json`
  (`--offline` without a network). If it fails, say so: the file is not
  ready for review. Its warnings are worth reading either way.
- Read the paper independently first: its experimental section and the
  tables and figures the file cites (`figure.py text`, then `figure.py
  render` each cited page and look at it). Know what the paper says before
  you read what the file claims.

## What to check

**Each new ligand.** Is it the compound the paper measured? Compare the
SMILES's formula and mass (`ligand.py formula`) with the paper's, and its
chain lengths, branching and substituents with the scheme and name. Is it
already in the database (`ligand.py lookup --smiles`, `--name`)? Does it
have at least one condition with two or more elements? A new ligand known
only from single-element series is out of scope (GOVERNANCE.md §1): say so
in its `evidence` and give it `Unsure`, for the maintainer to decide.

**Each condition.**

1. **Molecule**: the series belongs to the extractant the condition names.
2. **Measured here**: the paper measured it, not a cited paper.
3. **Conditions**: acid and its concentration, extractant concentration,
   diluent, temperature, contact time: as the paper states them, and none
   the paper does not state.
4. **Values**, at the locator:
   - a table: every element; agreement means within 0.02 dex of log10 of
     the printed D;
   - a figure: redraw the contributor's calibration
     (`figure.py read contributions/figures/<doi-slug>/<file>.json --overlay
     outputs/contributions/<name>.png`) and look whether each cross sits on
     its marker; without a calibration, make your own (`references/figures.md`).
     Agreement means within 0.15 dex on a log axis (15 % of D on a linear
     one); more than 0.3 dex is wrong; between is `Unsure` unless the
     reading is clean. Check at least the lightest, the heaviest and the
     extreme points, and every value that breaks its neighbours' trend;
   - bounds are bounds in the paper too, and the value_origin fits the
     source (`reported` table or text, `digitised` figure).

## Verdicts

One per condition, exactly one of:

| Verdict | When |
|---|---|
| `Correct` | molecule, conditions and every value checked agree |
| `Wrong value` | a value disagrees (give the paper's in `correction`, e.g. `{"log_D_Eu": 1.05}`) |
| `Wrong molecule` | the values belong to another compound (name it; its SMILES in `correction`) |
| `Wrong condition` | a stated condition disagrees, or one is stated that the paper does not give (`correction`: `{"Acid_Concentration_M": 5.0}`) |
| `Wrong acid` | the aqueous acid is another acid |
| `Not measured here` | the paper took these values from another paper |
| `Not in paper` | the values are nowhere in the paper or the SI you have |
| `Duplicate` | another condition of the file, or of the database, is the same measurement |
| `Unsure` | you cannot decide: say what is missing (always `low` confidence) |

Several things wrong: the most serious (molecule > acid > condition >
value), all listed in `correction`. For a new ligand: `Correct`,
`Wrong molecule`, `Already in the database` or `Unsure`.

Confidence: `high` when you saw the evidence directly and it is
unambiguous; `medium` when a reading was imprecise or only part could be
checked; `low` otherwise.

## The review file

Write it with a Python script (`json.dump`, `indent=1`):

```json
{
 "contribution": "contributions/<doi-slug>.json",
 "contribution_sha256": "<sha256 of that file, as reviewed>",
 "doi": "<the file's DOI>",
 "pages_viewed": [3, 7],
 "paper_notes": "what the paper measured, where, and anything inconsistent in it",
 "new_ligands": [
  {"Name": "<as in the file>", "verdict": "Correct", "confidence": "high",
   "evidence": "p.3 Scheme 1, compound 4b; Table 1 prints C30H52N2O4, MW 504.75 = the SMILES's"}
 ],
 "conditions": [
  {"id": "A0001", "verdict": "Correct", "confidence": "high",
   "elements_checked": ["La", "Eu", "Lu"],
   "evidence": "p.7 Table 2, row 3 M HNO3: D(La) 1.9, D(Eu) 12.4 ... D(Lu) 85; all match",
   "correction": {}, "note": ""}
 ]
}
```

`evidence`: page, table or figure, what you read. Short, specific, factual.
Never claim a reading you did not make.

Then run

```bash
python .agents/skills/reel-quest-contribute/scripts/review.py check outputs/contributions/<doi-slug>_review.json
```

and fix whatever it reports. It writes the summary
`outputs/contributions/<doi-slug>_review.md` to post in the pull request.
Tell the person the verdict counts, every condition that is not `Correct`,
and that a person must still spot-check and sign `reviewed_by`.
