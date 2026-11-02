---
name: reel-quest-contribute
description: Draft or review a REEL-QUEST data contribution. Draft turns the solvent-extraction measurements a paper reports (distribution ratios D of the lanthanides, Y and Sc, from its tables, figures or supporting information) into contributions/<doi-slug>.json that passes the repository's contribution check. Review is the second read of such a file against its paper, with page-level evidence, before a maintainer applies it. Use when someone wants to contribute, add, submit or ingest a paper's data into REEL-QUEST, digitise distribution ratios from a figure, or review, verify or check a contribution file or data pull request against the paper.
compatibility: Runs in a clone of the REEL-QUEST repository with Python 3.10+ and its requirements (rdkit, pypdf, pypdfium2, Pillow).
---

# Contributing a paper's data to REEL-QUEST

REEL-QUEST takes new measurements one paper at a time: one file,
`contributions/<doi-slug>.json`, in the format of `contributions/README.md`.
This skill does one of two jobs on such a file:

- **Draft** it from the paper: read `references/draft.md`.
- **Review** it against the paper (the second read GOVERNANCE.md requires
  before anything is applied): read `references/review.md`.

Either job reads values off figures with `references/figures.md`.

The skill does not decide what enters the database. That is done by
`scripts/check_contribution.py` (run on every pull request), the person who
reviews the file, and the maintainer who applies it. Your job is to make the
file, or the review, right and easy to check.

## Setup

Work from the repository root. Every command below is run from there.

```bash
python -c "import rdkit, pypdf, pypdfium2, PIL"
```

If that fails, ask the person before running `pip install -r requirements.txt`.
The skill's scripts live in `.agents/skills/reel-quest-contribute/scripts/`:

| Script | Does |
|---|---|
| `scripts/figure.py` | `text` (a PDF's text by page), `render` (a page as PNG), `crop` (an enlarged region with a pixel grid), `read` (values from a figure calibration) |
| `scripts/ligand.py` | `lookup` (is this molecule or name already in the database?), `formula` (formula and masses of a SMILES), `subtypes` (the classification vocabulary) |
| `scripts/review.py` | `check` (validate a review and write its summary) |

Call them as `python .agents/skills/reel-quest-contribute/scripts/<script> ...`;
`--help` shows each one's arguments.

## Reading a paper

The text layer (`figure.py text`) finds sections, names and captions fast,
but it garbles tables and cannot show figures. For every value, render the
page (`figure.py render <pdf> <page>`) and look at the image. Open the PNG
with whatever lets you view images; zoom with `figure.py crop`.

## Rules for both jobs

These are the errors found in the data this database was built from. The
check cannot catch them; you have to.

1. **The paper is the only authority.** Never take a value or condition from
   memory, a review article, a compilation (SAFE, Liu et al. 2022, IDEaL) or
   this database. Look at the page.
2. **Molecule first.** Work out which compound each data series belongs to
   from the paper's scheme, name and compound number, and check the SMILES
   against the paper's formula, molecular weight or mass spectrum
   (`ligand.py formula`). Chain length, branching and substituents matter:
   most of the old errors were values filed under a sibling compound.
3. **The DOI is the paper that measured the values.** A series the paper
   takes from earlier work ("as reported previously [17]") was measured
   there: it is not this paper's. Leave it out, or flag it in `note` if it
   cannot be told apart.
4. **One condition is one experiment**: one extractant at one concentration,
   one aqueous phase at one acid concentration, one diluent. Never average
   or pool experiments. In a plot of D against [HNO3] with one curve per
   element, a condition is the set of points at one x.
5. **State only what the paper states.** A field the paper does not give is
   left out: no assumed 25 °C, no guessed concentration. "Room temperature"
   is not a number. Units are converted (mM to M) and the conversion said
   in `why`.
6. **Bounds are not values.** "D < 0.01", a point on the axis floor, or an
   arrow is a bound: `censor` `<` or `>` with `log_bound`.
7. **Two extractants in one organic phase** go in `co_extractants`.
8. **Leave the database alone.** Never edit `data/`, never run
   `apply_contributions.py`, never write `reviewed_by` or `status`. A person
   signs and a maintainer applies.
9. **Write JSON with a script** (`json.dump`), never by hand, and keep what
   you are unsure about in the file's `note` fields, where the reviewer reads
   it.

## Files you leave

| Job | File |
|---|---|
| Draft | `contributions/<doi-slug>.json`; one calibration per figure read, `contributions/figures/<doi-slug>/<locator>.json` |
| Review | `outputs/contributions/<doi-slug>_review.json` and its summary `<doi-slug>_review.md` |

`<doi-slug>` is the DOI with every character other than letters, digits and
dots replaced by `_`.

When you finish, tell the person in a few lines: what the file holds, what
you were unsure of, and the next step (draft: open a pull request with the
file; review: post the summary in the pull request for the person who signs).
