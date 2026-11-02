# Claude Code skill

`skills/reel-quest/` is a skill for REEL-QUEST in the
[Agent Skills](https://agentskills.io) format, which Claude Code, Codex and
Gemini CLI all read. The MCP server (`docs/MCP_SERVER.md`) answers ten fixed
questions with the processing enforced. The skill gives an agent that runs
code the database and `rel_quest`'s canonical helpers, so it can answer new
questions without reinventing the processing. Examples: a separation factor
restricted to one acid range, a trend recomputed on a subset, or a per-element
summary.

| | MCP server | Skill |
|---|---|---|
| Works in | any MCP client, including chat without code execution | agents that run code (Claude Code, Codex, Gemini CLI) |
| Questions | the 10 tools | anything, built on the canonical helpers |
| Canonical processing | enforced: the agent never sees raw values | guided: SKILL.md says which helpers to use and the rules, but an agent can still go around them |
| Updates | redeploy once | each user rebuilds their copy |
| Context cost | the answer | ~200 tokens of metadata until a question calls for it; then the body (~2.3k) and the script output |

## Install

```bash
python scripts/build_skill.py --out ~/.claude/skills/reel-quest     # Claude Code
python scripts/build_skill.py --out ~/.agents/skills/reel-quest     # Codex
python scripts/build_skill.py --out ~/.gemini/skills/reel-quest     # Gemini CLI
```

Each is for you, in every project; the agent loads it in its next session.
Claude Code reads only `.claude/skills/` and Codex only `.agents/skills/`
(both checked 2026-10-23, Codex 0.146). The Gemini CLI path is from its
documentation and untested here: the installed 0.24.5 predates skills. It needs a Python 3.10+ on PATH with
pandas, numpy, scipy and plotly (`skills/reel-quest/requirements.txt`).
Rebuild it after a database update. To remove it, delete the directory.

## What's in it

The bundle that `scripts/build_skill.py` assembles:

- `SKILL.md` covers:
  - when to use the skill;
  - `query.py` (`skills/reel-quest/scripts/query.py`), the MCP server's ten tools as one command that
    prints JSON;
  - five recipes on `rel_quest`'s helpers: look up a ligand, paired
    separation factors under chosen conditions, citing the papers behind a
    result, one element across ligands, recomputing a trend on a subset;
  - the rules: bounds, mixtures, pairing within a condition, pooled
    medians, missing parameters, `tested_elements`, Subtype, DGA naming;
  - citing.
- `rel_quest/`, `vocab/` and `data/sx_ligands_merged_v14_flat.jsonl`, in the
  checkout's layout, so the helpers run wherever the skill is installed;
- `reference/data_dictionary.csv` (from `schema/`), read only when needed;
- the licence files and `CITATION.cff`.

The source is `skills/reel-quest/`. Everything else is copied in at build
time, so don't edit an installed copy.

## A finding the skill documents

`plotting.extract_element_data`, the dashboard's per-element frame, puts 0.0
for a missing acid concentration and 25 °C for a missing temperature. A
filter on those columns silently includes conditions that never reported
them. `extract_separation_factor_data` leaves them NaN. The skill tells
agents to filter on the latter or on the condition dict. The dashboard code
is unchanged; whether to change it is a separate decision.

## Testing

`tests/test_skill.py` builds the bundle into a temporary directory. It runs
everything from another working directory with `PYTHONPATH` cleared, so it
tests what an installed copy sees:

- **K1:** the bundle holds exactly what the skill reads.
- **K2:** the front matter is what Claude Code loads: name, description
  length and the words that trigger the skill.
- **K3, K4:** every python block in SKILL.md runs, and its numbers match the
  raw JSONL read without `rel_quest`. The paired-SF recipe has to report 3
  conditions for TODGA in HNO3 at 3 M or more, where pooling bounds would
  give 16.
- **K5:** `query.py` answers as the `agent_tools` function of the same name
  does, and it exits non-zero on errors.

The live check runs real Claude Code sessions in an empty project with the
skill installed, so it costs three `claude -p` runs:

```bash
RELQUEST_SKILL_LIVE=1 python -m pytest tests/test_skill.py -k InClaudeCode
```

It asks one question the fixed tools answer, one that needs custom code, and
one that does not name the skill. It checks that the skill loaded and the
numbers are right. The transcripts go to `outputs/skill_*_transcript.json`.

**First run (2026-10-13): all three passed** (Opus, $0.32, $1.05 and $0.74;
5 to 17 turns):
- the skill loaded each time, including for the question that does not name
  it;
- DBPhthA 4.30, with its single condition noted;
- 1.675 over 3 paired conditions, with the 9 bound-only pairs explained;
- all 4 yttrium ligands, with CyMe4-BTPhen flagged as limits only.

The cost is the agent's own context, re-read on every turn (126k to 509k
cached input tokens), not the data. The skill saves data tokens, not agent
turns.

Claude Code's shell check stopped some `python - <<'EOF'` scripts with
"Contains brace with quote character (expansion obfuscation)". Most recipes
have such braces (dict literals, f-strings). In a normal session that means
an approval prompt per script, which SKILL.md's `.py`-file route avoids. The
live check allows python commands but not Write, so that route was not
exercised there. The check used a project skill (`.claude/skills/` of the
temporary project); a user-level install uses the same discovery.

## The contribution skill

A second skill, `.agents/skills/reel-quest-contribute/`, ships inside the
repository rather than being installed: it drafts a contribution file from a
paper, or does the second read of one (`contributions/README.md`). Codex and
Gemini CLI find it in `.agents/skills/`; Claude Code through
`.claude/skills/reel-quest-contribute/SKILL.md`, a pointer that repeats its
front matter (a test keeps the two equal) and sends the agent to the shared
instructions. A symbolic link would also work, but not on a Windows checkout
without symlink support. `tests/test_contribute_skill.py` covers it; its live
part runs real agents on two papers whose verified data are stripped from a
public-export clone: `RELQUEST_CONTRIBUTE_LIVE=claude,codex`.

**First live run (2026-10-23)**, Claude Code (Opus) and Codex 0.146 (default
model), each in its own clone. Draft: Ustynyuk 2015 (Figure 2, one new
ligand) and Alyapyshev 2011 (Table 5, one new ligand, plus Eu-only figures).
Review: the verified files with defects planted. Both agents:
- found and used the skill unprompted;
- matched every verified value: within 0.15 dex on the figure, exactly on
  the table, with the right InChIKey for each new ligand;
- left out Table 5's first column, which the text calls earlier work
  (ref. 17); the curator had kept it with a provenance flag;
- drafted the Eu-only series of the Am/Eu figures too (Claude 67
  conditions, Codex 64), which the verified file had not covered and which
  no truth exists for here;
- called the planted defects Wrong value, Wrong condition and Wrong
  molecule, and Table 5's first column Not measured here.

Claude cost $18.61 over the four sessions (the Alyapyshev draft alone $8.92,
111 turns); Codex used 0.7-3.8 M input tokens per session, mostly cached.
Two test-side errors were fixed from this run and the saved files rescored:
conditions are matched on the locator as well, and Not measured here is an
accepted verdict for that column. One skill fix: Codex put the phase ratio
in `Solvent_Ratio`, which holds a mixed diluent's ratio.
