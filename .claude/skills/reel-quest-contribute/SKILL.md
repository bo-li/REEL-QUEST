---
name: reel-quest-contribute
description: Draft or review a REEL-QUEST data contribution. Draft turns the solvent-extraction measurements a paper reports (distribution ratios D of the lanthanides, Y and Sc, from its tables, figures or supporting information) into contributions/<doi-slug>.json that passes the repository's contribution check. Review is the second read of such a file against its paper, with page-level evidence, before a maintainer applies it. Use when someone wants to contribute, add, submit or ingest a paper's data into REEL-QUEST, digitise distribution ratios from a figure, or review, verify or check a contribution file or data pull request against the paper.
compatibility: Runs in a clone of the REEL-QUEST repository with Python 3.10+ and its requirements (rdkit, pypdf, pypdfium2, Pillow).
---

# Contributing a paper's data to REEL-QUEST

This skill is shared with other agents (Codex, Gemini CLI), which read
`.agents/skills/`. Its instructions live there, not here:

1. Read `.agents/skills/reel-quest-contribute/SKILL.md` and follow it.
2. Paths in it such as `references/draft.md` and `scripts/figure.py` are
   relative to `.agents/skills/reel-quest-contribute/`.
