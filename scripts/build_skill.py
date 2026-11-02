#!/usr/bin/env python3
"""Assemble the REEL-QUEST query skill (Agent Skills format: Claude Code, Codex, Gemini CLI).

The skill is skills/reel-quest/ (SKILL.md, scripts/query.py, requirements.txt)
beside copies of rel_quest/, vocab/, the database, the licence files and the
data dictionary, in the checkout's own layout, so the canonical helpers run
wherever the skill is installed. Like the Space bundle, it is generated: change
skills/reel-quest/, not an installed copy.

    python scripts/build_skill.py --out build/reel-quest               # to look at
    python scripts/build_skill.py --out ~/.claude/skills/reel-quest    # for you, in Claude Code
    python scripts/build_skill.py --out ~/.agents/skills/reel-quest    # for you, in Codex
    python scripts/build_skill.py --out ~/.gemini/skills/reel-quest    # for you, in Gemini CLI (from its docs; untested)

The directory's name must stay reel-quest, the skill's name. Tested by
tests/test_skill.py.
"""

import argparse
import shutil
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from build_space import ROOT, copy_database, prepare  # noqa: E402

SKILL = ROOT / "skills" / "reel-quest"
SKILL_FILES = ("SKILL.md", "requirements.txt", "scripts/query.py")
REFERENCE = ROOT / "schema" / "data_dictionary.csv"


def assemble(out: Path) -> Path:
    if out.name != "reel-quest":
        raise SystemExit(f"the skill directory must be named reel-quest, not {out.name}")
    prepare(out, ("SKILL.md", "rel_quest"), "skill bundle")
    for name in SKILL_FILES:
        (out / name).parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(SKILL / name, out / name)
    (out / "reference").mkdir()
    shutil.copyfile(REFERENCE, out / "reference" / REFERENCE.name)
    copy_database(out)
    return out


def main(argv=None) -> None:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--out", type=Path, required=True,
                        help="the skill directory to write, named reel-quest (replaced if it is one)")
    args = parser.parse_args(argv)
    out = assemble(args.out.expanduser().resolve())
    print(f"REEL-QUEST skill at {out}", file=sys.stderr)


if __name__ == "__main__":
    main()
