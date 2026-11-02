#!/usr/bin/env python3
"""Assemble the Hugging Face Space bundle, and optionally upload it.

The Space repository holds space/app.py, its card and requirements beside
copies of rel_quest/, vocab/, the database and the licence files, in this
repository's own layout (vocab/ is read beside the package). Nothing is
duplicated in git: this script is how the Space gets its code, and how it
gets a new database version.

    python scripts/build_space.py --out build/space
    python scripts/build_space.py --out build/space --upload bolivu/reel-quest

--upload runs `hf upload` (pip install huggingface_hub; hf auth login). The
bundle is tested before it ships: tests/test_mcp_server.py launches the app
from a bundle this script built (docs/MCP_SERVER.md, "Testing").
"""

import argparse
import shutil
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
DATABASE = "data/sx_ligands_merged_v14_flat.jsonl"
SPACE_FILES = ("app.py", "README.md", "requirements.txt")  # from space/
FILES = ("LICENSE", "LICENSE-DATA", "CITATION.cff", DATABASE)
PACKAGES = (("rel_quest", "*.py"), ("vocab", "*.json"))


def prepare(out: Path, markers, kind: str) -> None:
    """An empty `out`, replacing only a bundle this script made earlier (it
    holds every marker file): a mistyped --out must not delete an unrelated
    directory."""
    if out.exists():
        if any(out.iterdir()) and not all((out / name).exists() for name in markers):
            raise SystemExit(f"{out} exists and is not a {kind}; choose another --out")
        shutil.rmtree(out)
    out.mkdir(parents=True)


def copy_database(out: Path) -> None:
    """rel_quest/, vocab/, the database and the licence files, in the checkout's layout."""
    for folder, pattern in PACKAGES:
        (out / folder).mkdir()
        for source in sorted((ROOT / folder).glob(pattern)):
            shutil.copyfile(source, out / folder / source.name)
    for name in FILES:
        (out / name).parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(ROOT / name, out / name)


def assemble(out: Path) -> Path:
    prepare(out, SPACE_FILES, "Space bundle")
    for name in SPACE_FILES:
        shutil.copyfile(ROOT / "space" / name, out / name)
    copy_database(out)
    return out


def source_revision() -> str:
    """The commit the bundle was built from, marked when the tree had changes."""
    def git(*args):
        return subprocess.run(["git", *args], cwd=ROOT, capture_output=True, text=True).stdout.strip()

    revision = git("rev-parse", "--short", "HEAD") or "unknown"
    return revision + ("+uncommitted" if git("status", "--porcelain") else "")


def upload(out: Path, repo: str) -> None:
    # `hf upload` creates a missing repo, public unless told otherwise, and
    # without the Space's SDK, hardware or GRADIO_SSR_MODE: the database
    # must not land somewhere like that by a typo.
    if subprocess.run(["hf", "spaces", "info", repo], capture_output=True).returncode:
        raise SystemExit(f"No Space {repo} (or no access to it). Create it first: "
                         "docs/MCP_SERVER.md, Hugging Face Space, Deploy.")
    subprocess.run(["hf", "upload", repo, str(out), ".", "--repo-type", "space",
                    "--commit-message", f"Deploy REEL-QUEST from {source_revision()}"],
                   check=True)


def main(argv=None) -> None:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--out", type=Path, required=True,
                        help="directory to assemble the bundle in (replaced if it is one)")
    parser.add_argument("--upload", metavar="NAMESPACE/SPACE",
                        help="then upload the bundle to this Space with `hf upload`")
    args = parser.parse_args(argv)

    out = assemble(args.out.resolve())
    print(f"Space bundle at {out}", file=sys.stderr)
    if args.upload:
        upload(out, args.upload)


if __name__ == "__main__":
    main()
