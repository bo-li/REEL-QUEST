#!/usr/bin/env python3
"""Quick lookups in REEL-QUEST: the MCP server's ten tools, as one command.

    python query.py <tool> ['<JSON arguments>']
    python query.py vocab        # valid trends, subtypes, mechanisms, elements

Prints the answer as JSON on stdout. A ligand that is not found exits 1; an
unknown tool, bad JSON, a wrong argument or a value outside the vocabulary
exits 2. The reason goes to stderr.
"""

import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))  # the bundled rel_quest

from rel_quest import agent_tools  # noqa: E402
from rel_quest.core import LigandDatabase  # noqa: E402
from rel_quest.schema import TREND_CLASSIFICATIONS, TREND_DIRECTIONS  # noqa: E402

ELEMENTS = ("Sc", "Y") + tuple(agent_tools.LANTHANIDES)

TOOLS = {
    "get_database_statistics": agent_tools.get_database_statistics,
    "search_ligands_by_name": agent_tools.search_ligands_by_name,
    "get_ligand_details": agent_tools.get_ligand_details,
    "get_ligand_conditions": agent_tools.get_ligand_conditions,
    "search_by_trend": agent_tools.search_by_trend,
    "search_by_mechanism": agent_tools.search_by_mechanism,
    "search_by_subtype": lambda db, subtype, limit=10:
        agent_tools.search_by_extractant_type(db, subtype, limit),
    "get_best_separators": agent_tools.get_best_separators,
    "compare_ligands": agent_tools.compare_ligands,
    "get_ligands_for_element": agent_tools.get_ligands_for_element,
}

# Argument -> the vocabulary its value must come from.
CHECKED = {"trend": "trends", "mechanism": "mechanisms", "subtype": "subtypes",
           "element": "elements", "element1": "elements", "element2": "elements"}


def vocabulary(db):
    ligands = db.get_all_ligands()
    return {
        "trends": {t: TREND_DIRECTIONS[t] for t in TREND_CLASSIFICATIONS},
        "subtypes": sorted({l.subtype for l in ligands if l.subtype}),
        "mechanisms": sorted({l.mechanism for l in ligands if l.mechanism}),
        "elements": [e for e in ELEMENTS if any(e in l.tested_elements for l in ligands)],
    }


def fail(code: int, message: str):
    print(message, file=sys.stderr)
    sys.exit(code)


def main(argv):
    if not argv or argv[0] not in (*TOOLS, "vocab"):
        fail(2, f"usage: query.py <tool> ['<JSON arguments>']; tools: vocab, {', '.join(TOOLS)}")
    tool = argv[0]
    try:
        args = json.loads(argv[1]) if len(argv) > 1 else {}
    except json.JSONDecodeError as error:
        fail(2, f"arguments are not JSON: {error}")
    if not isinstance(args, dict):
        fail(2, "arguments must be a JSON object, e.g. '{\"name\": \"TODGA\"}'")

    db = LigandDatabase()
    vocab = vocabulary(db)
    if tool == "vocab":
        result = vocab
    else:
        for key, value in args.items():
            valid = vocab.get(CHECKED.get(key))
            if valid is not None and value is not None and value not in valid:
                fail(2, f"{key}={value!r} is not valid; use one of: {', '.join(valid)}")
        try:
            result = TOOLS[tool](db, **args)
        except TypeError as error:
            fail(2, f"{tool}: {error}")
        if isinstance(result, dict) and "error" in result:
            fail(1, result["error"])
    print(json.dumps(result, ensure_ascii=False, default=str))


if __name__ == "__main__":
    main(sys.argv[1:])
