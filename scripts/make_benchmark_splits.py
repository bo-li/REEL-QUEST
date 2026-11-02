"""Leakage-controlled cross-validation splits for machine-learning users.

Writes data/splits/benchmark_splits.json with two split definitions keyed by
ligand Name and InChIKey:

  grouped_5fold      5 folds in which no *group* of ligands straddles train
                     and test. A group is a ligand plus every ligand that
                     shares an identical condition with it (same conditions,
                     same log D for every element). Such groups exist because
                     the Liu 2022 re-ingest attributed each structural
                     family's measurements to all of its variants (see
                     docs/BUILD_LINEAGE.md); splitting by ligand alone would
                     put the same measurements on both sides.
  leave_subtype_out  one fold per chemical subtype (21), for extrapolation to
                     an unseen structural family.

Folds are balanced by number of D-values with a seeded greedy assignment.
Harish (ACS Omega 2026, 10.1021/acsomega.6c07482) showed that on the Liu
1,202-point set a random row split gives R^2 = 0.94 while grouped splits give
0.20 to -0.03; publishing the split makes that comparison reproducible.

Run:
  python scripts/make_benchmark_splits.py [--folds 5] [--seed 20260923]
"""

from __future__ import annotations

import argparse
import collections
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

DB_PATH = ROOT / "data" / "sx_ligands_merged_v14_flat.jsonl"
OUT = ROOT / "data" / "splits" / "benchmark_splits.json"


def condition_key(condition):
    return (condition.get("Acid_Name"), condition.get("Acid_Concentration_M"),
            condition.get("Extractant_Concentration_M"), condition.get("Solvent_Name"),
            condition.get("Temperature_C"),
            tuple(sorted((el, v.get("log_value"))
                         for el, v in condition["distribution_ratios"].items())))


def shared_measurement_groups(records):
    """Connected components of ligands linked by an identical condition."""
    parent = {r["ligand"]["Name"]: r["ligand"]["Name"] for r in records}

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    owners = collections.defaultdict(set)
    for r in records:
        for c in r["experimental_data"]:
            owners[condition_key(c)].add(r["ligand"]["Name"])
    for names in owners.values():
        names = sorted(names)
        for other in names[1:]:
            parent[find(other)] = find(names[0])
    groups = collections.defaultdict(list)
    for name in parent:
        groups[find(name)].append(name)
    return sorted((sorted(g) for g in groups.values()), key=lambda g: g[0])


def build(records, folds=5, seed=20260923):
    size = {r["ligand"]["Name"]: sum(len(c["distribution_ratios"]) for c in r["experimental_data"])
            for r in records}
    key = {r["ligand"]["Name"]: r["ligand"].get("InChI_Key") for r in records}
    groups = shared_measurement_groups(records)
    rng = np.random.default_rng(seed)
    order = rng.permutation(len(groups))
    ordered = sorted((groups[i] for i in order), key=lambda g: -sum(size[n] for n in g))
    load = [0] * folds
    assignment = {}
    for group in ordered:
        fold = int(np.argmin(load))
        load[fold] += sum(size[n] for n in group)
        for name in group:
            assignment[name] = fold

    subtype = {r["ligand"]["Name"]: r["ligand"]["Subtype"] for r in records}
    return {
        "description": ("Cross-validation splits for REEL-QUEST. Keys are ligand Names "
                        "(unique) with InChIKeys; see scripts/make_benchmark_splits.py."),
        "database": DB_PATH.name,
        "seed": seed,
        "grouped_5fold": {
            "n_folds": folds,
            "n_groups": len(groups),
            "n_multi_ligand_groups": sum(1 for g in groups if len(g) > 1),
            "d_values_per_fold": load,
            "folds": [[{"name": n, "inchi_key": key[n]} for n in sorted(assignment)
                       if assignment[n] == f] for f in range(folds)],
        },
        "leave_subtype_out": {
            s: [{"name": n, "inchi_key": key[n]} for n in sorted(subtype) if subtype[n] == s]
            for s in sorted(set(subtype.values()))
        },
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--folds", type=int, default=5)
    parser.add_argument("--seed", type=int, default=20260923)
    args = parser.parse_args()
    records = [json.loads(line) for line in DB_PATH.open(encoding="utf-8") if line.strip()]
    splits = build(records, args.folds, args.seed)
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(splits, indent=1) + "\n", encoding="utf-8")
    g = splits["grouped_5fold"]
    print(f"{g['n_groups']} groups ({g['n_multi_ligand_groups']} with shared measurements); "
          f"D-values per fold {g['d_values_per_fold']}")
    print(f"wrote {OUT.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
