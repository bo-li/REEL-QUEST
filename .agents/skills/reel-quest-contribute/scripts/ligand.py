#!/usr/bin/env python3
"""Is a paper's extractant already in REEL-QUEST, and what is its formula?

    python ligand.py lookup --smiles '<SMILES>'   records with this structure (InChIKey), its
                                                  stereoisomers (same connectivity), and names
    python ligand.py lookup --name '<name>'       records whose Name or Alternative_Names match
    python ligand.py formula '<SMILES>'           formula, average and monoisotopic mass, InChIKey
    python ligand.py subtypes                     the Subtype vocabulary, with the Mechanism and
                                                  Extractant_Type its records carry

Prints JSON. A molecule found under `same_structure` is not new: file its
conditions under that Name. One under `stereoisomers` is a different molecule
only if the paper states the other configuration.
"""

import argparse
import collections
import json
import sys
from pathlib import Path


def repository_root() -> Path:
    for parent in Path(__file__).resolve().parents:
        if (parent / "rel_quest" / "__init__.py").is_file():
            return parent
    sys.exit("ligand.py must run inside a REEL-QUEST checkout")


ROOT = repository_root()
sys.path.insert(0, str(ROOT))

from rel_quest.contributions import inchikey, load_records  # noqa: E402

DB = ROOT / "data" / "sx_ligands_merged_v14_flat.jsonl"


def summary(record) -> dict:
    lig = record["ligand"]
    return {"Name": lig["Name"], "SMILES": lig["SMILES"], "Subtype": lig.get("Subtype"),
            "Alternative_Names": lig.get("Alternative_Names") or [],
            "conditions": len(record.get("experimental_data", []))}


def cmd_lookup(args) -> int:
    records = load_records(args.database)
    out = {}
    if args.smiles:
        key = inchikey(args.smiles)
        if not key:
            sys.exit(f"RDKit cannot parse {args.smiles!r}")
        same, stereo = [], []
        for r in records:
            other = inchikey(r["ligand"]["SMILES"])
            if other == key:
                same.append(summary(r))
            elif other and other.split("-")[0] == key.split("-")[0]:
                stereo.append(summary(r))
        out.update(inchikey=key, same_structure=same, stereoisomers=stereo)
    if args.name:
        wanted = args.name.casefold().replace(" ", "")
        hits = []
        for r in records:
            names = [r["ligand"]["Name"], *(r["ligand"].get("Alternative_Names") or [])]
            if any(wanted == n.casefold().replace(" ", "") or wanted in n.casefold().replace(" ", "")
                   for n in names if n):
                hits.append(summary(r))
        out["by_name"] = hits
    if not out:
        sys.exit("give --smiles or --name")
    print(json.dumps(out, indent=1, ensure_ascii=False))
    return 0


def cmd_formula(args) -> int:
    from rdkit import Chem
    from rdkit.Chem import Descriptors, rdMolDescriptors
    mol = Chem.MolFromSmiles(args.smiles)
    if mol is None:
        sys.exit(f"RDKit cannot parse {args.smiles!r}")
    print(json.dumps({"formula": rdMolDescriptors.CalcMolFormula(mol),
                      "mw": round(Descriptors.MolWt(mol), 2),
                      "monoisotopic": round(Descriptors.ExactMolWt(mol), 4),
                      "inchikey": Chem.MolToInchiKey(mol),
                      "canonical_smiles": Chem.MolToSmiles(mol)}, indent=1))
    return 0


def cmd_subtypes(args) -> int:
    vocab = json.loads((ROOT / "vocab" / "chemical_families.json").read_text(encoding="utf-8"))["subtypes"]
    seen = collections.defaultdict(collections.Counter)
    for r in load_records(args.database):
        lig = r["ligand"]
        seen[lig.get("Subtype")][(lig.get("Mechanism"), lig.get("Extractant_Type"))] += 1
    out = []
    for name, entry in vocab.items():
        common = seen[name].most_common(1)
        out.append({"Subtype": name, "donor_set": entry.get("donor_set"), "description": entry.get("description"),
                    "representative": entry.get("representative"), "records": sum(seen[name].values()),
                    "Mechanism": common[0][0][0] if common else None,
                    "Extractant_Type": common[0][0][1] if common else None})
    print(json.dumps(out, indent=1, ensure_ascii=False))
    return 0


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0],
                                     formatter_class=argparse.RawDescriptionHelpFormatter, epilog=__doc__)
    parser.add_argument("--database", type=Path, default=DB)
    sub = parser.add_subparsers(dest="command", required=True)
    p = sub.add_parser("lookup", help="find a molecule in the database")
    p.add_argument("--smiles")
    p.add_argument("--name")
    p.set_defaults(func=cmd_lookup)
    p = sub.add_parser("formula", help="formula and masses of a SMILES")
    p.add_argument("smiles")
    p.set_defaults(func=cmd_formula)
    p = sub.add_parser("subtypes", help="the Subtype vocabulary")
    p.set_defaults(func=cmd_subtypes)
    args = parser.parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    raise SystemExit(main())
