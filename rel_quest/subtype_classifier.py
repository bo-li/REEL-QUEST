"""Structure-based Subtype classification for the v13 taxonomy.

The 21-subtype taxonomy was assigned by human review, and only the four
P-subtypes carry SMARTS in vocab/chemical_families.json. Ingesting new ligands
therefore needs a classifier that works from structure alone.

This derives the Subtype from SMARTS, following the two rules the taxonomy memo
documents (docs/LIGAND_FAMILY_CLASSIFICATION_MEMO.md):

  1. For P-containing ligands, P-bond topology decides — not the presence of an
     amide arm or a heterocyclic scaffold. This is what moved CMPO from v12's
     Amide_Phosphoryl to Phosphine_oxide.
  2. For everything else, a heterocyclic scaffold wins: a phenanthroline with
     amide arms is DAPhen/MAPhen (by arm count), not Amide.

`classify_subtype` returns (subtype, confidence, reason). Confidence is "high"
when a specific pattern matched, "low" when it fell through to a broad rule, and
None subtype means undecidable — callers must route those to human review
rather than guessing.

Accuracy against the 166 curated records is measured by
tests/test_subtype_classifier.py, which is the only reason to trust it.
"""

from __future__ import annotations

import json
import os
from functools import lru_cache

VOCAB_PATH = os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
    "vocab", "chemical_families.json",
)

# Ordered most-specific first; the first match wins.
# alpha-carbons may be substituted (Me-DODGA etc.), so CX4 not CH2.
DGA_SMARTS = "[NX3][CX3](=O)[CX4][OX2][CX4][CX3](=O)[NX3]"
DOODA_SMARTS = "[NX3][CX3](=O)[CX4][OX2][CX4][CX4][OX2][CX4][CX3](=O)[NX3]"

SUBTYPE_SMARTS = [
    # --- O-donor scaffolds, most specific first ---------------------------
    ("Polyether_Diamide", DOODA_SMARTS, "O4 polyether diamide"),
    ("DGA", DGA_SMARTS, "diglycolyl bisamide"),
    ("Malonamide", "[NX3][CX3](=O)[CX4][CX3](=O)[NX3]", "malonyl bisamide"),
    ("Butanediamide", "[NX3][CX3](=O)[CX4][CX4][CX3](=O)[NX3]", "succinyl bisamide"),
]

# Phenanthroline / pyridine / bipyridine scaffolds, resolved by arm count below.
PHENANTHROLINE = "c1cc2ccc3cccnc3c2nc1"
PYRIDINE_DIAMIDE = "[NX3][CX3](=O)c1cccc(n1)[CX3](=O)[NX3]"
# Two pyridine rings joined, each bearing an amide — written loosely because
# substitution patterns vary across the BPyDA series.
BIPYRIDYL = "c1ccncc1-c1ccncc1|c1cccc(n1)-c1ncccc1"
BIPYRIDINE_RING = "n1ccccc1-c1ccccn1"
AROMATIC_AMIDE = "c-[CX3](=O)[NX3]"
# Bis-lactam cores (BLPhen, BLBPy): the amide C=O and N sit in a ring fused to the
# heterocycle, next to its ring N. Bond orders are wildcards because an unsaturated
# lactam ring is perceived as aromatic (a pyridinone) and breaks the aromatic
# PHENANTHROLINE pattern.
RING_LACTAM = "[#6;R](=[#8])(~[#7;R])~[#6;R]~[#7;R]"
PHENANTHROLINE_ANY = "[#6]1~[#6]~[#6]2~[#6]~[#6]~[#6]3~[#6]~[#6]~[#6]~[#7]~[#6]3~[#6]2~[#7]~[#6]1"
BIPYRIDINE_ANY = "[#7]1~[#6]~[#6]~[#6]~[#6]~[#6]1-[#6]1~[#6]~[#6]~[#6]~[#6]~[#7]1"
TRIAZINE_OR_PYRAZOLE = "c1nnc(nn1)|c1cc[nH]n1"


@lru_cache(maxsize=1)
def _vocab():
    with open(VOCAB_PATH) as handle:
        return json.load(handle)


def _mol(smiles):
    from rdkit import Chem, RDLogger

    RDLogger.DisableLog("rdApp.*")
    return Chem.MolFromSmiles(smiles) if smiles else None


def _matches(mol, smarts):
    from rdkit import Chem

    patt = Chem.MolFromSmarts(smarts)
    return bool(patt) and mol.HasSubstructMatch(patt)


def _count(mol, smarts):
    from rdkit import Chem

    patt = Chem.MolFromSmarts(smarts)
    return len(mol.GetSubstructMatches(patt)) if patt else 0


def classify_phosphorus(mol):
    """Return the P-subtype by bond topology, or None if there is no P."""
    vocab = _vocab()
    patterns = vocab["p_functional_smarts"]
    for name in vocab["p_match_order"]:
        smarts = patterns.get(name)
        if smarts and _matches(mol, smarts):
            return name
    return None


def classify_subtype(smiles):
    """Classify a SMILES into a v13 Subtype.

    Returns:
        (subtype, confidence, reason). subtype is None when undecidable.
    """
    mol = _mol(smiles)
    if mol is None:
        return None, "none", "unparseable SMILES"

    # Rule 1: P-bond topology decides for anything containing phosphorus.
    if any(a.GetSymbol() == "P" for a in mol.GetAtoms()):
        p_subtype = classify_phosphorus(mol)
        if p_subtype:
            return p_subtype, "high", f"P topology -> {p_subtype}"
        return None, "none", "contains P but matched no P pattern"

    # Rule 2: heterocyclic scaffolds win over the plain amide classes.
    if _count(mol, RING_LACTAM) >= 2:
        if _matches(mol, PHENANTHROLINE_ANY):
            return "BLPhen", "high", "phenanthroline + 2 fused lactams"
        if _matches(mol, BIPYRIDINE_ANY):
            return "Bipydiamide", "high", "bipyridine + 2 fused lactams"
    if _matches(mol, PHENANTHROLINE):
        arms = _count(mol, AROMATIC_AMIDE)
        if arms >= 2:
            return "DAPhen", "high", "phenanthroline + 2 amide arms"
        if arms == 1:
            return "MAPhen", "high", "phenanthroline + 1 amide arm"
        return "Multident_N_donor", "low", "phenanthroline, no amide arm"

    if _matches(mol, BIPYRIDINE_RING) and _count(mol, AROMATIC_AMIDE) >= 2:
        return "Bipydiamide", "high", "bipyridyl dicarboxamide"
    if _matches(mol, PYRIDINE_DIAMIDE):
        return "Pyridinediamide", "high", "pyridine-2,6-dicarboxamide"

    # Tripodal-DGA: three diglycolamide arms on one scaffold.
    if _count(mol, DGA_SMARTS) >= 3:
        return "Tripodal-DGA", "high", "three DGA arms"

    # Crown ether: a macrocycle carrying several ring oxygens.
    ring_info = mol.GetRingInfo()
    for ring in ring_info.AtomRings():
        if len(ring) >= 12:
            oxygens = sum(1 for i in ring if mol.GetAtomWithIdx(i).GetSymbol() == "O")
            if oxygens >= 4:
                return "Crown_ether", "high", f"macrocycle with {oxygens} ring O"

    # Aminopolycarboxylate chelators (EDTA family) and N-centred tripods.
    if _matches(mol, "[NX3]([CX4][CX3](=O)[OX2H1])[CX4][CX3](=O)[OX2H1]"):
        return "EDTA", "high", "aminopolycarboxylic acid"
    if _count(mol, "[NX3]([CX4][CX3](=O)[NX3])([CX4][CX3](=O)[NX3])") >= 1:
        return "Nitrilotriacetamide", "high", "N-centred triamide tripod"

    # Specific O-donor scaffolds, most specific first.
    for subtype, smarts, reason in SUBTYPE_SMARTS:
        if _matches(mol, smarts):
            return subtype, "high", reason

    # All-nitrogen multidentates (BTP/BTBP/BTPhen family).
    if _matches(mol, "c1nnc(nn1)") or _count(mol, "n") >= 4:
        if not _matches(mol, "[CX3](=O)[NX3]"):
            return "Multident_N_donor", "low", "polyaza, no amide"

    # Simple amide, no recognised scaffold.
    if _matches(mol, "[CX3](=O)[NX3]"):
        return "Amide", "low", "monoamide, no recognised scaffold"

    # Phenolic O-donors.
    if _matches(mol, "c[OX2H]"):
        return "Polyphenol", "low", "phenolic OH"

    return None, "none", "no pattern matched"
