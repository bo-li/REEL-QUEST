# REEL-QUEST Ligand Classification — Flat 21-Subtype Taxonomy

**Memo for external presentation**
Date: 2026-07-27 · Database: REEL-QUEST v13 (166 ligands) · Vocabulary: `vocab/chemical_families.json` v2.1

> **Supersedes** the 2026-06-02 version of this memo, which described the two-layer (13 families × 26 subtypes) scheme used in v12. **v13 retired the two-layer taxonomy** in favour of a single flat layer of 21 subtypes. The v12 memo is preserved in `docs/history/` for reference.

---

## 1. Purpose

REEL-QUEST catalogs experimental rare-earth-element (REE) extraction data for 166 ligands compiled from three sources (Liu 2022 JACS Au, the SAFE archive, and the OECD-NEA IDEAL database). Every ligand carries exactly one **Subtype** label drawn from a controlled vocabulary of 21 categories. This is the canonical and only classification.

Each subtype additionally carries a **donor group** attribute (`O_donor` / `N_O_hybrid` / `N_donor` / `P_donor`) on the vocabulary side. Donor groups exist for documentation and figure-coloring only — they are **not stored on records** and are not a second classification layer.

---

## 2. Why a single layer

The v12 scheme nested 26 subtypes under 13 families. In practice the nesting bought little: several families contained exactly one subtype (Crown_ether, Polyphenol, Polyether_Diamide, Multident_N_donor), two subtypes were empty (BisDGA, Phosphonate-under-Organophosphorus), and the family boundary for P-containing ligands was drawn by scaffold (phenanthroline vs not) rather than by P chemistry — which meant chemically identical P=O donors landed in different families depending on what they were attached to.

v13 flattens to the layer researchers actually name, and re-cuts the P-containing categories by **P-bond topology**, which is what determines whether a ligand is neutral or acidic:

| Subtype | Rule |
|---|---|
| `Phosphine_oxide` | three P–C bonds (R₃P=O) |
| `Phosphonate` | ≥1 P–C bond **and** ≥1 P–O ester/acid bond |
| `Phosphate` | three P–OR bonds, no P–C, no P–OH |
| `Phosphoric_acid` | no P–C bond, ≥1 P–OH (acidic) |

This is applied by SMARTS matching in a defined precedence order (`p_match_order` in the vocabulary file), so the assignment is reproducible from structure alone.

The remaining design principles are unchanged from v12:

- **Group by donor chemistry, not historical naming.** The donor-atom set a ligand presents to the metal center is what drives selectivity.
- **Preserve established sub-family identities.** Subtype names match the labels the literature already uses (DGA, BTBP, CMPO-family, DAPhen), so researchers can query for what they know.
- **Encode chemistry in the palette.** Subtype colors are tinted within their donor group (blues = O, violets = N+O hybrid, amber = pure N, reds/oranges = P), so a pie chart still reads as donor chemistry at a glance.

---

## 3. The 21 subtypes

Counts are from REEL-QUEST v13 (166 ligands, sums to 166).

### O-donors (9 subtypes — blue/green palette)

| Subtype | Donor set | Description | Representative | n |
|---|---|---|---|---|
| **DGA** | O3 | Tridentate diglycolyl-bisamide (2 amide-O + 1 ether-O); monomeric TODGA-type and pre-organized BisDGA | TODGA | 75 |
| **Malonamide** | O2 | 1,3-Bisamide with malonyl (C1) bridge | DMDOHEMA | 11 |
| **Crown_ether** | O6 | Macrocyclic polyether; 18-crown-6 and dibenzo/dicyclohexyl/dipyridino derivatives | 18-crown-6 | 4 |
| **Tripodal-DGA** | O3 × 3 | Tripodal scaffold bearing three diglycolamide arms | TriEtTriDoDGA | 3 |
| **Amide** | O1 | Simple N,N-dialkyl alkanamide monoamide | DHOA | 2 |
| **Polyether_Diamide** | O4 | Tetradentate polyether-diamide; extended-glycol analogue of DGA | DOODA | 2 |
| **Phthalamic_acid** | O2 | Amide-carboxylate hybrid, chelating via amide-O + COOH-O | DBPhthA | 1 |
| **Butanediamide** | O2 | 1,4-Bisamide with succinyl (C2) bridge | TBBDA | 1 |
| **Polyphenol** | O2–O4 | Polyhydroxy-aromatic donor (calixarene-OH, catechol) | Quercetin | 1 |

### N,O-hybrid donors (7 subtypes — violet/magenta palette)

| Subtype | Donor set | Description | Representative | n |
|---|---|---|---|---|
| **DAPhen** | N2O2 | 1,10-Phenanthroline with diamide arms at 2,9-positions | Et-Tol-DAPhen | 11 |
| **Bipydiamide** | N2O2 | 6,6′-Bipyridyl-2,2′-dicarboxamide scaffold | EHPhBPDA | 6 |
| **Pyridinediamide** | N1O2 | Pyridine-2,6-dicarboxamide scaffold | PDA | 5 |
| **MAPhen** | N2O1 | 1,10-Phenanthroline with a single amide arm at the 2-position | DOPTA | 4 |
| **EDTA** | N2O4–N2O6 | Aminopolycarboxylic acid chelators (renamed from v12's "DTA") | HEDTA | 2 |
| **Nitrilotriacetamide** | N1O3 | N-centered tripod with three amide arms | HONTA | 2 |
| **BLPhen** | N4O2–N4O4 | Bi-loop phenanthroline diamide bridging two phen-amide cores | DHTPhBLPhen | 2 |

### N-donors (1 subtype — amber palette)

| Subtype | Donor set | Description | Representative | n |
|---|---|---|---|---|
| **Multident_N_donor** | N3–N4 | All-nitrogen multidentate donors with no amide or P=O arm: BTP, BTBP, BTPhen, BPyrPhen (consolidated from three v12 subtypes) | C5-BTBP | 8 |

### P-donors (4 subtypes — red/orange palette)

| Subtype | Donor set | Description | Representative | n |
|---|---|---|---|---|
| **Phosphine_oxide** | O1 | Three P–C bonds (R₃P=O); includes pure phosphine oxides, the CMPO-family amidophosphine oxides, and phen-PO | CMPO | 19 |
| **Phosphonate** | O1–O2 | ≥1 P–C bond and ≥1 P–O ester/acid bond; includes POPPhen and HEH[EHP] (PC88A) | DHDECMP | 5 |
| **Phosphate** | O1 | Three P–OR bonds, no P–C, no P–OH | TBP | 1 |
| **Phosphoric_acid** | O1–O2 acidic | No P–C bond, ≥1 P–OH | DHPA | 1 |

---

## 4. Worked examples

| Ligand | Subtype | Reasoning |
|---|---|---|
| **TODGA** | DGA | Monomeric N,N,N′,N′-tetraoctyl-diglycolamide; O3 tridentate |
| **DMDOHEMA** | Malonamide | 1,3-Bisamide on a malonyl backbone; O2 chelating |
| **CMPO** | Phosphine_oxide | Three P–C bonds. The amide arm does *not* move it out of the P category in v13 — this is the main change from v12, where it was `Amide_Phosphoryl` |
| **HDEHP** | Phosphoric_acid | No P–C bond, acidic P–OH |
| **TBP** | Phosphate | Trialkyl phosphate; (RO)₃P=O |
| **HEH[EHP]** (PC88A) | Phosphonate | One P–C bond plus a P–OH |
| **C5-BTBP** | Multident_N_donor | All-N N4 scaffold; no amide, no P=O |
| **Et-Tol-DAPhen** | DAPhen | Phenanthroline with two amide arms |

**Key v13 convention:** for P-containing ligands, **P-bond topology decides the subtype** — not the presence of an amide arm or a heterocyclic scaffold. A phosphorylated phenanthroline is `Phosphine_oxide` or `Phosphonate` according to its P bonds. For non-P ligands, a heterocyclic scaffold still wins: a phenanthroline with amide arms is `DAPhen`/`MAPhen`, not `Amide`.

---

## 5. Palette (for slides and figures)

Colors are tinted within donor group, so the donor chemistry is legible even without the legend.

```
O DONORS (blue / green)
  DGA                   #2E86AB      Tripodal-DGA        #4A9BC0
  Polyether_Diamide     #1F6F8B      Amide               #8DB4D2
  Phthalamic_acid       #A8C5DC      Malonamide          #5BB3D6
  Butanediamide         #7AC4E0      Crown_ether         #6A994E
  Polyphenol            #A7C957

N,O HYBRID DONORS (violet / magenta)
  Nitrilotriacetamide   #A23B72      EDTA                #B85A8E
  DAPhen                #5C4B8A      MAPhen              #7867A4
  BLPhen                #9785BD      Pyridinediamide     #9E8FC0
  Bipydiamide           #B5A8D0

N DONORS (amber)
  Multident_N_donor     #E9C46A

P DONORS (red / orange)
  Phosphine_oxide       #C73E1D      Phosphonate         #E76F51
  Phosphate             #F18F01      Phosphoric_acid     #9B2D11
```

The palette lives in `rel_quest/plotting.py` as `CHEMICAL_FAMILY_COLORS` and is shared by the dashboard and all figure-export scripts.

---

## 6. How the taxonomy is used

- **Dashboard overview.** The distribution pie shows subtype counts, colored by donor group.
- **Filtering.** The sidebar multiselect filters every visualization to one or more subtypes.
- **Box plots.** logD and logSF distributions can be grouped by subtype for cross-family comparison.
- **AI assistant.** Natural-language query results carry the subtype so the assistant speaks the same taxonomy as the rest of the dashboard.
- **Programmatic access.** `LigandDatabase.filter_by_subtype(name)` and `filter_by_chemical_family(name)` both accept a subtype label and return the same result set.

---

## 7. Governance

- **Canonical reference:** `vocab/chemical_families.json` (v2.1, versioned). It also carries `v12_to_v13_consolidations`, a machine-readable map of how each v12 family/subtype pair folded into a v13 subtype.
- **Record fields:** each ligand carries `Subtype` **and** `Chemical_Family`, both holding the *same* 21-label value. The duplication is deliberate — it keeps pre-v13 consumer code working. Do not read them as two layers.
- **Migration audit:** `metadata.prior_classification` on every record preserves its v12 family/subtype pair, and `metadata.taxonomy_version` reads `v13_flat`.
- **Coverage:** 166 / 166 ligands classified. There is no `Unknown` subtype in v13.
- **Classification accuracy:** mechanism inference (Solvation vs Cation_Exchange) and acidity classification reach 97.2% agreement against SMARTS pattern matching.

---

*Prepared for external review of the REEL-QUEST ligand classification scheme. Questions or proposed revisions: contact the database maintainer.*
