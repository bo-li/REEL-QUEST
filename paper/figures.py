"""Manuscript figures, generated from the database (no manual editing).

Main text:
  fig1_pipeline        what the database holds: inputs -> REEL-QUEST -> capabilities
  fig2_platform        the platform in use: dashboard, AI assistant, an MCP tool call
                       (composed from paper/capture_platform.py's captures)
  fig3_overview        ligands by subtype (with each subtype's core), trend classes, coverage
  fig4_trend_gallery   one representative ligand per trend class, spread across conditions
  fig5_eu_la_showcase  (a) log D(Eu) and (b) SF(Eu/La) per condition, by subtype
  fig6_uncertainty     one ligand whose trend label and SF sign its uncertainty can flip
  fig7_trend_by_subtype  share of each subtype's ligands in each trend class
  fig8_data_gaps       per subtype, where more experiments would change the answer
Supporting information:
  figS1_ci_width       slope CI width vs conditions behind the fit
  figS2_cross_source   log D where Liu 2022 and IDEaL transcribed the same condition
  figS3_sf_viewer      the dashboard's separation-factor viewer set to one pair (captured)

Writes PDF (vector, for submission) and PNG (300 dpi) to paper/figures/, plus
CAPTIONS.md, whose numbers come from paper/stats.py rather than being typed.

Every choice is data-driven (which ligand illustrates which class, which
subtypes are shown, the order of rows), so the set regenerates unchanged in
code when the database changes; `--db` draws it from another snapshot, e.g.
a candidate, without touching paper/figures/ unless `--out` says so.

Style: black frame on every plot and box, no grid. Trend classes use slots
1-5 of a categorical palette validated for colour-vision deficiency
(adjacent-pair CVD dE >= 9.1, normal-vision dE >= 19.6 on a light surface);
Mixed and Insufficient_Data are neutral greys, because they are the absence of
a shape. `Reverse` is shown as "Descending" (display only: the database, API
and dashboard still store `Reverse`). Text is never drawn in a series colour.

Requires matplotlib (the `paper` extra); RDKit for the structures in Figure 1
(without it the structures are left out and the rest is drawn).

Run:
  python paper/figures.py [--only NAME] [--no-export] [--db PATH --out DIR]
"""

from __future__ import annotations

import argparse
import collections
import csv
import io
import json
import string
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

import matplotlib  # noqa: E402

matplotlib.use("Agg")
import matplotlib.image as mpimg  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.patches import FancyArrow, FancyBboxPatch, Rectangle, Ellipse  # noqa: E402
from matplotlib.ticker import FormatStrFormatter  # noqa: E402
from scipy.stats import spearmanr  # noqa: E402

from paper.stats import (ADJACENT, adjacent_sf, complete, compute_stats, pooled_conditions,  # noqa: E402
                         measured_log, showcase_groups, log_d_by_subtype,
                         has_source_labels, load_stored,
                         duplicated_conditions, label_stability,
                         perturbed_draws, sharp_extremum)
from rel_quest.trends import classify_selectivity_trend  # noqa: E402
from rel_quest.schema import TREND_DISPLAY_ORDER  # noqa: E402
from rel_quest.uncertainty import (condition_trend_points, is_censored,  # noqa: E402
                                   per_condition_slopes)
from rel_quest.validation import cross_source_pairs  # noqa: E402

DB_PATH = ROOT / "data" / "sx_ligands_merged_v14_flat.jsonl"
OUT = ROOT / "paper" / "figures"
FAMILIES = ROOT / "vocab" / "chemical_families.json"

INK = "#000000"
INK_2 = "#3d3c39"
MUTED = "#898781"
LIGHT = "#c3c2b7"
SURFACE = "#ffffff"
SERIES = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4"]
PANEL_FILL = "#eef4fb"
DB_FILL = "#dce9f7"

TREND_FIGURE_COLORS = {
    "Ascending": SERIES[0],
    "Peaked": SERIES[1],
    "Flat": SERIES[2],
    "U_Shaped": SERIES[3],
    "Reverse": SERIES[4],
    "Mixed": MUTED,
    "Insufficient_Data": "#d8d6cf",
}
TREND_LABELS = {"U_Shaped": "U-shaped", "Insufficient_Data": "Insufficient data",
                "Reverse": "Descending"}
NOT_A_SHAPE = {"Mixed", "Insufficient_Data", "Unknown"}
CURVED = {"Peaked", "U_Shaped"}
# linear classes, then curved, then Mixed; derived so the vocabulary stays in schema
GALLERY_ORDER = tuple(sorted((t for t in TREND_DISPLAY_ORDER
                              if t not in NOT_A_SHAPE or t == "Mixed"),
                             key=lambda t: (t == "Mixed", t in CURVED,
                                            TREND_DISPLAY_ORDER.index(t))))
MIN_GALLERY_SPAN_DEX = 1.5  # a narrower axis makes noise look like a trend
LANTHANIDES = ["La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy",
               "Ho", "Er", "Tm", "Yb", "Lu"]
COVERAGE_ELEMENTS = LANTHANIDES + ["Y", "Sc"]
Z = {el: 57 + i for i, el in enumerate(LANTHANIDES)}
SYMBOL = {v: k for k, v in Z.items()}
DONOR_GROUP_LABELS = {"O_donor": "O donors", "N_O_hybrid": "N,O hybrids",
                      "N_donor": "N donors", "P_donor": "P donors"}
SUBTYPE_LABELS = {"Multident_N_donor": "Multidentate N-donor",
                  "Phosphine_oxide": "Phosphine oxide", "Phosphoric_acid": "Phosphoric acid",
                  "Phthalamic_acid": "Phthalamic acid", "Crown_ether": "Crown ether",
                  "Polyether_Diamide": "Polyether diamide"}
SUBTYPE_CORES = {  # generic core per subtype; * = R (alkyl or aryl), [1*] = a whole DGA arm
    "DGA": "O=C(N(*)*)COCC(=O)N(*)*",
    "Phosphine_oxide": "O=C(N(*)*)CP(=O)(*)c1ccccc1",  # CMPO, 12 of the 18
    "DAPhen": "O=C(N(*)*)c1ccc2ccc3ccc(C(=O)N(*)*)nc3c2n1",
    "Malonamide": "O=C(N(*)*)C(*)C(=O)N(*)*",
    "Multident_N_donor": "*c1nnc(-c2cccc(-c3nnc(*)c(*)n3)n2)nc1*",  # BTP
    "Bipydiamide": "O=C(N(*)*)c1cccc(-c2cccc(C(=O)N(*)*)n2)n1",
    "Pyridinediamide": "O=C(N(*)*)c1cccc(C(=O)N(*)*)n1",
    "MAPhen": "O=C(N(*)*)c1ccc2ccc3cccnc3c2n1",
    "Phosphonate": "*P(=O)(O*)O*",
    "Crown_ether": "C1COCCOCCOCCOCCOCCO1",
    "BLPhen": "O=C1N(*)C(*)=C(*)c2cc3ccc4cc5C(*)=C(*)N(*)C(=O)c5nc4c3nc21",
    "Tripodal-DGA": "[1*]NCCc1cc(CCN[1*])cc(CCN[1*])c1",
    "Polyether_Diamide": "O=C(N(*)*)COCCOCC(=O)N(*)*",
    "Nitrilotriacetamide": "O=C(N(*)*)CN(CC(=O)N(*)*)CC(=O)N(*)*",
    "Amide": "*C(=O)N(*)*",
    "Polyphenol": "OC(C=C1)=C(O)C=C1C2=C(O)C(C3=C(O)C=C(O)C=C3O2)=O",  # quercetin
    "Phthalamic_acid": "O=C(N(*)*)c1ccccc1C(=O)O",
    "Phosphoric_acid": "O=P(O)(O*)O*",
    "Phosphate": "O=P(O*)(O*)O*",
    "Butanediamide": "O=C(N(*)*)CCC(=O)N(*)*",
    "EDTA": "O=C(O)CN(CC(=O)O)CCN(CC(=O)O)CC(=O)O",
}
GENERIC_CORES = tuple((name, SUBTYPE_CORES[subtype]) for name, subtype in (  # Figure 1
    ("DGA", "DGA"), ("Malonamide", "Malonamide"), ("CMPO", "Phosphine_oxide"),
    ("BTP", "Multident_N_donor"), ("DAPhen", "DAPhen")))


def style():
    plt.rcParams.update({
        "font.family": "sans-serif",
        "font.size": 8,
        "axes.edgecolor": INK,
        "axes.labelcolor": INK,
        "axes.linewidth": 0.8,
        "axes.spines.top": True,
        "axes.spines.right": True,
        "axes.grid": False,
        "xtick.color": INK,
        "ytick.color": INK,
        "xtick.direction": "out",
        "ytick.direction": "out",
        "legend.frameon": False,
        "figure.facecolor": SURFACE,
        "axes.facecolor": SURFACE,
        "savefig.facecolor": SURFACE,
        "pdf.fonttype": 42,
    })


def load_records(path=None):
    path = Path(path) if path else DB_PATH
    return [json.loads(line) for line in path.open(encoding="utf-8") if line.strip()]


def label(trend):
    return TREND_LABELS.get(trend, trend)


def subtype_label(subtype):
    return SUBTYPE_LABELS.get(subtype, subtype.replace("_", " "))


def panel_label(ax, letter, x=-0.02, y=1.0):
    ax.text(x, y, letter, transform=ax.transAxes, fontsize=10, fontweight="bold",
            color=INK, ha="right", va="bottom")


def donor_groups():
    return json.loads(FAMILIES.read_text(encoding="utf-8"))["donor_groups"]


def shared_ligands(records):
    return {name for group in duplicated_conditions(records) for name, _ in group}


def n_fitted_conditions(record):
    return len(per_condition_slopes(record["experimental_data"]))


class Context:
    """Everything a figure needs, computed once per build."""

    def __init__(self, records):
        self.records = records
        self.stability = label_stability(records)  # one Monte Carlo run serves all
        self.stats = complete(compute_stats(records, stability=self.stability), load_stored())


# --- Figure 1 ------------------------------------------------------------------

def core_image(smiles):
    """RGBA drawing of a generic core at a fixed bond length; None without RDKit."""
    try:
        from rdkit import Chem
        from rdkit.Chem.Draw import rdMolDraw2D
    except ImportError:  # pragma: no cover - RDKit is optional here
        return None
    mol = Chem.MolFromSmiles(smiles)
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 0:
            atom.SetProp("atomLabel", "DGA" if atom.GetIsotope() == 1 else "R")
    drawer = rdMolDraw2D.MolDraw2DCairo(-1, -1)
    options = drawer.drawOptions()
    options.clearBackground = False
    options.fixedBondLength = 34
    options.bondLineWidth = 2
    options.fixedFontSize = 22
    options.padding = 0.08
    rdMolDraw2D.PrepareAndDrawMolecule(drawer, mol)
    drawer.FinishDrawing()
    return mpimg.imread(io.BytesIO(drawer.GetDrawingText()))


def structure_labels(ax, subtypes, width_pt, height_pt, max_zoom=0.13):
    """Draw each subtype's core just left of its y tick; the tick text moves left of it.

    Structures share one scale (max_zoom points per pixel, i.e. one bond
    length) unless that would overflow the cell, when they shrink to fit.
    Without RDKit the axis keeps its plain labels.
    """
    from matplotlib.offsetbox import AnnotationBbox, OffsetImage
    images = [core_image(SUBTYPE_CORES[s]) if s in SUBTYPE_CORES else None for s in subtypes]
    if any(image is None for image in images):
        return
    for y, image in enumerate(images):
        h, w = image.shape[:2]
        zoom = min(max_zoom, width_pt / w, height_pt / h)
        ax.add_artist(AnnotationBbox(OffsetImage(image, zoom=zoom, interpolation="lanczos"),
                                     (-0.01, y), xycoords=ax.get_yaxis_transform(),
                                     box_alignment=(1, 0.5), frameon=False, pad=0,
                                     annotation_clip=False))
    ax.tick_params(axis="y", length=0, pad=width_pt + 6)


CHAIN_LABELS = (  # n-alkyl and 2-ethylhexyl chains drawn as labels, so cores stay legible
    ("Dodec", "*CCCCCCCCCCCC"), ("Dec", "*CCCCCCCCCC"), ("Oct", "*CCCCCCCC"),
    ("Hex", "*CCCCCC"), ("EH", "*CC(CC)CCCC"))


def layout_clearance(mol):
    """Closest approach of two non-bonded atoms, in median bond lengths."""
    pos = mol.GetConformer().GetPositions()[:, :2]
    bonds = [np.linalg.norm(pos[b.GetBeginAtomIdx()] - pos[b.GetEndAtomIdx()])
             for b in mol.GetBonds()]
    scale = float(np.median(bonds)) if bonds else 1.0
    gaps = [np.linalg.norm(pos[i] - pos[j]) / scale
            for i in range(mol.GetNumAtoms()) for j in range(i + 1, mol.GetNumAtoms())
            if mol.GetBondBetweenAtoms(i, j) is None]
    return min(gaps) if gaps else float("inf")


def ligand_image(smiles):
    """RGBA drawing of a real ligand with its alkyl chains condensed to labels.

    Both RDKit layouts are tried and the one whose atoms stand furthest apart
    is drawn (each crowds a different molecule). None without RDKit.
    """
    try:
        from rdkit import Chem
        from rdkit.Chem import rdAbbreviations, rdDepictor
        from rdkit.Chem.Draw import rdMolDraw2D
    except ImportError:  # pragma: no cover - RDKit is optional here
        return None
    abbreviations = rdAbbreviations.ParseAbbreviations(
        "".join(f"{label} {smarts}\n" for label, smarts in CHAIN_LABELS))
    layouts = []
    for coordgen in (False, True):
        rdDepictor.SetPreferCoordGen(coordgen)
        mol = rdAbbreviations.CondenseMolAbbreviations(Chem.MolFromSmiles(smiles),
                                                       abbreviations, maxCoverage=0.9)
        rdDepictor.Compute2DCoords(mol)
        layouts.append((layout_clearance(mol), not coordgen, mol))
    rdDepictor.SetPreferCoordGen(False)
    mol = max(layouts, key=lambda x: x[:2])[2]
    drawer = rdMolDraw2D.MolDraw2DCairo(-1, -1)
    options = drawer.drawOptions()
    options.clearBackground = False
    options.fixedBondLength = 40
    options.bondLineWidth = 2
    options.fixedFontSize = 20
    options.padding = 0.04
    drawer.DrawMolecule(mol)
    drawer.FinishDrawing()
    return mpimg.imread(io.BytesIO(drawer.GetDrawingText()))


def rounded(ax, x, y, w, h, fill=SURFACE, edge=INK, lw=1.0, radius=2.0, z=1):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle=f"round,pad=0,rounding_size={radius}",
                                facecolor=fill, edgecolor=edge, linewidth=lw, zorder=z))


def block_arrow(ax, x, y, length=4.6):
    ax.add_patch(FancyArrow(x, y, length, 0, width=3.4, head_width=7.0, head_length=2.6,
                            length_includes_head=True, color=SERIES[1], zorder=3,
                            linewidth=0))


def database_icon(ax, cx, cy, w=9.0, h=11.0):
    ew = 2.4
    ax.add_patch(Rectangle((cx - w / 2, cy - h / 2), w, h, facecolor=SURFACE, edgecolor="none",
                           zorder=3))
    for y in (cy - h / 2, cy - h / 6, cy + h / 6):
        ax.add_patch(Ellipse((cx, y), w, ew, facecolor=SURFACE, edgecolor=INK, linewidth=0.9,
                             zorder=4))
    ax.add_patch(Ellipse((cx, cy + h / 2), w, ew, facecolor=DB_FILL, edgecolor=INK,
                         linewidth=0.9, zorder=5))
    ax.plot([cx - w / 2] * 2, [cy - h / 2, cy + h / 2], color=INK, linewidth=0.9, zorder=5)
    ax.plot([cx + w / 2] * 2, [cy - h / 2, cy + h / 2], color=INK, linewidth=0.9, zorder=5)
    ax.add_patch(Rectangle((cx - w / 2 + 0.05, cy - h / 2), w - 0.1, h, facecolor=SURFACE,
                           edgecolor="none", zorder=3.5))


def document_icon(ax, x, y, w=8.5, h=11.0, text="JSONL"):
    fold = 2.4
    xs = [x, x + w - fold, x + w, x + w, x, x]
    ys = [y + h, y + h, y + h - fold, y, y, y + h]
    ax.fill(xs, ys, facecolor=SURFACE, edgecolor=INK, linewidth=0.9, zorder=4)
    ax.plot([x + w - fold, x + w - fold, x + w], [y + h, y + h - fold, y + h - fold],
            color=INK, linewidth=0.9, zorder=5)
    for k in range(3):
        ax.plot([x + 1.4, x + w - 1.6], [y + h - 4.2 - 1.6 * k] * 2, color=LIGHT,
                linewidth=0.8, zorder=5)
    ax.text(x + w / 2, y + 2.0, text, ha="center", va="center", fontsize=5.8,
            fontweight="bold", color=INK, zorder=6)


def glyph_axes(ax, x, y, w, h):
    inset = ax.inset_axes((x, y, w, h), transform=ax.transData)
    inset.set_xticks([])
    inset.set_yticks([])
    for side in ("top", "right"):
        inset.spines[side].set_visible(False)
    inset.patch.set_alpha(0)
    return inset


def fig1_pipeline(ctx):
    stats = ctx.stats
    fig = plt.figure(figsize=(7.1, 3.35))
    ax = fig.add_axes((0, 0, 1, 1))
    ax.set_xlim(0, 142)
    ax.set_ylim(0, 67)
    ax.axis("off")
    panels = ((1, 44, "Data inputs"), (50.5, 45, "REEL-QUEST"),
              (102, 39, "Analytical capabilities"))
    for x, w, title in panels:
        rounded(ax, x, 1, w, 65, fill=PANEL_FILL, radius=2.5)
        ax.text(x + w / 2, 62.8, title, ha="center", va="top", fontsize=9.5,
                fontweight="bold", color=INK, zorder=4)
    block_arrow(ax, 45.2, 36)
    block_arrow(ax, 96.7, 36)

    # Data inputs: generic cores of the largest families, then the counts
    box_x, box_y, box_w, box_h = 3, 22.5, 40, 36
    rounded(ax, box_x, box_y, box_w, box_h, fill=SURFACE, radius=2.0, z=2)
    images = {name: core_image(smiles) for name, smiles in GENERIC_CORES}
    if all(image is not None for image in images.values()):
        # rows of structures at one common scale, as large as the box allows
        names = [name for name, _ in GENERIC_CORES]
        rows = [names[0:2], names[2:4], names[4:]]
        gap, label_h, pad = 2.0, 2.4, 1.2
        width_scale = min((box_w - 2 * pad - gap * (len(r) - 1))
                          / sum(images[n].shape[1] for n in r) for r in rows)
        height_scale = ((box_h - 2 * pad - label_h * len(rows))
                        / sum(max(images[n].shape[0] for n in r) for r in rows))
        scale = min(width_scale, height_scale)
        y = box_y + box_h - pad
        for r in rows:
            row_h = max(images[n].shape[0] for n in r) * scale
            widths = [images[n].shape[1] * scale for n in r]
            x = box_x + (box_w - sum(widths) - gap * (len(r) - 1)) / 2
            for n, w in zip(r, widths):
                h = images[n].shape[0] * scale
                cy = y - row_h / 2
                ax.imshow(images[n], extent=(x, x + w, cy - h / 2, cy + h / 2), zorder=3,
                          interpolation="lanczos")
                ax.text(x + w / 2, y - row_h - 0.2, n, ha="center", va="top", fontsize=5.8,
                        fontweight="bold", color=INK_2, zorder=4)
                x += w + gap
            y -= row_h + label_h
    ax.set_xlim(0, 142)
    ax.set_ylim(0, 67)
    ax.text(23, 19.8, f"{stats['n_ligands']} extractants", ha="center", va="top",
            fontsize=8.5, fontweight="bold", color=INK, zorder=4)
    ax.text(23, 15.4, f"{stats['n_d_values']:,} distribution ratios\n"
            f"{stats['n_conditions']:,} experimental conditions\n"
            f"from {stats['n_distinct_primary_dois']} original papers",
            ha="center", va="top", fontsize=6.8, color=INK, linespacing=1.45, zorder=4)

    # REEL-QUEST: storage, what one record holds, interfaces
    database_icon(ax, 63.5, 50.5)
    document_icon(ax, 77, 45)
    ax.text(73, 41.5, "one record per ligand", ha="center", va="top", fontsize=6.6,
            fontstyle="italic", color=INK, zorder=4)
    record = (("Ligand", "name, SMILES, subtype"),
              ("Conditions", "acid, diluent, T, conc."),
              ("D values", "La–Lu, Y, Sc, with σ"),
              ("Trend", "class, slope, 95% CI"),
              ("Provenance", "DOI, figure or table"))
    rounded(ax, 53, 17.5, 40, 21, fill=SURFACE, radius=1.5, z=2)
    for i, (key, value) in enumerate(record):
        y = 36.6 - 3.9 * i
        ax.text(55, y, key, fontsize=6.3, fontweight="bold", color=INK, va="top", zorder=4)
        ax.text(67.5, y, value, fontsize=6.2, color=INK, va="top", zorder=4)
    # how people and AI agents reach it: interactive routes on top
    access = (("Dashboard", "AI assistant", "MCP server"), ("Agent skill", "REST API", "Python · CLI"))
    for row, names in enumerate(access):
        for i, name in enumerate(names):
            x, y = 53 + i * 13.6, 10 - row * 6.8
            rounded(ax, x, y, 12.8, 5.8, fill=SURFACE, radius=1.2, z=2)
            ax.text(x + 6.4, y + 2.9, name, ha="center", va="center", fontsize=5.9,
                    fontweight="bold" if row == 0 else "normal", color=INK, zorder=4)

    # Capabilities: each with a glyph drawn from the paper's own palette
    rows = ((46.5, "Selectivity-\ntrend\nclassification"),
            (27.5, "Uncertainty\nquantification"),
            (8.5, "Separation\nfactor for\nany Ln pair"))
    g = np.linspace(0, 1, 20)
    for i, (y, text) in enumerate(rows):
        inset = glyph_axes(ax, 104.5, y, 14, 12)
        if i == 0:
            for trend, curve in (("Ascending", 0.15 + 0.7 * g), ("Flat", 0.5 + 0 * g),
                                 ("Reverse", 0.8 - 0.65 * g),
                                 ("Peaked", 0.3 + 0.55 * np.exp(-((g - 0.5) / 0.2) ** 2))):
                inset.plot(g, curve, color=TREND_FIGURE_COLORS[trend], linewidth=1.4)
            inset.set_ylim(0, 1.05)
        elif i == 1:
            xs = np.linspace(0.08, 0.92, 7)
            ys = 0.25 + 0.55 * xs + np.array([0.05, -0.06, 0.04, 0.07, -0.05, 0.03, -0.04])
            band = 0.07 + 0.35 * (xs - 0.5) ** 2
            inset.fill_between(g, 0.25 + 0.55 * g - (0.07 + 0.35 * (g - 0.5) ** 2),
                               0.25 + 0.55 * g + (0.07 + 0.35 * (g - 0.5) ** 2),
                               color=SERIES[0], alpha=0.2, linewidth=0)
            inset.plot(g, 0.25 + 0.55 * g, color=SERIES[0], linewidth=1.3)
            inset.errorbar(xs, ys, yerr=band * 0.8, fmt="o", markersize=2.2, color=INK,
                           elinewidth=0.7, capsize=1.2)
            inset.set_ylim(0, 1.05)
            inset.set_xlim(0, 1)
        else:
            inset.bar([0.3, 0.7], [0.35, 0.85], width=0.28, color=[SERIES[0], SERIES[1]])
            inset.annotate("", xy=(0.7, 0.93), xytext=(0.3, 0.43),
                           arrowprops={"arrowstyle": "-|>", "color": INK, "linewidth": 0.8})
            inset.text(0.3, 0.05, "Ln₁", ha="center", va="bottom", fontsize=5.5, color=SURFACE)
            inset.text(0.7, 0.05, "Ln₂", ha="center", va="bottom", fontsize=5.5, color=SURFACE)
            inset.set_xlim(0, 1)
            inset.set_ylim(0, 1.05)
        ax.text(121, y + 6, text, fontsize=7, fontweight="bold", color=INK, va="center",
                linespacing=1.2, zorder=4)
    return fig


# --- Figure 2 ------------------------------------------------------------------

PLATFORM = OUT / "platform"  # written by paper/capture_platform.py from the running software
MCP_FIELDS = ("name", "ligand_family", "median_log_SF", "n_paired_conditions")


def mcp_lines(exchange, width=62):
    """The recorded MCP exchange as wrapped JSON lines, each result entry abridged to MCP_FIELDS."""
    result = dict(exchange["result"])
    top = result.pop("top_separators")

    def wrapped(items, first, indent, close):
        lines, line = [], first
        for i, item in enumerate(items):
            piece = item + ("," if i < len(items) - 1 else close)
            if len(line) + len(piece) + 1 > width and line.strip() not in ("", "{"):
                lines.append(line.rstrip())
                line = indent
            line += piece + " "
        return lines + [line.rstrip()]

    lines = [f"→ tools/call  {exchange['call']}", "  " + json.dumps(exchange["arguments"]), "← result"]
    def short(value):  # a long string (the server's note) to its first clause
        if isinstance(value, str) and len(value) > width - 14:
            clause = value.split(";")[0]
            if len(clause) > width - 18:
                clause = clause[:width - 18].rsplit(" ", 1)[0]
            return json.dumps(clause + " …", ensure_ascii=False)
        return json.dumps(value, ensure_ascii=False)

    lines += wrapped([f"{json.dumps(k)}: {short(v)}" for k, v in result.items()] + ['"top_separators": ['],
                     "  {", "   ", "")
    for i, entry in enumerate(top):
        items = [f"{json.dumps(k)}: {json.dumps(entry[k])}" for k in MCP_FIELDS if k in entry] + ["…"]
        lines += wrapped(items, "    {", "     ", "}" + ("]}" if i == len(top) - 1 else ","))
    return lines


def trimmed(path, pad=6, white=0.97):
    """A captured PNG without the page's blank margins (pixels lighter than `white`)."""
    img = mpimg.imread(path)[..., :3]
    ink = (img < white).any(axis=2)
    rows, cols = np.flatnonzero(ink.any(axis=1)), np.flatnonzero(ink.any(axis=0))
    r0, r1 = max(rows[0] - pad, 0), min(rows[-1] + pad + 1, img.shape[0])
    c0, c1 = max(cols[0] - pad, 0), min(cols[-1] + pad + 1, img.shape[1])
    return img[r0:r1, c0:c1]


def fig2_platform(ctx):
    """The platform in use: ligand browser, AI assistant, an MCP tool call (captured, not drawn)."""
    names = ("browser_table.png", "browser_detail.png", "assistant.png", "mcp_exchange.json")
    missing = [n for n in names if not (PLATFORM / n).exists()]
    if missing:
        raise FileNotFoundError(f"{', '.join(missing)} missing in {PLATFORM}: run paper/capture_platform.py")
    table, detail, assistant = (trimmed(PLATFORM / n) for n in names[:3])
    exchange = json.loads((PLATFORM / names[3]).read_text(encoding="utf-8"))

    width, gap, title = 7.1, 0.14, 0.2  # inches; `title` is the strip above each panel
    half = (width - gap) / 2

    def tall(img, w):
        return w * img.shape[0] / img.shape[1]

    h_table, h_detail, h_b = tall(table, width), tall(detail, width), tall(assistant, half)
    height = title + h_table + 0.04 + h_detail + gap + title + h_b
    fig = plt.figure(figsize=(width, height))

    def axes(x, top, w, h):  # x and top in inches from the top-left corner
        ax = fig.add_axes([x / width, (height - top - h) / height, w / width, h / height])
        ax.set_xticks([])
        ax.set_yticks([])
        for spine in ax.spines.values():
            spine.set_color(INK)
            spine.set_linewidth(0.8)
        return ax

    def heading(x, top, letter, text):
        fig.text(x / width, 1 - (top + title * 0.8) / height, letter, fontsize=10, fontweight="bold",
                 color=INK, ha="left", va="baseline")
        fig.text((x + 0.17) / width, 1 - (top + title * 0.8) / height, text, fontsize=7.5, color=INK,
                 ha="left", va="baseline")

    heading(0, 0, "a", "Dashboard: ligand browser (search result and the selected ligand)")
    a1 = axes(0, title, width, h_table)
    a1.imshow(table, interpolation="lanczos")
    a2 = axes(0, title + h_table + 0.04, width, h_detail)
    a2.imshow(detail, interpolation="lanczos")

    row = title + h_table + 0.04 + h_detail + gap
    heading(0, row, "b", "Dashboard: AI assistant")
    b = axes(0, row + title, half, h_b)
    b.imshow(assistant, interpolation="lanczos")
    heading(half + gap, row, "c", "MCP server: one tool call from an agent")
    c = axes(half + gap, row + title, half, h_b)
    c.set_facecolor("#f4f4f1")
    c.set_xlim(0, 1)
    c.set_ylim(0, 1)
    tools = exchange["tools"]
    c.text(0.04, 0.95, f"MCP client  ⇄  {exchange['server']} server (stdio)", fontsize=6.8,
           fontweight="bold", color=INK, va="top")
    c.text(0.04, 0.86, f"tools/list → {len(tools)} tools, e.g. {tools[0]}",
           fontsize=5.6, family="monospace", color=INK_2, va="top")
    c.text(0.04, 0.77, "\n".join(mcp_lines(exchange, width=60)), fontsize=5.6, family="monospace",
           color=INK, va="top", linespacing=1.35)
    return fig


# --- Figure 3 ------------------------------------------------------------------

def donut(ax, values, labels, colors, centre):
    wedges, _ = ax.pie(values, colors=colors, startangle=90, counterclock=False,
                       wedgeprops={"width": 0.38, "edgecolor": SURFACE, "linewidth": 1.0})
    ax.text(0, 0, centre, ha="center", va="center", fontsize=7.5, fontweight="bold",
            color=INK)
    ax.set_aspect("equal")
    return wedges


def fig3_overview(ctx):
    stats = ctx.stats
    groups = donor_groups()
    group_of = {s: g for g, spec in groups.items() for s in spec["subtypes"]}
    fig = plt.figure(figsize=(7.1, 7.0))

    # (a) ligands per subtype, coloured by donor chemistry, each with its core
    ax = fig.add_axes((0.395, 0.065, 0.235, 0.91))
    subtypes = sorted(stats["subtypes"].items(), key=lambda kv: (-kv[1], kv[0]))
    names = [s for s, _ in subtypes]
    counts = [n for _, n in subtypes]
    colors = [groups.get(group_of.get(s), {}).get("color", MUTED) for s in names]
    ax.barh(range(len(names)), counts, height=0.66, color=colors, edgecolor=INK,
            linewidth=0.5)
    for i, n in enumerate(counts):
        ax.text(n + 1, i, str(n), va="center", fontsize=6.6, color=INK)
    ax.set_yticks(range(len(names)))
    ax.set_yticklabels([subtype_label(s) for s in names], fontsize=7)
    ax.set_ylim(len(names) - 0.5, -0.5)
    ax.set_xlim(0, max(counts) * 1.18)
    ax.set_xlabel("Ligands")
    row_pt = 0.91 * 7.0 * 72 / len(names)
    structure_labels(ax, names, width_pt=104, height_pt=row_pt * 0.94, max_zoom=0.19)
    handles = [Rectangle((0, 0), 1, 1, facecolor=spec["color"], edgecolor=INK, linewidth=0.5)
               for spec in groups.values()]
    ax.legend(handles, [DONOR_GROUP_LABELS.get(g, g) for g in groups], loc="lower right",
              fontsize=6.6, handlelength=1.0, title="Donor chemistry", title_fontsize=6.8)
    ax.text(0.005, 0.985, "a", transform=fig.transFigure, fontsize=10, fontweight="bold",
            color=INK, ha="left", va="top")

    # (b) selectivity-trend classes
    ax = fig.add_axes((0.69, 0.66, 0.26, 0.29))
    trends = [t for t in TREND_DISPLAY_ORDER if stats["trend_classes"].get(t)]
    donut(ax, [stats["trend_classes"][t] for t in trends], trends,
          [TREND_FIGURE_COLORS[t] for t in trends], f"{stats['n_ligands']}\nligands")
    ax.legend([Rectangle((0, 0), 1, 1, facecolor=TREND_FIGURE_COLORS[t]) for t in trends],
              [f"{label(t)} ({stats['trend_classes'][t]})" for t in trends],
              loc="upper center", bbox_to_anchor=(0.5, 0.0), ncol=2, fontsize=6.2,
              handlelength=0.9, columnspacing=0.8, handletextpad=0.4)
    ax.set_title("Selectivity trend", fontsize=8, color=INK)
    ax.text(0.655, 0.985, "b", transform=fig.transFigure, fontsize=10, fontweight="bold",
            color=INK, ha="left", va="top")

    # (c) lanthanide coverage
    ax = fig.add_axes((0.745, 0.065, 0.22, 0.42))
    measured = [stats["ligands_measuring_element"].get(el, 0) for el in COVERAGE_ELEMENTS]
    ax.barh(range(len(COVERAGE_ELEMENTS)), measured, height=0.7, color=SERIES[0],
            edgecolor=INK, linewidth=0.5)
    for i, n in enumerate(measured):
        ax.text(n + 2, i, str(n), va="center", ha="left", fontsize=5.8, color=INK)
    ax.set_yticks(range(len(COVERAGE_ELEMENTS)))
    ax.set_yticklabels(COVERAGE_ELEMENTS, fontsize=6.6)
    ax.set_ylim(len(COVERAGE_ELEMENTS) - 0.4, -0.6)
    ax.set_xlim(0, max(measured) * 1.22)
    ax.set_xlabel("Ligands measured")
    ax.set_title("Lanthanide coverage", fontsize=8, color=INK)
    ax.text(0.655, 0.53, "c", transform=fig.transFigure, fontsize=10, fontweight="bold",
            color=INK, ha="left", va="top")
    return fig


# --- Figure 4 ------------------------------------------------------------------

def label_supported(record):
    """Does the ligand's own slope CI agree with its label?

    Ascending needs the whole interval above zero, Reverse below, Flat an
    interval containing zero. Curved and Mixed labels are not statements
    about the slope, so they always pass.
    """
    sel = record["selectivity"]
    trend, lo, hi = sel["trend_classification"], sel.get("slope_CI_lower"), sel.get("slope_CI_upper")
    if lo is None:
        return False
    return {"Ascending": lo > 0, "Reverse": hi < 0, "Flat": lo <= 0 <= hi}.get(trend, True)


def gallery_examples(records):
    """One representative per trend class, chosen by rule, never by name.

    Candidates: the class's ligands with a fitted slope. Ranked by: the
    extremum, if any, is smooth (paper/stats.py sharp_extremum is False); the
    label agrees with the ligand's own slope CI; its measurements are not
    shared with another ligand (duplicated_conditions); most conditions
    behind the fit; quality score; name. A curved class whose every member
    has a sharp extremum is left out: its label rests on one element off the
    curve (in v13, the U-shaped minima sit on log D = -4.0 placeholder
    cells), and the caption says so.
    """
    shared = shared_ligands(records)
    picks = {}
    for trend in GALLERY_ORDER:
        candidates = [r for r in records if r["selectivity"]["trend_classification"] == trend
                      and r["selectivity"].get("slope") is not None]
        if trend in CURVED and all(sharp_extremum(r) for r in candidates):
            continue
        if candidates:
            picks[trend] = max(candidates, key=lambda r: (
                not sharp_extremum(r), label_supported(r), r["ligand"]["Name"] not in shared,
                n_fitted_conditions(r), r["data_quality"]["overall_quality_score"],
                r["ligand"]["Name"]))
    return picks


def fitted_profiles(record):
    """Sorted (Z, log D) arrays of the conditions the trend fit uses (>= 4 lanthanides)."""
    profiles = [np.array(sorted(condition_trend_points(c)), dtype=float)
                for c in record["experimental_data"]]
    return [p for p in profiles if len(p) >= 4]


def aligned(profiles, iterations=50):
    """Each profile shifted by its own offset so that all share one shape.

    Fits log D(condition j, element i) = a_j + f_i by alternating means: a
    condition's acidity or extractant concentration lifts all of its values
    together, and conditions that cover different elements still line up
    (centring each on its own mean would not: a condition measuring only
    La-Nd under a rising trend would sit too high). Offsets average zero, so
    the result stays on the ligand's own log D scale. Returns the shifted
    profiles in the same order.
    """
    if len(profiles) < 2:
        return [p.copy() for p in profiles]
    offsets = np.array([p[:, 1].mean() for p in profiles])
    for _ in range(iterations):
        by_z = collections.defaultdict(list)
        for p, a in zip(profiles, offsets):
            for z, y in p:
                by_z[z].append(y - a)
        shape = {z: np.mean(v) for z, v in by_z.items()}
        offsets = np.array([np.mean([y - shape[z] for z, y in p]) for p in profiles])
        offsets -= offsets.mean()
    return [np.column_stack([p[:, 0], p[:, 1] - a]) for p, a in zip(profiles, offsets)]


def spread_plot(ax, profiles, color, light_alpha=0.18, dark_alpha=0.38):
    """Median line over the aligned profiles, with the range and IQR across them.

    One profile: its points joined by a line. Returns {Z: median}.
    """
    by_z = collections.defaultdict(list)
    for p in profiles:
        for z, y in p:
            by_z[z].append(y)
    zs = sorted(by_z)
    med = np.array([np.median(by_z[z]) for z in zs])
    if len(profiles) > 1:
        lo = [min(by_z[z]) for z in zs]
        hi = [max(by_z[z]) for z in zs]
        q1 = [np.percentile(by_z[z], 25) for z in zs]
        q3 = [np.percentile(by_z[z], 75) for z in zs]
        ax.fill_between(zs, lo, hi, color=color, alpha=light_alpha, linewidth=0, zorder=1)
        ax.fill_between(zs, q1, q3, color=color, alpha=dark_alpha, linewidth=0, zorder=2)
    ax.plot(zs, med, color=color, linewidth=1.6, zorder=3, marker="o", markersize=2.6,
            markeredgewidth=0)
    return dict(zip(zs, med))


def acid_range(conditions):
    """'0.1-0.2 M HNO3' for a group of conditions (subscripts as Unicode)."""
    acids = {str(c.get("Acid_Name") or "").upper().replace("HCL", "HCl") for c in conditions}
    acid = acids.pop() if len(acids) == 1 else "acid"
    acid = acid.replace("HNO3", "HNO₃").replace("H2SO4", "H₂SO₄")
    values = sorted(c["Acid_Concentration_M"] for c in conditions
                    if isinstance(c.get("Acid_Concentration_M"), (int, float)))
    if not values:
        return acid
    span = f"{values[0]:g}" if values[0] == values[-1] else f"{values[0]:g}–{values[-1]:g}"
    return f"{span} M {acid}"


def condition_votes(record):
    votes = []
    for c in record["experimental_data"]:
        if len(condition_trend_points(c)) >= 4:
            votes.append(classify_selectivity_trend([c])["trend_classification"])
    return votes


def same_table_extrema(records, record):
    """Other ligands in the table a one-condition curved example came from, and how
    many of them show the same kind of extremum within one element of it.
    Returns (k, n, extremum element).

    A dip (or peak) that every ligand of one table shares at the same element
    points at that table's method, not at one ligand. None when the example has
    several conditions or no table siblings.
    """
    sel = record["selectivity"]
    conditions = [c for c in record["experimental_data"] if len(condition_trend_points(c)) >= 4]
    if sel["trend_classification"] not in CURVED or len(conditions) != 1:
        return None
    ref = conditions[0]["reference"]
    key = (ref.get("DOI"), ref.get("Locator"))
    sign = 1 if sel["trend_classification"] == "Peaked" else -1
    near = {Z[sel["peak_element"]] + d for d in (-1, 0, 1)}
    hits, n = set(), 0
    for r in records:
        if r is record:
            continue
        for c in r["experimental_data"]:
            if (c["reference"].get("DOI"), c["reference"].get("Locator")) != key:
                continue
            pts = sorted(condition_trend_points(c))
            n += 1
            local = [z for (z0, y0), (z, y), (z1, y1) in zip(pts, pts[1:], pts[2:])
                     if z in near and sign * (y - y0) > 0 and sign * (y - y1) > 0]
            if local:
                hits.add(r["ligand"]["Name"])
            break
    if not n:
        return None
    return len(hits), n, sel["peak_element"]


STRUCTURE_BOX = (0.5, 0.36)  # the largest drawing fills this much of a panel


def data_vertices(ax):
    """Every plotted vertex (lines, bands, markers) in axes fractions."""
    to_axes = ax.transData + ax.transAxes.inverted()
    points = [line.get_xydata() for line in ax.lines]
    points += [path.vertices for coll in ax.collections for path in coll.get_paths()]
    points = [p for p in points if len(p)]
    return to_axes.transform(np.vstack(points)) if points else np.empty((0, 2))


def place_structure(ax, image, zoom, size, reserved):
    """Draw a ligand in the panel corner the data leave empty.

    Corners are scored by the plotted vertices they would cover; a corner
    overlapping a note or legend (reserved, axes fractions x0, y0, x1, y1) is
    out. If every corner covers data, the y-range grows at the top until the
    drawing fits above it.
    """
    from matplotlib.offsetbox import AnnotationBbox, OffsetImage
    w, h = size
    pad = 0.02

    def corners():
        for x0 in (pad, 1 - pad - w):
            for y0 in (1 - pad - h, 0.22, pad):
                yield x0, y0

    def score(x0, y0, pts):
        if any(x0 < rx1 and x0 + w > rx0 and y0 < ry1 and y0 + h > ry0
               for rx0, ry0, rx1, ry1 in reserved):
            return None
        inside = ((pts[:, 0] > x0 - pad) & (pts[:, 0] < x0 + w + pad)
                  & (pts[:, 1] > y0 - pad) & (pts[:, 1] < y0 + h + pad))
        return int(inside.sum())

    pts = data_vertices(ax)
    scored = [(s, -y0, x0, y0) for x0, y0 in corners() if (s := score(x0, y0, pts)) is not None]
    best = min(scored) if scored else None
    if best is None or best[0] > 0:
        low, high = ax.get_ylim()
        ax.set_ylim(low, high + (h + 2 * pad) / (1 - h - 2 * pad) * (high - low))
        pts = data_vertices(ax)
        top = [(score(x0, 1 - pad - h, pts), x0) for x0 in (pad, 1 - pad - w)]
        x0 = min((s, x) for s, x in top if s is not None)[1]
        best = (0, 0, x0, 1 - pad - h)
    _, _, x0, y0 = best
    ax.add_artist(AnnotationBbox(OffsetImage(image, zoom=zoom, interpolation="lanczos"),
                                 (x0 + w / 2, y0 + h / 2), xycoords="axes fraction",
                                 box_alignment=(0.5, 0.5), frameon=False, pad=0, zorder=6))


def fig4_trend_gallery(ctx, cols=3, figsize=None, only=None):
    """only: a subset of trend classes (slides split the gallery); letters stay global."""
    all_picks = gallery_examples(ctx.records)
    letters = {trend: string.ascii_lowercase[i] for i, trend in enumerate(all_picks)}
    picks = {t: r for t, r in all_picks.items() if only is None or t in only}
    rows = (len(picks) + cols - 1) // cols
    fig, axes = plt.subplots(rows, cols, figsize=figsize or (7.1, 2.6 * rows), sharex=True,
                             squeeze=False)
    fig.subplots_adjust(left=0.08, right=0.99, top=0.93, bottom=0.07, wspace=0.3,
                        hspace=0.26)
    # every ligand at one scale (one bond length), as large as the biggest allows
    # inside a box of STRUCTURE_BOX (axes fractions) in a corner of its panel
    images = {trend: ligand_image(r["ligand"]["SMILES"]) for trend, r in picks.items()}
    box = axes.flat[0].get_position()
    box_w = STRUCTURE_BOX[0] * box.width * fig.get_figwidth() * 72
    box_h = STRUCTURE_BOX[1] * box.height * fig.get_figheight() * 72
    zoom = min((min(box_w / im.shape[1], box_h / im.shape[0])
                for im in images.values() if im is not None), default=None)
    for ax in axes.flat[len(picks):]:
        ax.set_visible(False)
    for i in range(len(picks), rows * cols):
        axes.flat[i - cols].xaxis.set_tick_params(labelbottom=True)
    for k, (ax, (trend, record)) in enumerate(zip(axes.flat, picks.items())):
        sel = record["selectivity"]
        color = TREND_FIGURE_COLORS[trend]
        n = n_fitted_conditions(record)
        conditions = f"{n} condition{'s' if n != 1 else ''}"
        if n > 1:
            conditions += " (offsets removed)"
        if trend == "Mixed":
            # one band per vote: which conditions agree, and what sets them apart
            groups = collections.defaultdict(list)
            for c in record["experimental_data"]:
                if len(condition_trend_points(c)) >= 4:
                    groups[classify_selectivity_trend([c])["trend_classification"]].append(c)
            order = sorted(groups, key=lambda t: (-len(groups[t]), TREND_DISPLAY_ORDER.index(t)))
            handles = []
            for vote in order:
                profiles = [np.array(sorted(condition_trend_points(c)), dtype=float)
                            for c in groups[vote]]
                vote_color = TREND_FIGURE_COLORS.get(vote, MUTED)
                spread_plot(ax, aligned(profiles), vote_color, light_alpha=0.15,
                            dark_alpha=0.3)
                handles.append((plt.Line2D([], [], color=vote_color, linewidth=1.6),
                                f"{label(vote)} ({len(groups[vote])}): "
                                f"{acid_range(groups[vote])}"))
            ax.legend(*zip(*handles), loc="lower left", fontsize=5.6, frameon=True,
                      framealpha=1, edgecolor=INK, handlelength=1.2,
                      title="condition votes", title_fontsize=5.6)
            note = None
        else:
            medians = spread_plot(ax, aligned(fitted_profiles(record)), color)
            if trend in CURVED:
                element = sel["peak_element"]
                ax.scatter([Z[element]], [medians[Z[element]]], s=60, facecolors="none",
                           edgecolors=INK, linewidths=1.0, zorder=4)
                note = f"{'maximum' if trend == 'Peaked' else 'minimum'} at {element} " \
                       f"(ringed)\n{conditions}"
                siblings = same_table_extrema(ctx.records, record)
                if siblings:
                    hits, others, element = siblings
                    note += (f"\nsame table: {'peak' if trend == 'Peaked' else 'dip'} near "
                             f"{element} in {hits}/{others} others")
            else:
                note = (f"slope {sel['slope']:+.3f} dex/Z, 95% CI [{sel['slope_CI_lower']:+.3f}, "
                        f"{sel['slope_CI_upper']:+.3f}]\n{conditions}")
        image = images.get(trend)
        ax.set_title(f"{label(trend)}: {record['ligand']['Name']}", fontsize=8, color=INK,
                     loc="left")
        if note:
            ax.text(0.03, 0.04, note, transform=ax.transAxes, fontsize=5.5, color=INK,
                    bbox={"facecolor": SURFACE, "edgecolor": "none", "pad": 0.8}, zorder=5)
        ax.margins(y=0.12)
        low, high = ax.get_ylim()
        if high - low < MIN_GALLERY_SPAN_DEX:
            middle = (low + high) / 2
            low, high = middle - MIN_GALLERY_SPAN_DEX / 2, middle + MIN_GALLERY_SPAN_DEX / 2
        # room below the data for the note or legend
        ax.set_ylim(low - (0.55 if trend == "Mixed" else 0.22) * (high - low), high)
        ax.set_xlim(56.3, 71.7)
        ax.set_xticks([57, 59, 61, 63, 65, 67, 69, 71])
        ax.set_xticklabels([SYMBOL[z] for z in [57, 59, 61, 63, 65, 67, 69, 71]])
        if image is not None:
            size = (image.shape[1] * zoom / (box.width * fig.get_figwidth() * 72),
                    image.shape[0] * zoom / (box.height * fig.get_figheight() * 72))
            reserved = [(0.0, 0.0, 0.85, 0.42 if trend == "Mixed" else 0.2)]  # note / legend
            place_structure(ax, image, zoom, size, reserved)
        panel_label(ax, letters[trend], x=-0.1, y=1.02)
    for ax in axes[:, 0]:
        ax.set_ylabel("log D")
    return fig


# --- Figure 5 ------------------------------------------------------------------

def natural_ticks(axis, low, high, steps=(1,), every=1):
    """Label a log10 axis in natural units: ticks at 10^k x steps inside [low, high]."""
    from matplotlib.ticker import FixedLocator, FuncFormatter
    ticks = [np.log10(m) + k for k in range(int(np.floor(low)) - 1, int(np.ceil(high)) + 2)
             for m in steps if low <= np.log10(m) + k <= high and k % every == 0]
    axis.set_major_locator(FixedLocator(ticks))
    axis.set_major_formatter(FuncFormatter(lambda v, _: f"{round(10.0 ** v, 6):g}"))
    minor = [np.log10(m) + k for k in range(int(np.floor(low)) - 1, int(np.ceil(high)) + 2)
             for m in range(2, 10) if low <= np.log10(m) + k <= high]
    if steps == (1,):
        axis.set_minor_locator(FixedLocator(minor))


def violin_box(ax, groups, order, color, font=6.6):
    for y, subtype in enumerate(order):
        values = np.array(groups[subtype])
        if len(set(values)) > 1:
            parts = ax.violinplot(values, positions=[y], vert=False, widths=0.82,
                                  showextrema=False, bw_method=0.5)
            for body in parts["bodies"]:
                body.set_facecolor(color)
                body.set_edgecolor(color)
                body.set_alpha(0.35)
                body.set_linewidth(0.6)
        ax.boxplot(values, positions=[y], vert=False, widths=0.24, showfliers=True,
                   patch_artist=True,
                   boxprops={"facecolor": SURFACE, "edgecolor": INK, "linewidth": 0.8},
                   medianprops={"color": INK, "linewidth": 1.5},
                   whiskerprops={"color": INK, "linewidth": 0.8},
                   capprops={"color": INK, "linewidth": 0.8},
                   flierprops={"marker": "o", "markersize": 2, "markerfacecolor": INK,
                               "markeredgewidth": 0})
        ax.text(1.02, y, f"{len(values)}", transform=ax.get_yaxis_transform(), va="center",
                ha="left", fontsize=font, color=INK)
    ax.text(1.02, -0.7, "n", transform=ax.get_yaxis_transform(), va="bottom", ha="left",
            fontsize=font, color=INK, fontstyle="italic")
    ax.axvline(0, color=MUTED, linewidth=0.8, linestyle=(0, (3, 2)), zorder=0)
    ax.tick_params(axis="y", length=0)



def fig5_eu_la_showcase(ctx, figsize=None, font=7.0):
    """(a) log D(Eu), (b) log D(La) on the same axis, (c) SF(Eu/La), by subtype."""
    kept, log_eu, log_sf = showcase_groups(ctx.records)
    log_la = log_d_by_subtype(ctx.records, "La")
    figsize = figsize or (7.1, 0.44 * len(kept) + 1.2)
    fig, (ax_a, ax_b, ax_c) = plt.subplots(1, 3, figsize=figsize, sharey=True,
                                           gridspec_kw={"wspace": 0.28})
    small = font * 0.9
    violin_box(ax_a, log_eu, kept, SERIES[0], small)
    violin_box(ax_b, log_la, kept, SERIES[2], small)
    violin_box(ax_c, log_sf, kept, SERIES[1], small)  # density in log SF; ticks read as SF
    low = min(min(v[s]) for v in (log_eu, log_la) for s in kept if v.get(s))
    high = max(max(v[s]) for v in (log_eu, log_la) for s in kept if v.get(s))
    pad = 0.05 * (high - low)
    for ax in (ax_a, ax_b):
        ax.set_xlim(low - pad, high + pad)
        ax.tick_params(axis="x", labelsize=font)
    ax_c.tick_params(axis="x", labelsize=font)
    ax_a.set_yticks(range(len(kept)))
    ax_a.set_yticklabels([subtype_label(s) for s in kept], fontsize=font)
    ax_a.set_ylim(len(kept) - 0.5, -0.7)
    fig.subplots_adjust(left=0.3, right=0.96, top=0.95,
                        bottom=font * 4.6 / 72 / figsize[1])
    row_pt = 0.95 * figsize[1] * 72 / (len(kept) + 0.2)
    structure_labels(ax_a, kept, width_pt=figsize[0] * 72 * 0.11, height_pt=row_pt * 0.8,
                     max_zoom=0.17 * font / 7)
    low_sf, high_sf = ax_c.get_xlim()
    natural_ticks(ax_c.xaxis, low_sf, high_sf, every=2 if high_sf - low_sf > 4 else 1)
    ax_a.set_xlabel("log D(Eu)", fontsize=font * 1.05)
    ax_b.set_xlabel("log D(La)", fontsize=font * 1.05)
    ax_c.set_xlabel("SF(Eu/La)", fontsize=font * 1.05)
    for ax, letter in zip((ax_a, ax_b, ax_c), "abc"):
        panel_label(ax, letter, x=-0.02)
    return fig


# --- Figure 6 ------------------------------------------------------------------

CASE_MAX_CONDITIONS = 3


def uncertainty_case(ctx):
    """The shape-labelled ligand whose label its own sigma most often overturns.

    Rule: stored label is a shape (not Mixed); measurements not shared with
    another ligand; at most CASE_MAX_CONDITIONS fitted conditions, so every
    point stays legible; lowest share of perturbed draws keeping the label;
    then more conditions; then name.
    """
    shared = shared_ligands(ctx.records)
    by_name = {r["ligand"]["Name"]: r for r in ctx.records}
    candidates = [(share, -n_fitted_conditions(by_name[name]), name)
                  for name, (stored, share, _) in ctx.stability.items()
                  if stored not in NOT_A_SHAPE and name not in shared
                  and n_fitted_conditions(by_name[name]) <= CASE_MAX_CONDITIONS]
    return by_name[min(candidates)[2]] if candidates else None


def fig6_uncertainty(ctx):
    record = uncertainty_case(ctx)
    fig, (ax_a, ax_b, ax_c) = plt.subplots(1, 3, figsize=(7.1, 2.75),
                                           gridspec_kw={"width_ratios": [1.45, 0.8, 1.25],
                                                        "wspace": 0.42})
    if record is None:
        return fig
    stored, share, counts = ctx.stability[record["ligand"]["Name"]]
    fitted = [c for c in record["experimental_data"] if len(condition_trend_points(c)) >= 4]

    # (a) perturbed profiles, coloured by the label each draw receives
    for k, (conditions, draw_label) in enumerate(perturbed_draws(record)):
        if k >= 60:
            break
        for c in conditions:
            pts = sorted(condition_trend_points(c))
            if len(pts) >= 4:
                x, y = zip(*pts)
                ax_a.plot(x, y, color=TREND_FIGURE_COLORS.get(draw_label, MUTED),
                          linewidth=0.5, alpha=0.35, zorder=1)
    for c in fitted:
        cells = [(Z[el], v["log_value"], v.get("sigma_log_D") or 0.0)
                 for el, v in c["distribution_ratios"].items()
                 if el in Z and v.get("log_value") is not None and not is_censored(v)]
        x, y, s = zip(*sorted(cells))
        ax_a.errorbar(x, y, yerr=2 * np.array(s), fmt="o", markersize=2.8, color=INK,
                      elinewidth=0.8, capsize=1.5, zorder=3)
    if "Peaked" in (stored, *counts):
        # the rule the draws cross: an interior extremum 0.3 dex beyond both ends
        pts = sorted(condition_trend_points(fitted[0]))
        top = max(p[1] for p in pts)
        ax_a.axhline(top - 0.3, color=INK, linewidth=0.7, linestyle=(0, (3, 2)), zorder=2)
        ax_a.text(0.5, (top - 0.3) - 0.02, "maximum − 0.3 dex",
                  transform=ax_a.get_yaxis_transform(), ha="center", va="top", fontsize=5.8,
                  color=INK, bbox={"facecolor": SURFACE, "edgecolor": "none", "pad": 0.5})
    ax_a.set_xticks([57, 59, 61, 63, 65, 67, 69, 71])
    ax_a.set_xticklabels([SYMBOL[z] for z in [57, 59, 61, 63, 65, 67, 69, 71]])
    ax_a.set_xlim(56.3, 71.7)
    ax_a.set_ylabel("log D")
    ax_a.set_title(f"{record['ligand']['Name']}: stored {label(stored)}", fontsize=8,
                   loc="left", color=INK)
    panel_label(ax_a, "a", x=-0.2)

    # (b) share of draws per label
    labels = [t for t in TREND_DISPLAY_ORDER if counts.get(t)]
    total = sum(counts.values())
    shares = [counts[t] / total for t in labels]
    ax_b.barh(range(len(labels)), shares, color=[TREND_FIGURE_COLORS[t] for t in labels],
              edgecolor=INK, linewidth=0.5, height=0.6)
    for i, (t, v) in enumerate(zip(labels, shares)):
        ax_b.text(v + 0.03, i, f"{v:.0%}" + (" (stored)" if t == stored else ""),
                  va="center", fontsize=6.3, color=INK)
    ax_b.set_yticks(range(len(labels)))
    ax_b.set_yticklabels([label(t) for t in labels])
    ax_b.set_ylim(len(labels) - 0.4, -0.6)
    ax_b.set_xlim(0, 1.6)
    ax_b.set_xticks([0, 0.5, 1])
    ax_b.set_xticklabels(["0", "50", "100"])
    ax_b.set_xlabel(f"Draws (%), n = {total}")
    panel_label(ax_b, "b", x=-0.45)

    # (c) adjacent-pair separation factors with 2 sigma
    order = {f"{h}/{l}": i for i, (h, l) in enumerate(ADJACENT)}
    per_condition = [adjacent_sf(c) for c in fitted]
    pairs = sorted({r[0] for rows_ in per_condition for r in rows_}, key=order.get)
    xpos = {pair: i for i, pair in enumerate(pairs)}
    rows = [r for rows_ in per_condition for r in rows_]
    for k, rows_ in enumerate(per_condition):
        shift = (k - (len(per_condition) - 1) / 2) * 0.22
        for pair, value, sigma in rows_:
            settled = abs(value) >= 2 * sigma
            ax_c.errorbar(xpos[pair] + shift, value, yerr=2 * sigma, fmt="o", markersize=3.2,
                          color=INK, markerfacecolor=INK if settled else SURFACE,
                          elinewidth=0.8, capsize=1.5)
    ax_c.axhline(0, color=MUTED, linewidth=0.8, linestyle=(0, (3, 2)))
    ax_c.set_xticks(range(len(pairs)))
    ax_c.set_xticklabels(pairs, rotation=90, fontsize=6)
    natural_ticks(ax_c.yaxis, *ax_c.get_ylim(), steps=(0.5, 0.7, 1, 1.5, 2, 3, 5))
    ax_c.set_ylabel("SF, adjacent pair (heavier/lighter)")
    ax_c.set_xlim(-0.7, len(pairs) - 0.3)
    undetermined = sum(1 for _, v, s in rows if abs(v) < 2 * s)
    ax_c.set_title(f"open: SF within 2σ of 1 ({undetermined} of {len(rows)})", fontsize=7,
                   loc="left", color=INK)
    panel_label(ax_c, "c", x=-0.25)
    fig.subplots_adjust(left=0.08, right=0.99, top=0.9, bottom=0.2)
    return fig


# --- Figure 7 ------------------------------------------------------------------

STACK_ORDER = GALLERY_ORDER + ("Insufficient_Data",)


def ink_on(color):
    """Black or white, whichever reads on a fill colour."""
    r, g, b = matplotlib.colors.to_rgb(color)
    return INK if 0.2126 * r + 0.7152 * g + 0.0722 * b > 0.5 else SURFACE


def fig7_trend_by_subtype(ctx, figsize=(7.1, 6.4), font=7.0, structures=True):
    """Share of each subtype's ligands in each trend class (100 % bars, Fig 3a order)."""
    gaps = ctx.stats["by_subtype"]
    names = list(gaps)
    fig = plt.figure(figsize=figsize)
    bottom = font * 4.8 / 72 / figsize[1]       # tick labels and x label
    top = 1 - font * 4.4 / 72 / figsize[1]      # two-row legend
    ax = fig.add_axes((0.40, bottom, 0.53, top - bottom))
    for y, subtype in enumerate(names):
        n, left = gaps[subtype]["n_ligands"], 0.0
        for trend in STACK_ORDER:
            k = gaps[subtype]["trend_classes"].get(trend, 0)
            if not k:
                continue
            color = TREND_FIGURE_COLORS[trend]
            ax.barh(y, k / n, left=left, height=0.7, color=color, edgecolor=SURFACE,
                    linewidth=0.6)
            if k / n >= 0.09:
                ax.text(left + k / n / 2, y, str(k), ha="center", va="center",
                        fontsize=font * 0.85, color=ink_on(color))
            left += k / n
        ax.text(1.015, y, str(n), transform=ax.get_yaxis_transform(), va="center",
                fontsize=font * 0.9, color=INK)
    ax.text(1.015, -0.75, "n", transform=ax.get_yaxis_transform(), va="bottom",
            fontsize=font * 0.9, color=INK, fontstyle="italic")
    ax.set_yticks(range(len(names)))
    ax.set_yticklabels([subtype_label(s) for s in names], fontsize=font)
    ax.set_ylim(len(names) - 0.5, -0.5)
    ax.set_xlim(0, 1)
    ax.set_xticks([0, 0.25, 0.5, 0.75, 1])
    ax.set_xticklabels(["0", "25", "50", "75", "100"], fontsize=font)
    ax.set_xlabel("Ligands in each trend class (%)", fontsize=font * 1.1)
    if structures:
        row_pt = (top - bottom) * figsize[1] * 72 / len(names)
        structure_labels(ax, names, width_pt=figsize[0] * 72 * 0.14, height_pt=row_pt * 0.92,
                         max_zoom=0.19 * font / 7)
    else:
        ax.tick_params(axis="y", length=0)
    present = [t for t in STACK_ORDER if any(g["trend_classes"].get(t) for g in gaps.values())]
    ax.legend([Rectangle((0, 0), 1, 1, facecolor=TREND_FIGURE_COLORS[t]) for t in present],
              [label(t) for t in present], loc="lower center", bbox_to_anchor=(0.5, 1.005),
              ncol=4, fontsize=font * 0.9, handlelength=1.0, columnspacing=1.0,
              handletextpad=0.4)
    return fig


# --- Figure 8 ------------------------------------------------------------------

def gap_columns(gaps):
    """(header, [(value for colour, label, denominator)] per subtype, larger_is_worse)."""
    def share(k, n):
        return (k / n if n else None, f"{k}/{n}" if n else "–", n)

    cols = [
        ("No trend class\n(ligands)",
         [share(g["no_trend_class"], g["n_ligands"]) for g in gaps.values()], True),
        ("Lanthanides\nmeasured\n(median of 14)",
         [(g["median_lanthanides_measured"], f"{g['median_lanthanides_measured']:g}",
           g["n_ligands"]) for g in gaps.values()], False),
        ("Conditions\nper ligand\n(median)",
         [(g["median_conditions"], f"{g['median_conditions']:g}", g["n_ligands"])
          for g in gaps.values()], False),
        ("Class kept in\n< 90% of σ draws\n(fitted ligands)",
         [share(g["fitted_class_kept_under_90pct"], g["n_fitted"]) for g in gaps.values()], True),
        ("Adjacent SF\nwithin 2σ of 1",
         [((g["adjacent_sf_within_2_sigma"] / g["n_adjacent_sf"]) if g["n_adjacent_sf"] else None,
           f"{g['adjacent_sf_within_2_sigma'] / g['n_adjacent_sf']:.0%}" if g["n_adjacent_sf"]
           else "–", g["n_ligands_with_adjacent_sf"]) for g in gaps.values()], True),
    ]
    return cols


GAP_MIN_LIGANDS = 3  # a cell resting on fewer ligands is drawn hollow


def fig8_data_gaps(ctx, figsize=(7.1, 6.0), font=7.0):
    """Where more experiments would change the answer, per subtype (Fig 3a order).

    Each column is shaded on its own scale, darker = larger gap; a cell
    resting on fewer than GAP_MIN_LIGANDS ligands is left white with grey text.
    """
    from matplotlib.colors import LinearSegmentedColormap
    gaps = ctx.stats["by_subtype"]
    names = list(gaps)
    cols = gap_columns(gaps)
    ramp = LinearSegmentedColormap.from_list("gap", [SURFACE, SERIES[1]])
    fig = plt.figure(figsize=figsize)
    ax = fig.add_axes((0.25, 0.02, 0.73, 0.84))
    for j, (_, cells, worse) in enumerate(cols):
        values = [v for v, _, n in cells if v is not None and n >= GAP_MIN_LIGANDS]
        lo, hi = (min(values), max(values)) if values else (0, 1)
        for i, (v, text, n) in enumerate(cells):
            solid = v is not None and n >= GAP_MIN_LIGANDS
            if solid:
                need = (v - lo) / (hi - lo) if hi > lo else 0.0
                fill = ramp(0.08 + 0.92 * (need if worse else 1 - need))
            else:
                fill = SURFACE
            ax.add_patch(Rectangle((j, i), 1, 1, facecolor=fill, edgecolor=LIGHT,
                                   linewidth=0.5))
            ax.text(j + 0.5, i + 0.5, text, ha="center", va="center", fontsize=font,
                    color=INK if solid else MUTED)
    ax.set_xlim(0, len(cols))
    ax.set_ylim(len(names), 0)
    ax.set_xticks([j + 0.5 for j in range(len(cols))])
    ax.set_xticklabels([h for h, _, _ in cols], fontsize=font)
    ax.xaxis.tick_top()
    ax.tick_params(axis="x", length=0, pad=3)
    ax.set_yticks([i + 0.5 for i in range(len(names))])
    ax.set_yticklabels([f"{subtype_label(s)} ({gaps[s]['n_ligands']})" for s in names],
                       fontsize=font)
    ax.tick_params(axis="y", length=0)
    for side in ax.spines.values():
        side.set_visible(True)
    return fig


# --- Supporting information ---------------------------------------------------

def ci_width_points(records):
    xs, ys = [], []
    for r in records:
        sel = r["selectivity"]
        if sel.get("slope_CI_lower") is None:
            continue
        xs.append(n_fitted_conditions(r))
        ys.append(sel["slope_CI_upper"] - sel["slope_CI_lower"])
    return xs, ys


def figS1_ci_width(ctx):
    xs, ys = ci_width_points(ctx.records)
    rho, _ = spearmanr(xs, ys)
    jitter = np.exp(np.random.default_rng(3).uniform(-0.06, 0.06, len(xs)))
    fig, ax = plt.subplots(figsize=(3.5, 3.0))
    ax.scatter(np.array(xs) * jitter, ys, s=14, color=SERIES[0], alpha=0.8,
               edgecolors=SURFACE, linewidths=0.5)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xticks([1, 2, 5, 10, 20, 50])
    ax.set_yticks([0.01, 0.03, 0.1, 0.3, 1])
    ax.xaxis.set_major_formatter(FormatStrFormatter("%g"))
    ax.yaxis.set_major_formatter(FormatStrFormatter("%g"))
    ax.minorticks_off()
    ax.set_xlabel("Conditions behind the fit")
    ax.set_ylabel("95% CI width of slope (dex per Z)")
    ax.text(0.97, 0.95, f"Spearman ρ = {rho:.2f} (n = {len(xs)})", transform=ax.transAxes,
            ha="right", va="top", fontsize=7, color=INK)
    fig.tight_layout()
    return fig


CROSS_SOURCE_CSV = ROOT / "outputs" / "validation" / "cross_source_agreement.csv"


def cross_source(records):
    """(relaxed pairs, number of strict pairs) of cells two sources transcribed for one condition.

    Pairing needs the source labels; the public release has none and reads the
    pairs scripts/validation/cross_source_agreement.py wrote here instead.
    """
    if has_source_labels(records):
        strict = cross_source_pairs(records, temp_tol=2.0)
        return cross_source_pairs(records, temp_tol=5.0), len(strict)
    rows = list(csv.DictReader(CROSS_SOURCE_CSV.open(encoding="utf-8")))
    pairs = [{**row, "same_paper": row["same_paper"] == "True",
              "log_D_a": float(row["log_D_a"]), "log_D_b": float(row["log_D_b"])} for row in rows]
    return pairs, sum(1 for row in rows if row["tier"] == "strict")


def figS2_cross_source(ctx):
    pairs, n_strict = cross_source(ctx.records)
    fig, ax = plt.subplots(figsize=(3.4, 3.4))
    for same, color, text in ((True, SERIES[0], "same paper"),
                              (False, SERIES[1], "different papers")):
        subset = [p for p in pairs if p["same_paper"] == same]
        if subset:
            ax.scatter([p["log_D_a"] for p in subset], [p["log_D_b"] for p in subset], s=16,
                       color=color, edgecolors=SURFACE, linewidths=0.6,
                       label=f"{text}, n = {len(subset)}", zorder=3)
    values = [v for p in pairs for v in (p["log_D_a"], p["log_D_b"])]
    lo, hi = (min(values) - 0.3, max(values) + 0.3) if values else (-2, 2)
    ax.plot([lo, hi], [lo, hi], color=MUTED, linewidth=1, zorder=1)
    ax.set_xlim(lo, hi)
    ax.set_ylim(lo, hi)
    ax.set_aspect("equal")
    ax.set_xlabel("log D, one source")
    ax.set_ylabel("log D, the other source")
    ax.legend(loc="upper left", fontsize=6.8, handletextpad=0.2)
    ax.text(0.98, 0.02, f"{n_strict} strict / {len(pairs)} relaxed pairs",
            transform=ax.transAxes, ha="right", fontsize=6.6, color=INK)
    fig.tight_layout()
    return fig



def figS3_sf_viewer(ctx):
    """The separation-factor viewer: its element selectors and the plot they drive (captured)."""
    names = ("sf_selectors.png", "sf_plot.png")
    missing = [n for n in names if not (PLATFORM / n).exists()]
    if missing:
        raise FileNotFoundError(f"{', '.join(missing)} missing in {PLATFORM}: "
                                "run paper/capture_platform.py --only sf_viewer")
    selectors, plot = (trimmed(PLATFORM / n, pad=30) for n in names)
    width, gap, title = 7.1, 0.12, 0.2
    h_a, h_b = (width * img.shape[0] / img.shape[1] for img in (selectors, plot))
    height = title + h_a + gap + title + h_b
    fig = plt.figure(figsize=(width, height))
    for letter, text, top, img, h in (("a", "Element selectors", 0, selectors, h_a),
                                      ("b", "log SF by ligand subtype", title + h_a + gap, plot, h_b)):
        fig.text(0, 1 - (top + title * 0.8) / height, letter, fontsize=10, fontweight="bold",
                 color=INK, ha="left", va="baseline")
        fig.text(0.17 / width, 1 - (top + title * 0.8) / height, text, fontsize=7.5, color=INK,
                 ha="left", va="baseline")
        ax = fig.add_axes([0, (height - top - title - h) / height, 1, h / height])
        ax.set_xticks([])
        ax.set_yticks([])
        for spine in ax.spines.values():
            spine.set_color(INK)
            spine.set_linewidth(0.8)
        ax.imshow(img, interpolation="lanczos")
    return fig


FIGURES = {
    "fig1_pipeline": fig1_pipeline,
    "fig2_platform": fig2_platform,
    "fig3_overview": fig3_overview,
    "fig4_trend_gallery": fig4_trend_gallery,
    "fig5_eu_la_showcase": fig5_eu_la_showcase,
    "fig6_uncertainty": fig6_uncertainty,
    "fig7_trend_by_subtype": fig7_trend_by_subtype,
    "fig8_data_gaps": fig8_data_gaps,
    "figS1_ci_width": figS1_ci_width,
    "figS2_cross_source": figS2_cross_source,
    "figS3_sf_viewer": figS3_sf_viewer,
}


# --- Captions ------------------------------------------------------------------

def captions(ctx):
    """Captions, kept short: method detail is in the paper's Methods and ESI, findings
    in its text. Every number is computed from the records or stats.json."""
    records, stats = ctx.records, ctx.stats
    kept, log_d, log_sf = showcase_groups(records)
    stab = stats["trend_label_stability"]
    case = uncertainty_case(ctx)
    xs, ys = ci_width_points(records)
    rho, _ = spearmanr(xs, ys)
    pairs, _ = cross_source(records)
    diff = np.abs([p["log_D_a"] - p["log_D_b"] for p in pairs])

    case_text, case_sf = "No case", "No case"
    if case is not None:
        stored, share, counts = ctx.stability[case["ligand"]["Name"]]
        total = sum(counts.values())
        rest = ", ".join(f"{label(t)} {v / total:.0%}" for t, v in counts.most_common() if t != stored)
        fitted = [c for c in case["experimental_data"] if len(condition_trend_points(c)) >= 4]
        rows = [x for c in fitted for x in adjacent_sf(c)]
        open_pairs = sum(1 for _, v, s in rows if abs(v) < 2 * s)
        case_text = (f"{case['ligand']['Name']} ({subtype_label(case['ligand']['Subtype'])}), stored as "
                     f"{label(stored)}, keeps that class in {share:.0%} of draws ({rest})")
        if "Peaked" in (stored, *counts):
            case_text += "; dashed line, the Peaked threshold"
        case_sf = f"{open_pairs} of {len(rows)}"
    capture = json.loads((PLATFORM / "capture.json").read_text(encoding="utf-8"))
    sf_a, sf_b = capture["sf_viewer"]["pair"].split("/")
    return f"""# Figure captions (draft)

Generated by `python paper/figures.py` from {len(records)} ligand records; do
not edit by hand. Trend class "Descending" is stored as `Reverse`.

## Figure 1 — fig1_pipeline

**REEL-QUEST.** Distribution ratios for {stats["n_ligands"]} extractants ({stats["n_d_values"]:,} values,
{stats["n_conditions"]:,} conditions, {stats["n_distinct_primary_dois"]} original papers; generic cores of the example
ligand families) are stored as one record per ligand and served to people and
AI agents through a dashboard with an AI assistant, an MCP server, an agent skill, a REST
API and a Python package.

## Figure 2 — fig2_platform

**The platform in use**, captured from the running Streamlit dashboard. (a) Dashboard ligand browser: a name
search for {capture["browser"]["ligand"]} lists the matching ligands; selecting a row opens its
structure and selectivity-trend fit (median log D per element, with the 95% interval of the
slope); the chemical identifiers, conditions and references below are not shown. (b) The dashboard's AI assistant answering a question from database query results.
(c) An MCP client calling the server's
{capture["mcp_call"]["tool"]} tool; JSON abridged (…).

## Figure 3 — fig3_overview

**Database overview.** (a) Ligands per structural subtype, colored by donor chemistry,
with each subtype's generic core. (b) Selectivity-trend class of all {stats["n_ligands"]} ligands.
(c) Ligands with at least one measured D for each element.

## Figure 4 — fig4_trend_gallery

**One representative ligand per trend class.** Conditions measuring at
least 4 lanthanides, each shifted by its fitted offset; line, median; dark and light
bands, interquartile and full range across conditions. Slope and 95% interval are the
stored values; for curved classes the stored extremum is ringed. Mixed: one band per
group of conditions sharing a class, labeled with their acid concentrations.

## Figure 5 — fig5_eu_la_showcase

**Extraction strength and Eu/La selectivity by subtype.** (a) log D(Eu), (b) log D(La)
and (c) SF(Eu/La) per condition, for the {len(kept)} subtypes with at least 5 conditions measuring
both, ordered by median SF. Violins, density; boxes, quartiles and median; n, conditions.
Measurements stored on several structural variants count once; conditions measured with
another extractant present are left out.

## Figure 6 — fig6_uncertainty

**Measurement uncertainty can change a ligand's class.** Every log D is perturbed by its
σ and the classifier re-run ({stab["draws"]} seeded draws). (a) {case_text}; points ± 2σ,
thin lines, draws colored by class. (b) Share of draws per class. (c) Adjacent-pair
separation factors ± 2σ; {case_sf} (open circles) lie within 2σ of SF = 1.

## Figure 7 — fig7_trend_by_subtype

**Selectivity trend by ligand family.** Share of each subtype's ligands in each trend
class; numbers in bars, ligands; n, subtype total.

## Figure 8 — fig8_data_gaps

**Where more experiments would change the answer.** Per subtype (rows as in Fig. 3a):
ligands without a trend class, median lanthanides and conditions per ligand, ligands
whose class survives fewer than 90% of the draws of Fig. 6, and adjacent-pair SF values
within 2σ of 1. Each column is shaded on its own scale (darker, larger gap); cells
resting on fewer than {GAP_MIN_LIGANDS} ligands are white.

## Figure S1 — figS1_ci_width

Width of the 95% bootstrap interval of the slope against the number of conditions, for
the {len(xs)} ligands with a fitted slope (Spearman ρ = {rho:.2f}); x jittered.

## Figure S2 — figS2_cross_source

log D for the {len(pairs)} values that two compilations report for the same ligand, element and
condition (criteria in Section S5), colored by whether both cite the same paper; line,
y = x. Median |Δ log D| {np.median(diff) if len(diff) else float("nan"):.2f} dex.

## Figure S3 — figS3_sf_viewer

**The dashboard's separation-factor viewer**, captured from the running software
({capture["sf_viewer"]["captured"][:10]}). (a) The two periodic-table selectors, set to {sf_a} and {sf_b}.
(b) The plot they drive: log SF({sf_a}/{sf_b}) for each condition measuring both elements, by
ligand subtype, ordered by median. Violins, density; boxes, quartiles and median; dashed
diamonds, mean ± s.d.; filled circles, outliers; open triangles, values stored as bounds;
dashed line, no selectivity. The options between the two panels (x-axis grouping, solvent
filter) are not shown.
"""


def build(names=None, export=True, records=None, out=None):
    style()
    ctx = Context(records if records is not None else load_records())
    out = Path(out) if out else OUT
    made = {}
    for name, make in FIGURES.items():
        if names and name not in names:
            continue
        fig = make(ctx)
        if export:
            out.mkdir(parents=True, exist_ok=True)
            fig.savefig(out / f"{name}.pdf", metadata={"CreationDate": None, "ModDate": None})
            fig.savefig(out / f"{name}.png", dpi=300)
        made[name] = fig
        plt.close(fig)
    text = captions(ctx)  # always built, so CI exercises it
    if export and not names:
        (out / "CAPTIONS.md").write_text(text, encoding="utf-8")
    return made


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--only", action="append", choices=sorted(FIGURES))
    parser.add_argument("--no-export", action="store_true",
                        help="build the figures without writing files (for CI)")
    parser.add_argument("--db", type=Path, help=f"database to draw from (default {DB_PATH})")
    parser.add_argument("--out", type=Path, help=f"output directory (default {OUT})")
    args = parser.parse_args()
    if args.db and not args.out and not args.no_export:
        parser.error("--db draws a different snapshot; give --out so paper/figures/ keeps "
                     "the canonical set")
    made = build(args.only, export=not args.no_export, records=load_records(args.db),
                 out=args.out)
    for name in made:
        print(f"{'built' if args.no_export else 'wrote'} {name}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
