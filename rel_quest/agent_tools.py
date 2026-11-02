"""
REEL-QUEST Agent Tools
Database query functions exposed as LLM-callable tools.
"""

import json
import math
from typing import Dict, List, Optional, Any

from .core import LigandDatabase
from .schema import TREND_CLASSIFICATIONS
from . import plotting
from .uncertainty import is_censored, is_mixture

# Lanthanide series in order
LANTHANIDES = ['La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd',
               'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu']


def _family(lig) -> str:
    """Canonical type label for a ligand — Layer-1 Chemical_Family with legacy fallback.

    Matches the dashboard convention so the AI assistant surfaces the same taxonomy
    as the rest of the UI (see vocab/chemical_families.json).
    """
    return getattr(lig, "chemical_family", None) or lig.extractant_type or "Unknown"

# Tool schemas for Ollama / OpenAI-compatible API
TOOL_SCHEMAS = [
    {
        "type": "function",
        "function": {
            "name": "get_database_statistics",
            "description": (
                "Get overall statistics about the REEL-QUEST database: "
                "total ligands, breakdown by ligand family (Layer-1), mechanism, "
                "selectivity trend, and data source."
            ),
            "parameters": {"type": "object", "properties": {}, "required": []},
        },
    },
    {
        "type": "function",
        "function": {
            "name": "search_ligands_by_name",
            "description": (
                "Search for ligands by name (partial match, case-insensitive). "
                "Returns matching ligand names with basic info (mechanism, trend, type)."
            ),
            "parameters": {
                "type": "object",
                "properties": {
                    "name": {
                        "type": "string",
                        "description": "Name or partial name to search (e.g. 'TODGA', 'DGA', 'CMPO')",
                    },
                    "limit": {
                        "type": "integer",
                        "description": "Max results to return (default 15)",
                    },
                },
                "required": ["name"],
            },
        },
    },
    {
        "type": "function",
        "function": {
            "name": "get_ligand_details",
            "description": (
                "Get detailed information about a specific ligand: SMILES, mechanism, "
                "extractant type, selectivity trend, D-values for all measured lanthanides, "
                "experimental conditions (acid, concentration, solvent), and data source."
            ),
            "parameters": {
                "type": "object",
                "properties": {
                    "name": {
                        "type": "string",
                        "description": "Ligand name (exact or partial)",
                    }
                },
                "required": ["name"],
            },
        },
    },
    {
        "type": "function",
        "function": {
            "name": "search_by_trend",
            "description": (
                "Find ligands with a specific selectivity trend across the lanthanide series. "
                "Ascending: La→Lu D-values increase (heavy REE preference). "
                "Reverse: Lu→La D-values increase (light REE preference). "
                "Flat: no clear trend. Mixed: conflicting trends across conditions. "
                "Peaked: D-values peak in the middle of the series. "
                "U_Shaped: D-values dip in the middle."
            ),
            "parameters": {
                "type": "object",
                "properties": {
                    "trend": {
                        "type": "string",
                        "enum": list(TREND_CLASSIFICATIONS),
                        "description": "Selectivity trend classification",
                    },
                    "limit": {
                        "type": "integer",
                        "description": "Max results to return (default 10)",
                    },
                },
                "required": ["trend"],
            },
        },
    },
    {
        "type": "function",
        "function": {
            "name": "search_by_mechanism",
            "description": "Find ligands by extraction mechanism.",
            "parameters": {
                "type": "object",
                "properties": {
                    "mechanism": {
                        "type": "string",
                        "enum": ["Solvation", "Cation_Exchange", "Anion_Exchange",
                                 "Ion_Pair", "Unknown"],
                        "description": "Extraction mechanism",
                    },
                    "limit": {
                        "type": "integer",
                        "description": "Max results to return (default 10)",
                    },
                },
                "required": ["mechanism"],
            },
        },
    },
    {
        "type": "function",
        "function": {
            "name": "search_by_extractant_type",
            "description": (
                "Find ligands by Subtype (single-layer taxonomy, 21 categories). "
                "See vocab/chemical_families.json for descriptions."
            ),
            "parameters": {
                "type": "object",
                "properties": {
                    "extractant_type": {
                        "type": "string",
                        "enum": [
                            "DGA", "Tripodal-DGA", "Polyether_Diamide", "Amide",
                            "Phthalamic_acid", "Malonamide", "Butanediamide",
                            "Crown_ether", "Polyphenol", "Nitrilotriacetamide",
                            "EDTA", "DAPhen", "MAPhen", "BLPhen",
                            "Pyridinediamide", "Bipydiamide", "Multident_N_donor",
                            "Phosphine_oxide", "Phosphonate", "Phosphate", "Phosphoric_acid",
                        ],
                        "description": "Ligand Subtype label",
                    },
                    "limit": {
                        "type": "integer",
                        "description": "Max results to return (default 10)",
                    },
                },
                "required": ["extractant_type"],
            },
        },
    },
    {
        "type": "function",
        "function": {
            "name": "get_best_separators",
            "description": (
                "Find the best ligands for separating two specific lanthanide elements, "
                "ranked by separation factor SF = D(elem1)/D(elem2). "
                "Returns ligands with the highest separation factors."
            ),
            "parameters": {
                "type": "object",
                "properties": {
                    "element1": {
                        "type": "string",
                        "description": "First element symbol (e.g. 'Eu', 'La', 'Nd')",
                    },
                    "element2": {
                        "type": "string",
                        "description": "Second element symbol (e.g. 'La', 'Pr', 'Sm')",
                    },
                    "limit": {
                        "type": "integer",
                        "description": "Max results to return (default 10)",
                    },
                },
                "required": ["element1", "element2"],
            },
        },
    },
    {
        "type": "function",
        "function": {
            "name": "compare_ligands",
            "description": (
                "Compare multiple ligands side by side: mechanism, extractant type, "
                "selectivity trend, median D-values for common elements, and best "
                "separation pairs."
            ),
            "parameters": {
                "type": "object",
                "properties": {
                    "names": {
                        "type": "array",
                        "items": {"type": "string"},
                        "description": "List of ligand names to compare (2–5 ligands)",
                    }
                },
                "required": ["names"],
            },
        },
    },
    {
        "type": "function",
        "function": {
            "name": "get_ligands_for_element",
            "description": (
                "Find ligands that have experimental distribution-ratio data "
                "for a specific lanthanide element."
            ),
            "parameters": {
                "type": "object",
                "properties": {
                    "element": {
                        "type": "string",
                        "description": "Element symbol (e.g. 'Eu', 'La', 'Nd', 'Tb')",
                    },
                    "limit": {
                        "type": "integer",
                        "description": "Max results to return (default 10)",
                    },
                },
                "required": ["element"],
            },
        },
    },
    # ── Plot tools ───────────────────────────────────────────────────────────
    {
        "type": "function",
        "function": {
            "name": "plot_selectivity_trend",
            "description": (
                "Generate a selectivity trend plot for a ligand: logD vs atomic number "
                "with linear fit line and confidence bands. Use when the user asks to "
                "'show', 'plot', or 'visualize' a ligand's selectivity trend."
            ),
            "parameters": {
                "type": "object",
                "properties": {
                    "ligand_name": {
                        "type": "string",
                        "description": "Ligand name (e.g. 'TODGA', 'CMPO')",
                    },
                },
                "required": ["ligand_name"],
            },
        },
    },
    {
        "type": "function",
        "function": {
            "name": "plot_element_distribution",
            "description": (
                "Generate a box+violin plot showing the distribution of logD values "
                "for a specific element, grouped by category. Use when the user asks "
                "to plot or visualize element distribution data."
            ),
            "parameters": {
                "type": "object",
                "properties": {
                    "element": {
                        "type": "string",
                        "description": "Element symbol (e.g. 'Eu', 'La', 'Nd')",
                    },
                    "group_by": {
                        "type": "string",
                        "enum": ["ligand_type", "solvent", "solvent_category", "ligand_acidity"],
                        "description": "Category to group by (default: ligand_type)",
                    },
                },
                "required": ["element"],
            },
        },
    },
    {
        "type": "function",
        "function": {
            "name": "plot_separation_factor",
            "description": (
                "Generate a box+violin plot of logSF (separation factor) between two "
                "elements, grouped by category. logSF = logD[elem1] - logD[elem2]. "
                "Use when the user asks to plot separation factors."
            ),
            "parameters": {
                "type": "object",
                "properties": {
                    "element1": {
                        "type": "string",
                        "description": "First element symbol (e.g. 'Eu')",
                    },
                    "element2": {
                        "type": "string",
                        "description": "Second element symbol (e.g. 'La')",
                    },
                    "group_by": {
                        "type": "string",
                        "enum": ["ligand_type", "solvent", "solvent_category", "ligand_acidity"],
                        "description": "Category to group by (default: ligand_type)",
                    },
                },
                "required": ["element1", "element2"],
            },
        },
    },
    {
        "type": "function",
        "function": {
            "name": "plot_database_overview",
            "description": (
                "Generate a pie/donut chart showing database composition. "
                "Use when the user asks to 'show' or 'visualize' the database overview."
            ),
            "parameters": {
                "type": "object",
                "properties": {
                    "chart_type": {
                        "type": "string",
                        "enum": ["extractant_types", "selectivity_trends", "element_coverage"],
                        "description": "Type of overview chart to generate",
                    },
                },
                "required": ["chart_type"],
            },
        },
    },
    {
        "type": "function",
        "function": {
            "name": "plot_ligand_comparison",
            "description": (
                "Generate an overlay plot comparing selectivity trends of multiple ligands. "
                "Use when the user asks to compare, overlay, or visualize multiple ligands together."
            ),
            "parameters": {
                "type": "object",
                "properties": {
                    "ligand_names": {
                        "type": "array",
                        "items": {"type": "string"},
                        "description": "List of ligand names to compare (2-8 ligands)",
                    },
                },
                "required": ["ligand_names"],
            },
        },
    },
]


# ---------------------------------------------------------------------------
# Tool implementations
# ---------------------------------------------------------------------------

def get_database_statistics(db: LigandDatabase) -> Dict[str, Any]:
    """Return high-level database statistics."""
    all_ligands = db.get_all_ligands()

    family_counts: Dict[str, int] = {}
    mechanism_counts: Dict[str, int] = {}
    trend_counts: Dict[str, int] = {}
    source_counts: Dict[str, int] = {}
    element_counts: Dict[str, int] = {}

    for lig in all_ligands:
        fam = _family(lig)
        family_counts[fam] = family_counts.get(fam, 0) + 1
        mechanism_counts[lig.mechanism] = mechanism_counts.get(lig.mechanism, 0) + 1
        trend_counts[lig.selectivity_trend] = trend_counts.get(lig.selectivity_trend, 0) + 1

        for src in getattr(lig, "data_sources", []):
            source_counts[src] = source_counts.get(src, 0) + 1

        # tested_elements: coverage counts only fitted La-Lu points,
        # so an element-frequency tally built on it omits Y and Sc.
        for elem in lig.tested_elements:
            element_counts[elem] = element_counts.get(elem, 0) + 1

    top_elements = sorted(element_counts.items(), key=lambda x: x[1], reverse=True)[:5]

    return {
        "total_ligands": len(all_ligands),
        "by_ligand_family": dict(sorted(family_counts.items(), key=lambda x: x[1], reverse=True)),
        "by_mechanism": dict(sorted(mechanism_counts.items(), key=lambda x: x[1], reverse=True)),
        "by_selectivity_trend": dict(sorted(trend_counts.items(), key=lambda x: x[1], reverse=True)),
        **({"by_data_source": dict(sorted(source_counts.items(), key=lambda x: x[1], reverse=True))}
           if source_counts else {}),
        "most_measured_elements": [f"{e}({c})" for e, c in top_elements],
        "lanthanide_coverage": "All 15 (La–Lu including Pm)",
    }


def _sources(lig, key: str) -> Dict:
    """{key: the record's source labels}, or {} where it has none (the public release)."""
    sources = getattr(lig, "data_sources", None) or []
    return {key: sources} if sources else {}


def search_ligands_by_name(db: LigandDatabase, name: str, limit: int = 15) -> List[Dict]:
    """Fuzzy name search returning basic ligand info."""
    name_lower = name.lower()
    results = []
    for lig in db.get_all_ligands():
        if name_lower in lig.name.lower():
            results.append({
                "name": lig.name,
                "ligand_family": _family(lig),
                "mechanism": lig.mechanism,
                "selectivity_trend": lig.selectivity_trend,
                "elements_measured": len(lig.tested_elements),
                **_sources(lig, "data_source"),
            })
    results = results[:limit]
    return results if results else [{"message": f"No ligands found matching '{name}'"}]


def _rank_matches(all_ligands, name_lower):
    """Return ligands matching name_lower, exact matches first."""
    exact = [l for l in all_ligands if l.name.lower() == name_lower]
    partial = [l for l in all_ligands if name_lower in l.name.lower() and l.name.lower() != name_lower]
    return exact + partial


def _condition_parameters(cond: Dict) -> Dict[str, Any]:
    """The experimental parameters of one condition, under the stored keys."""
    return {
        "condition_id": cond.get("condition_id"),
        "acid": cond.get("Acid_Name"),
        "acid_conc_M": cond.get("Acid_Concentration_M"),
        "extractant_conc_M": cond.get("Extractant_Concentration_M"),
        "solvent": cond.get("Solvent_Name"),
        "temperature_C": cond.get("Temperature_C"),
    }


def get_ligand_details(db: LigandDatabase, name: str) -> Dict[str, Any]:
    """Full details for a single ligand."""
    name_lower = name.lower()
    matches = _rank_matches(db.get_all_ligands(), name_lower)
    if not matches:
        return {"error": f"No ligand found matching '{name}'"}

    lig = matches[0]  # best match (first/exact)

    # Aggregate D-values across conditions
    agg = lig.get_aggregated_d_values(aggregation="median", lanthanides_only=True)

    d_values_summary = {}
    for elem in LANTHANIDES:
        if elem in agg:
            d_values_summary[elem] = {
                "median_logD": round(agg[elem]["log_value"], 3) if agg[elem]["log_value"] is not None else None,
                "n_conditions": agg[elem]["n_conditions"],
            }
            if agg[elem].get("n_bounds"):   # bounds (D > x, D < x) are not in the median
                d_values_summary[elem]["n_bounds_not_in_median"] = agg[elem]["n_bounds"]

    # Conditions summary (first 3)
    cond_summaries = []
    for cond in lig.experimental_conditions[:3]:
        cs = _condition_parameters(cond)
        cs["elements"] = list(cond.get("distribution_ratios", {}).keys())
        cs["DOI"] = (cond.get("reference") or {}).get("DOI")
        cond_summaries.append(cs)

    return {
        "name": lig.name,
        "smiles": lig.smiles,
        "ligand_family": _family(lig),
        "subtype": getattr(lig, "subtype", None),
        "mechanism": lig.mechanism,
        **({"structural_class": lig.structural_class} if lig.structural_class else {}),
        "selectivity_trend": lig.selectivity_trend,
        "trend_confidence": lig.trend_confidence,
        "elements_measured": lig.tested_elements,
        "median_logD_by_element": d_values_summary,
        "total_conditions": len(lig.experimental_conditions),
        "sample_conditions": cond_summaries,
        **_sources(lig, "data_sources"),
        "other_matches": [l.name for l in matches[1:5]] if len(matches) > 1 else [],
    }


def search_by_trend(db: LigandDatabase, trend: str, limit: int = 10) -> List[Dict]:
    """List ligands with a given selectivity trend."""
    results = []
    for lig in db.get_all_ligands():
        if lig.selectivity_trend == trend:
            results.append({
                "name": lig.name,
                "ligand_family": _family(lig),
                "mechanism": lig.mechanism,
                "elements_measured": len(lig.tested_elements),
                **_sources(lig, "data_source"),
            })
    total = len(results)
    return {
        "trend": trend,
        "total_found": total,
        "results": results[:limit],
        "showing": min(limit, total),
    }


def search_by_mechanism(db: LigandDatabase, mechanism: str, limit: int = 10) -> List[Dict]:
    """List ligands by extraction mechanism."""
    results = []
    for lig in db.get_all_ligands():
        if lig.mechanism == mechanism:
            results.append({
                "name": lig.name,
                "ligand_family": _family(lig),
                "selectivity_trend": lig.selectivity_trend,
                "elements_measured": len(lig.tested_elements),
            })
    total = len(results)
    return {
        "mechanism": mechanism,
        "total_found": total,
        "results": results[:limit],
        "showing": min(limit, total),
    }


def search_by_extractant_type(db: LigandDatabase, extractant_type: str, limit: int = 10) -> List[Dict]:
    """List ligands by Layer-1 Chemical_Family (legacy parameter name retained for tool-schema stability)."""
    results = []
    for lig in db.get_all_ligands():
        if _family(lig) == extractant_type:
            results.append({
                "name": lig.name,
                "subtype": getattr(lig, "subtype", None),
                "mechanism": lig.mechanism,
                "selectivity_trend": lig.selectivity_trend,
                "elements_measured": len(lig.tested_elements),
            })
    total = len(results)
    return {
        "ligand_family": extractant_type,
        "total_found": total,
        "results": results[:limit],
        "showing": min(limit, total),
    }


def get_best_separators(db: LigandDatabase, element1: str, element2: str, limit: int = 10) -> Dict:
    """Rank ligands by separation factor for two elements.

    log SF is taken within each condition that measured both elements
    (plotting.extract_separation_factor_data, the logSF viewer's data), then
    summarised per ligand by its median. Dividing per-element medians instead
    would pair a D from one paper with a D from another. Bounds and
    mixed-extractant conditions stay out of the median.
    """
    elem1 = element1.capitalize()
    elem2 = element2.capitalize()

    pairs = plotting.extract_separation_factor_data(db, elem1, elem2)
    ranked = []
    if not pairs.empty:
        for name, group in pairs.groupby("ligand_name", sort=False):
            values = group[group["censor"].isna()]
            if values.empty:
                continue
            median = float(values["log_sf"].median())
            ranked.append({
                "name": name,
                "ligand_family": group["chemical_family"].iloc[0],
                "median_log_SF": round(median, 3),
                "SF": round(10 ** median, 3),
                "min_log_SF": round(float(values["log_sf"].min()), 3),
                "max_log_SF": round(float(values["log_sf"].max()), 3),
                "n_paired_conditions": len(values),
                # pairs where log SF is only a bound; left out of the median
                "n_bound_pairs": len(group) - len(values),
                "prefers": elem1 if median >= 0 else elem2,
            })

    ranked.sort(key=lambda x: abs(x["median_log_SF"]), reverse=True)
    total = len(ranked)

    return {
        "element1": elem1,
        "element2": elem2,
        "total_with_both_elements": total,
        "top_separators": ranked[:limit],
        "note": (
            f"log SF = log D({elem1}) - log D({elem2}) within one condition; the median "
            f"across a ligand's paired conditions is shown. SF > 1 means {elem1} is "
            f"preferentially extracted, SF < 1 means {elem2}. Ranked by |log SF|."
        ),
    }


def get_ligand_conditions(db: LigandDatabase, name: str, element: Optional[str] = None,
                          limit: int = 50) -> Dict[str, Any]:
    """Every experimental condition of one ligand, with its D values and source."""
    matches = _rank_matches(db.get_all_ligands(), name.lower())
    if not matches:
        return {"error": f"No ligand found matching '{name}'"}
    lig = matches[0]

    conditions = lig.experimental_conditions
    if element:
        elem = element.capitalize()
        conditions = [c for c in conditions if elem in c.get("distribution_ratios", {})]

    rows = []
    for cond in conditions[:limit]:
        row = _condition_parameters(cond)
        row.update({
            "phase_modifier": cond.get("Phase_Modifier_Name"),
            "phase_modifier_conc_M": cond.get("Phase_Modifier_Concentration_M"),
            "additives": cond.get("Additives"),
            "co_extractants": cond.get("Co_Extractants"),
            "is_mixture": is_mixture(cond),
            "D": {
                elem: {
                    "D": cell.get("value"),
                    "log_D": cell.get("log_value"),
                    "censor": cell["censor"] if is_censored(cell) else None,
                    "sigma_log_D": cell.get("sigma_log_D"),
                }
                for elem, cell in cond.get("distribution_ratios", {}).items()
            },
        })
        ref = cond.get("reference") or {}
        row["reference"] = {
            "DOI": ref.get("DOI"),
            "URL": ref.get("URL"),
            "Locator": ref.get("Locator"),
        }
        if "Provenance_Level" in ref:   # absent in the public release
            row["reference"]["Provenance_Level"] = ref["Provenance_Level"]
            row["reference"]["Aggregator"] = (ref.get("Aggregator") or {}).get("name")
        rows.append(row)

    return {
        "name": lig.name,
        "smiles": lig.smiles,
        "total_conditions": len(conditions),
        "showing": len(rows),
        "conditions": rows,
        "other_matches": [l.name for l in matches[1:5]],
        "note": (
            "censor '>' or '<' marks a bound (D > x, D < x), not a measured value. "
            "is_mixture marks a condition measured with other extractants present. "
            "Cite reference.DOI (the paper that measured the value)."
        ),
    }


def compare_ligands(db: LigandDatabase, names: List[str]) -> Dict:
    """Side-by-side comparison of multiple ligands."""
    comparison = []
    not_found = []

    for name in names:
        name_lower = name.lower()
        matches = _rank_matches(db.get_all_ligands(), name_lower)
        if not matches:
            not_found.append(name)
            continue

        lig = matches[0]
        agg = lig.get_aggregated_d_values(aggregation="median", lanthanides_only=True)

        # Summarize median logD for key elements
        key_elements = ["La", "Nd", "Eu", "Gd", "Ho", "Er", "Lu"]
        logd_summary = {}
        for elem in key_elements:
            if elem in agg and agg[elem]["log_value"] is not None:
                logd_summary[elem] = round(agg[elem]["log_value"], 2)

        comparison.append({
            "name": lig.name,
            "ligand_family": _family(lig),
            "mechanism": lig.mechanism,
            "selectivity_trend": lig.selectivity_trend,
            "trend_confidence": lig.trend_confidence,
            "elements_measured": len(lig.tested_elements),
            "median_logD_key_elements": logd_summary,
            **_sources(lig, "data_source"),
        })

    return {
        "comparison": comparison,
        "not_found": not_found,
        "note": "logD values are medians across all experimental conditions.",
    }


def get_ligands_for_element(db: LigandDatabase, element: str, limit: int = 10) -> Dict:
    """Find ligands with data for a specific element."""
    elem = element.capitalize()
    results = []
    for lig in db.get_all_ligands():
        # tested_elements so Y/Sc questions return their real answers.
        if elem in lig.tested_elements:
            # Get median logD for this element
            agg = lig.get_aggregated_d_values(aggregation="median")
            logd = None
            if elem in agg and agg[elem]["log_value"] is not None:
                logd = round(agg[elem]["log_value"], 2)

            results.append({
                "name": lig.name,
                "ligand_family": _family(lig),
                "selectivity_trend": lig.selectivity_trend,
                f"median_logD_{elem}": logd,
            })

    # Sort by logD descending (best extractors first)
    results.sort(
        key=lambda x: x.get(f"median_logD_{elem}") or -999,
        reverse=True,
    )
    total = len(results)
    return {
        "element": elem,
        "total_ligands_with_data": total,
        "results": results[:limit],
        "showing": min(limit, total),
    }


# ---------------------------------------------------------------------------
# Plot tool implementations
# ---------------------------------------------------------------------------

PLOT_TOOL_NAMES = frozenset([
    "plot_selectivity_trend", "plot_element_distribution",
    "plot_separation_factor", "plot_database_overview",
    "plot_ligand_comparison",
])


def _dispatch_plot_tool(db: LigandDatabase, tool_name: str, tool_args: Dict) -> Optional[Dict]:
    """
    Execute a plot tool and return a dict with keys:
      - "type": "figure"
      - "figure": plotly Figure object
      - "summary": text description for the LLM

    Returns None if the tool_name is not a plot tool.
    """
    if tool_name not in PLOT_TOOL_NAMES:
        return None

    fig = None
    summary = ""

    if tool_name == "plot_selectivity_trend":
        ligand_name = tool_args.get("ligand_name", "")
        fig = plotting.plot_selectivity_trend(db, ligand_name)
        summary = plotting.get_selectivity_trend_summary(db, ligand_name)

    elif tool_name == "plot_element_distribution":
        element = tool_args.get("element", "Eu")
        group_by = tool_args.get("group_by", "ligand_type")
        fig = plotting.plot_element_distribution(db, element, group_by=group_by)
        df = plotting.extract_element_data(db, element)
        values, bounds = plotting.split_bounds(df)
        n_lig = df['ligand_name'].nunique() if len(df) > 0 else 0
        summary = (
            f"Generated box+violin plot of logD[{element}] grouped by {group_by}. "
            f"{len(values)} measured values across {n_lig} ligands"
            + (f", plus {len(bounds)} bounds (D > x or D < x, drawn as triangles, not in the statistics)."
               if len(bounds) else ".")
        )

    elif tool_name == "plot_separation_factor":
        e1 = tool_args.get("element1", "Eu")
        e2 = tool_args.get("element2", "La")
        group_by = tool_args.get("group_by", "ligand_type")
        fig = plotting.plot_separation_factor(db, e1, e2, group_by=group_by)
        df = plotting.extract_separation_factor_data(db, e1, e2)
        values, bounds = plotting.split_bounds(df)
        n = len(values)
        median_sf = round(values['log_sf'].median(), 3) if n > 0 else 0
        summary = (
            f"Generated logSF[{e1}/{e2}] distribution grouped by {group_by}. "
            f"{n} paired measurements, median logSF = {median_sf}"
            + (f"; {len(bounds)} more pairs give only a bound on logSF (not in the median)."
               if len(bounds) else ".")
        )

    elif tool_name == "plot_database_overview":
        chart_type = tool_args.get("chart_type", "extractant_types")
        fig = plotting.plot_database_overview(db, chart_type)
        summary = f"Generated {chart_type} overview chart for the database."

    elif tool_name == "plot_ligand_comparison":
        names = tool_args.get("ligand_names", [])
        fig = plotting.plot_ligand_comparison(db, names)
        summary = f"Generated comparison plot for {', '.join(names[:5])}."

    if fig is None:
        return {"type": "figure", "figure": None,
                "summary": f"Could not generate plot: no data found for the given parameters."}

    return {"type": "figure", "figure": fig, "summary": summary}


# ---------------------------------------------------------------------------
# Dispatcher: call a tool by name
# ---------------------------------------------------------------------------

def dispatch_tool(db: LigandDatabase, tool_name: str, tool_args: Dict):
    """
    Execute a tool by name.

    For data-query tools: returns a JSON string.
    For plot tools: returns a dict with "type": "figure", "figure", and "summary".
    """
    # Check plot tools first
    plot_result = _dispatch_plot_tool(db, tool_name, tool_args)
    if plot_result is not None:
        return plot_result

    try:
        if tool_name == "get_database_statistics":
            result = get_database_statistics(db)
        elif tool_name == "search_ligands_by_name":
            result = search_ligands_by_name(db, **tool_args)
        elif tool_name == "get_ligand_details":
            result = get_ligand_details(db, **tool_args)
        elif tool_name == "search_by_trend":
            result = search_by_trend(db, **tool_args)
        elif tool_name == "search_by_mechanism":
            result = search_by_mechanism(db, **tool_args)
        elif tool_name == "search_by_extractant_type":
            result = search_by_extractant_type(db, **tool_args)
        elif tool_name == "get_best_separators":
            result = get_best_separators(db, **tool_args)
        elif tool_name == "compare_ligands":
            result = compare_ligands(db, **tool_args)
        elif tool_name == "get_ligands_for_element":
            result = get_ligands_for_element(db, **tool_args)
        else:
            result = {"error": f"Unknown tool: {tool_name}"}
    except Exception as exc:
        result = {"error": f"Tool '{tool_name}' raised an exception: {exc}"}

    return json.dumps(result, ensure_ascii=False)
