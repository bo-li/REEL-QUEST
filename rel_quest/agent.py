"""
REEL-QUEST AI Agent
Conversational agent that answers questions about the ligand database.

Two modes (chosen automatically):
  1. Tool-calling mode  – for models that support the `tools` parameter
     (qwen2.5:7b, llama3.1:8b, llama3.2:3b, …)
  2. Context-injection mode – fallback for small/cloud models
     Pre-runs relevant database queries and injects results into the prompt.

OLLAMA_BASE_URL env var overrides the default localhost address.
"""

import json
import os
import re
from collections import Counter
from typing import Any, Dict, Generator, List, Optional

import requests

from .core import LigandDatabase
from .schema import TREND_CLASSIFICATIONS, subtype_aliases
from .agent_tools import (
    TOOL_SCHEMAS, dispatch_tool, PLOT_TOOL_NAMES,
    get_database_statistics, get_best_separators,
    search_ligands_by_name, get_ligand_details,
    search_by_trend, search_by_mechanism,
    search_by_extractant_type, get_ligands_for_element,
    compare_ligands, LANTHANIDES,
)
from . import plotting
from .llm_backends import LLMBackend, OllamaBackend, OpenAICompatibleBackend

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

OLLAMA_BASE_URL = os.environ.get("OLLAMA_BASE_URL", "http://localhost:11434")
DEFAULT_MODEL = os.environ.get("RELQUEST_LLM_MODEL", "qwen3.5:0.8b")
MAX_TOOL_ROUNDS = 8
REQUEST_TIMEOUT = 120

# Models known to support tool-calling via Ollama
TOOL_CAPABLE_PATTERNS = [
    "qwen2.5", "llama3.1", "llama3.2", "llama3.3",
    "mistral-nemo", "command-r", "firefunction",
    "nemotron-mini", "granite3",
]

# The prompt below is a TEMPLATE: every count is filled in from the live
# database by build_system_prompt(). Do not paste literal numbers into it.
# The previous version hardcoded them and drifted — it told the model
# Phosphonate (4) and Phosphoric_acid (2) against a real 5 and 1, and quoted
# data-source counts summing to 217 for a 166-ligand database.
SYSTEM_PROMPT_TEMPLATE = """\
You are REEL-QUEST AI, an expert assistant for the REEL-QUEST database \
(Rare-Earth Extraction Ligands – Quantitative Understanding of Elemental Selectivity Trends). \
You help researchers explore rare earth element (REE) liquid-liquid extraction data.

## Database Overview
- {n_ligands} unique extractant ligands with experimental data
- Distribution ratios (D-values) across the lanthanide series, plus {extras}
- Selectivity trend classifications: {trends}
- Extraction mechanisms: {mechanisms}
- Classification (single-layer taxonomy): {n_subtypes} Subtypes. Both Chemical_Family \
and Subtype fields carry the same Subtype label on each record. \
Subtypes (with counts): {subtypes}.
- Data sources: {sources}

## Key Domain Concepts
- **D-value**: Distribution ratio = [M]_organic / [M]_aqueous. logD > 0 → good extraction.
- **Separation factor SF**: SF = D(M1)/D(M2). Higher SF = better selectivity.
- **Ascending trend**: D increases La→Lu (heavy REE preference, typical DGA ligands).
- **Reverse trend**: D increases Lu→La (light REE preference).
- **Lanthanide contraction**: ionic radius decreases La→Lu; adjacent pairs (Nd/Pr, Tb/Dy) \
are hardest to separate.
- **TODGA**: benchmark DGA extractant with ascending trend.
- **CMPO**: key organophosphorus extractant.

## Instructions
- Provide concise, scientifically accurate answers grounded in the data.
- When reporting D-values use logD (base-10 logarithm) unless raw D is requested.
- If a ligand is not in the database, say so clearly.
- Keep rankings in the order the data gives them, and do not add experimental \
conditions, protocols or version numbers that the data do not state.
"""


def build_system_prompt(db: LigandDatabase) -> str:
    """Fill the prompt template from the database that is actually loaded.

    Stating counts to the model is only useful if they are true: a stale
    figure is a confident wrong answer the model has no way to check. So
    everything countable is counted here rather than transcribed.
    """
    ligands = db.get_all_ligands()

    def tally(values):
        counts = Counter(v for v in values if v)
        return ", ".join(f"{name} ({n})"
                         for name, n in counts.most_common())

    lanthanides = set(LANTHANIDES)
    extra_elements = sorted(
        {e for lig in ligands for e in lig.tested_elements} - lanthanides)

    # Primary attribution, one source per ligand, so the counts sum to the
    # ligand total. A flat tally of LigandQuery.data_sources would not: 40
    # records list two contributing sources and would be counted twice,
    # summing to 206 for 166 ligands — the arithmetic that made the old
    # hardcoded line ("102 / 84 / 31") impossible to sanity-check.
    sources = tally(lig.metadata.get("data_source") for lig in ligands)
    multi_source = sum(1 for lig in ligands if len(lig.data_sources) > 1)
    if multi_source:
        sources += (f"; {multi_source} of these ligands additionally carry "
                    "measurements merged from a second source")
    if not sources:   # the public release carries no source labels
        sources = "every condition cites the paper that measured it (DOI and table/figure)"

    return SYSTEM_PROMPT_TEMPLATE.format(
        n_ligands=len(ligands),
        extras=(", ".join(extra_elements) if extra_elements
                else "no non-lanthanide elements"),
        trends=", ".join(TREND_CLASSIFICATIONS),
        mechanisms=tally(lig.mechanism for lig in ligands),
        n_subtypes=len({lig.subtype for lig in ligands if lig.subtype}),
        subtypes=tally(lig.subtype for lig in ligands),
        sources=sources,
    )


# ---------------------------------------------------------------------------
# Ollama module-level helpers (kept for backwards compat with app.py imports)
# ---------------------------------------------------------------------------

def check_ollama_connection() -> bool:
    return OllamaBackend(OLLAMA_BASE_URL).is_available()


def list_available_models() -> List[str]:
    return OllamaBackend(OLLAMA_BASE_URL).list_models()


def pull_model(model: str) -> bool:
    return OllamaBackend(OLLAMA_BASE_URL).pull_model(model)


def _supports_tools(model: str) -> bool:
    """Heuristic: return True if model name suggests tool-calling support."""
    m = model.lower()
    return any(p in m for p in TOOL_CAPABLE_PATTERNS)


def _describe_backend_error(exc: Exception) -> str:
    """
    Turn a provider exception into advice the reader can act on.

    Matched on the exception's class name and text rather than by importing
    openai, so this stays correct whether or not the SDK is installed.
    """
    name = type(exc).__name__
    text = str(exc)

    if "RateLimit" in name or "429" in text:
        return (
            "\n\n⚠️ **Rate limit reached.** The free models run on a shared daily "
            "quota for all visitors of this app. Wait a while and retry, or open "
            "⚙️ Settings and switch to **Your own API key** — free OpenRouter keys "
            "are available at https://openrouter.ai/keys"
        )
    if "Authentication" in name or "PermissionDenied" in name or "401" in text:
        return (
            "\n\n⚠️ **The API key was rejected.** Check the key and Base URL in "
            "⚙️ Settings — the key must belong to the provider that Base URL points at."
        )
    if "NotFound" in name or "404" in text:
        return (
            "\n\n⚠️ **Model not found.** Check the model ID in ⚙️ Settings; the "
            "provider may have retired it."
        )
    return f"\n\n⚠️ **Error**: {exc}"


# ---------------------------------------------------------------------------
# Context-injection fallback (for models without tool-calling)
# ---------------------------------------------------------------------------

def _extract_elements(text: str) -> List[str]:
    """Find lanthanide symbols mentioned in the text, in the order they are mentioned.

    Order matters: "separate Eu from La" asks for SF(Eu/La), so the first
    element named becomes element1. Sorting by atomic number instead turned it
    into La/Eu and reported every Eu-selective ligand with a negative log SF.
    """
    found = []
    for elem in LANTHANIDES:
        # Match whole-word element symbols
        match = re.search(r'\b' + elem + r'\b', text)
        if match:
            found.append((match.start(), elem))
    return [elem for _, elem in sorted(found)]


def _extract_ligand_names(db: LigandDatabase, text: str) -> List[str]:
    """Find ligand names from the database mentioned in the text."""
    text_lower = text.lower()
    found = []
    for lig in db.get_all_ligands():
        if lig.name.lower() in text_lower:
            found.append(lig.name)
    return found[:5]


def _build_context(db: LigandDatabase, user_message: str) -> str:
    """
    Pre-run relevant database queries based on question keywords and
    return a formatted context block to inject into the prompt.
    """
    context_parts = []
    msg_lower = user_message.lower()

    # Always include brief stats
    stats = get_database_statistics(db)
    context_parts.append(
        f"DATABASE STATS: {stats['total_ligands']} ligands | "
        f"Ligand families: {', '.join(list(stats['by_ligand_family'].keys())[:5])} | "
        f"Trends: {', '.join(list(stats['by_selectivity_trend'].keys())[:5])}"
    )

    elements = _extract_elements(user_message)
    ligand_names = _extract_ligand_names(db, user_message)

    # Best separators query
    if len(elements) >= 2 and any(
        kw in msg_lower for kw in ["separat", "best", "top", "highest", "factor"]
    ):
        result = get_best_separators(db, elements[0], elements[1], limit=10)
        context_parts.append(
            f"\nBEST SEPARATORS for {elements[0]}/{elements[1]} "
            f"(ranked by |log SF|, total {result['total_with_both_elements']} ligands):\n"
            + json.dumps(result["top_separators"], indent=2)
        )

    # Single element query
    elif len(elements) == 1 and any(
        kw in msg_lower for kw in ["ligand", "extract", "data", "available", "which"]
    ):
        result = get_ligands_for_element(db, elements[0], limit=10)
        context_parts.append(
            f"\nLIGANDS WITH DATA FOR {elements[0]} "
            f"(top {result['showing']} of {result['total_ligands_with_data']}):\n"
            + json.dumps(result["results"], indent=2)
        )

    # Specific ligand details
    if ligand_names:
        if len(ligand_names) >= 2 and any(
            kw in msg_lower for kw in ["compar", "vs", "versus", "difference"]
        ):
            result = compare_ligands(db, ligand_names)
            context_parts.append("\nLIGAND COMPARISON:\n" + json.dumps(result, indent=2))
        else:
            for name in ligand_names[:2]:
                result = get_ligand_details(db, name)
                context_parts.append(f"\nLIGAND DETAILS for {name}:\n" + json.dumps(result, indent=2))

    # Trend query. Iterate the shared vocabulary rather than a local list —
    # this one had been left at six values, so no phrasing could reach
    # Insufficient_Data.
    for trend in TREND_CLASSIFICATIONS:
        if trend.lower() in msg_lower:
            result = search_by_trend(db, trend, limit=8)
            context_parts.append(
                f"\nLIGANDS WITH {trend.upper()} TREND "
                f"({result['total_found']} total, showing {result['showing']}):\n"
                + json.dumps(result["results"], indent=2)
            )
            break

    # Mechanism query. Iterate the mechanisms the database actually holds,
    # so this cannot inject a confident "(0 total)" for a vocabulary value
    # no record uses.
    for mech in sorted({lig.mechanism for lig in db.get_all_ligands() if lig.mechanism}):
        if mech.lower().replace("_", " ") in msg_lower or mech.lower() in msg_lower:
            result = search_by_mechanism(db, mech, limit=8)
            context_parts.append(
                f"\nLIGANDS WITH {mech} MECHANISM "
                f"({result['total_found']} total):\n"
                + json.dumps(result["results"], indent=2)
            )
            break

    # Subtype query. The aliases come from vocab/chemical_families.json via
    # schema.subtype_aliases(), so they name v13 subtypes — which is what
    # search_by_extractant_type matches on. This used to be a hardcoded map of
    # v12 two-layer family labels ("Diglycolamide", "Organophosphorus", ...);
    # ten of its thirteen entries matched nothing, so a question about
    # phosphine oxides was answered with "(0 total)" against 19 records.
    #
    # Longest alias first, so "phenanthroline diamide" is not swallowed by a
    # shorter alias of a different subtype.
    matches = [
        (alias, subtype)
        for subtype, aliases in subtype_aliases().items()
        for alias in aliases
        if alias in msg_lower
    ]
    if matches:
        _, subtype = max(matches, key=lambda pair: len(pair[0]))
        result = search_by_extractant_type(db, subtype, limit=8)
        context_parts.append(
            f"\n{subtype.upper()} LIGANDS "
            f"({result['total_found']} total):\n"
            + json.dumps(result["results"], indent=2)
        )

    # General name search
    search_terms = re.findall(r'\b[A-Z][A-Z0-9\-]{2,}\b', user_message)
    if search_terms and not ligand_names:
        for term in search_terms[:2]:
            result = search_ligands_by_name(db, term, limit=5)
            if result and not (len(result) == 1 and "message" in result[0]):
                context_parts.append(
                    f"\nSEARCH RESULTS for '{term}':\n" + json.dumps(result, indent=2)
                )

    return "\n".join(context_parts)


# ---------------------------------------------------------------------------
# Agent class
# ---------------------------------------------------------------------------

class RELQuestAgentStreaming:
    """
    Conversational agent for REL-QUEST.

    Automatically switches between:
    - Tool-calling mode (larger models: qwen2.5, llama3.1, …)
    - Context-injection mode (small models: qwen3.5:0.8b, etc.)
    """

    def __init__(
        self,
        db: LigandDatabase,
        backend: Optional[LLMBackend] = None,
        model: str = DEFAULT_MODEL,
        ollama_url: str = OLLAMA_BASE_URL,
    ):
        self.db = db
        self.model = model
        self.ollama_url = ollama_url
        self.backend: LLMBackend = backend or OllamaBackend(ollama_url, model)
        self.history: List[Dict] = []
        self._use_tools: Optional[bool] = None  # auto-detect on first call
        self.last_figures: list = []  # Plotly figures from most recent chat()

    def reset(self):
        self.history = []
        self._use_tools = None
        self.last_figures = []

    def _build_messages(self, user_message: str, extra_context: str = "") -> List[Dict]:
        system = build_system_prompt(self.db)
        if extra_context:
            system += f"\n\n## Relevant Database Data\n{extra_context}"
        messages = [{"role": "system", "content": system}]
        messages.extend(self.history)
        messages.append({"role": "user", "content": user_message})
        return messages

    # ── Tool-calling path ────────────────────────────────────────────────────

    def _run_tool_loop(self, user_message: str) -> str:
        messages = self._build_messages(user_message)
        tool_rounds = 0

        while tool_rounds < MAX_TOOL_ROUNDS:
            response = self.backend.chat_sync(messages, tools=TOOL_SCHEMAS)
            msg = response.get("message", {})
            tool_calls = msg.get("tool_calls", [])

            if not tool_calls:
                return msg.get("content", "")

            messages.append({
                "role": "assistant",
                "content": msg.get("content", ""),
                "tool_calls": tool_calls,
            })

            for tc in tool_calls:
                fn = tc.get("function", {})
                tool_name = fn.get("name", "")
                raw_args = fn.get("arguments", {})
                if isinstance(raw_args, str):
                    try:
                        tool_args = json.loads(raw_args)
                    except json.JSONDecodeError:
                        tool_args = {}
                else:
                    tool_args = raw_args

                tool_result = dispatch_tool(self.db, tool_name, tool_args)

                # Handle plot tool results (dict with figure)
                if isinstance(tool_result, dict) and tool_result.get("type") == "figure":
                    fig = tool_result.get("figure")
                    if fig is not None:
                        self.last_figures.append(fig)
                    # Feed only the text summary back to the LLM
                    messages.append({"role": "tool", "content": tool_result.get("summary", "")})
                else:
                    messages.append({"role": "tool", "content": tool_result})

            tool_rounds += 1

        # Fallback: ask for a final answer without tools
        messages.append({
            "role": "user",
            "content": "Provide your final answer based on the tool results above.",
        })
        response = self.backend.chat_sync(messages, tools=None)
        return response.get("message", {}).get("content", "")

    # ── Context-injection path ───────────────────────────────────────────────

    def _run_context_injection(self, user_message: str) -> Generator[str, None, None]:
        context = _build_context(self.db, user_message)

        # Pre-generate plots for context-injection mode (no tool-calling)
        self._generate_plots_from_keywords(user_message)

        messages = self._build_messages(user_message, extra_context=context)
        yield from self.backend.chat_stream(messages)

    # ── Plot keyword detection (context-injection fallback) ──────────────────

    def _generate_plots_from_keywords(self, user_message: str):
        """Detect plot-related keywords and pre-generate relevant figures."""
        msg_lower = user_message.lower()
        plot_keywords = ["plot", "show", "visualize", "chart", "graph", "display", "compare"]
        if not any(kw in msg_lower for kw in plot_keywords):
            return

        elements = _extract_elements(user_message)
        ligand_names = _extract_ligand_names(self.db, user_message)

        # Comparison keywords + multiple ligands → comparison plot (check first)
        compare_kw = ["compare", "comparison", "overlay", "vs", "versus"]
        if len(ligand_names) >= 2 and any(kw in msg_lower for kw in compare_kw):
            fig = plotting.plot_ligand_comparison(self.db, ligand_names)
            if fig:
                self.last_figures.append(fig)
            return

        # Trend / selectivity keywords + specific ligand → selectivity trend plot
        trend_kw = ["trend", "selectivity", "logd", "log d", "log(d)"]
        if ligand_names and any(kw in msg_lower for kw in trend_kw):
            for name in ligand_names[:2]:
                fig = plotting.plot_selectivity_trend(self.db, name)
                if fig:
                    self.last_figures.append(fig)
            return

        # Separation factor keywords + two elements
        sf_kw = ["separation", "sf", "logsf", "factor"]
        if len(elements) >= 2 and any(kw in msg_lower for kw in sf_kw):
            group_by = "ligand_type"
            if "solvent" in msg_lower:
                group_by = "solvent_category" if "categor" in msg_lower else "solvent"
            fig = plotting.plot_separation_factor(self.db, elements[0], elements[1], group_by=group_by)
            if fig:
                self.last_figures.append(fig)
            return

        # Distribution keywords + single element
        dist_kw = ["distribution", "box", "violin"]
        if elements and any(kw in msg_lower for kw in dist_kw):
            group_by = "ligand_type"
            if "solvent" in msg_lower:
                group_by = "solvent_category" if "categor" in msg_lower else "solvent"
            elif "acid" in msg_lower:
                group_by = "ligand_acidity"
            fig = plotting.plot_element_distribution(self.db, elements[0], group_by=group_by)
            if fig:
                self.last_figures.append(fig)
            return

        # Single element mentioned with plot keyword → element distribution
        if elements:
            group_by = "ligand_type"
            if "solvent" in msg_lower:
                group_by = "solvent_category" if "categor" in msg_lower else "solvent"
            fig = plotting.plot_element_distribution(self.db, elements[0], group_by=group_by)
            if fig:
                self.last_figures.append(fig)
            return

        # Single ligand mentioned with plot keyword → selectivity trend
        if ligand_names:
            fig = plotting.plot_selectivity_trend(self.db, ligand_names[0])
            if fig:
                self.last_figures.append(fig)
            return

        # Overview keywords
        overview_kw = ["overview", "database", "summary", "pie", "donut"]
        if any(kw in msg_lower for kw in overview_kw):
            fig = plotting.plot_database_overview(self.db, "extractant_types")
            if fig:
                self.last_figures.append(fig)
            return

    # ── Main chat entry point ────────────────────────────────────────────────

    def chat(self, user_message: str) -> Generator[str, None, None]:
        self.last_figures = []  # Reset figures for this response
        try:
            # Auto-detect tool support on first call (or after reset).
            # Cloud backends (OpenAICompatibleBackend) always use context-injection.
            if self._use_tools is None:
                self._use_tools = (
                    isinstance(self.backend, OllamaBackend)
                    and _supports_tools(self.model)
                )

            if self._use_tools:
                try:
                    answer = self._run_tool_loop(user_message)
                    self.history.append({"role": "user", "content": user_message})
                    self.history.append({"role": "assistant", "content": answer})
                    yield answer
                    return
                except requests.exceptions.HTTPError as e:
                    if e.response is not None and e.response.status_code in (400, 500):
                        # Model doesn't support tools — fall back permanently
                        self._use_tools = False
                    else:
                        raise

            # Context-injection fallback
            full_answer = []
            for chunk in self._run_context_injection(user_message):
                full_answer.append(chunk)
                yield chunk

            self.history.append({"role": "user", "content": user_message})
            self.history.append({"role": "assistant", "content": "".join(full_answer)})

        except requests.exceptions.ConnectionError:
            if isinstance(self.backend, OllamaBackend):
                yield (
                    "\n\n⚠️ **Cannot connect to Ollama.** Make sure it is running:\n"
                    "```bash\nollama serve\n# or, with Docker:\ndocker compose up -d ollama\n```"
                )
            else:
                yield (
                    "\n\n⚠️ **Cannot reach the LLM provider.** Check your network "
                    "connection and the Base URL in ⚙️ Settings."
                )
        except requests.exceptions.Timeout:
            if isinstance(self.backend, OllamaBackend):
                yield "\n\n⚠️ **Ollama request timed out.** The model may still be loading — try again."
            else:
                yield "\n\n⚠️ **The LLM provider timed out.** Try again, or pick a smaller model."
        except requests.exceptions.HTTPError as e:
            yield f"\n\n⚠️ **LLM backend error**: {e}"
        except Exception as e:
            yield _describe_backend_error(e)

    def chat_sync(self, user_message: str) -> str:
        return "".join(self.chat(user_message))


# Keep old name as alias for compatibility
RELQuestAgent = RELQuestAgentStreaming
