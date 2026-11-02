"""REST API for the REEL-QUEST ligand database.

A read-only FastAPI layer over LigandDatabase, so the database is usable
without cloning the repo. The dashboard (app.py) remains the interactive UI;
this exists for programmatic access and deep links.

Run locally:
    uvicorn rel_quest.api:app --reload

Endpoints (all GET, JSON):
    /api/v1/health                 liveness + record counts
    /api/v1/stats                  database-level statistics
    /api/v1/vocab                  filterable vocabularies with counts
    /api/v1/ligands                filtered, paginated ligand summaries
    /api/v1/ligands/{name}         full record (exact, alternative, or fuzzy name)

Interactive docs: /api/docs (Swagger) and /api/redoc (ReDoc); the OpenAPI
spec is served at /api/v1/openapi.json.

Design notes:
  - The element filter reads condition-level distribution_ratios directly
    instead of selectivity.element_coverage, which is lanthanide-only and
    stale on three merged records (PROJECT_REVIEW_2026-07-31.md, finding 2).
    This makes Y and Sc queryable, unlike core.filter_by_elements().
  - The trend filter accepts the full 7-value vocabulary including Peaked
    and U_Shaped (which scripts/query_ligands.py still rejects).
  - Enum-ish query params are validated against the live database, not a
    hardcoded list, so a rebuild cannot silently desynchronize them.
"""

from __future__ import annotations

from typing import Any, Dict, List, Optional

from fastapi import FastAPI, HTTPException, Query, Request
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import RedirectResponse
from pydantic import BaseModel, Field

from . import __version__
from .core import LigandDatabase, LigandQuery
from .plotting import ATOMIC_NUMBERS
from .uncertainty import is_censored

API_PREFIX = "/api/v1"  # the public prefix, standalone or mounted (create_app)


# ---------------------------------------------------------------------------
# Response models
# ---------------------------------------------------------------------------

class LigandSummary(BaseModel):
    """The list-view projection of a ligand record."""
    name: str
    smiles: Optional[str] = None
    iupac_name: Optional[str] = None
    community_name: Optional[str] = None
    subtype: Optional[str] = None
    chemical_family: Optional[str] = None
    mechanism: Optional[str] = None
    selectivity_trend: Optional[str] = None
    trend_confidence: Optional[str] = None
    data_sources: List[str] = Field(default_factory=list)
    n_conditions: int = 0
    elements_measured: List[str] = Field(
        default_factory=list,
        description="Derived from condition-level measurements (periodic order), "
                    "not the stale element_coverage summary field.",
    )
    cas_number: Optional[str] = None
    pubchem_cid: Optional[int] = None


class LigandDetail(LigandSummary):
    """The full record, plus how the requested name was resolved."""
    matched_by: str = Field(
        description="'name' (exact), 'alternative_name', or 'fuzzy_name'")
    alternative_names: List[str] = Field(default_factory=list)
    inchi: Optional[str] = None
    inchi_key: Optional[str] = None
    selectivity: Optional[Dict[str, Any]] = None
    performance_summary: Optional[Dict[str, Any]] = None
    data_quality: Optional[Dict[str, Any]] = None
    metadata: Optional[Dict[str, Any]] = None
    experimental_data: List[Dict[str, Any]] = Field(default_factory=list)


class LigandListResponse(BaseModel):
    total: int = Field(description="Ligands matching the filters, before paging")
    offset: int
    limit: int
    count: int = Field(description="Ligands in this page")
    results: List[LigandSummary]


# ---------------------------------------------------------------------------
# Projection helpers
# ---------------------------------------------------------------------------

# plotting.ATOMIC_NUMBERS covers the lanthanides only; Y and Sc are real
# database elements too. Kept local — plot code iterates the shared dict.
_NON_LANTHANIDE_Z = {"Sc": 21, "Y": 39}


def _element_sort_key(symbol: str):
    z = ATOMIC_NUMBERS.get(symbol) or _NON_LANTHANIDE_Z.get(symbol, 999)
    return (z, symbol)


def elements_measured(ligand: LigandQuery) -> List[str]:
    """Every element with an actual measurement, in periodic order.

    Reads LigandQuery.tested_elements (the stored `elements_measured`
    field), not selectivity.element_coverage — see the module docstring.
    """
    return sorted(ligand.tested_elements, key=_element_sort_key)


def _summary(ligand: LigandQuery) -> Dict[str, Any]:
    return {
        "name": ligand.name,
        "smiles": ligand.smiles,
        "iupac_name": ligand.iupac_name,
        "community_name": ligand.community_name,
        "subtype": ligand.subtype,
        "chemical_family": ligand.chemical_family,
        "mechanism": ligand.mechanism,
        "selectivity_trend": ligand.selectivity_trend,
        "trend_confidence": ligand.trend_confidence,
        "data_sources": ligand.data_sources,
        "n_conditions": len(ligand.experimental_conditions),
        "elements_measured": elements_measured(ligand),
        "cas_number": ligand.cas_number,
        "pubchem_cid": ligand.pubchem_cid,
    }


def _detail(ligand: LigandQuery, matched_by: str) -> Dict[str, Any]:
    payload = _summary(ligand)
    payload.update({
        "matched_by": matched_by,
        "alternative_names": ligand.alternative_names or [],
        "inchi": ligand.inchi,
        "inchi_key": ligand.inchi_key,
        "selectivity": ligand.selectivity,
        "performance_summary": ligand.performance_summary,
        "data_quality": ligand.data_quality,
        "metadata": ligand.metadata,
        "experimental_data": ligand.experimental_conditions,
    })
    return payload


# ---------------------------------------------------------------------------
# Filtering (mirrors app.py::_filter_ligands semantics)
# ---------------------------------------------------------------------------

def apply_filters(
    ligands: List[LigandQuery],
    q: Optional[str] = None,
    subtype: Optional[str] = None,
    mechanism: Optional[str] = None,
    trend: Optional[str] = None,
    source: Optional[str] = None,
    elements: Optional[List[str]] = None,
    logd_min: Optional[float] = None,
    logd_max: Optional[float] = None,
) -> List[LigandQuery]:
    """Apply every active filter; None means 'not filtered'."""
    if q:
        needle = q.lower()
        ligands = [l for l in ligands if needle in l.name.lower()]
    if subtype:
        ligands = [l for l in ligands if l.subtype == subtype]
    if mechanism:
        ligands = [l for l in ligands if l.mechanism == mechanism]
    if trend:
        ligands = [l for l in ligands if l.selectivity_trend == trend]
    if source:
        ligands = [l for l in ligands if source in l.data_sources]
    if elements:
        required = set(elements)
        ligands = [l for l in ligands
                   if required.issubset(elements_measured(l))]
    if logd_min is not None or logd_max is not None:
        lo = logd_min if logd_min is not None else float("-inf")
        hi = logd_max if logd_max is not None else float("inf")

        def in_range(ligand: LigandQuery) -> bool:
            for condition in ligand.experimental_conditions:
                for data in (condition.get("distribution_ratios") or {}).values():
                    log_value = data.get("log_value") if isinstance(data, dict) else None
                    # a bound (D > x, D < x) is not a value in range
                    if log_value is not None and not is_censored(data) and lo <= log_value <= hi:
                        return True
            return False

        ligands = [l for l in ligands if in_range(l)]
    return ligands


# ---------------------------------------------------------------------------
# Application factory
# ---------------------------------------------------------------------------

def create_app(db_path: Optional[str] = None, db: Optional[LigandDatabase] = None,
               mount: str = "/api") -> FastAPI:
    """The API app. `db` shares an already loaded database (the Space's);
    `mount` is the path the API lives under in this app: "/api" when it runs
    on its own, "" when another app mounts it at /api (space/app.py), so the
    public paths are /api/v1/... either way."""
    db = db if db is not None else LigandDatabase(db_path)
    prefix = f"{mount}/v1"
    all_ligands = sorted(db.get_all_ligands(), key=lambda l: l.name.lower())

    # Live vocabularies for query-param validation and /vocab.
    def counted(values) -> Dict[str, int]:
        counts: Dict[str, int] = {}
        for value in values:
            if value:
                counts[value] = counts.get(value, 0) + 1
        return dict(sorted(counts.items(), key=lambda kv: (-kv[1], kv[0])))

    vocab = {
        "subtypes": counted(l.subtype for l in all_ligands),
        "mechanisms": counted(l.mechanism for l in all_ligands),
        "trends": counted(l.selectivity_trend for l in all_ligands),
        "sources": counted(s for l in all_ligands for s in l.data_sources),
    }
    known_elements = sorted(
        {e for l in all_ligands for e in elements_measured(l)},
        key=_element_sort_key,
    )
    total_measurements = sum(len(l.experimental_conditions) for l in all_ligands)

    app = FastAPI(
        title="REEL-QUEST API",
        version=__version__,
        description=(
            "Read-only REST access to the REEL-QUEST rare-earth-element "
            "ligand database: curated solvent-extraction distribution ratios "
            "with mechanism classification and selectivity-trend analysis."
        ),
        openapi_url=f"{prefix}/openapi.json",
        docs_url=f"{mount}/docs",
        redoc_url=f"{mount}/redoc",
    )
    # Read-only public data; let browser-based tools call it cross-origin.
    app.add_middleware(
        CORSMiddleware,
        allow_origins=["*"], allow_methods=["GET"], allow_headers=["*"],
    )
    app.state.db = db

    def validated(param: str, value: Optional[str], allowed: Dict[str, int]) -> Optional[str]:
        if value is not None and value not in allowed:
            raise HTTPException(status_code=422, detail={
                "param": param,
                "given": value,
                "allowed": sorted(allowed),
            })
        return value

    @app.get("/", include_in_schema=False)
    def root(request: Request) -> RedirectResponse:
        # root_path is the mount point another app put this one under
        return RedirectResponse(url=f"{request.scope.get('root_path', '')}{mount}/docs")

    @app.get(f"{prefix}/health")
    def health() -> Dict[str, Any]:
        return {
            "status": "ok",
            "version": __version__,
            "ligands": len(all_ligands),
            "measurements": total_measurements,
        }

    @app.get(f"{prefix}/stats")
    def stats() -> Dict[str, Any]:
        return db.get_statistics()

    @app.get(f"{prefix}/vocab")
    def vocabularies() -> Dict[str, Any]:
        """Every filterable value with its ligand count."""
        return {**vocab, "elements": known_elements}

    @app.get(f"{prefix}/ligands", response_model=LigandListResponse)
    def list_ligands(
        q: Optional[str] = Query(None, description="Case-insensitive name substring"),
        subtype: Optional[str] = Query(None, description="Subtype, e.g. DGA"),
        mechanism: Optional[str] = Query(None, description="e.g. Solvation"),
        trend: Optional[str] = Query(
            None, description="Selectivity trend (7-value vocabulary, "
                              "including Peaked and U_Shaped)"),
        source: Optional[str] = Query(None, description="Source compilation, where records "
                                                        "carry one (listed in /vocab)"),
        elements: Optional[str] = Query(
            None, description="Comma-separated symbols the ligand must have "
                              "measurements for, e.g. 'Eu,La'. Y and Sc work."),
        logd_min: Optional[float] = Query(None, description="Keep ligands with at "
                                          "least one log D >= this"),
        logd_max: Optional[float] = Query(None, description="Keep ligands with at "
                                          "least one log D <= this"),
        limit: int = Query(50, ge=1, le=200),
        offset: int = Query(0, ge=0),
    ) -> Dict[str, Any]:
        validated("subtype", subtype, vocab["subtypes"])
        validated("mechanism", mechanism, vocab["mechanisms"])
        validated("trend", trend, vocab["trends"])
        validated("source", source, vocab["sources"])

        required_elements = None
        if elements:
            required_elements = [e.strip().capitalize()
                                 for e in elements.split(",") if e.strip()]
            unknown = [e for e in required_elements if e not in known_elements]
            if unknown:
                raise HTTPException(status_code=422, detail={
                    "param": "elements",
                    "given": unknown,
                    "allowed": known_elements,
                })

        matches = apply_filters(
            all_ligands, q=q, subtype=subtype, mechanism=mechanism, trend=trend,
            source=source, elements=required_elements,
            logd_min=logd_min, logd_max=logd_max,
        )
        page = matches[offset:offset + limit]
        return {
            "total": len(matches),
            "offset": offset,
            "limit": limit,
            "count": len(page),
            "results": [_summary(l) for l in page],
        }

    @app.get(prefix + "/ligands/{name}", response_model=LigandDetail)
    def ligand_detail(name: str) -> Dict[str, Any]:
        # 1. Exact primary name (Name is unique across the database).
        exact = db.query_by_name(name)
        if exact:
            return _detail(exact[0], matched_by="name")

        # 2. Exact alternative name — renamed/merged records keep their old
        #    names there, so pre-2026-07-28 references still resolve.
        wanted = name.lower()
        alternates = [l for l in all_ligands
                      if any(a.lower() == wanted for a in (l.alternative_names or []))]
        if len(alternates) == 1:
            return _detail(alternates[0], matched_by="alternative_name")

        # 3. Unambiguous substring match as a convenience.
        fuzzy = db.query_by_name(name, fuzzy=True)
        if len(fuzzy) == 1:
            return _detail(fuzzy[0], matched_by="fuzzy_name")

        suggestions = sorted({l.name for l in (alternates + fuzzy)})[:10]
        raise HTTPException(status_code=404, detail={
            "error": f"No ligand named {name!r}",
            "suggestions": suggestions,
        })

    return app


app = create_app()
