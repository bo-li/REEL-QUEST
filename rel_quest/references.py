"""Reference provenance for experimental conditions (schema v0.3).

Every condition says where its numbers came from. Before v0.3 that was one
free-text `reference.DOI` string, and it did not parse as a DOI on a third of
the database: 214 SAFE-derived conditions stored two URLs joined by a comma,

    "https://doi.org/10.1080/07366299.2015.1087209 , https://doi.org/10.1021/jacs.5c19738"

and the rest used six different prefix styles (`doi:`, `DOI: `, `dx.doi.org/`,
`https://doi.org/`, bare). The first URL is the paper that measured the
values. The second is the SAFE-MolGen paper (Zhang et al., J. Am. Chem. Soc.
2026, 148, 5520), which SAFE appends to every row curated for that work. It
is not a measurement source, and the dashboard previously misattributed it as
"the Liu 2025 review".

v0.3 splits that string into fields that each answer one question:

    DOI               bare DOI of the paper that reported the measurement
    URL               link to the source when it has no DOI (6 IDEAL rows)
    Locator           where in that paper: "Figure 3", "Table S2", ...
    Provenance_Level  "primary"    — transcribed by us from the paper itself
                      "aggregator" — taken from a compilation that cites it
    Aggregator        for aggregator rows: which compilation, its own DOI/URL,
                      the row's identifier inside it (SAFE exp_id, Liu
                      reference number, IDEAL reference index), and, where the
                      compilation tagged it, the publication the row was
                      curated for (`curation_DOI`)

`upgrade_reference` is the single code path that produces this shape. The
one-off repair `scripts/fix_reference_provenance.py` applies it to the stored
database; the ingest writers call it so newly ingested records arrive in the
same shape.
"""

from __future__ import annotations

import json
import re
from functools import lru_cache
from pathlib import Path
from typing import Dict, List, Optional, Tuple

from .schema import PROFILE, PROVENANCE_LEVELS

# A DOI is "10." + registrant code + "/" + suffix (DOI Handbook §2.2). The
# suffix may hold almost any printable character, but never whitespace; we
# also exclude the comma, which is how two DOIs ended up in one field.
DOI_PATTERN = re.compile(r"^10\.\d{4,9}/[^\s,]+$")

LIU_2022_DOI = "10.1021/jacsau.2c00122"
SAFE_DATABASE_DOI = "10.1080/07366299.2025.2564381"
SAFE_MOLGEN_DOI = "10.1021/jacs.5c19738"

# DOIs that describe a compilation or a design study, never a measurement.
# A condition must not cite one of these as the paper its D-values came from.
NON_MEASUREMENT_DOIS = {
    SAFE_DATABASE_DOI: "SAFE database announcement (Leite et al., Solvent Extr. "
                       "Ion Exch. 2025)",
    SAFE_MOLGEN_DOI: "SAFE-MolGen LLM ligand-design study (Zhang et al., J. Am. "
                     "Chem. Soc. 2026); SAFE tags rows curated for it",
}

# Compilations the database draws on, keyed by the condition_id prefix that
# every ingest writer already stamps. `name` matches metadata.data_source.
AGGREGATORS = {
    "Liu2022": {
        "name": "Liu_2022_JACSAu",
        "DOI": LIU_2022_DOI,
        "URL": None,
    },
    "SAFE": {
        "name": "SAFE_01262026",
        "DOI": SAFE_DATABASE_DOI,
        "URL": "https://safe.lanl.gov/",
    },
    "IDEAL": {
        "name": "OECD_NEA_IDEAL",
        "DOI": None,
        "URL": "https://www.oecd-nea.org/ideal/",
    },
}

# Condition-id prefixes for data transcribed by this project directly from the
# measuring paper, with no compilation in between. "Lit" marks conditions added
# from outputs/verification/additions/ (scripts/apply_verification.py).
PRIMARY_PREFIXES = ("Pramanik2024", "Lit")

_DOI_DECORATION = re.compile(
    r"^\s*(?:https?://)?(?:dx\.)?(?:doi\.org/)?(?:doi\s*:?\s*)?", re.IGNORECASE)

# "Figure 3", "fig1 b", "figS4-A", "Table SI1", "SI Figure 2", "TableS2",
# "Graph 5". The panel letter is taken only when it stands alone ("fig1 b",
# "S4A"), so "fig1 graphreader3" and "figure2 Ctodga" do not grow a panel.
_LOCATOR = re.compile(
    r"(?<![A-Za-z])(SI\s+)?(fig(?:ure)?|table|graph|scheme)\.?\s*(SI|S)?\s*(\d+)"
    r"(?:\s*-?\s*([A-Za-z])(?![A-Za-z]))?",
    re.IGNORECASE,
)
_LOCATOR_NAMES = {"fig": "Figure", "figure": "Figure", "table": "Table",
                  "graph": "Graph", "scheme": "Scheme"}
_DATA_LOCATION = re.compile(r"Data Location:\s*([^;]*)")
_IDEAL_REF = re.compile(r"Condition refs:\s*\[(\d+)\]")
_IDEAL_NAME = re.compile(r"references for\s+(.+?)\.\s")


def normalize_doi(raw) -> Optional[str]:
    """Return the bare, lower-cased DOI in `raw`, or None if it is not a DOI.

    DOIs are case-insensitive (DOI Handbook §3.4), so storing one case makes
    identical papers compare equal: "10.1039/D4GC01146E" and
    "10.1039/d4gc01146e" were counted as two sources before.
    """
    if raw is None:
        return None
    text = _DOI_DECORATION.sub("", str(raw)).strip().lstrip("/").rstrip(".").strip()
    if not DOI_PATTERN.match(text):
        return None
    return text.lower()


def split_reference_dois(raw) -> Tuple[List[str], List[str], List[str]]:
    """Split a legacy DOI field into (dois, urls, unparsed) tokens.

    Tokens are separated by commas or semicolons. A token that normalises to a
    DOI goes in `dois`, a remaining http(s) link in `urls`, anything else in
    `unparsed` (so the caller can report it rather than lose it).
    """
    dois, urls, unparsed = [], [], []
    if raw is None:
        return dois, urls, unparsed
    for token in re.split(r"[,;]", str(raw)):
        token = token.strip()
        if not token:
            continue
        doi = normalize_doi(token)
        if doi:
            if doi not in dois:
                dois.append(doi)
        elif re.match(r"https?://", token, re.IGNORECASE):
            urls.append(token)
        else:
            unparsed.append(token)
    return dois, urls, unparsed


def parse_locator(text) -> Optional[str]:
    """Canonical figure/table locator(s) in `text`, or None.

    >>> parse_locator("Data Location: Table SI1; Additional Comments: nan")
    'Table S1'
    >>> parse_locator("fig1 b nan nan extraction:nan")
    'Figure 1b'
    >>> parse_locator("Data Location: Figure 10 and Table 1; ...")
    'Figure 10; Table 1'
    """
    if not text:
        return None
    text = str(text)
    located = _DATA_LOCATION.search(text)
    if located:
        text = located.group(1)
    found = []
    for si_prefix, kind, supplementary, number, panel in _LOCATOR.findall(text):
        label = _LOCATOR_NAMES[kind.lower()]
        s = "S" if (si_prefix or supplementary) else ""
        entry = f"{label} {s}{number}{panel or ''}"
        if entry not in found:
            found.append(entry)
    return "; ".join(found) or None


@lru_cache(maxsize=1)
def liu_reference_dois() -> Dict[str, str]:
    """Liu 2022 SI reference number -> DOI of the paper it cites.

    The Liu ingest writers stamp the Liu paper's own DOI on every condition;
    the SI's reference number is what identifies the measuring paper. This is
    the table scripts/add_references_to_database.py applied historically.
    """
    path = Path(__file__).resolve().parents[1] / "data" / "liu_reference_database_with_dois.json"
    try:
        table = json.loads(path.read_text(encoding="utf-8"))
    except FileNotFoundError:
        return {}
    return {str(v["reference_number"]): normalize_doi(v["doi"])
            for v in table.values() if v.get("doi") and normalize_doi(v["doi"])}


def _aggregator_entry(prefix: str, condition_id: str, reference: Dict) -> Optional[str]:
    """The row's identifier inside its compilation, in that compilation's terms."""
    if prefix == "SAFE":
        tail = condition_id.rsplit("_", 1)[-1]
        return tail if tail.isdigit() else None
    if prefix == "Liu2022":
        number = reference.get("Reference_Number")
        return str(number) if number is not None else None
    if prefix == "IDEAL":
        notes = str(reference.get("Notes") or "")
        ref = _IDEAL_REF.search(notes)
        name = _IDEAL_NAME.search(notes)
        if ref:
            return f"{name.group(1) if name else '?'} ref [{ref.group(1)}]"
    return None


def upgrade_reference(reference: Dict, condition_id: str,
                      aggregator_name: Optional[str] = None) -> Dict:
    """Return `reference` in the v0.3 shape. Idempotent.

    Legacy keys (Source_Type, Data_Used, Reference_Number, Notes) are kept
    verbatim; `DOI` is replaced by the bare primary DOI. For Liu 2022 rows
    that cite only the Liu paper itself, the measuring paper is resolved from
    the SI reference number (`liu_reference_dois`).

    `aggregator_name` overrides the compilation label, for a newer snapshot
    of the same compilation (e.g. "SAFE_02102026").

    Raises:
        ValueError: if the condition_id prefix names no known source, or the
            DOI field holds text that is neither a DOI nor a URL. Both mean
            the record needs a curator, not a guess.
    """
    reference = dict(reference or {})
    prefix = str(condition_id).split("_", 1)[0]

    dois, urls, unparsed = split_reference_dois(reference.get("DOI"))
    if unparsed:
        raise ValueError(
            f"{condition_id}: DOI field holds unparseable text {unparsed!r}")

    curation = [d for d in dois if d in NON_MEASUREMENT_DOIS]
    measured = [d for d in dois if d not in NON_MEASUREMENT_DOIS]
    previous = reference.get("Aggregator") or {}
    if previous.get("curation_DOI") and previous["curation_DOI"] not in curation:
        curation.insert(0, previous["curation_DOI"])

    if prefix == "Liu2022" and measured in ([], [LIU_2022_DOI]):
        cited = liu_reference_dois().get(str(reference.get("Reference_Number")))
        if cited:
            measured = [cited]

    if prefix in AGGREGATORS:
        level = "aggregator"
        aggregator = dict(AGGREGATORS[prefix])
        if aggregator_name:
            aggregator["name"] = aggregator_name
        elif previous.get("name"):
            aggregator["name"] = previous["name"]
        # A split part's ID or a rewritten note no longer yields the entry; the row it came from is kept.
        aggregator["entry"] = _aggregator_entry(prefix, condition_id, reference) or previous.get("entry")
        aggregator["curation_DOI"] = curation[0] if curation else None
    elif prefix in PRIMARY_PREFIXES:
        level = "primary"
        aggregator = None
    else:
        raise ValueError(
            f"{condition_id}: unknown condition_id prefix {prefix!r}; add it to "
            f"rel_quest.references.AGGREGATORS or PRIMARY_PREFIXES")

    notes, data_used = reference.get("Notes"), reference.get("Data_Used")
    if reference.get("Locator"):
        # A structured locator (curator- or verification-set) wins over re-parsing the notes
        locator = reference["Locator"]
    elif prefix == "SAFE":
        locator = parse_locator(notes)
    else:
        locator = parse_locator(data_used)

    upgraded = {
        "DOI": measured[0] if measured else None,
    }
    url = urls[0] if urls else reference.get("URL")
    if upgraded["DOI"] is None and url:
        upgraded["URL"] = url
    upgraded["Locator"] = locator
    upgraded["Provenance_Level"] = level
    upgraded["Aggregator"] = aggregator
    for key in ("Source_Type", "Data_Used", "Reference_Number", "Notes"):
        if key in reference:
            upgraded[key] = reference[key]
    return upgraded


def check_reference(reference: Dict, require_provenance: bool = PROFILE == "dev") -> List[str]:
    """Problems with a v0.3 reference block (empty list if it is sound).

    The public release keeps only DOI/URL and Locator (schema.PUBLIC_FIELDS), so
    there the compilation fields are checked only when present.
    """
    problems = []
    if not isinstance(reference, dict):
        return ["reference is not an object"]
    doi = reference.get("DOI")
    if doi is not None:
        if not isinstance(doi, str) or not DOI_PATTERN.match(doi):
            problems.append(f"reference.DOI is not a bare DOI: {doi!r}")
        elif doi != doi.lower():
            problems.append(f"reference.DOI is not lower-case: {doi!r}")
        elif doi in NON_MEASUREMENT_DOIS:
            problems.append(
                f"reference.DOI {doi} is {NON_MEASUREMENT_DOIS[doi]}, not the "
                f"paper that measured these values")
    elif not reference.get("URL"):
        problems.append("reference has neither DOI nor URL")

    level = reference.get("Provenance_Level")
    if level is None and "Aggregator" not in reference and not require_provenance:
        pass
    elif level not in PROVENANCE_LEVELS:
        problems.append(
            f"reference.Provenance_Level {level!r} not in {list(PROVENANCE_LEVELS)}")
    aggregator = reference.get("Aggregator")
    if level == "aggregator":
        if not isinstance(aggregator, dict) or not aggregator.get("name"):
            problems.append("aggregator reference without Aggregator.name")
        else:
            curation = aggregator.get("curation_DOI")
            if curation is not None and not DOI_PATTERN.match(str(curation)):
                problems.append(f"Aggregator.curation_DOI is not a bare DOI: {curation!r}")
    elif level == "primary" and aggregator is not None:
        problems.append("primary reference must have Aggregator = null")

    locator = reference.get("Locator")
    if locator is not None and not isinstance(locator, str):
        problems.append(f"reference.Locator is not a string: {locator!r}")
    return problems


def primary_doi(reference: Dict) -> Optional[str]:
    """The measuring paper's DOI for any record generation (v0.1–v0.3).

    Reads the v0.3 field directly; for older records, splits the legacy
    string and skips compilation DOIs. Consumers that may be handed an
    archived generation should call this rather than read `DOI` raw.
    """
    if not isinstance(reference, dict):
        return None
    dois, _, _ = split_reference_dois(reference.get("DOI") or reference.get("doi"))
    for doi in dois:
        if doi not in NON_MEASUREMENT_DOIS:
            return doi
    return None


def doi_url(doi: Optional[str]) -> Optional[str]:
    """Resolver link for a bare DOI."""
    return f"https://doi.org/{doi}" if doi else None
