"""
REEL-QUEST Query System
Provides classes and methods for querying the ligand database
"""

import json
import os
import sys
import warnings
from dataclasses import dataclass
from typing import List, Dict, Optional

from .schema import TREND_DIRECTIONS
from .uncertainty import is_censored, is_mixture
from .utils import get_lanthanide_order


class DatabaseLoadError(Exception):
    """A record in the database file could not be loaded.

    Raised only when LigandDatabase(strict=True); the default is to skip the
    record, record it in `load_errors`, and warn.
    """


@dataclass
class LigandQuery:
    """
    Query result object representing a single ligand
    """
    smiles: Optional[str]
    name: str
    mechanism: str
    structural_class: str
    extractant_type: str
    selectivity_trend: str
    trend_confidence: str
    element_coverage: List[str]
    experimental_conditions: List[Dict]
    performance_summary: Dict
    metadata: Dict
    selectivity: Dict = None
    data_quality: Dict = None
    iupac_name: Optional[str] = None
    inchi: Optional[str] = None
    inchi_key: Optional[str] = None
    pubchem_cid: Optional[int] = None
    cas_number: Optional[str] = None
    name_source: Optional[str] = None
    name_confidence: Optional[str] = None
    community_name: Optional[str] = None
    chemical_family: Optional[str] = None  # Layer 1 (13 cats, see vocab/chemical_families.json)
    subtype: Optional[str] = None  # Layer 2 (~27 cats)
    alternative_names: Optional[List[str]] = None
    # Every element with an actual measurement, Y and Sc included. Distinct
    # from element_coverage, which lists only the La-Lu points the trend fit
    # used. Ask this one "was it measured?"; ask that one "what was fitted?".
    elements_measured: Optional[List[str]] = None

    def get_d_values(self, condition_id: Optional[str] = None) -> Dict[str, float]:
        """
        Get raw D-values for elements

        Args:
            condition_id: Specific condition ID, or None for first condition

        Returns:
            Dictionary of {element: D_value}
        """
        if not self.experimental_conditions:
            return {}

        if condition_id:
            cond = next((c for c in self.experimental_conditions
                        if c['condition_id'] == condition_id), None)
            if not cond:
                raise ValueError(f"Condition ID '{condition_id}' not found")
        else:
            cond = self.experimental_conditions[0]

        return {elem: data['value']
                for elem, data in cond['distribution_ratios'].items()
                if data['value'] is not None}

    def get_log_d_values(self, condition_id: Optional[str] = None) -> Dict[str, float]:
        """
        Get log(D) values for elements

        Args:
            condition_id: Specific condition ID, or None for first condition

        Returns:
            Dictionary of {element: log_D_value}
        """
        if not self.experimental_conditions:
            return {}

        if condition_id:
            cond = next((c for c in self.experimental_conditions
                        if c['condition_id'] == condition_id), None)
            if not cond:
                raise ValueError(f"Condition ID '{condition_id}' not found")
        else:
            cond = self.experimental_conditions[0]

        return {elem: data['log_value']
                for elem, data in cond['distribution_ratios'].items()
                if data['log_value'] is not None}

    def get_trend_summary(self) -> Dict:
        """
        Get selectivity trend classification summary

        Returns:
            Dictionary with classification, direction, and elements tested
        """
        return {
            'classification': self.selectivity_trend,
            'confidence': self.trend_confidence,
            'direction': TREND_DIRECTIONS.get(self.selectivity_trend, 'Unknown'),
            'elements_tested': self.element_coverage,
            'coverage_count': len(self.element_coverage)
        }

    def get_statistics(self, condition_id: Optional[str] = None) -> Dict:
        """
        Get statistical summary of D-values

        Args:
            condition_id: Specific condition ID, or None for first condition

        Returns:
            The condition's stored `statistics` block — since 2026-08-25 the
            one canonical shape: mean_log_d, std_log_d, min_log_d,
            max_log_d, n_elements (rel_quest.utils.condition_statistics).
        """
        if not self.experimental_conditions:
            return {}

        if condition_id:
            cond = next((c for c in self.experimental_conditions
                        if c['condition_id'] == condition_id), None)
            if not cond:
                raise ValueError(f"Condition ID '{condition_id}' not found")
        else:
            cond = self.experimental_conditions[0]

        return cond['statistics']

    def get_condition_list(self) -> List[str]:
        """Get list of available condition IDs"""
        return [exp['condition_id'] for exp in self.experimental_conditions]

    def get_separation_factor(self, elem1: str, elem2: str,
                              condition_id: Optional[str] = None) -> Optional[float]:
        """
        Calculate separation factor between two elements

        SF = D(elem1) / D(elem2) if D(elem1) > D(elem2), else D(elem2) / D(elem1)

        Args:
            elem1: First element symbol
            elem2: Second element symbol
            condition_id: Specific condition ID, or None for first condition

        Returns:
            Separation factor (always >= 1), or None if data not available
            or either D is censored (a bound, not a value)
        """
        d_values = self.get_d_values(condition_id)

        if elem1 not in d_values or elem2 not in d_values:
            return None
        # A censored D (stored as its bound, e.g. D > 100) gives no ratio.
        cond = (next(c for c in self.experimental_conditions if c['condition_id'] == condition_id)
                if condition_id else self.experimental_conditions[0])
        cells = cond.get('distribution_ratios') or {}
        if any((cells.get(e) or {}).get('censor') in ('>', '<') for e in (elem1, elem2)):
            return None

        d1 = d_values[elem1]
        d2 = d_values[elem2]

        if d1 == 0 or d2 == 0:
            return None

        return max(d1, d2) / min(d1, d2)

    def get_all_elements_tested(self, lanthanides_only: bool = False) -> List[str]:
        """
        Get all elements actually tested across all conditions

        Args:
            lanthanides_only: If True, exclude Y and Sc (return only La-Lu)

        Returns:
            Sorted list of element symbols tested

        Note: element_coverage only includes lanthanides (La-Lu) used for trend
        analysis. This method can return either all rare earths (including Y, Sc)
        or just lanthanides.
        """
        from .utils import get_lanthanide_order

        all_elements = set()
        for cond in self.experimental_conditions:
            all_elements.update(cond['distribution_ratios'].keys())

        if lanthanides_only:
            lanthanides = set(get_lanthanide_order())
            all_elements = all_elements & lanthanides

        return sorted(all_elements)

    def get_aggregated_d_values(self, aggregation: str = 'median',
                               lanthanides_only: bool = False) -> Dict[str, Dict]:
        """
        Get aggregated D-values across all experimental conditions

        Useful for ligands with multiple conditions (especially SAFE data)
        where elements are measured separately across many conditions.

        Args:
            aggregation: Method to aggregate ('median', 'mean', 'min', 'max')
            lanthanides_only: If True, exclude Y and Sc (return only La-Lu)

        Returns:
            Dictionary of {element: {'value': aggregated_D, 'log_value': aggregated_log_D,
                                     'n_conditions': count, 'std': std_deviation,
                                     'n_bounds': censored cells left out}}

        A censored cell (D > x, D < x) is a bound, not a value: it is left
        out of the aggregate and counted in 'n_bounds'. An element measured
        only as bounds has no aggregate and is absent. A mixed-extractant
        condition (uncertainty.is_mixture) is left out entirely.
        """
        import numpy as np
        from .utils import get_lanthanide_order

        # Get lanthanides set if filtering
        lanthanides = set(get_lanthanide_order()) if lanthanides_only else None

        # Collect all D-values by element
        element_data = {}
        for cond in self.experimental_conditions:
            if is_mixture(cond):
                continue
            for elem, data in cond['distribution_ratios'].items():
                # Skip non-lanthanides if filtering
                if lanthanides_only and elem not in lanthanides:
                    continue

                if elem not in element_data:
                    element_data[elem] = {'d_values': [], 'log_d_values': [], 'n_bounds': 0}

                if is_censored(data):
                    element_data[elem]['n_bounds'] += 1
                    continue
                if data['value'] is not None:
                    element_data[elem]['d_values'].append(data['value'])
                if data['log_value'] is not None:
                    element_data[elem]['log_d_values'].append(data['log_value'])

        # Aggregate
        aggregated = {}
        for elem, values in element_data.items():
            d_vals = values['d_values']
            log_d_vals = values['log_d_values']

            if not d_vals:
                continue

            # Choose aggregation method
            if aggregation == 'median':
                agg_d = float(np.median(d_vals))
                agg_log_d = float(np.median(log_d_vals)) if log_d_vals else None
            elif aggregation == 'mean':
                agg_d = float(np.mean(d_vals))
                agg_log_d = float(np.mean(log_d_vals)) if log_d_vals else None
            elif aggregation == 'min':
                agg_d = float(np.min(d_vals))
                agg_log_d = float(np.min(log_d_vals)) if log_d_vals else None
            elif aggregation == 'max':
                agg_d = float(np.max(d_vals))
                agg_log_d = float(np.max(log_d_vals)) if log_d_vals else None
            else:
                raise ValueError(f"Unknown aggregation method: {aggregation}")

            aggregated[elem] = {
                'value': agg_d,
                'log_value': agg_log_d,
                'n_conditions': len(d_vals),
                'std': float(np.std(d_vals)) if len(d_vals) > 1 else 0.0,
                'min': float(np.min(d_vals)),
                'max': float(np.max(d_vals)),
                'n_bounds': values['n_bounds'],
            }

        return aggregated

    def get_conditions_by_element(self, element: str) -> List[Dict]:
        """
        Get all experimental conditions containing a specific element

        Useful for analyzing how D-values vary with experimental conditions
        for a particular element.

        Args:
            element: Element symbol (e.g., 'Eu', 'Nd', 'Y')

        Returns:
            List of condition dictionaries containing the specified element
        """
        conditions = []
        for cond in self.experimental_conditions:
            if element in cond['distribution_ratios']:
                conditions.append({
                    'condition_id': cond['condition_id'],
                    'acid_name': cond.get('Acid_Name', 'Unknown'),
                    'acid_concentration': cond.get('Acid_Concentration_M', 'Unknown'),
                    'D_value': cond['distribution_ratios'][element]['value'],
                    'log_D_value': cond['distribution_ratios'][element]['log_value'],
                    'n_replicates': cond['distribution_ratios'][element].get('n_replicates', 1)
                })
        return conditions

    def get_conditions_by_source(self, source: str) -> List[Dict]:
        """
        Get experimental conditions from a specific data source

        For merged ligands with data from multiple sources, this filters
        conditions to only those from the specified source.

        Args:
            source: Data source name (e.g., 'SAFE_Database', 'Liu_2022_JACSAu', 'OECD_NEA_IDEAL')

        Returns:
            List of experimental condition dictionaries from that source
        """
        # Check if this is a single-source ligand
        if not isinstance(self.metadata.get('data_sources'), list):
            # Single source ligand
            if self.metadata.get('data_source') == source:
                return self.experimental_conditions
            else:
                return []

        # Merged ligand - filter by condition.data_source
        return [
            cond for cond in self.experimental_conditions
            if cond.get('data_source') == source
        ]

    def get_trend_by_source(self, source: str) -> Dict:
        """
        Get selectivity trend from a specific data source

        For merged ligands, this returns the trend classification specific
        to one data source rather than the consensus trend.

        Args:
            source: Data source name (e.g., 'SAFE_Database')

        Returns:
            Dictionary with trend_classification and other trend info for that source,
            or {'trend_classification': 'Unknown'} if source not found
        """
        # Get trends_by_source from selectivity
        # This is populated by the merge process for merged ligands
        trends_by_source = {}

        # Check if metadata has selectivity info
        if hasattr(self, 'selectivity') and isinstance(self.selectivity, dict):
            trends_by_source = self.selectivity.get('trends_by_source', {})

        if source in trends_by_source:
            return trends_by_source[source]

        # If not merged or source not found, check if this is single-source
        if not isinstance(self.metadata.get('data_sources'), list):
            if self.metadata.get('data_source') == source:
                # Return overall trend for single-source ligand
                return {
                    'trend_classification': self.selectivity_trend,
                    'trend_confidence': self.trend_confidence,
                    'element_coverage': self.element_coverage
                }

        return {'trend_classification': 'Unknown'}

    @property
    def tested_elements(self) -> List[str]:
        """Elements this ligand has measurements for, Y and Sc included.

        Prefer this over `element_coverage` for "was it measured?" questions.
        Coverage lists only the La-Lu points the trend fit used, so filtering
        on it drops Y and Sc entirely. Falls back to deriving from the
        conditions for records written before `elements_measured` existed.
        """
        if self.elements_measured is not None:
            return list(self.elements_measured)

        found = set()
        for condition in self.experimental_conditions:
            for element, data in (condition.get('distribution_ratios') or {}).items():
                if isinstance(data, dict) and (
                    data.get('value') is not None or data.get('log_value') is not None
                ):
                    found.add(element)
        return sorted(found)

    @property
    def data_sources(self):
        """
        Get list of data sources for this ligand.

        Returns:
            list: Data sources (e.g., ['Liu_2022_JACSAu', 'SAFE_01262026'])
        """
        # Check for plural form first (merged ligands)
        if isinstance(self.metadata.get('data_sources'), list):
            return self.metadata['data_sources']

        # Check for singular form (single-source ligands)
        elif self.metadata.get('data_source'):
            return [self.metadata['data_source']]

        # No source labels: the public release (schema.PUBLIC_FIELDS) carries none
        return []


class LigandDatabase:
    """
    Main database interface for querying ligand data
    """

    def __init__(self, db_path: Optional[str] = None, strict: bool = False):
        """
        Initialize database and build indexes

        Args:
            db_path: Path to JSONL database file. If None, uses default path.
            strict: Raise DatabaseLoadError on the first unusable record
                instead of skipping it. Skipped records are reported either
                way — see `load_errors`.
        """
        if db_path is None:
            # Default path relative to this file. The canonical file comes
            # first; archived generations are a last-resort fallback that is
            # loudly warned about — silently serving stale data to a research
            # workflow is worse than failing.
            db_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
            candidates = [
                os.path.join('data', 'sx_ligands_merged_v14_flat.jsonl'),                  # canonical (verified, 162 ligs)
                os.path.join('data', 'archive', 'sx_ligands_merged_v13_flat.jsonl'),
                os.path.join('data', 'archive', 'sx_ligands_merged_v12_tddga_fixed.jsonl'),
                os.path.join('data', 'archive', 'sx_ligands_merged_v11_deduplicated.jsonl'),
            ]
            for idx, fname in enumerate(candidates):
                p = os.path.join(db_dir, fname)
                if os.path.exists(p):
                    if idx > 0:
                        warnings.warn(
                            f"Canonical database '{candidates[0]}' not found; "
                            f"falling back to older generation '{fname}'. "
                            "Results may be outdated.",
                            RuntimeWarning,
                            stacklevel=2,
                        )
                    db_path = p
                    break
            if db_path is None:
                raise FileNotFoundError(
                    f"No database file found. Expected '{candidates[0]}' "
                    f"relative to {db_dir}."
                )

        self.db_path = db_path

        if not os.path.exists(db_path):
            raise FileNotFoundError(f"Database file not found: {db_path}")

        self._records = []
        # Records that could not be loaded, as {'line', 'reason'} dicts. Empty
        # on a clean load; callers that must not operate on partial data can
        # assert it. Populated by _build_index.
        self.load_errors: List[Dict] = []
        self._index = self._build_index(strict=strict)

    @staticmethod
    def _required_fields(record) -> Dict:
        """Pull the fields every index needs, or raise ValueError saying which.

        Separated from indexing so a record is fully checked before anything
        is mutated: the loader must never half-admit a record.
        """
        if not isinstance(record, dict):
            raise ValueError(
                f"record is {type(record).__name__}, expected a JSON object")

        missing = [key for key in ('ligand', 'selectivity') if key not in record]
        if missing:
            raise ValueError(f"missing required container(s): {missing}")
        for key in ('ligand', 'selectivity'):
            if not isinstance(record[key], dict):
                raise ValueError(
                    f"'{key}' is {type(record[key]).__name__}, expected an object")

        ligand, selectivity = record['ligand'], record['selectivity']
        for container, key in (('ligand', 'Name'), ('ligand', 'Mechanism'),
                               ('selectivity', 'trend_classification')):
            if not record[container].get(key):
                raise ValueError(f"missing required field: {container}.{key}")
        if not isinstance(ligand['Name'], str):
            raise ValueError(
                f"ligand.Name is {type(ligand['Name']).__name__}, expected a string")

        return {
            'smiles': ligand.get('SMILES'),
            'name': ligand['Name'],
            'mechanism': ligand['Mechanism'],
            'trend': selectivity['trend_classification'],
            'structural_class': ligand.get('Structural_Class', ''),
            'community_name': ligand.get('Community_Name'),
            'chem_family': ligand.get('Chemical_Family'),
            'subtype': ligand.get('Subtype'),
        }

    def _build_index(self, strict: bool = False) -> Dict:
        """
        Build in-memory indexes for fast lookup

        A record is validated first, then appended to `self._records` and
        indexed together. Previously the append came first and the field
        access after it, inside one broad `except` — so a record missing a
        required field was warned about, left in `_records`, and never
        indexed. `get_statistics()` then counted a ligand that no query could
        return, and `get_all_ligands()` raised on it. A research database must
        not serve a record set it cannot explain.

        Args:
            strict: raise DatabaseLoadError on the first unusable record
                rather than skipping it.

        Returns:
            Dictionary with indexes: by_smiles, by_name, by_mechanism, by_trend
        """
        index = {
            'by_smiles': {},
            'by_name': {},
            'by_mechanism': {},
            'by_trend': {},
            'by_structural_class': {},
            'by_community_name': {},
            'by_chemical_family': {},
            'by_subtype': {},
        }

        with open(self.db_path, 'r') as f:
            for line_num, line in enumerate(f, 1):
                if not line.strip():
                    continue

                try:
                    record = json.loads(line)
                    fields = self._required_fields(record)
                except json.JSONDecodeError as e:
                    reason = f"invalid JSON: {e}"
                except ValueError as e:
                    reason = str(e)
                else:
                    # Only past this point does the record exist to the
                    # database — append and index in one uninterrupted block.
                    self._records.append(record)

                    if fields['smiles']:
                        index['by_smiles'][fields['smiles']] = record
                    index['by_name'][fields['name'].lower()] = record
                    index['by_mechanism'].setdefault(
                        fields['mechanism'], []).append(record)
                    index['by_trend'].setdefault(fields['trend'], []).append(record)
                    if fields['structural_class']:
                        index['by_structural_class'].setdefault(
                            fields['structural_class'], []).append(record)
                    if fields['community_name']:
                        index['by_community_name'].setdefault(
                            fields['community_name'], []).append(record)
                    # Two-layer classification fields (added in v10/v12)
                    if fields['chem_family']:
                        index['by_chemical_family'].setdefault(
                            fields['chem_family'], []).append(record)
                    if fields['subtype']:
                        index['by_subtype'].setdefault(
                            fields['subtype'], []).append(record)
                    continue

                # Unusable record: recorded, never half-admitted.
                self.load_errors.append({'line': line_num, 'reason': reason})
                if strict:
                    raise DatabaseLoadError(
                        f"{self.db_path} line {line_num}: {reason}")

        if self.load_errors:
            warnings.warn(
                f"{len(self.load_errors)} record(s) in {self.db_path} could not "
                f"be loaded and were skipped: "
                + "; ".join(f"line {e['line']}: {e['reason']}"
                            for e in self.load_errors[:5])
                + ("; ..." if len(self.load_errors) > 5 else "")
                + ". Inspect LigandDatabase.load_errors, or pass strict=True "
                  "to fail instead of skipping.",
                RuntimeWarning,
                stacklevel=3,
            )

        # stderr, not stdout: the MCP server's stdout is its protocol channel.
        print(f"Loaded {len(self._records)} ligands from {self.db_path}", file=sys.stderr)
        return index

    def _record_to_query(self, record: Dict) -> LigandQuery:
        """Convert JSON record to LigandQuery object"""
        return LigandQuery(
            smiles=record['ligand'].get('SMILES'),
            name=record['ligand']['Name'],
            mechanism=record['ligand']['Mechanism'],
            structural_class=record['ligand'].get('Structural_Class', ''),
            extractant_type=record['ligand']['Extractant_Type'],
            selectivity_trend=record['selectivity']['trend_classification'],
            trend_confidence=record['selectivity']['trend_confidence'],
            element_coverage=record['selectivity']['element_coverage'],
            experimental_conditions=record['experimental_data'],
            performance_summary=record['performance_summary'],
            metadata=record['metadata'],
            selectivity=record.get('selectivity'),
            data_quality=record.get('data_quality'),
            iupac_name=record['ligand'].get('IUPAC_Name'),
            inchi=record['ligand'].get('InChI'),
            inchi_key=record['ligand'].get('InChI_Key'),
            pubchem_cid=record['ligand'].get('PubChem_CID'),
            cas_number=record['ligand'].get('CAS_Number'),
            name_source=record['ligand'].get('Name_Source'),
            name_confidence=record['ligand'].get('Name_Confidence'),
            community_name=record['ligand'].get('Community_Name'),
            chemical_family=record['ligand'].get('Chemical_Family'),
            subtype=record['ligand'].get('Subtype'),
            alternative_names=record['ligand'].get('Alternative_Names'),
            elements_measured=record.get('elements_measured'),
        )

    def query_by_smiles(self, smiles: str) -> Optional[LigandQuery]:
        """
        Query by exact SMILES match

        Args:
            smiles: SMILES string

        Returns:
            LigandQuery object, or None if not found
        """
        record = self._index['by_smiles'].get(smiles)
        return self._record_to_query(record) if record else None

    def query_by_name(self, name: str, fuzzy: bool = False) -> List[LigandQuery]:
        """
        Query by ligand name

        Args:
            name: Ligand name (case-insensitive)
            fuzzy: If True, perform substring matching

        Returns:
            List of LigandQuery objects (may be empty)
        """
        name_lower = name.lower()

        if not fuzzy:
            record = self._index['by_name'].get(name_lower)
            return [self._record_to_query(record)] if record else []
        else:
            # Fuzzy matching: partial substring
            matches = [r for n, r in self._index['by_name'].items()
                      if name_lower in n]
            return [self._record_to_query(r) for r in matches]

    def filter_by_mechanism(self, mechanism: str) -> List[LigandQuery]:
        """
        Filter ligands by extraction mechanism

        Args:
            mechanism: Mechanism type (e.g., 'Solvation', 'Cation_Exchange')

        Returns:
            List of LigandQuery objects
        """
        records = self._index['by_mechanism'].get(mechanism, [])
        return [self._record_to_query(r) for r in records]

    def filter_by_trend(self, trend: str) -> List[LigandQuery]:
        """
        Filter ligands by selectivity trend

        Args:
            trend: Trend classification (Ascending/Reverse/Flat/Mixed/Insufficient_Data)

        Returns:
            List of LigandQuery objects
        """
        records = self._index['by_trend'].get(trend, [])
        return [self._record_to_query(r) for r in records]

    def filter_by_structural_class(self, structural_class: str) -> List[LigandQuery]:
        """
        Filter ligands by structural class

        Args:
            structural_class: Structural class from Liu 2022

        Returns:
            List of LigandQuery objects
        """
        records = self._index['by_structural_class'].get(structural_class, [])
        return [self._record_to_query(r) for r in records]

    def filter_by_community_name(self, community_name: str) -> List[LigandQuery]:
        """
        Filter ligands by community-recognized name (DAPhen, TODGA, BTBP, ...).

        Args:
            community_name: Community-form Type label (matches Community_Name field).

        Returns:
            List of LigandQuery objects
        """
        records = self._index['by_community_name'].get(community_name, [])
        return [self._record_to_query(r) for r in records]

    def filter_by_chemical_family(self, chemical_family: str) -> List[LigandQuery]:
        """
        Filter ligands by Layer-1 Chemical_Family (Diglycolamide, Phenanthroline_Amide, ...).

        See vocab/chemical_families.json for the 13 valid family names.

        Args:
            chemical_family: Layer-1 Chemical_Family label.

        Returns:
            List of LigandQuery objects
        """
        records = self._index['by_chemical_family'].get(chemical_family, [])
        return [self._record_to_query(r) for r in records]

    def filter_by_subtype(self, subtype: str) -> List[LigandQuery]:
        """
        Filter ligands by Layer-2 Subtype (DGA, BisDGA, Malonamide, DAPhen, BTP, ...).

        See vocab/chemical_families.json for the ~27 valid subtype labels.

        Args:
            subtype: Layer-2 Subtype label.

        Returns:
            List of LigandQuery objects
        """
        records = self._index['by_subtype'].get(subtype, [])
        return [self._record_to_query(r) for r in records]

    def filter_by_elements(self, required_elements: List[str]) -> List[LigandQuery]:
        """
        Find ligands tested on all specified elements

        Reads `elements_measured`, which records what was actually measured.
        It previously read `selectivity.element_coverage`, but that field
        lists only the La-Lu points the trend fit consumed, so Y and Sc
        queries returned nothing however much data existed. Records written
        before elements_measured existed fall back to the old field.

        Args:
            required_elements: List of element symbols (e.g., ['La', 'Nd', 'Dy'])

        Returns:
            List of LigandQuery objects that have data for all specified elements
        """
        results = []
        required_set = set(required_elements)

        for record in self._records:
            measured = record.get('elements_measured')
            if measured is None:
                measured = record['selectivity']['element_coverage']
            if required_set.issubset(set(measured)):
                results.append(self._record_to_query(record))

        return results

    def get_all_ligands(self) -> List[LigandQuery]:
        """Get all ligands in database"""
        return [self._record_to_query(r) for r in self._records]

    def has_source_labels(self) -> bool:
        """True when records name the compilation they came from (metadata.data_source).

        The developer database does; the public release does not
        (schema.PUBLIC_FIELDS), and every view that filters or labels by
        source hides itself there instead of showing 'Unknown'.
        """
        return any((r.get('metadata') or {}).get('data_source') for r in self._records)

    def get_statistics(self) -> Dict:
        """
        Get database-level statistics

        Returns:
            Dictionary with counts and distributions
        """
        mechanism_counts = {k: len(v) for k, v in self._index['by_mechanism'].items()}
        trend_counts = {k: len(v) for k, v in self._index['by_trend'].items()}
        structural_class_counts = {k: len(v) for k, v in self._index['by_structural_class'].items()}

        return {
            'total_ligands': len(self._records),
            'mechanisms': mechanism_counts,
            'trends': trend_counts,
            'structural_classes': structural_class_counts
        }

    def search(self, query: str, search_fields: List[str] = None) -> List[LigandQuery]:
        """
        Search across multiple fields

        Args:
            query: Search term (case-insensitive)
            search_fields: Fields to search (default: ['name', 'smiles', 'subtype',
                'structural_class']; the legacy structural_class is absent from the public release)

        Returns:
            List of LigandQuery objects matching the query
        """
        if search_fields is None:
            search_fields = ['name', 'smiles', 'subtype', 'structural_class']

        query_lower = query.lower()
        results = []

        for record in self._records:
            match = False

            if 'name' in search_fields:
                if query_lower in record['ligand']['Name'].lower():
                    match = True

            if 'smiles' in search_fields:
                smiles = record['ligand'].get('SMILES')
                if smiles and query_lower in smiles.lower():
                    match = True

            if 'subtype' in search_fields:
                subtype = record['ligand'].get('Subtype') or ''
                if query_lower in subtype.lower():
                    match = True

            if 'structural_class' in search_fields:
                structural_class = record['ligand'].get('Structural_Class') or ''
                if query_lower in structural_class.lower():
                    match = True

            if match:
                results.append(self._record_to_query(record))

        return results
