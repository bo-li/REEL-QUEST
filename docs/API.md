# REEL-QUEST Python API Documentation

## Installation

The REEL-QUEST Python library is located in the `rel_quest/` directory. Add the parent directory to your Python path:

```python
import sys
sys.path.append('/path/to/Database')

from rel_quest import LigandDatabase, LigandQuery
```

## Quick Start

```python
# Initialize database
db = LigandDatabase()

# Query by name
results = db.query_by_name('TODGA')
if results:
    ligand = results[0]
    print(f"Mechanism: {ligand.mechanism}")
    print(f"Trend: {ligand.selectivity_trend}")

# Get D-values
d_values = ligand.get_d_values()
print(d_values)  # {'La': 0.50, 'Ce': 0.85, ...}

# Get trend summary
trend = ligand.get_trend_summary()
print(trend['direction'])  # 'La→Lu increasing (light to heavy REE preference)'
```

## Core Classes

### `LigandDatabase`

Main database interface for querying ligand data.

#### Constructor

```python
LigandDatabase(db_path: Optional[str] = None, strict: bool = False)
```

**Parameters:**
- `db_path` (str, optional): Path to the JSONL database. If `None`, resolves
  the canonical file (`data/sx_ligands_merged_v14_flat.jsonl`) relative to the
  installed package, so it works from any working directory. Falling back to
  an archived generation emits a `RuntimeWarning`.
- `strict` (bool, optional): raise `DatabaseLoadError` on the first record that
  cannot be loaded, instead of skipping it. Default `False`; skipped records
  are always listed in `LigandDatabase.load_errors` and warned about.

**Raises:**
- `FileNotFoundError`: If database file not found.

**Example:**
```python
# Use default database
db = LigandDatabase()

# Use custom path
db = LigandDatabase('/path/to/custom_database.jsonl')
```

#### Methods

##### `query_by_smiles(smiles: str) -> Optional[LigandQuery]`

Query by exact SMILES match.

**Parameters:**
- `smiles` (str): SMILES string

**Returns:**
- `LigandQuery` object if found, `None` otherwise

**Example:**
```python
ligand = db.query_by_smiles('CCCCCCCCN(CCCCCCCC)C(=O)COCC(=O)N(CCCCCCCC)CCCCCCCC')
if ligand:
    print(f"Found: {ligand.name}")
```

##### `query_by_name(name: str, fuzzy: bool = False) -> List[LigandQuery]`

Query by ligand name (case-insensitive).

**Parameters:**
- `name` (str): Ligand name
- `fuzzy` (bool, optional): Enable substring matching. Default: `False`

**Returns:**
- List of `LigandQuery` objects (may be empty)

**Example:**
```python
# Exact match
results = db.query_by_name('TODGA')

# Fuzzy match (finds TODGA, DMDODGA, etc.)
results = db.query_by_name('DGA', fuzzy=True)
for ligand in results:
    print(ligand.name)
```

##### `filter_by_mechanism(mechanism: str) -> List[LigandQuery]`

Filter ligands by extraction mechanism.

**Parameters:**
- `mechanism` (str): Mechanism type

**Valid Mechanisms:**
- `'Cation_Exchange'`
- `'Solvation'`
- `'Anion_Exchange'`
- `'Ion_Pair'`
- `'Supramolecular'`
- `'Hybrid'`

**Returns:**
- List of `LigandQuery` objects

**Example:**
```python
solvation_ligands = db.filter_by_mechanism('Solvation')
print(f"Found {len(solvation_ligands)} solvation ligands")
```

##### `filter_by_trend(trend: str) -> List[LigandQuery]`

Filter ligands by selectivity trend.

**Parameters:**
- `trend` (str): Trend classification

**Valid Trends:**
- `'Ascending'` - La→Lu increasing
- `'Reverse'` - Lu→La increasing
- `'Flat'` - No clear trend
- `'Mixed'` - Conflicting trends
- `'Insufficient_Data'` - Too few data points

**Returns:**
- List of `LigandQuery` objects

**Example:**
```python
ascending = db.filter_by_trend('Ascending')
for ligand in ascending:
    print(f"{ligand.name}: {ligand.trend_confidence} confidence")
```

##### `filter_by_structural_class(structural_class: str) -> List[LigandQuery]`

Filter by Liu 2022 structural class.

**Parameters:**
- `structural_class` (str): Structural class name

**Valid Classes:**
- `'Diglycolamide (DGA)'`
- `'Aromatic ligand'`
- `'Phosphine oxide'`
- `'Other ligand'`

**Returns:**
- List of `LigandQuery` objects

**Example:**
```python
dga_ligands = db.filter_by_structural_class('Diglycolamide (DGA)')
```

##### `filter_by_elements(required_elements: List[str]) -> List[LigandQuery]`

Find ligands tested on all specified elements.

**Parameters:**
- `required_elements` (List[str]): List of element symbols

**Returns:**
- List of `LigandQuery` objects with data for all specified elements

**Example:**
```python
# Find ligands tested on Nd and Dy
ligands = db.filter_by_elements(['Nd', 'Dy'])
```

##### `search(query: str, search_fields: List[str] = None) -> List[LigandQuery]`

Search across multiple fields.

**Parameters:**
- `query` (str): Search term (case-insensitive)
- `search_fields` (List[str], optional): Fields to search. Default: `['name', 'smiles', 'structural_class']`

**Returns:**
- List of matching `LigandQuery` objects

**Example:**
```python
# Search for "BTP" in all fields
results = db.search('BTP')

# Search only in name field
results = db.search('BTP', search_fields=['name'])
```

##### `get_all_ligands() -> List[LigandQuery]`

Get all ligands in database.

**Returns:**
- List of all `LigandQuery` objects

**Example:**
```python
all_ligands = db.get_all_ligands()
print(f"Total ligands: {len(all_ligands)}")
```

##### `get_statistics() -> Dict`

Get database-level statistics.

**Returns:**
- Dictionary with counts and distributions

**Example:**
```python
stats = db.get_statistics()
print(f"Total ligands: {stats['total_ligands']}")
print(f"Mechanisms: {stats['mechanisms']}")
print(f"Trends: {stats['trends']}")
```

---

### `LigandQuery`

Query result object representing a single ligand.

#### Attributes

| Attribute | Type | Description |
|-----------|------|-------------|
| `smiles` | str | SMILES molecular structure |
| `name` | str | Ligand name |
| `mechanism` | str | Extraction mechanism |
| `structural_class` | str | Liu 2022 classification |
| `extractant_type` | str | Extractant functional type |
| `selectivity_trend` | str | Trend classification |
| `trend_confidence` | str | Confidence level |
| `element_coverage` | List[str] | Elements tested |
| `experimental_conditions` | List[Dict] | Experimental data |
| `performance_summary` | Dict | Performance metrics |
| `metadata` | Dict | Record metadata |

#### Methods

##### `get_d_values(condition_id: Optional[str] = None) -> Dict[str, float]`

Get raw distribution coefficients.

**Parameters:**
- `condition_id` (str, optional): Specific condition ID. If `None`, uses first condition.

**Returns:**
- Dictionary of `{element: D_value}`

**Example:**
```python
d_values = ligand.get_d_values()
print(f"La: {d_values['La']:.2f}")

# Get D-values for specific condition
d_values = ligand.get_d_values(condition_id='Liu2023_TODGA_HNO3_4.3M_295K')
```

##### `get_log_d_values(condition_id: Optional[str] = None) -> Dict[str, float]`

Get log₁₀(D) values.

**Parameters:**
- `condition_id` (str, optional): Specific condition ID

**Returns:**
- Dictionary of `{element: log_D_value}`

**Example:**
```python
log_d = ligand.get_log_d_values()
```

##### `get_trend_summary() -> Dict`

Get selectivity trend summary.

**Returns:**
- Dictionary with keys:
  - `'classification'` (str): Trend type
  - `'confidence'` (str): Confidence level
  - `'direction'` (str): Human-readable description
  - `'elements_tested'` (List[str]): Element symbols
  - `'coverage_count'` (int): Number of elements

**Example:**
```python
trend = ligand.get_trend_summary()
print(trend['classification'])  # 'Ascending'
print(trend['direction'])  # 'La→Lu increasing (light to heavy REE preference)'
```

##### `get_statistics(condition_id: Optional[str] = None) -> Dict`

Get statistical summary of D-values.

**Parameters:**
- `condition_id` (str, optional): Specific condition ID

**Returns:**
- The condition's `statistics` block in the canonical shape (2026-08-25): `mean_log_d`, `std_log_d`, `min_log_d`, `max_log_d`, `n_elements`

**Example:**
```python
stats = ligand.get_statistics()
print(f"Mean log D: {stats['mean_log_d']:.2f}")
print(f"Tested {stats['n_elements']} elements")
```

##### `get_condition_list() -> List[str]`

Get list of available condition IDs.

**Returns:**
- List of condition ID strings

**Example:**
```python
conditions = ligand.get_condition_list()
for cond_id in conditions:
    print(cond_id)
```

##### `get_separation_factor(elem1: str, elem2: str, condition_id: Optional[str] = None) -> Optional[float]`

Calculate separation factor between two elements.

**Parameters:**
- `elem1` (str): First element symbol
- `elem2` (str): Second element symbol
- `condition_id` (str, optional): Specific condition ID

**Returns:**
- Separation factor (≥ 1.0), or `None` if data unavailable

**Formula:**
- SF = max(D₁, D₂) / min(D₁, D₂)

**Example:**
```python
sf = ligand.get_separation_factor('Lu', 'La')
print(f"SF(Lu/La) = {sf:.1f}")
```

## Utilities

### `load_vocabulary(vocab_name: str) -> Dict`

Load controlled vocabulary file.

**Parameters:**
- `vocab_name` (str): Vocabulary name (`'mechanisms'`, `'elements'`, `'acids'`, `'solvents'`)

**Returns:**
- Dictionary with vocabulary data

**Example:**
```python
from rel_quest.utils import load_vocabulary

elements = load_vocabulary('elements')
lanthanides = elements['lanthanides']
```

### `canonicalize_smiles(smiles: str) -> Optional[str]`

Convert SMILES to canonical form using RDKit.

**Parameters:**
- `smiles` (str): SMILES string

**Returns:**
- Canonical SMILES, or `None` if invalid

**Requires:** RDKit (optional dependency)

**Example:**
```python
from rel_quest.utils import canonicalize_smiles

canonical = canonicalize_smiles('C(C)CC')
print(canonical)  # 'CCC'
```

### `validate_smiles(smiles: str) -> bool`

Check if SMILES string is valid.

**Parameters:**
- `smiles` (str): SMILES string

**Returns:**
- `True` if valid, `False` otherwise

**Example:**
```python
from rel_quest.utils import validate_smiles

is_valid = validate_smiles('CCO')  # True
is_valid = validate_smiles('invalid')  # False
```

## Complete Example

```python
from rel_quest import LigandDatabase

# Initialize database
db = LigandDatabase()

# Find all DGA ligands with ascending trends
dga_ligands = db.filter_by_structural_class('Diglycolamide (DGA)')
ascending_dga = [l for l in dga_ligands if l.selectivity_trend == 'Ascending']

print(f"Found {len(ascending_dga)} ascending DGA ligands:")
for ligand in ascending_dga:
    # Get trend info
    trend = ligand.get_trend_summary()

    # Get D-values
    d_values = ligand.get_d_values()

    # Calculate Lu/La separation factor
    if 'Lu' in d_values and 'La' in d_values:
        sf = ligand.get_separation_factor('Lu', 'La')
        print(f"  {ligand.name}: SF(Lu/La) = {sf:.1f}, confidence = {trend['confidence']}")
```

## Dependencies

**Required:**
- Python 3.7+
- pandas (for data ingestion)
- numpy (for statistics)
- scipy (for linear regression)

**Optional:**
- rdkit (for SMILES validation and canonicalization)

## Error Handling

All methods handle errors gracefully:
- Invalid queries return empty lists or `None`
- Missing data returns `None` for numeric values
- File I/O errors raise standard Python exceptions (`FileNotFoundError`, `IOError`)
