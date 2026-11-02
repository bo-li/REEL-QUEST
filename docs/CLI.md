# REEL-QUEST Command-Line Interface

## Overview

The REEL-QUEST CLI provides command-line tools for querying and validating the rare earth ligand database.

## Scripts

### 2. `query_ligands.py` - Query Database

Command-line interface for querying the ligand database.

#### Usage

```bash
python scripts/query_ligands.py [query_type] [options]
```

#### Query Types (Mutually Exclusive)

| Query | Description | Example |
|-------|-------------|---------|
| `--name NAME` | Query by ligand name | `--name TODGA` |
| `--smiles SMILES` | Query by SMILES string | `--smiles "CCC..."` |
| `--mechanism MECH` | Filter by mechanism | `--mechanism Solvation` |
| `--trend TREND` | Filter by selectivity trend | `--trend Ascending` |
| `--search TERM` | Search across fields | `--search DGA` |
| `--list-all` | List all ligands | `--list-all` |
| `--stats` | Show database statistics | `--stats` |

#### Output Formats

| Option | Values | Description |
|--------|--------|-------------|
| `--format` | `raw`, `trend`, `stats`, `full` | What information to display |
| `--output` | `json`, `table`, `csv` | Output style |

**Format Types:**
- `raw` - Distribution ratios (D-values) only
- `trend` - Selectivity trend information
- `stats` - Statistical summary (min/max/mean/std)
- `full` - Complete information (default)

**Output Styles:**
- `json` - Pretty JSON format (default)
- `table` - ASCII table format
- `csv` - CSV format

#### Additional Options

| Option | Description |
|--------|-------------|
| `--condition COND_ID` | Specify experimental condition ID |
| `--fuzzy` | Enable fuzzy name matching |

#### Examples

##### Basic Queries

```bash
# Query by name with full output
python scripts/query_ligands.py --name TODGA --format full

# Get raw D-values only
python scripts/query_ligands.py --name TODGA --format raw

# Query with fuzzy matching
python scripts/query_ligands.py --name DGA --fuzzy
```

##### Filtering

```bash
# Find all ascending trend ligands
python scripts/query_ligands.py --trend Ascending --format trend

# Find all solvation mechanism ligands
python scripts/query_ligands.py --mechanism Solvation

# Search for BTP-containing ligands
python scripts/query_ligands.py --search BTP
```

##### Output Formats

```bash
# JSON output (default)
python scripts/query_ligands.py --name DMDODGA --format full --output json

# Table output
python scripts/query_ligands.py --name DMDODGA --format trend --output table

# CSV output (good for spreadsheets)
python scripts/query_ligands.py --trend Ascending --format stats --output csv
```

##### Database Statistics

```bash
# Show database overview
python scripts/query_ligands.py --stats
```

**Output:**
```
Database Statistics:
==================================================
Total ligands: 169

Mechanisms:
  Cation_Exchange: 49
  Solvation: 120

Selectivity Trends:
  Ascending: 40
  Flat: 11
  Insufficient_Data: 64
  Mixed: 15
  Peaked: 8
  Reverse: 27
  U_Shaped: 4

Structural Classes:
  Aromatic ligand: 53
  Unknown: 40
  Diglycolamide (DGA): 31
  Other ligand: 12
  Phosphine oxide: 9
  N,N-dialkylamide: 1
  ...
```

---

### 3. `validate_database.py` - Validate Database

Checks database integrity and reports validation errors.

#### Usage

```bash
python scripts/validate_database.py [options]
```

#### Options

| Option | Description | Default |
|--------|-------------|---------|
| `--database`, `-d` | Path to JSONL database | canonical v14 file |
| `--max-records`, `-m` | Maximum records to validate | All |
| `--verbose`, `-v` | Show all error details | False |

#### Examples

```bash
# Validate full database
python scripts/validate_database.py

# Validate first 10 records only
python scripts/validate_database.py --max-records 10

# Show detailed errors
python scripts/validate_database.py --verbose
```

#### Validation Checks

The validator checks:
- Schema version present
- Required fields populated
- Valid SMILES structures
- Mechanism in canonical list
- Non-empty experimental_data
- Valid references with DOI
- Valid trend classifications
- JSON format correctness

**Success output:**
```
============================================================
REEL-QUEST Validation
============================================================
Database: data/sx_ligands_merged_v14_flat.jsonl
============================================================

Validation Results:
  Total records: 169
  Valid: 169
  Invalid: 0

Validation PASSED - All records are valid
```

**Failure output:**
```
Validation Results:
  Total records: 170
  Valid: 166
  Invalid: 4

Validation FAILED - 4 error(s) found

Error Details:
  - Line 5: Missing required ligand field: SMILES
  - Line 12: Invalid mechanism 'Unknown_Mechanism'
  ...
```

---

## Complete Workflows

### Initial Database Setup

```bash
# Validate database
python scripts/validate_database.py

# Check database statistics
python scripts/query_ligands.py --stats
```

### Finding Ligands for Nd/Dy Separation

```bash
# 1. Find ligands with ascending trends (prefer heavy REEs)
python scripts/query_ligands.py --trend Ascending --format csv --output csv > ascending_ligands.csv

# 2. Query specific ligand
python scripts/query_ligands.py --name DMDODGA --format full --output json

# 3. Search for all DGA ligands
python scripts/query_ligands.py --search DGA --format trend
```

### Analyzing Mechanisms

```bash
# Get all solvation ligands
python scripts/query_ligands.py --mechanism Solvation --output csv > solvation_ligands.csv

# Get all cation exchange ligands
python scripts/query_ligands.py --mechanism Cation_Exchange --output csv > cation_exchange.csv
```

## Output Formats

### JSON Output

```json
{
  "ligand": "DMDODGA",
  "mechanism": "Solvation",
  "structural_class": "Diglycolamide (DGA)",
  "trend": {
    "classification": "Ascending",
    "confidence": "High",
    "direction": "La→Lu increasing (light to heavy REE preference)",
    "elements_tested": ["La", "Ce", "Pr", "Nd", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu"],
    "coverage_count": 14
  },
  "statistics": {
    "D_min": 0.868,
    "D_max": 0.868,
    "D_mean": 0.868,
    "D_std": 0.0,
    "element_count": 1
  },
  "D_values": {
    "La": 0.868,
    ...
  }
}
```

### Table Output

```
======================================================================
Ligand: DMDODGA
======================================================================
Mechanism: Solvation
Structural Class: Diglycolamide (DGA)

Selectivity Trend: Ascending (High confidence)
Direction: La→Lu increasing (light to heavy REE preference)
Elements tested: La, Ce, Pr, Nd, Sm, Eu, Gd, Tb, Dy, Ho, Er, Tm, Yb, Lu (14 total)
======================================================================
```

### CSV Output

```csv
Ligand,Mechanism,Trend,Confidence,D_mean,D_std,Elements_Tested
DMDODGA,Solvation,Ascending,High,0.868,0.0,14
```

---

## Exit Codes

All scripts use standard exit codes:
- `0` - Success
- `1` - Error (with error message printed to stderr)

---

## Tips

1. **Piping output**: Use `--output csv` with `>` to save results to file
   ```bash
   python scripts/query_ligands.py --trend Ascending --output csv > results.csv
   ```

2. **Fuzzy search**: Use `--fuzzy` with `--name` to find similar ligands
   ```bash
   python scripts/query_ligands.py --name BTP --fuzzy
   ```

3. **JSON processing**: Pipe JSON output to `jq` for filtering
   ```bash
   python scripts/query_ligands.py --name TODGA | jq '.D_values'
   ```

4. **Validation after any change**: validate the records and check that every derived field reproduces
   ```bash
   python scripts/validate_database.py && python scripts/recompute_derived_fields.py --check
   ```

---

## Troubleshooting

### Database file not found

```
ERROR: Database file not found: data/sx_ligands_merged_v14_flat.jsonl
```

**Solution**: Run from a checkout that contains `data/sx_ligands_merged_v14_flat.jsonl`
(the path is resolved relative to the `rel_quest` package, so install it editable:
`pip install -e .`).

### No ligands found

```
No ligands found matching query
```

**Solutions**:
- Check spelling (names are case-insensitive)
- Try fuzzy matching with `--fuzzy`
- Use `--search` for broader matching
- Use `--stats` to see available ligands

## Environment

**Tested on:**
- Python 3.8+
- Linux (WSL2), macOS, Windows
- Dependencies: pandas, numpy, scipy, openpyxl

**Performance:**
- Query latency: <1 second (in-memory indexing)
- Memory usage: ~30 MB for the 170-ligand canonical database
