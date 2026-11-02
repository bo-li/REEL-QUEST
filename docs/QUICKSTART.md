# REEL-QUEST Quick Start Guide

Get started with the Rare Earth eXtraction Database in 5 minutes.

## Installation

```bash
# Clone or navigate to the database directory
cd Database

# Install dependencies
pip install pandas rdkit scipy numpy

# Optional: for visualization
pip install matplotlib
```

## Quick Examples

### 1. Load the Database

```python
from rel_quest import LigandDatabase

# Initialize — resolves the canonical v13 database (170 ligands)
db = LigandDatabase()
```

### 2. Query a Ligand by Name

```python
# Query TODGA (popular diglycolamide extractant)
results = db.query_by_name('TODGA')
if results:
    todga = results[0]
    print(f"Name: {todga.name}")
    print(f"Mechanism: {todga.mechanism}")
    print(f"Trend: {todga.selectivity_trend}")
```

**Output:**
```
Name: TODGA
Mechanism: Solvation
Trend: Flat
```

### 3. Get D-values for All Elements

```python
# Get distribution ratios
d_values = todga.get_d_values()

# Print first few elements
for elem, d_val in list(d_values.items())[:5]:
    print(f"{elem}: D = {d_val:.3f}")
```

**Output:**
```
La: D = 4.560
Ce: D = 5.840
Pr: D = 7.650
Nd: D = 8.410
Sm: D = 11.800
```

### 4. Get Trend Summary

```python
trend = todga.get_trend_summary()
print(f"Classification: {trend['classification']}")
print(f"Direction: {trend['direction']}")
print(f"Elements tested: {trend['coverage_count']}")
```

**Output:**
```
Classification: Flat
Direction: No clear selectivity trend across lanthanide series
Elements tested: 14
```

### 5. Calculate Separation Factor

```python
# Calculate SF between La and Lu
sf = todga.get_separation_factor('La', 'Lu')
print(f"SF(La/Lu) = {sf:.2f}")
```

**Output:**
```
SF(La/Lu) = 2.77
```

---

## Common Queries

### Find Ligands by Mechanism

```python
# Get all solvation extractants
solvation = db.filter_by_mechanism('Solvation')
print(f"Found {len(solvation)} solvation ligands")

# Get all cation exchange ligands
cation_exchange = db.filter_by_mechanism('Cation_Exchange')
print(f"Found {len(cation_exchange)} cation exchange ligands")
```

### Find Ligands by Selectivity Trend

```python
# Find ligands with ascending trend (La→Lu increasing)
ascending = db.filter_by_trend('Ascending')
print(f"Ascending trend ligands: {len(ascending)}")

# Print first 5
print("\nTop 5 ascending trend ligands:")
for lig in ascending[:5]:
    print(f"  - {lig.name} ({lig.mechanism})")
```

**Output:**
```
Ascending trend ligands: 26

Top 5 ascending trend ligands:
  - CA-BTP (Cation_Exchange)
  - n-propyl-BTP (Cation_Exchange)
  - C5-BTBP (Cation_Exchange)
  - 4 (Solvation)
  - DMDODGA (Solvation)
```

### Search by Keyword

```python
# Search for all DGA ligands
dga_ligands = db.search('DGA')
print(f"Found {len(dga_ligands)} ligands containing 'DGA'")

for lig in dga_ligands[:3]:
    print(f"  - {lig.name}")
```

**Output:**
```
Found 11 ligands containing 'DGA'
  - TODGA
  - DMDODGA
  - TEHDGA
```

---

## Using the Command-Line Interface

### Get Database Statistics

```bash
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
```

### Query a Ligand

```bash
# Table format (human-readable)
python scripts/query_ligands.py --name TODGA --format full --output table
```

**Output:**
```
======================================================================
Ligand: TODGA
======================================================================
Mechanism: Solvation
Structural Class: diglycolamide

Selectivity Trend: Flat (High confidence)
Direction: No clear selectivity trend across lanthanide series
Elements tested: La, Ce, Pr, Nd, Sm, Eu, Gd, Tb, Dy, Ho, Er, Tm, Yb, Lu (14 total)

Statistics:
  D_min: 4.56
  D_max: 12.64
  D_mean: 8.37
  D_std: 2.43
  ...
```

```bash
# JSON format (for scripts)
python scripts/query_ligands.py --name TBP --format full --output json
```

```bash
# CSV format (for Excel)
python scripts/query_ligands.py --name TBP --format raw --output csv > tbp_data.csv
```

### Filter Ligands

```bash
# Get all ascending trend ligands
python scripts/query_ligands.py --trend Ascending --format trend

# Get all solvation ligands
python scripts/query_ligands.py --mechanism Solvation --format trend
```

---

## Real-World Example: Find Best Ligand for Nd/Dy Separation

Neodymium (Nd) and dysprosium (Dy) separation is critical for recycling rare earth magnets (NdFeB). Let's find the best ligand for this separation.

```python
from rel_quest import LigandDatabase

# Load database
db = LigandDatabase()

# Find ligands with data for both Nd and Dy
candidates = db.filter_by_elements(['Nd', 'Dy'])

print(f"Found {len(candidates)} ligands with Nd and Dy data\n")

# Calculate separation factors and rank
sf_data = []
for ligand in candidates:
    sf = ligand.get_separation_factor('Nd', 'Dy')
    if sf:
        sf_data.append({
            'name': ligand.name,
            'SF': sf,
            'mechanism': ligand.mechanism,
            'trend': ligand.selectivity_trend
        })

# Sort by SF (descending)
sf_data.sort(key=lambda x: x['SF'], reverse=True)

# Print top 10
print("Top 10 ligands for Nd/Dy separation:")
print("-" * 70)
for i, data in enumerate(sf_data[:10], 1):
    print(f"{i:2d}. {data['name']:20s}  SF = {data['SF']:6.2f}  ({data['mechanism']})")
```

**Example Output:**
```
Found 31 ligands with Nd and Dy data

Top 10 ligands for Nd/Dy separation:
----------------------------------------------------------------------
 1. CA-BTP               SF =   3.45  (Cation_Exchange)
 2. n-propyl-BTP         SF =   2.89  (Cation_Exchange)
 3. C5-BTBP              SF =   2.54  (Cation_Exchange)
 4. DMDODGA              SF =   2.12  (Solvation)
 5. TODGA                SF =   1.84  (Solvation)
 ...
```

---

## Data Sources

The database was compiled from the Liu et al. 2022 JACS Au Supporting
Information, the SAFE database, the OECD-NEA IDEaL database and papers read
directly, and every condition was checked against the paper that measured it.
Each condition cites that paper and where in it the values are:

```python
ligand = db.query_by_name('TBP')[0]
for condition in ligand.experimental_conditions[:3]:
    ref = condition['reference']
    print(ref['DOI'], ref.get('Locator'))
```

---

## Next Steps

- **Explore more examples**: See [USAGE_EXAMPLES.md](USAGE_EXAMPLES.md) for advanced use cases
- **Learn the API**: Read [API.md](API.md) for complete API documentation
- **Use the CLI**: Check [CLI.md](CLI.md) for command-line usage
- **Understand the data**: Review [DATABASE_STRUCTURE.md](DATABASE_STRUCTURE.md) for schema details

---

## Getting Help

- **Documentation**: Check the `docs/` directory
- **Examples**: See `docs/USAGE_EXAMPLES.md`
- **Tests**: Run `python tests/test_query_system.py` to verify installation
- **Issues**: See `CONTRIBUTING.md` for how to report issues

---

## Quick Reference Card

### Python API
```python
# Load database
from rel_quest import LigandDatabase
db = LigandDatabase()

# Query methods
db.query_by_name(name, fuzzy=False)          # Query by name
db.query_by_smiles(smiles)                   # Query by SMILES
db.filter_by_mechanism(mechanism)            # Filter by mechanism
db.filter_by_trend(trend)                    # Filter by trend
db.filter_by_elements(element_list)          # Filter by element coverage
db.search(query, search_fields=None)         # Multi-field search
db.get_statistics()                          # Database statistics

# LigandQuery methods
ligand.get_d_values(condition_id=None)       # Get D-values
ligand.get_log_d_values(condition_id=None)   # Get log(D) values
ligand.get_trend_summary()                   # Get trend classification
ligand.get_statistics(condition_id=None)     # Get statistics
ligand.get_separation_factor(elem1, elem2)   # Calculate SF
ligand.get_condition_list()                  # List available conditions
```

### CLI
```bash
# Query commands
python scripts/query_ligands.py --name <name>
python scripts/query_ligands.py --smiles <smiles>
python scripts/query_ligands.py --mechanism <mechanism>
python scripts/query_ligands.py --trend <trend>
python scripts/query_ligands.py --search <keyword>
python scripts/query_ligands.py --stats

# Output formats
--format raw|trend|stats|full
--output json|table|csv
```
