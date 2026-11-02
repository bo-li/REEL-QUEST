# REEL-QUEST Usage Examples

This guide demonstrates how to use the REEL-QUEST ligand database for querying rare earth extraction data.

## Table of Contents
1. [Python API Examples](#python-api-examples)
2. [CLI Examples](#cli-examples)
3. [Common Use Cases](#common-use-cases)
4. [Data Analysis Examples](#data-analysis-examples)

---

## Python API Examples

### Basic Setup

```python
from rel_quest import LigandDatabase

# Initialize database (uses merged database by default)
db = LigandDatabase()

# Or specify database path explicitly
db = LigandDatabase('data/sx_ligands_merged.jsonl')
```

### Query by Name

```python
# Exact name match
results = db.query_by_name('TODGA')
if results:
    todga = results[0]
    print(f"Mechanism: {todga.mechanism}")
    print(f"Trend: {todga.selectivity_trend}")

# Fuzzy search (substring matching)
dga_ligands = db.query_by_name('DGA', fuzzy=True)
print(f"Found {len(dga_ligands)} ligands containing 'DGA'")
```

### Query by SMILES

```python
smiles = "CCCCCCCCN(CCCCCCCC)C(=O)COCC(=O)N(CCCCCCCC)CCCCCCCC"
result = db.query_by_smiles(smiles)
if result:
    print(f"Found: {result.name}")
```

### Filter by Mechanism

```python
# Get all solvation extractants
solvation_ligands = db.filter_by_mechanism('Solvation')
print(f"Solvation extractants: {len(solvation_ligands)}")

# Get all cation exchange ligands
cation_ligands = db.filter_by_mechanism('Cation_Exchange')
print(f"Cation exchange ligands: {len(cation_ligands)}")
```

### Filter by Selectivity Trend

```python
# Find ligands with ascending trend (La→Lu increasing)
ascending = db.filter_by_trend('Ascending')
print(f"Ascending trend ligands: {len(ascending)}")

# Find ligands with reverse trend (Lu→La increasing)
reverse = db.filter_by_trend('Reverse')
print(f"Reverse trend ligands: {len(reverse)}")
```

### Filter by Element Coverage

```python
# Find ligands tested on specific elements
required_elements = ['La', 'Nd', 'Eu', 'Lu']
results = db.filter_by_elements(required_elements)
print(f"Ligands with data for {required_elements}: {len(results)}")
```

### Search Across Multiple Fields

```python
# Search by keyword in name, SMILES, or structural class
results = db.search('phosph')
print(f"Found {len(results)} ligands containing 'phosph'")

# Search specific fields only
results = db.search('BTP', search_fields=['name', 'structural_class'])
```

---

## Working with Query Results

### Get Distribution Ratios (D-values)

```python
ligand = db.query_by_name('TBP')[0]

# Get D-values as dictionary
d_values = ligand.get_d_values()
for element, d_val in sorted(d_values.items()):
    print(f"{element}: {d_val:.3f}")

# Get log(D) values
log_d_values = ligand.get_log_d_values()
```

### Get Selectivity Trend Summary

```python
ligand = db.query_by_name('TODGA')[0]

trend = ligand.get_trend_summary()
print(f"Classification: {trend['classification']}")
print(f"Confidence: {trend['confidence']}")
print(f"Direction: {trend['direction']}")
print(f"Elements tested: {', '.join(trend['elements_tested'])}")
```

### Get Statistical Summary

```python
ligand = db.query_by_name('TBP')[0]

stats = ligand.get_statistics()
print(f"D_min: {stats['D_min']:.3f}")
print(f"D_max: {stats['D_max']:.3f}")
print(f"D_mean: {stats['D_mean']:.3f}")
print(f"D_std: {stats['D_std']:.3f}")
print(f"Elements tested: {stats['element_count']}")
```

### Calculate Separation Factors

```python
ligand = db.query_by_name('TBP')[0]

# Calculate SF between two elements
sf_la_lu = ligand.get_separation_factor('La', 'Lu')
print(f"SF(La/Lu) = {sf_la_lu:.2f}")

# Calculate multiple separation factors
elements = ['La', 'Nd', 'Eu', 'Dy', 'Lu']
d_values = ligand.get_d_values()

print("Separation factors from La:")
for elem in elements[1:]:
    sf = ligand.get_separation_factor('La', elem)
    if sf:
        print(f"  SF(La/{elem}) = {sf:.2f}")
```

### Get Database Statistics

```python
stats = db.get_statistics()

print(f"Total ligands: {stats['total_ligands']}")
print("\nMechanisms:")
for mech, count in stats['mechanisms'].items():
    print(f"  {mech}: {count}")

print("\nSelectivity Trends:")
for trend, count in stats['trends'].items():
    print(f"  {trend}: {count}")
```

---

## CLI Examples

### Basic Queries

```bash
# Query by name with full output (JSON)
python scripts/query_ligands.py --name TODGA --format full

# Query with table output
python scripts/query_ligands.py --name TODGA --format full --output table

# Get raw D-values only
python scripts/query_ligands.py --name TODGA --format raw
```

### Filter Operations

```bash
# Find all ascending trend ligands
python scripts/query_ligands.py --trend Ascending --format trend

# Find all solvation extractants
python scripts/query_ligands.py --mechanism Solvation --format trend

# Search for ligands containing keyword
python scripts/query_ligands.py --search DGA
```

### Output Formats

```bash
# JSON output (default)
python scripts/query_ligands.py --name TBP --format full --output json

# Table output (human-readable)
python scripts/query_ligands.py --name TBP --format full --output table

# CSV output (for Excel/plotting)
python scripts/query_ligands.py --name TBP --format raw --output csv
```

### Database Statistics

```bash
# Show database statistics
python scripts/query_ligands.py --stats

# List all ligands
python scripts/query_ligands.py --list-all --format trend --output table
```

---

## Common Use Cases

### Use Case 1: Find Best Ligand for La/Lu Separation

```python
from rel_quest import LigandDatabase

db = LigandDatabase()

# Get all ligands with data for La and Lu
candidates = db.filter_by_elements(['La', 'Lu'])

# Calculate separation factors
sf_data = []
for ligand in candidates:
    sf = ligand.get_separation_factor('La', 'Lu')
    if sf:
        sf_data.append({
            'name': ligand.name,
            'mechanism': ligand.mechanism,
            'SF_La_Lu': sf,
            'trend': ligand.selectivity_trend
        })

# Sort by separation factor
sf_data.sort(key=lambda x: x['SF_La_Lu'], reverse=True)

# Print top 10
print("Top 10 ligands for La/Lu separation:")
for i, data in enumerate(sf_data[:10], 1):
    print(f"{i}. {data['name']} (SF={data['SF_La_Lu']:.2f}, {data['mechanism']})")
```

### Use Case 2: Compare DGA Ligands

```python
from rel_quest import LigandDatabase

db = LigandDatabase()

# Get all DGA ligands
dga_ligands = db.search('DGA')

print(f"Found {len(dga_ligands)} DGA ligands\n")

for ligand in dga_ligands:
    trend = ligand.get_trend_summary()
    stats = ligand.get_statistics()

    print(f"{ligand.name}:")
    print(f"  Trend: {trend['classification']} ({trend['confidence']} confidence)")
    print(f"  Elements tested: {trend['coverage_count']}")
    print(f"  D_mean: {stats.get('D_mean', 'N/A'):.3f}")
    print()
```

### Use Case 3: Find Ligands for Specific Element Pair

```python
from rel_quest import LigandDatabase

db = LigandDatabase()

# Target: Nd/Dy separation (important for NdFeB magnets)
target_elements = ['Nd', 'Dy']

# Get ligands with data for both elements
candidates = db.filter_by_elements(target_elements)

print(f"Ligands with Nd and Dy data: {len(candidates)}\n")

# Analyze separation factors
for ligand in candidates[:10]:  # Show first 10
    sf = ligand.get_separation_factor('Nd', 'Dy')
    d_values = ligand.get_d_values()

    print(f"{ligand.name}:")
    print(f"  D(Nd) = {d_values['Nd']:.3f}")
    print(f"  D(Dy) = {d_values['Dy']:.3f}")
    print(f"  SF(Nd/Dy) = {sf:.2f}")
    print(f"  Mechanism: {ligand.mechanism}")
    print()
```

### Use Case 4: Generate Dataset for Machine Learning

```python
from rel_quest import LigandDatabase
import pandas as pd

db = LigandDatabase()

# Extract features for all ligands
data = []
for ligand in db.get_all_ligands():
    # Get statistical features
    stats = ligand.get_statistics()

    # Build feature vector
    row = {
        'name': ligand.name,
        'mechanism': ligand.mechanism,
        'extractant_type': ligand.extractant_type,
        'trend': ligand.selectivity_trend,
        'D_mean': stats.get('D_mean'),
        'D_std': stats.get('D_std'),
        'log_D_mean': stats.get('log_D_mean'),
        'log_D_std': stats.get('log_D_std'),
        'element_count': stats.get('element_count'),
        'smiles': ligand.smiles
    }
    data.append(row)

# Create DataFrame
df = pd.DataFrame(data)

# Save to CSV
df.to_csv('ligand_features.csv', index=False)
print(f"Exported {len(df)} ligands to ligand_features.csv")
```

---

## Data Analysis Examples

### Example 1: Plot D-value Trends

```python
import matplotlib.pyplot as plt
from rel_quest import LigandDatabase

db = LigandDatabase()

# Get ligand
ligand = db.query_by_name('TBP')[0]

# Get atomic numbers and D-values
lanthanides = ['La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu']
atomic_numbers = list(range(57, 72))

d_values = ligand.get_d_values()
log_d_values = ligand.get_log_d_values()

# Extract data for plotting
elements = []
atomic_nums = []
log_d_vals = []

for elem, atomic_num in zip(lanthanides, atomic_numbers):
    if elem in log_d_values:
        elements.append(elem)
        atomic_nums.append(atomic_num)
        log_d_vals.append(log_d_values[elem])

# Plot
plt.figure(figsize=(10, 6))
plt.plot(atomic_nums, log_d_vals, 'o-', linewidth=2, markersize=8)
plt.xlabel('Atomic Number', fontsize=12)
plt.ylabel('log(D)', fontsize=12)
plt.title(f'{ligand.name} - Selectivity Trend: {ligand.selectivity_trend}', fontsize=14)
plt.grid(True, alpha=0.3)
plt.xticks(atomic_nums, elements, rotation=45)
plt.tight_layout()
plt.savefig('tbp_selectivity_trend.png', dpi=300)
print("Plot saved to tbp_selectivity_trend.png")
```

### Example 2: Compare Multiple Ligands

```python
import matplotlib.pyplot as plt
from rel_quest import LigandDatabase

db = LigandDatabase()

# Get DGA ligands with ascending trend
ligands_to_compare = ['TODGA', 'DMDODGA', 'TEHDGA']

plt.figure(figsize=(12, 6))

for ligand_name in ligands_to_compare:
    results = db.query_by_name(ligand_name)
    if not results:
        continue

    ligand = results[0]
    log_d_values = ligand.get_log_d_values()

    lanthanides = ['La', 'Ce', 'Pr', 'Nd', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu']
    atomic_numbers = list(range(57, 72))

    # Extract data
    atomic_nums = []
    log_d_vals = []
    for elem, atomic_num in zip(lanthanides, atomic_numbers):
        if elem in log_d_values:
            atomic_nums.append(atomic_num)
            log_d_vals.append(log_d_values[elem])

    plt.plot(atomic_nums, log_d_vals, 'o-', label=ligand_name, linewidth=2, markersize=6)

plt.xlabel('Atomic Number', fontsize=12)
plt.ylabel('log(D)', fontsize=12)
plt.title('Comparison of DGA Ligands', fontsize=14)
plt.legend(fontsize=10)
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig('dga_comparison.png', dpi=300)
print("Plot saved to dga_comparison.png")
```

### Example 3: Mechanism Distribution Analysis

```python
import matplotlib.pyplot as plt
from rel_quest import LigandDatabase

db = LigandDatabase()

stats = db.get_statistics()

# Create pie chart
mechanisms = stats['mechanisms']
labels = list(mechanisms.keys())
sizes = list(mechanisms.values())

plt.figure(figsize=(10, 6))
plt.pie(sizes, labels=labels, autopct='%1.1f%%', startangle=90)
plt.title(f'Extraction Mechanisms (Total: {sum(sizes)} ligands)', fontsize=14)
plt.axis('equal')
plt.tight_layout()
plt.savefig('mechanism_distribution.png', dpi=300)
print("Plot saved to mechanism_distribution.png")
```

---

## Advanced Queries

### Find Ligands with Specific Properties

```python
from rel_quest import LigandDatabase

db = LigandDatabase()

# Find solvation ligands with ascending trend
results = []
for ligand in db.filter_by_mechanism('Solvation'):
    if ligand.selectivity_trend == 'Ascending':
        results.append(ligand)

print(f"Solvation ligands with ascending trend: {len(results)}")
```

### Get All Experimental Conditions

```python
ligand = db.query_by_name('TODGA')[0]

conditions = ligand.get_condition_list()
print(f"Available conditions for {ligand.name}:")
for cond_id in conditions:
    d_vals = ligand.get_d_values(cond_id)
    stats = ligand.get_statistics(cond_id)
    print(f"  {cond_id}:")
    print(f"    Elements: {len(d_vals)}")
    print(f"    D_mean: {stats['D_mean']:.3f}")
```

---

## Tips and Best Practices

1. **Use fuzzy search for exploratory analysis**: `query_by_name('DGA', fuzzy=True)` helps find related ligands

2. **Filter by elements for specific applications**: If you need separation of Nd/Dy for magnets, use `filter_by_elements(['Nd', 'Dy'])`

3. **Check where a value comes from**: every condition cites the paper that measured it, `condition['reference']['DOI']`, and where in the paper the values are, `condition['reference']['Locator']`

4. **SMILES**: every record carries one; the SMILES, not the name, is the structure

5. **Use statistics for quick comparisons**: `get_statistics()` provides D_mean, D_std for rapid screening

6. **Combine filters**: Chain multiple filters for targeted searches
   ```python
   # Find all solvation ligands with ascending trend and La/Lu data
   results = []
   for lig in db.filter_by_mechanism('Solvation'):
       if lig.selectivity_trend == 'Ascending' and 'La' in lig.element_coverage and 'Lu' in lig.element_coverage:
           results.append(lig)
   ```

---

For more information, see:
- [API Documentation](API.md)
- [CLI Documentation](CLI.md)
- [Database Structure](DATABASE_STRUCTURE.md)
