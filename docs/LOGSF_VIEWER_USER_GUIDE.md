# logSF Viewer - Quick User Guide

## What is logSF?

**logSF** (log Separation Factor) quantifies how well a ligand separates two lanthanides:

```
logSF[Ln1/Ln2] = logD[Ln1] - logD[Ln2]
```

- **Positive value**: Ligand prefers Ln1 (extracts more Ln1 to organic phase)
- **Negative value**: Ligand prefers Ln2 (extracts more Ln2 to organic phase)
- **Zero value**: No selectivity (equal extraction of both elements)

## How to Access

1. Start dashboard: `streamlit run app.py`
2. Navigate to: **⚖️ logSF (Separation Factor) Viewer** section
3. Located below "🔬 Single Element Distribution Viewer"

## Interface Layout

```
⚖️ logSF (Separation Factor) Viewer
────────────────────────────────────────────────────────────

[Brief explanation of logSF]

Select First Element (Ln1):
┌────┬────┬────┬────┬────┬────┬────┬────┬────┬────┬────┬────┬────┬────┬────┐
│ 57 │ 58 │ 59 │ 60 │ 61 │ 62 │ 63 │ 64 │ 65 │ 66 │ 67 │ 68 │ 69 │ 70 │ 71 │
│ La │ Ce │ Pr │ Nd │ Pm │ Sm │ Eu │ Gd │ Tb │ Dy │ Ho │ Er │ Tm │ Yb │ Lu │
└────┴────┴────┴────┴────┴────┴────┴────┴────┴────┴────┴────┴────┴────┴────┘
Ln1: Europium (Eu)

Select Second Element (Ln2):
┌────┬────┬────┬────┬────┬────┬────┬────┬────┬────┬────┬────┬────┬────┬────┐
│ 57 │ 58 │ 59 │ 60 │ 61 │ 62 │ 63 │ 64 │ 65 │ 66 │ 67 │ 68 │ 69 │ 70 │ 71 │
│ La │ Ce │ Pr │ Nd │ Pm │ Sm │ Eu │ Gd │ Tb │ Dy │ Ho │ Er │ Tm │ Yb │ Lu │
└────┴────┴────┴────┴────┴────┴────┴────┴────┴────┴────┴────┴────┴────┴────┘
Ln2: Lanthanum (La)

Current Separation: logSF[Eu/La] = logD[Eu] - logD[La]

Visualization Options:
X-Axis: ○ Ligand Types  ○ Solvents  ○ Solvent Categories  ○ Ligand Acidity
        ○ Acid Concentration  ○ Extractant Concentration  ○ Temperature

Filters:
▼ 🧪 Filter by Solvent

[Interactive Plot Displayed Here]
- Y-axis: logSF[Eu/La]
- X-axis: [Selected option]
- Gray dashed line at y=0 (no selectivity)

▼ ℹ️ How to Use This Viewer
```

## Step-by-Step Workflow

### Example: Find Best Ligands for Eu/La Separation

**Goal:** Identify which ligand types show highest Eu selectivity

**Steps:**

1. **Select Elements**
   - Click **Eu** in first row (Ln1) → Button turns blue
   - Click **La** in second row (Ln2) → Button turns blue
   - Verify: "logSF[Eu/La] = logD[Eu] - logD[La]"

2. **Choose Visualization**
   - Click **Ligand Types** radio button
   - Plot updates automatically

3. **Interpret Results**
   - Look for categories on **LEFT side** (highest median logSF)
   - **Positive logSF = Eu preference**
   - Hover over boxes to see:
     - Median, mean, quartiles
     - Number of measurements
     - Top 3 ligands in that category

4. **Explore Details**
   - Hover over **outlier points** (individual measurements)
   - See: ligand name, logD[Eu], logD[La], solvent, acid conc, temp
   - Note down promising ligands

5. **Filter by Solvent** (optional)
   - Expand "🧪 Filter by Solvent"
   - Select only "n-dodecane"
   - See how performance changes in standard diluent

**Expected Result:**
- **Diglycolamides** (blue) show highest Eu selectivity
- Median logSF ≈ +0.5 to +1.5
- Top ligands: TODGA, TEHDGA, T2EHDGA

### Example: Examine Solvent Effects on Nd/Pr

**Goal:** See how organic diluent affects Nd/Pr separation

**Steps:**

1. Select **Nd** (Ln1) and **Pr** (Ln2)
2. Choose **Solvent Categories** x-axis
3. Compare median logSF across:
   - aliphatic (n-dodecane)
   - aromatic (toluene)
   - chlorinated_hydrocarbon (chloroform)
   - alcohol (1-octanol)

4. Expand "ℹ️ Solvent Category Definitions" to understand classifications

**Expected Result:**
- Different solvent polarity affects selectivity
- Aliphatic solvents: moderate selectivity
- Polar solvents: enhanced or reversed selectivity

### Example: Optimize Acid Concentration

**Goal:** Find optimal HNO₃ molarity for Tb/Dy separation

**Steps:**

1. Select **Tb** (Ln1) and **Dy** (Ln2)
2. Choose **Acid Concentration** x-axis
3. Observe scatter plot trend
4. Hover over points to see:
   - Ligand name
   - Exact acid concentration
   - logSF value

**Expected Result:**
- Scatter plot may show trend: logSF vs [HNO₃]
- Identify optimal range (e.g., 1-3 M)
- Some ligands perform better at high/low acidity

## Interpretation Guide

### Reading Box + Violin Plots

```
        ╭─────╮
   3    │  ●  │        ← Outlier (exceptional performance)
        │     │
   2    ├─────┤        ← 75th percentile (Q3)
        │█████│
   1    │█████│        ← Median (typical performance)
        │█████│
   0 ───┼─────┼─── ─ ─  ← No selectivity (logSF = 0)
        │█████│
  -1    ├─────┤        ← 25th percentile (Q1)
        │     │
  -2    ╰─────╯
     Diglycolamide
```

- **Violin (colored background)**: Probability density
  - Wide = many measurements at that logSF
  - Narrow = few measurements at that logSF
- **Box (white overlay)**: Statistical summary
  - Thick line = median
  - Box edges = Q1 and Q3 (interquartile range)
  - Whiskers = min/max (excluding outliers)
  - Dots = outliers (>1.5 IQR from box)

### Reading Scatter Plots

```
   3    ●
        ●  ●
   2       ●  ●
           ● ●  ●
   1         ●    ●
              ●  ●  ●
   0 ─ ─ ─ ─ ─ ─ ─ ─ ─  ← No selectivity
                ●
  -1              ●

     0   1   2   3   4   5
     Acid Concentration (M)
```

- Each dot = one experimental condition
- Color = solvent type (if ≤10 unique solvents)
- Trend line (if present) = linear regression

### logSF Magnitude Guide

| logSF Range | Selectivity | Example |
|-------------|-------------|---------|
| > +2.0 | Excellent Ln1 preference | α > 100 |
| +1.0 to +2.0 | Good Ln1 preference | α = 10-100 |
| +0.3 to +1.0 | Moderate Ln1 preference | α = 2-10 |
| -0.3 to +0.3 | Poor selectivity | α ≈ 1 |
| -1.0 to -0.3 | Moderate Ln2 preference | α = 0.5-0.1 |
| -2.0 to -1.0 | Good Ln2 preference | α = 0.1-0.01 |
| < -2.0 | Excellent Ln2 preference | α < 0.01 |

Where α = SF = D[Ln1]/D[Ln2] = 10^logSF

## Common Use Cases

### 1. Ligand Discovery
**Task:** Find new ligand candidates for specific separation

**Workflow:**
- Choose target element pair
- Use "Ligand Types" x-axis
- Identify top-performing categories
- Expand "How to Use" and note ligand names
- Look up structures in Ligand Browser section

### 2. Process Optimization
**Task:** Optimize extraction conditions for known ligand

**Workflow:**
- Choose element pair
- Try "Acid Concentration" x-axis → Find optimal [HNO₃]
- Try "Extractant Concentration" → Find optimal [L]org
- Try "Temperature" → Check thermal stability
- Filter by specific solvent of interest

### 3. Solvent Selection
**Task:** Choose best diluent for extraction process

**Workflow:**
- Choose element pair
- Use "Solvent Categories" x-axis
- Compare median logSF across categories
- Filter by category (expand, select specific solvents)
- Consider cost, safety, environmental impact

### 4. Mechanism Investigation
**Task:** Understand how ligand acidity affects selectivity

**Workflow:**
- Choose element pair (e.g., Eu/Am for actinide/lanthanide)
- Use "Ligand Acidity" x-axis
- Compare acidic (cation exchange) vs neutral (solvation)
- Expand "ℹ️ Ligand Acidity Classification" for details

### 5. Comparative Analysis
**Task:** Compare selectivity for different element pairs

**Workflow:**
1. Eu/La → Note median logSF for top ligand type
2. Switch to Nd/Pr → Compare same ligand type
3. Switch to Tb/Dy → Compare same ligand type
4. Identify ligands with consistent high performance

## Tips and Tricks

### Navigation
- ⬆️ Scroll to see all lanthanide elements in periodic table
- 🔄 Click element again to verify selection (button stays blue)
- 🔍 Hover slowly over plots to see detailed tooltips

### Filtering
- 🧪 Start with all solvents, then narrow down
- 🔬 Try filtering to single solvent (e.g., "n-dodecane") to compare ligands in standard conditions
- ⚗️ Reset filter by selecting all solvents again

### Data Quality
- ✅ Look for **N > 10** in hover tooltips (reliable statistics)
- ⚠️ Be cautious with **N < 5** (limited data)
- 📊 Box plots sorted left-to-right by median (best on left)

### Performance
- ⚡ First plot load: ~200-500ms
- ⚡ Element switching: <100ms
- ⚡ Filter updates: ~50-200ms
- 💾 Database cached in memory (fast subsequent queries)

### Troubleshooting
- **"No paired measurements found"**: Try different element pair
- **Same element warning**: Select two different elements
- **Empty plot**: Adjust solvent filter or choose different x-axis
- **Slow rendering**: Close other browser tabs, restart Streamlit

## Advanced Features

### Comparing Heavy vs Light REE
```
Light REE preference (Reverse trend):
- Select La (Ln1) and Lu (Ln2)
- Positive logSF = Light REE selective ligands

Heavy REE preference (Ascending trend):
- Select Lu (Ln1) and La (Ln2)
- Positive logSF = Heavy REE selective ligands
```

### Identifying Universal Ligands
```
Test multiple pairs with same ligand type:
1. Eu/La (light vs middle)
2. Dy/Tb (middle vs middle)
3. Yb/Er (heavy vs heavy)

Look for ligand types that show:
- Consistent positive logSF → Heavy REE selective
- Consistent negative logSF → Light REE selective
- Near-zero logSF → Non-selective (good for group separation)
```

### Statistical Analysis
```
Hover over box plot → See:
- Median (typical performance)
- Mean (average including outliers)
- Q1/Q3 (spread of data)
- Min/Max (range)
- Top ligands (best performers)

Use this to assess:
- Reliability (narrow IQR = consistent performance)
- Potential (outliers show best-case results)
```

## Frequently Asked Questions

**Q: What does "paired measurements" mean?**
A: Only conditions where BOTH elements were measured in the same experiment are included. This ensures fair comparison.

**Q: Why are some element pairs empty?**
A: Not all combinations have been studied. Try common pairs: Eu/La, Nd/Pr, Tb/Dy, Er/Ho, Yb/Lu.

**Q: Can I see which specific ligands are in a category?**
A: Hover over box plot → "Top ligands" shows best 3. For full list, go to Ligand Browser section and filter by extractant type.

**Q: What's the difference between "Solvents" and "Solvent Categories"?**
A: "Solvents" shows all ~40 individual diluents. "Solvent Categories" groups them into ~10 chemical classes for easier pattern recognition.

**Q: How do I export the data?**
A: Currently, use screenshot tool or Plotly download button (camera icon). CSV export feature planned for future release.

**Q: Can I compare more than two elements at once?**
A: Not currently. logSF is a pairwise metric. For multi-element analysis, use Single Element Distribution Viewer to compare logD values side-by-side.

## Related Sections

- **🔬 Single Element Distribution Viewer**: View logD for individual elements
- **📊 Eu/La Separation Analysis**: Historical Eu/La-specific analysis (aggregated data)
- **🔍 Ligand Browser**: Search and filter ligands, see detailed properties
- **📈 Database Overview**: Statistics and pie charts

## Getting Help

- **In-dashboard help**: Expand "ℹ️ How to Use This Viewer"
- **Technical documentation**: See `LOGSF_VIEWER_IMPLEMENTATION_SUMMARY.md`
- **Test examples**: Run `python test_logsf_viewer.py`

---

**Last updated:** August 13, 2026
**Database version:** v13 (166 ligands, 332 Eu/La paired measurements)
