# logSF Viewer - Quick Start Guide

## 🚀 Launch Dashboard
```bash
cd <your checkout>   # the repository root
streamlit run app.py
```

Then navigate to: http://localhost:8501

## 📍 Location
Scroll down to: **⚖️ logSF (Separation Factor) Viewer**
(Below "🔬 Single Element Distribution Viewer")

## ⚡ Quick Start (30 seconds)

### 1. Select Elements
```
First row:  Click [Eu] (63)  → Ln1
Second row: Click [La] (57)  → Ln2
```

### 2. Choose View
```
X-Axis: ○ Ligand Types  [Select this first]
```

### 3. Read Plot
```
Left side  → Best performing categories (highest logSF)
Right side → Worst performing (lowest logSF)
y=0 line   → No selectivity

Positive values → Eu preferred
Negative values → La preferred
```

### 4. Get Details
```
Hover over box     → See statistics + top ligands
Hover over points  → See individual measurements
```

## 📊 What You See

### Default View (Eu/La by Ligand Types)
```
Expected results:
- 332 paired measurements
- 119 unique ligands
- Median logSF ≈ +0.46 (slight Eu preference)
- Top category: Diglycolamides (blue)
```

## 🎯 Key Actions

| Action | Result |
|--------|--------|
| Click element button | Change Ln1 or Ln2 |
| Click x-axis radio | Change plot type |
| Expand "Filter by Solvent" | Show/hide specific diluents |
| Hover over box | See category statistics |
| Hover over point | See measurement details |
| Expand "How to Use" | See full documentation |

## 💡 Common Tasks

### Find Best Ligands
```
1. Select your element pair
2. Choose "Ligand Types"
3. Look at leftmost category (highest median)
4. Hover over box → Note "Top ligands"
```

### Compare Solvents
```
1. Select your element pair
2. Choose "Solvent Categories"
3. Compare median values
4. Note: aliphatic vs aromatic differences
```

### Optimize Acid
```
1. Select your element pair
2. Choose "Acid Concentration"
3. Look for trend in scatter plot
4. Hover on peak points → Note [HNO₃]
```

## 📈 Interpreting logSF

```
logSF[Ln1/Ln2] = logD[Ln1] - logD[Ln2]

+2.0  █████████  Excellent Ln1 preference
+1.0  ██████     Good Ln1 preference
 0.0  ═══════    No selectivity
-1.0  ██████     Good Ln2 preference
-2.0  █████████  Excellent Ln2 preference
```

## 🔍 Example Queries

### "Which ligand type is best for Eu/La?"
```
Elements: Eu, La
X-Axis: Ligand Types
Look for: Highest median (leftmost)
Answer: Diglycolamides (median logSF ≈ +0.8)
```

### "How does n-dodecane compare to toluene?"
```
Elements: [Your choice]
X-Axis: Solvent Categories
Compare: aliphatic vs aromatic
Filter: Select only those two
```

### "What acid concentration is best for Tb/Dy?"
```
Elements: Tb, Dy
X-Axis: Acid Concentration
Look for: Peak in scatter plot
Hover: Get exact [HNO₃] values
```

## ⚠️ Troubleshooting

| Issue | Solution |
|-------|----------|
| "No paired measurements" | Try different element pair (Eu/La, Nd/Pr, Tb/Dy work well) |
| "Same element selected" | Click different element in second row |
| Plot is empty | Adjust solvent filter or choose different x-axis |
| Slow loading | First load takes ~500ms, subsequent <100ms |

## 📚 More Information

- Full guide: `LOGSF_VIEWER_USER_GUIDE.md`
- Implementation details: `LOGSF_VIEWER_IMPLEMENTATION_SUMMARY.md`
- Run tests: `python test_logsf_viewer.py`

## 🧪 Test Data Summary

Database contains **332 paired Eu/La measurements**:
- 119 unique ligands
- 41 solvents
- logSF range: -3.31 to +5.56
- Most common solvent: n-dodecane
- Best ligand type: Diglycolamides

## 🎓 Remember

1. **Positive logSF** = First element preferred
2. **Hover for details** = Most informative way to explore
3. **Box plots sorted** = Best on left, worst on right
4. **Gray line at y=0** = No selectivity baseline

---

**Ready to explore?** → `streamlit run app.py` 🚀
