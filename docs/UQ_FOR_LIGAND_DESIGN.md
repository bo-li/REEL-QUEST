# How Uncertainty Quantification Enables Rational Ligand Design

**Date**: March 4, 2026
**Key Question**: How does knowing uncertainty help us design better extractants?

---

## The Problem with Point Estimates (Current Approach)

### Scenario: Designing a New Diglycolamide

**Current approach WITHOUT uncertainty quantification:**

```
Literature search finds:
- TODGA:    slope = 0.085  → Ascending (heavy REE preference)
- TEHDGA:   slope = 0.092  → Ascending
- THDGA:    slope = 0.078  → Ascending
- TBDGA:    slope = 0.015  → Flat

Conclusion: "Alkyl chain length affects selectivity"
Design decision: "Make longer chains for better ascending trend"

Problem: Is TBDGA really Flat, or is it just poorly studied?
```

**With uncertainty quantification:**

```
- TODGA:    slope = 0.085 ± 0.006  (15 elements, 12 conditions) → HIGH confidence
- TEHDGA:   slope = 0.092 ± 0.012  (12 elements, 8 conditions)  → MEDIUM confidence
- THDGA:    slope = 0.078 ± 0.008  (14 elements, 10 conditions) → HIGH confidence
- TBDGA:    slope = 0.015 ± 0.067  (4 elements, 2 conditions)   → LOW confidence! ⚠️

New insight: TBDGA confidence interval [−0.052, +0.082] overlaps with TODGA!
  → We can't conclude alkyl chains matter
  → Need more data on TBDGA before making design decisions
```

**Outcome**: UQ **prevents false conclusions** that would waste synthesis effort!

---

## Five Ways UQ Enables Ligand Design

### 1. Identify Reliable Structure-Property Relationships

### 2. Guide Chemical Modifications with Confidence

### 3. Prioritize Synthesis Candidates (Risk-Aware Screening)

### 4. Optimize Process Conditions with Effect Size Quantification

### 5. Enable Active Learning for Efficient Exploration

Let's explore each with concrete chemistry examples:

---

## 1. Identify Reliable Structure-Property Relationships

### Question: "Do electron-withdrawing groups enhance ascending trends?"

#### Without UQ (Dangerous!)
```
Ligand Family: Aromatic heterocycles with substituents

Parent (no substituent):        slope = 0.045
With -CF₃ (electron-withdrawing): slope = 0.089
With -OCH₃ (electron-donating):  slope = 0.038

Conclusion: "EWGs enhance ascending selectivity by 2×"
Design decision: Add -CF₃ to all future designs

Chemistry reasoning:
  - EWGs reduce electron density on N-donors
  - Favors harder Lewis acids (heavy REEs have more covalent character)
  - Makes sense mechanistically! ✓
```

#### With UQ (Reveals the Truth!)
```
Parent:       slope = 0.045 ± 0.052  (3 elements measured, LOW confidence)
With -CF₃:    slope = 0.089 ± 0.018  (11 elements measured, HIGH confidence)
With -OCH₃:   slope = 0.038 ± 0.061  (4 elements measured, LOW confidence)

Statistical test (ANOVA with uncertainty):
  p-value = 0.38 → NOT statistically significant!

Why? Parent and -OCH₃ have huge uncertainties that overlap with -CF₃

New conclusion:
  "We can't conclude EWGs matter based on current data
   Need to measure more elements for Parent and -OCH₃ variants"

Revised design decision:
  1. First: Measure La, Pr, Eu, Dy, Lu for Parent and -OCH₃
  2. Then: If trend confirmed, synthesize -CF₃ variants
  3. Saves: 2-3 months synthesis time, $10K reagents
```

**Impact**: UQ prevents chasing **statistical noise** thinking it's real chemistry!

---

## 2. Guide Chemical Modifications with Confidence

### Question: "Where should I modify the ligand structure?"

#### Concrete Example: Optimizing a Diglycolamide

```
Current ligand: R-N-CO-CH₂-O-CH₂-CO-N-R
                │                       │
             Position A              Position B

Possible modifications:
  A: Change alkyl chain (R = octyl, hexyl, butyl)
  B: Add methyl on bridging carbon (-CH₂- → -CHMe-)
  C: Replace oxygen with sulfur (-O- → -S-)
  D: Add aromatic rings to R groups
```

#### Structure-Property Correlation with Confidence Intervals

**Step 1: Analyze existing data**

```python
# For each structural feature, calculate correlation with slope

Feature: Alkyl chain length
  Correlation: r = −0.15 ± 0.28  (p = 0.42)
  Confidence interval: [−0.43, +0.13]
  → NOT significant (crosses zero)

Feature: Aromatic vs aliphatic R groups
  Δslope = +0.047 ± 0.011  (p = 0.001)
  Confidence: HIGH
  → Aromatic groups increase ascending trend significantly!

Feature: Bridging atom (O vs S)
  Only 2 ligands with S (insufficient data)
  → Cannot conclude

Feature: Methyl substitution on bridge
  Δslope = +0.022 ± 0.031  (p = 0.18)
  → Marginally significant, needs more data
```

**Step 2: Design recommendations ranked by confidence**

```
Priority 1 (HIGH confidence): Add aromatic groups to R positions
  Expected effect: +0.047 ± 0.011 slope increase
  Synthesis difficulty: Medium
  Recommendation: DO THIS FIRST ✓

Priority 2 (MEDIUM confidence): Try bridging methylation
  Expected effect: +0.022 ± 0.031 (could be +0.05 or 0)
  Synthesis difficulty: Easy
  Recommendation: Test as secondary variant

Priority 3 (LOW confidence): Sulfur substitution
  Expected effect: Unknown (insufficient data)
  Synthesis difficulty: Hard
  Recommendation: SKIP for now, or measure existing S-variants first

Priority 4 (NOT significant): Alkyl chain length
  Expected effect: Essentially zero within uncertainty
  Recommendation: Keep chain length at C8 (good solubility), don't bother varying
```

**Outcome**: Focus synthesis effort on **modifications that actually work** (aromatic groups), avoid wasting time on ineffective changes (chain length).

---

## 3. Prioritize Synthesis Candidates (Risk-Aware Screening)

### Question: "I can synthesize 3 new ligands - which ones?"

#### Computational Screening with Uncertainty

You've designed 20 candidate ligands in silico. GPR model predicts slopes:

```
Candidate ligands (sorted by predicted slope):

Top predictions WITHOUT uncertainty:
1. Ligand_A:  predicted slope = 0.125  → Looks amazing! Synthesize?
2. Ligand_B:  predicted slope = 0.118
3. Ligand_C:  predicted slope = 0.112
4. Ligand_D:  predicted slope = 0.095
5. Ligand_E:  predicted slope = 0.089

Decision: Synthesize A, B, C (highest predicted slopes)
```

**BUT with prediction uncertainty:**

```
Top predictions WITH uncertainty (95% prediction intervals):

1. Ligand_A:  0.125 ± 0.072  [0.053, 0.197]  ← WIDE uncertainty!
   Why? Far from training data (novel scaffold)
   Risk: Could be anywhere from Flat to Strongly Ascending

2. Ligand_B:  0.118 ± 0.038  [0.080, 0.156]  ← Moderate uncertainty
   Why? Moderate distance from training data
   Risk: Medium

3. Ligand_C:  0.112 ± 0.021  [0.091, 0.133]  ← Narrow uncertainty
   Why? Close to training data (similar to TODGA variants)
   Risk: Low - likely to work as predicted

4. Ligand_D:  0.095 ± 0.015  [0.080, 0.110]  ← Very narrow!
   Why? Very similar to existing ligands
   Risk: Very low

5. Ligand_E:  0.089 ± 0.019  [0.070, 0.108]
```

#### Decision Matrix: Expected Value vs. Risk

```
Ligand | Best Case | Worst Case | Expected | Risk  | Synthesis | Decision
-------|-----------|------------|----------|-------|-----------|----------
A      | 0.197 (!) | 0.053      | 0.125    | HIGH  | Hard ($)  | RISKY
B      | 0.156     | 0.080      | 0.118    | MED   | Medium    | MAYBE
C      | 0.133     | 0.091      | 0.112    | LOW   | Easy      | SAFE BET ✓
D      | 0.110     | 0.080      | 0.095    | VLOW  | Easy      | SAFE ✓
E      | 0.108     | 0.070      | 0.089    | LOW   | Easy      | CONTROL ✓
```

**Smart synthesis strategy:**

```
Portfolio approach (hedge bets):

1. Synthesize C (safe bet, likely ~0.11 slope)
   → Guaranteed improvement over baseline

2. Synthesize D (very safe, likely ~0.095)
   → Backup in case C fails

3. Synthesize A (high risk, high reward)
   → If it works (0.197), huge win!
   → If it fails (0.053), we learn that scaffold doesn't work
   → Either way, reduces uncertainty for future designs

Do NOT synthesize B:
  → Medium risk, medium reward (worst of both worlds)
  → Better to take the safe bet (C) and the moonshot (A)
```

**Outcome**: **3 synthetic chemists work efficiently** - one makes safe improvements, one explores novel space, waste avoided.

---

## 4. Optimize Process Conditions with Effect Size Quantification

### Question: "Is it worth using 3M HNO₃ instead of 1M?"

#### Without UQ
```
TODGA at 1M HNO₃:  Average logSF[Eu/La] = 0.45
TODGA at 3M HNO₃:  Average logSF[Eu/La] = 0.62

Conclusion: "3M is better by 0.17 log units"
Decision: Use 3M HNO₃ in process
```

#### With UQ
```
TODGA at 1M HNO₃:  logSF[Eu/La] = 0.45 ± 0.08  (12 measurements)
TODGA at 3M HNO₃:  logSF[Eu/La] = 0.62 ± 0.15  (4 measurements)

Effect size: Δ = +0.17 ± 0.17
95% CI for difference: [0.00, +0.34]

Statistical test: p = 0.055 (marginally significant)

Interpretation:
  - Effect could be as large as +0.34 (worth it!)
  - Or as small as 0 (no benefit)
  - Need more measurements at 3M to be sure
```

**Cost-benefit analysis:**

```
Using 3M HNO₃ vs 1M:
  Benefits:
    - Possibly +0.17 better separation (if real)

  Costs:
    - 3× more acid consumption ($15K/year for pilot plant)
    - Higher corrosion (need expensive alloys)
    - More neutralization waste ($8K/year)
    - Total added cost: ~$25K/year

Decision:
  IF effect is +0.34 → Worth it (better separation justifies cost)
  IF effect is 0 → NOT worth it (pure waste)

Action:
  1. Measure logSF at 3M HNO₃ more carefully (6 replicates)
  2. Narrow uncertainty to ±0.05
  3. THEN decide if 3M is economically justified

Outcome: Don't commit to expensive process change on marginal data!
```

---

## 5. Enable Active Learning for Efficient Exploration

### Question: "I have budget for 10 experiments - which ligands and elements?"

This is where UQ becomes incredibly powerful!

#### Traditional Approach (No UQ)
```
Random selection:
  - Measure new ligand X at all 15 elements
  - Measure new ligand Y at all 15 elements
  - Total: 30 measurements (over budget)

OR:

Complete random coverage:
  - Pick 10 random (ligand, element) combinations
  - Hope for the best
```

#### Active Learning Approach (With UQ)
```
Current database state:
  - 176 ligands
  - Uncertainty quantified for each

Step 1: Identify high-impact measurements
  For each unmeasured (ligand, element) pair:
    Calculate: Expected Information Gain (EIG)

  EIG = How much would this measurement reduce overall uncertainty?

Step 2: Rank all possible measurements

Top 10 by EIG:

1. Ligand_047, Dy (EIG = 0.089)
   Why? Only has La, Ce, Eu, Lu → huge gap at Z=66
   Measuring Dy fills the gap, reduces slope uncertainty by 45%

2. Ligand_103, Pr (EIG = 0.076)
   Why? Has La, Nd, Eu, Lu → missing Pr
   Pr is critical for Nd/Pr separation (industrial target)

3. Ligand_047, Er (EIG = 0.071)
   Why? Same ligand as #1, filling another gap

4. Ligand_128, Tm (EIG = 0.063)
   Why? Novel scaffold, only 3 elements measured
   High uncertainty, Tm would help classify trend

5. Ligand_103, Dy (EIG = 0.061)
   Why? Pairs with #2, enables Tb/Dy separation analysis

... (6-10 continue similarly)
```

**Comparison:**

```
Random 10 measurements:
  Expected uncertainty reduction: ~8%
  Expected new knowledge: Low (random coverage)

Active learning 10 measurements:
  Expected uncertainty reduction: ~34%
  Expected new knowledge: High (fills critical gaps)

Result: 4× more efficient use of experimental resources!
```

**Practical outcome:**

```
After measuring the top 10:

Before:
  - 47 ligands with "Insufficient_Data"
  - 89 ligands with "Low" or "Medium" confidence

After:
  - 35 ligands with "Insufficient_Data" (−12)
  - 73 ligands with "Low/Medium" confidence (−16)
  - 28 new "High" confidence ligands (+16)

Now we can:
  ✓ Confidently identify structural trends (more high-quality data)
  ✓ Train better ML models (higher quality training set)
  ✓ Make safer process recommendations (less uncertainty)
```

---

## Real-World Ligand Design Workflow with UQ

### Scenario: Design a ligand for Nd/Pr separation (magnets application)

**Goal**: Maximize logSF[Nd/Pr] with confidence

#### Step 1: Query Database with Uncertainty Filter
```python
from request import LigandDatabase

db = LigandDatabase()

# Find ligands with HIGH confidence Ascending trends
candidates = db.query(
    trend="Ascending",
    confidence="High",
    min_R2=0.90,
    has_Nd_and_Pr=True
)

Results:
  23 ligands match criteria

Top 5 by logSF[Nd/Pr]:
1. CyMe4BTBP:  logSF = 0.28 ± 0.04  (HIGH confidence)
2. TODGA:      logSF = 0.24 ± 0.03  (HIGH confidence)
3. TEHDGA:     logSF = 0.26 ± 0.07  (MEDIUM confidence)
4. Ligand_X:   logSF = 0.31 ± 0.12  (LOW confidence - only 2 measurements!)
5. THDGA:      logSF = 0.22 ± 0.03  (HIGH confidence)
```

#### Step 2: Decide What to Do

```
Option A: Use CyMe4BTBP (safe bet)
  - logSF = 0.28 ± 0.04 is reliable
  - Well-studied (15 elements, 20 conditions)
  - Can scale up with confidence
  - Conservative choice ✓

Option B: Investigate Ligand_X (high risk/reward)
  - logSF could be 0.19 or 0.43 (huge range!)
  - If 0.43 → 50% better than CyMe4BTBP!
  - But only 2 measurements (not reliable)
  - Action: Measure 5 more elements first
  - THEN decide if worth scaling up

Option C: Modify CyMe4BTBP (guided optimization)
  - Start from HIGH confidence baseline (0.28)
  - Make systematic modifications
  - Each variant: predict slope ± uncertainty with GPR
  - Synthesize only high-confidence improvements
```

#### Step 3: Structure-Property Analysis Guides Modifications

```python
# Analyze what structural features correlate with Nd/Pr selectivity

features = db.analyze_structural_correlations(
    target="logSF[Nd/Pr]",
    min_confidence="Medium"
)

Results (with 95% CI):

Feature: Aromatic heterocycle (BTP, BTBP)
  Effect: +0.15 ± 0.05 logSF increase
  p-value: 0.002
  Confidence: HIGH ✓✓✓
  → Keep aromatic heterocycle core

Feature: Methyl substitution on rings
  Effect: +0.08 ± 0.09 logSF increase
  p-value: 0.18
  Confidence: LOW
  → Not conclusive, but worth trying

Feature: Alkyl chain length (C4 vs C8)
  Effect: −0.02 ± 0.11 logSF
  p-value: 0.76
  Confidence: NONE
  → No effect, don't waste time varying

Feature: Phenanthroline vs bipyridine core
  Effect: +0.06 ± 0.12
  p-value: 0.34
  Confidence: LOW
  → Unclear, needs more data
```

#### Step 4: Design Targets with Expected Performance

```
Based on high-confidence features:

Target 1: CyMe4BTBP-derivative with extra methyl groups
  Predicted logSF[Nd/Pr]: 0.28 + 0.08 = 0.36
  Uncertainty: ±0.10 (cumulative from model)
  95% CI: [0.26, 0.46]
  Risk: Medium
  Synthesis: 2 weeks, $2K

Target 2: Novel BTP scaffold (active learning suggested this)
  Predicted logSF[Nd/Pr]: 0.41
  Uncertainty: ±0.18 (far from training data)
  95% CI: [0.23, 0.59]
  Risk: High
  Synthesis: 6 weeks, $8K

Target 3: BTBP with fluorinated arms (EWG effect)
  Predicted logSF[Nd/Pr]: 0.33
  Uncertainty: ±0.07
  95% CI: [0.26, 0.40]
  Risk: Low-Medium
  Synthesis: 3 weeks, $3K

Decision matrix:

Make Target 1 + Target 3:
  → Target 1: Fast, cheap, medium upside (0.36)
  → Target 3: Moderate time/cost, likely improvement (0.33)
  → Skip Target 2: Too expensive for the uncertainty
  → After measuring 1 and 3: Use data to improve model
     Then revisit Target 2-type designs with lower uncertainty
```

---

## Summary: How UQ Transforms Ligand Design

### Without UQ (Old Way)
❌ "Ligand A has slope 0.12, Ligand B has 0.08 → A is better"
❌ "Add -CF₃ groups because one paper showed improvement"
❌ "Let's try these 20 random variants and see what happens"
❌ "3M acid seems better, change the whole process"
❌ Waste 6 months synthesizing ligands that don't work

### With UQ (RE-QUEST Way)
✅ "Ligand A: 0.12 ± 0.02 vs B: 0.08 ± 0.09 → A is significantly better (p<0.001)"
✅ "EWG effect: +0.04 ± 0.03 (p=0.04) → Real but small, calculate if worth synthesis cost"
✅ "Active learning says: measure these 10 specific (ligand, element) pairs for max info gain"
✅ "3M acid: +0.17 ± 0.17 → Uncertain, measure 6 replicates before committing to expensive process change"
✅ Focus synthesis on high-confidence predictions, explore high-uncertainty regions systematically

---

## Quantitative Impact

### Example: 3-Year PhD Student Project

**Scenario**: Design better extractant for Tb/Dy separation

**Without UQ:**
```
Year 1: Literature review, synthesize 10 variants based on intuition
  → 3 show improvement, 7 fail
  → Don't know why failures occurred
  → No systematic understanding

Year 2: Try 10 more variants (different strategy)
  → 4 show improvement
  → Still no clear structure-property relationships
  → Chasing random success

Year 3: Scale up the "best" one
  → Works in lab but fails at pilot (conditions weren't optimal)
  → Scramble to fix

Result: 1 marginal improvement, no predictive model
```

**With UQ (RE-QUEST-guided):**
```
Year 1:
  - Query RE-QUEST: Find 15 high-confidence ligands for Tb/Dy
  - Identify structural features with p<0.05
  - Active learning: Measure 20 strategic (ligand, element) pairs
  - Build GPR model with 95% prediction intervals

  → Clear understanding of what works and why
  → Confidence intervals guide synthesis priorities

Year 2:
  - Synthesize 5 high-confidence predictions (narrow CI)
  - Synthesize 2 high-risk/reward predictions (wide CI)
  - 4/5 high-confidence work as predicted ✓
  - 1/2 high-risk succeeds → New best performer!

  → Systematic improvement
  → Model validated and refined

Year 3:
  - Optimize process conditions with uncertainty quantification
  - Know exactly how acid/temp affect logSF[Tb/Dy] ± CI
  - Scale up with confidence (no surprises)

  → Pilot plant works first try
  → Publishable structure-property relationships
  → Predictive model for future designs

Result: 5 working improvements, validated predictive model, successful scale-up
```

**Efficiency Gain**: ~3-4× faster, ~50% lower reagent costs, higher success rate

---

## Bottom Line

### Uncertainty Quantification enables rational ligand design by:

1. **Preventing false positives**: Don't chase noise thinking it's real chemistry
2. **Guiding modifications**: Focus on structural changes that reliably work
3. **Risk-aware screening**: Synthesize the right mix of safe bets and moonshots
4. **Process optimization**: Don't commit to expensive changes on weak evidence
5. **Active learning**: Measure the most informative experiments, not random ones

### The workflow changes from:
**"Try stuff and hope" → "Design with confidence, explore with purpose"**

### This is why RE-QUEST's UQ capabilities are transformative:
Not just a database - it's a **decision engine** that tells you:
- What we know (with confidence intervals)
- What we don't know (where uncertainty is high)
- What to do next (active learning recommendations)

**Result**: Faster discovery, fewer failed syntheses, predictable scale-up, publishable science.

---

Would you like me to create:
1. A **worked example** showing the full workflow for designing a Nd/Pr extractant?
2. A **decision tree** for "given my data quality, what should I do next?"
3. A **cost-benefit calculator** for synthesis vs. measurement decisions?
