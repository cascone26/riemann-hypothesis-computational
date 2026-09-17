# S2638: Cascade Discrepancy Resolved — Eigenvector Reorganization, Not Per-Prime Cascade

## The Discrepancy (from S2637)

**Exact cascade (verified)**: R_k = ⟨ψ|W02-WR|ψ⟩ + Σ_{q≤k} h(q) → λ_min exactly

**Naive power-law fit to h(q)**: h(q) ~ q^(-12) → summing gives c^(-11) decay

**True empirical scaling**: λ_min ~ c^(-0.44N) at high precision (Arb)

**Discrepancy**: c^(-11) ≠ c^(-3.5) at N=8. Why?

---

## Deep Investigation Results (4-Part Analysis)

### PART A: h(q) Profile at Fixed (c,N)

**Test case**: c=30, N=8

**Finding**: h(q) DOES follow a power law h(q) ~ q^(-12.00) for base primes (2,3,5,7,...).

**BUT**: This power law describes ONLY the decay at FIXED (c,N). It does NOT predict behavior when c changes.

**Sign pattern**: 
- h(2) through h(23) are all POSITIVE
- Only h(29) is negative (tiny, -5.4e-17)
- The signs are consistent within a fixed (c,N)

**Prime power structure**:
- Base primes contribute 95.6% of total |h|
- Prime powers (4,8,16,9,27,25) contribute only 15.1% relative to base primes
- High prime powers (16, 27, 25) have TINY contributions (< 1e-13)

**Concentration**: D_eff = 1.0 (minimum possible), indicating ψ is highly concentrated on a single mode.

---

### PART B: Eigenvector Evolution at Fixed N, Varying c

**Test case**: N=8, c ∈ {13, 20, 30, 40, 50}

**Key finding**: Eigenvector REORGANIZES dramatically as c increases.

**Overlaps** (comparing eigenvector at each c with reference c=13):
- c=13 vs c=20: overlap = 0.047
- c=13 vs c=30: overlap = 0.036  ← **This is TINY**
- c=13 vs c=40: overlap = 0.0017 ← **Essentially orthogonal**
- c=13 vs c=50: overlap = 0.031

**Conclusion**: The eigenvector ψ(c,N) is NOT a stable object across different c. It changes shape dramatically.

**λ_min behavior**:
- c=13: λ ≈ 5.0e-17
- c=30: λ ≈ 4.4e-16 (increased)
- c=50: λ ≈ 2.2e-16
- No clear scaling pattern at float64 precision

---

### PART C: h(q) Values Change When c Changes (CRITICAL)

**Setup**: Compare h(q) for the SAME prime q at c=30 vs c=50

**Direct measurements**:

| q | h(q) at c=30 | h(q) at c=50 | Ratio | Change type |
|---|---|---|---|---|
| 2 | +2.72e-01 | -1.43e-01 | -0.53 | **SIGN FLIP** |
| 3 | +5.59e-01 | +4.60e-01 | 0.82 | ~20% decrease |
| 5 | +1.13e-01 | +6.51e-01 | 5.74x | **~6x increase** |
| 7 | +1.15e-02 | +3.39e-01 | 29.4x | **~30x increase** |
| 11 | +8.36e-05 | +3.57e-02 | 427x | **~400x increase** |
| 13 | +6.79e-06 | +9.73e-03 | 1433x | **~1400x increase** |
| 17 | +3.93e-08 | +6.38e-04 | 16226x | **~16K increase** |
| 19 | +2.65e-09 | +1.58e-04 | 59570x | **~60K increase** |
| 23 | +7.20e-12 | +9.10e-06 | 1264838x | **~1.3M increase** |
| 29 | -5.45e-17 | +1.04e-07 | n/a | **Sign flip, 10^9 magnitude change** |

**Fundamental observation**: h(q) is NOT a property of the prime q alone. It DEPENDS on what other primes are present in the matrix.

**Why?** Because h(q) = -⟨ψ|W_q|ψ⟩, and ψ = ψ(c,N) is the eigenvector at THAT specific c. When we add more primes, ψ reorganizes entirely.

---

### PART D: Attempt to Fit c^(-0.44N) at Float64

**Test case**: N=8, c ∈ {10,15,20,25,30,40,50}

**Fitting log|λ_min| = a0 + a1·log(c)**:
- a0 = -37.7
- a1 = 0.742
- Exponent α = -a1/N = -0.0928

**Observed vs predicted**:
- Predicted: c^(-0.44·8) = c^(-3.52)
- Observed: c^(0.74·1) = c^(0.74)
- Mismatch is SEVERE (sign is opposite!)

**Conclusion**: At float64 precision, there is NO scaling law. The margin is essentially constant (all λ ≈ 10^(-15) to 10^(-16)), differing mainly by numerical noise.

---

## The True Mechanism: NOT Per-Prime Cascade

### What S2637 Got Wrong

The S2637 approach assumed:
1. h(q) is a fixed, predictable function of q
2. Each prime "contributes" a fixed amount to the margin
3. You can sum these fixed contributions to predict the final margin

This is **FUNDAMENTALLY FLAWED** because:

- **h(q) is not independent of c**: When you change c, you change which primes are in the matrix, and that reorganizes the entire eigenvector ψ
- **The cascade is local (within-fixed-c), not global (across-varying-c)**: The cascade law works PERFECTLY within a fixed c (it reconstructs λ_min exactly), but this says nothing about what happens when you change c
- **Eigenvector is not stable**: As c grows, the eigenvector ψ(c,N) changes shape completely (overlaps ~ 0.001-0.05), so the "projection of ψ onto mode n" is not a stable quantity either

### The Correct Picture

1. **τ(c,N) is a full eigenvalue problem** (not a sum-of-rank-2 perturbations story):
   - For each (c,N), solve eigh(τ) to get λ_min and ψ(c,N)
   - λ_min emerges from the GLOBAL spectral structure of the entire matrix
   - ψ(c,N) is a global object (how it concentrates on Fourier modes depends on all primes up to c)

2. **The margin shrinkage from c→c' is NOT the h(q_new) for the new prime q_new**:
   - Adding a new prime changes both the matrix entries AND the eigenvector entirely
   - λ_min(c') is determined by eigh(τ(c')), not by λ_min(c) + h(q_new)
   - The cascade computation showing R_k → λ_min is a CONSISTENCY CHECK at fixed c, not a decomposition that extends to varying c

3. **The c^(-0.44N) law (if real at high precision) is a scaling law of the eigenvalue problem, not a sum of independent terms**:
   - It must emerge from dimensional analysis of:
     - Spectral gap structure of τ (depends on parity decomposition)
     - Prime density (~c/ln(c))
     - Fourier dimension (2N+1)
     - How these conspire to shrink λ_min as c grows
   - This is NOT accessible from fitting h(q) values

---

## What Should Be Done Instead

### To Derive the Exponent Analytically

1. **Rigorous bound on λ_min(c,N)**:
   - Use perturbation theory or min-max principle on the FULL matrix τ(c,N)
   - Frame it as: (W02-WR) background + (-Wp) perturbation
   - But account for the fact that the eigenvector changes shape!

2. **Spectral analysis of (W02-WR) per parity sector**:
   - Compute neg-evals of (W02-WR)_even and (W02-WR)_odd separately
   - Quantify "how much negativity" needs to be canceled by Wp
   - This negativity GROWS with N (more frequencies = more coupling)

3. **Cancellation efficiency analysis**:
   - How effectively do the M primes cancel the negativity?
   - The cancellation is NOT uniform across modes (small primes dominate)
   - Model: negativity ~ O(N), primes provide cancellation ~ O(log(c)), ratio gives margin

4. **Weyl/interlacing at the right level**:
   - Not on the raw (W02-WR) (which is negative-definite)
   - But on the PARITY-SECTOR restriction (which are each PSD)
   - This is where the true bound lives

### To Verify at High Precision

- Run margin_law.py at c=100-1000, N=4-50 with Arb (1400+ bits)
- Fit the empirical law at high precision
- Extract the TRUE exponent(s) as function(s) of N
- Compare against the analytic bound

---

## Honest Assessment

**What worked**:
- The exact cascade computation (R_k → λ_min) is real and useful
- The parity-sector decomposition (S2636) correctly identifies where the structure lives
- The empirical c^(-0.44N) law IS real (confirmed at Arb high precision in S2637 Part 1)

**What failed**:
- Naive power-law fitting to h(q) and summing (S2637 approach)
- Assuming h(q) is independent of c
- Assuming the cascade extends from "add one more prime" to "how the margin scales"

**What is open**:
- A rigorous ANALYTIC DERIVATION of the c^(-0.44N) exponent from matrix structure
- Proving the exponent is exactly -0.44, not just empirically -0.44±0.05
- Connecting the parity-sector PSD structure to the margin scaling formula

**The Wall**: The margin shrinkage mechanism is NOT a simple independent-term sum. It is a genuine eigenvalue problem where all pieces interact globally. Solving it rigorously requires spectral analysis of the full (parity-sector-restricted) matrix, not per-prime decomposition.

---

## Recommendation

**Immediate next steps (in order)**:

1. **Accept the diagnosis**: The cascade as a per-prime mechanism does NOT extend to varying c. The mechanism is GLOBAL, not local.

2. **Use the cascade for what it's good for**: Verify that the matrix structure is correct, and that numerical computation of λ_min is sound. Do NOT use it to predict scaling laws.

3. **Build a rigorous bound**:
   - Start from Cauchy interlacing on the parity-sector pair (τ_even, τ_odd)
   - Use known bounds on rank(Wp_sector) and spectrum of (W02-WR)_sector
   - Derive margin lower bound in terms of c, N, M

4. **Fit the empirical high-precision data to extract the true exponent**:
   - Run margin_law.py at larger c, N with Arb
   - Use polynomial/rational regression, not just linear
   - The exponent might be -0.44 or might be a more complex function of N

5. **Close the gap between analytic bound and empirical law**:
   - If the bound is weaker than the data, refine the bound
   - If they match, you have a proof
