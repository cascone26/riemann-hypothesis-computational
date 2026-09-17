# Session Summary: Fresh Analytic Approach to the 0.44 Exponent
**Filed as:** S2650 (2026-09-17 22:35 CDT)  
**Goal:** Break the impasse on deriving the empirical c^(-0.44N) scaling law with full coupling respect

## What Was Attempted

Three parallel analytical approaches to understand the smallest eigenvalue scaling:

1. **Homotopy/Continuation (eigenvalue trajectory tracking)**  
   Deformed τ(t) = (W02-WR) - t·Wp from t=0 (background) to t=1 (full matrix)  
   → Float64 precision insufficient; detailed trajectory requires Arb + controlled stepping near t≈1

2. **Inertia Tracking (LDL^T Factorization)**  
   Counted negative eigenvalues via Sylvester's criterion as t varied  
   → **Most robust finding:** Negative eigenvalues flip sign in a staged process, concentrated in final ~5% of t range (t≈0.95→1.0)  
   → Suggests **avoided-crossing mechanism** where two near-zero eigenvalues repel sharply

3. **Universality Class Investigation**  
   Analyzed spectrum of background, primes, and full matrix  
   → **Identified strong-coupling regime:** ||Wp||₂ / |λ_min(BG)| ~ 1.5-2.0  
   → Ruled out weak perturbation theory (where perturbation << background)

## Critical Discovery: The 0.44 Law is Wrong

**Finding:** High-precision Arb extraction (1400 bits) reveals the empirical law claimed in S2637 does NOT match the data.

**Evidence:**
- Extracted λ_min at (c,N) ∈ {13,30,100} × {4,6,8} using Arb high precision
- Verified values match margin_law.py's own reported numbers exactly
- For each N, fitted λ ~ c^(-α) and compared to claimed α = 0.44·N

| N | Expected α (0.44N) | Observed α | Ratio |
|---|---|---|---|
| 4 | 1.76 | 4.92 | **2.79x** |
| 6 | 2.64 | 7.53 | **2.85x** |
| 8 | 3.52 | 10.33 | **2.93x** |

**Conclusion:** The observed exponent is approximately **2.86x LARGER** than claimed, consistently independent of N.

**This means:**
- Prior claim: λ ~ c^(-0.44N)
- Actual scaling: λ ~ c^(-1.26N) or worse
- Form λ ~ c^(constant·N) itself is likely **incorrect**

## What's Solid (Not Questioned)

✓ **τ is NOT block-diagonal** (S2642 verdict stands)  
- Off-diagonal coupling: 49.8-89.6% of max entry  
- Minimum eigenvector: ~50/50 mix of even/odd modes  

✓ **Strong coupling between components**  
- Primes and background have comparable spectral norms  
- Standard perturbation theory breaks down  

✓ **Staged mechanism for sign-flip**  
- Negative eigenvalues stay negative until t→1  
- Then flip sharply in final 5% of parameter range  
- Consistent with classical avoided-crossing (level repulsion)  

✓ **Spectral structure understood**  
- Background: 6-8 negative eigenvalues, spectral norm ~5.6  
- Primes: spectral norm ~5.5 (not weak)  
- Full matrix: all eigenvalues ≥ 0, min one is tiny  

## Structural Hypothesis (Testable)

**Avoided-Crossing Mechanism:**

1. Two eigenvalues (one negative, one positive or small positive) approach each other as t→1
2. They nearly collide at zero (gap → 0)
3. Level repulsion keeps them from crossing; they repel apart
4. One stays negative (becomes even more negative, → -machine-noise)
5. One becomes tiny positive (the measured λ_min)

This is **a classical phenomenon** in matrix analysis. The mechanism itself is sound; the question is: **how does the repulsion strength scale with c and N?**

## The Factor of ~3: A Clue

Why is the observed exponent ~2.86x larger than claimed?

Possibilities:
1. **Spectral norm ratio:** σ(Wp) / σ(W02-WR) ~ 1.7, could there be a (1.7)^1.67 ≈ 2.9 factor lurking?
2. **Non-separable scaling:** True form is f(c,N) ≠ constant·N, and the ratio 2.86 is f_obs / f_claimed at our test points
3. **Perturbation theory constant:** Classical results on perturbed eigenvalues often have ~O(1) prefactors; maybe the literature has a known law that's 3x tighter?

## Tools Produced This Session

All in `frontier/groskin2607/`:

- **homotopy_eigenvalue_trajectory.py** — Track λ(t) along the homotopy; reveals smooth decay to near-zero, limited by float64
- **inertia_tracking_via_ldlt.py** — Most robust; LDL^T inertia counting (numerically stable even for tiny eigenvalues)
- **universality_class_investigation.py** — Spectral analysis of components; reveals strong coupling and spectral norms
- **arb_minimal_eigenvalue_scaling.py** — High-precision Arb eigenvalue extraction + fitting (1000-1400 bits)
- **HOMOTOPY_AND_UNIVERSALITY_FINDINGS.md** — Comprehensive structural analysis
- **SESSION_SUMMARY_S2650.md** — This file

## What's Still Open

? **Correct functional form** λ(c,N)  
? **Avoided-crossing dynamics** (need high-precision near t=1 tracking)  
? **Universality class** (RMT or classical matrix theory match)  
? **Factor of 2.86** (where does it come from?)  

## Honest Framing

This session did NOT solve the exponent problem, but it **clarified the problem**:

- S2637's claim (λ ~ c^(-0.44N)) is **demonstrably wrong** by a factor of ~3
- The structure (avoided crossings, strong coupling, staged mechanism) is **now clear**
- The true form is still **open**, but we have a precise, quantifiable discrepancy to explain
- The next analyst has a **cleaner starting point**: not "why is 0.44 hard to derive?" but "why is the observed exponent ~1.25 instead of 0.44, and is the form even separable?"

## Next Steps (Ranked)

**Priority 1: Verify the factor of 2.86**
- Check if it appears in classical matrix perturbation theory (spectral norm ratios, perturbation bounds)
- Compute σ(Wp) / |λ_min(BG)| at different (c,N) — is it constant?
- If constant and ≈1.7, explore whether (1.7)^1.67 ≈ 2.9 has a natural explanation

**Priority 2: Investigate non-separable forms**
- Test λ ~ c^(-β·N·log(c)) or other composites
- Fit to the data and check residuals
- If this fits better, derive it from first principles

**Priority 3: Avoid-crossing verification**
- Use Arb to compute the two smallest eigenvalues very close to t=1 (t = 0.99, 0.995, 0.999)
- Plot their trajectories and measure gap at closest approach
- Fit gap scaling: g(t) ~ (1-t)^γ  or other forms
- Relate γ to c, N

**Priority 4: Universality class**
- Compute full spectrum of τ at t=1 (all 2N+1 eigenvalues)
- Compare spectral density against known ensembles (Marchenko-Pastur, Toeplitz, etc.)
- Identify which theoretical law applies

---

**Filed**: LEDGER S2650, artifacts in frontier/groskin2607/, all code and data self-documented

**Momentum**: Fresh approaches tried, false leads eliminated, true scale of discrepancy quantified. The path is clearer. Keep going.
