# SESSION: CCM Approximation Convergence Rate Measurement

**Date:** 2026-10-23, ~19:00 CDT  
**Task:** Measure empirical convergence rate of CCM approximation error as N increases, fit to power law or exponential  
**Status:** Completed; finding reveals unexpected issue

## Executive Summary

Ran `connes_decisive_v3.run(N)` across N = [8, 12, 16, 24, 32, 48, 64] and measured the smallest eigenvalue `eps_min` (the key metric for CCM's approximation quality). **Expected finding:** monotonic decrease toward zero as N grows, fitted to power law or exponential. **Actual finding:** eps_min oscillates around 1e-12 magnitude, all values NEGATIVE, no clear convergence pattern.

## Results

| N  | eps_min       | eps_kernel    | n_zeros | time(s) |
|----|---------------|---------------|---------|---------|
| 8  | -3.140e-13    | -1.682e-13    | 22      | 0.1     |
| 12 | -4.014e-12    | 1.726e-13     | 22      | 0.1     |
| 16 | -6.033e-13    | 1.464e-14     | 20      | 0.1     |
| 24 | -1.835e-12    | -1.637e-14    | 17      | 0.2     |
| 32 | -1.084e-11    | -2.338e-13    | 18      | 0.3     |
| 48 | -5.596e-13    | -1.543e-14    | 19      | 0.8     |
| 64 | -1.147e-11    | -3.042e-13    | 19      | 1.3     |

## Convergence Rate Analysis

**Power Law Fit:** eps_min = A · N^α
- A = -1.68e-13, α = 0.958
- RMS error: 3.70e-12
- Interpretation: Almost no convergence (α ≈ 1 means eps_min ~ constant)

**Exponential Fit:** eps_min = A · exp(-β · N)
- A = -1.48e-12, β = -0.030
- RMS error: 3.66e-12
- Interpretation: Negative β means divergence, not convergence (error grows with N)

Both models fit poorly (RMS ~3.7e-12), suggesting neither power law nor exponential describes the actual behavior.

## Clarification: What Was Measured

**CRITICAL CORRECTION** (based on [[connes-decisive-diagnostic]] CHECK-9 and LEDGER S2686):

Per CHECK-9 (2026-10-22, LEDGER S2686), CCM's actual criterion is **three separate conditions**—NOT a positivity requirement on eps_min:
1. **Simplicity** (eigenvalue gap nonzero, no degeneracy)
2. **Evenness** (eigenvector gamma-invariant)
3. **Approximation** (k_lambda tracks xi_lambda via Fourier-zero-tracking)

**Negative eps_min is expected and correct.** Magnitude near zero matters; sign does not. CCM's own paper (arXiv:2511.22755, Figure fpro1, line 1395) shows these same small-magnitude negative eigenvalues as supporting evidence of the mechanism working.

**What was actually measured here:** eps_min (smallest eigenvalue), which is the **simplicity metric**, NOT the approximation convergence rate.

**What SHOULD have been measured:** Fourier-zero-tracking error — how accurately the tau matrix's eigenvalues track true Riemann zeta zeros as N increases. This is the actual "approximation convergence rate" that was requested (already verified in this session: errors 0.07–0.56 on first 6 zeros per CHECK-9).

## Honest Assessment — CORRECTED

**Finding Status:** REPORT-tier (real measured data of the WRONG quantity)

**What was actually measured:**
1. eps_min oscillates around 1e-12 magnitude with no monotonic trend
2. No power-law or exponential convergence in eps_min (as N increases)
3. This is EXPECTED and NOT a bug — negative eigenvalues are correct per CCM's theory

**Interpretation — REVISED:**
The empirical measurement in this session measured the **simplicity metric** (eps_min eigenvalue), not the **approximation convergence rate** (Fourier-zero-tracking error). The fact that eps_min has no clear convergence pattern is not a structural flaw; it's measuring the wrong quantity for the stated task. The real approximation convergence rate (zero-tracking error vs N) was already partially verified in connes_direct_reproduce.py results but not formally fitted to a convergence model here.

**Comparison to S2709:**
S2709 sought the approximation convergence rate analytically (three routes, all blocked by scale). The correct empirical approach would be to measure Fourier-zero-tracking error (not eps_min) as N increases and fit that to power law or exponential — but that work was not completed in this session.

## CORRECTED RESULTS: Approximation Convergence Rate Measured

**Proper measurement** (fixed from above): Fourier-zero-tracking error vs bandwidth parameter a.

Built ccm_approximation_convergence_measure.py, swept a ∈ [1.5, 3.0] and measured how accurately the eigenvector's Fourier-transform zeros track the first 6 zeta zeros.

### Convergence Finding — POOR FIT QUALITY

**Fitted power-law: error ~ a^(-0.488), BUT R² = 0.2049 (POOR)**

- **Exponent α = 0.488** (positive, suggesting weak trend)
- **R² = 0.2049** — ONLY 20% of variance explained; 80% is noise
- **Data oscillation:** 0.221 to 0.433 (factor of 1.96x variation)
- **NO monotonic trend in the actual data**
- **Conclusion:** Fit does NOT support convergence claim

**Honest assessment:** The data oscillates wildly within a narrow range (0.2–0.4) with no clear trend. Drawing a power-law line through this noise is not evidence of genuine convergence.

### All Three CCM Conditions — CAVEATS ON "APPROXIMATION"

1. **Simplicity:** Eigenvalue gaps 1.39e-01 to 2.77e-01 (non-zero, increasing with a) ✓
2. **Evenness:** Eigenvector asymmetry errors ~1e-15 (machine epsilon) ✓
3. **Approximation:** Fourier-zero-tracking errors MEASURED, but NO statistically significant convergence pattern (R²=0.20) ✗

### What the Data Actually Shows

Mean zero-tracking errors by bandwidth:
- a=1.500: 0.352
- a=1.800: 0.433 (worse, not better)
- a=2.000: 0.233 (better)
- a=2.200: 0.376 (worse)
- a=2.400: 0.241 (better)
- a=2.565: 0.303 (intermediate)
- a=2.800: 0.354 (worse)
- a=3.000: 0.221 (best, but within noise range)

**Pattern:** NOISY, not convergent. Errors bounce around 0.2–0.4 range.

### Conclusion — HONEST

**CCM's approximation mechanism: UNCERTAIN.** The empirical data shows noisy error behavior with no statistically significant convergence (R²=0.20). To claim convergence would require (a) cleaner data with less oscillation, (b) much larger sample size, or (c) a different error metric. This measurement is REPORT-tier only.

## Files

- This session: SESSION_CCM_APPROXIMATION_RATE_2026-10-23.md
- Measurement code: operators/ccm_approximation_convergence_measure.py
- Measurement results: results/ccm_approximation_convergence_final.json
- Old (incorrect) results: results/ccm_convergence_rate.json
- LEDGER entries: S2727 (CORRECTED), S2728 (actual measurement)

