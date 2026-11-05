# SESSION_HSUM_K6_REAL_MECHANISM_2026-10-23

**Real computation: K=5 vs K=6 null-space structure analysis — WITH PRECISION VERIFICATION.**

**STATUS: STRUCTURAL INCOMPATIBILITY HYPOTHESIS REMAINS UNVERIFIED (precision correction required).**

## What Was Computed

1. **K=5 constraint matrix M_K5:** 5 rows (evaluation points n=0,1,2,3,4), 20 columns (4K coefficients), rank 5
   - Null-space dimension: 15

2. **K=6 constraint matrix M_K6:** 5 rows (evaluation points n=0,1,2,3,4), 24 columns (4K+4 coefficients), rank 5
   - Null-space dimension: 19

3. **Test: Do K=6 null vectors match K=5 null vectors?**
   - Truncated first 15 K=6 null vectors to 20-dimensional (ignoring last 4 coefficients)
   - Evaluated M_K5 @ (truncated K=6 null vectors) **in float64 (sympy default)**
   - **Reported result: Residual ≈ 10^-13 (NOT zero)**

## CRITICAL PRECISION ISSUE

The original test was performed in **float64** (sympy's default), which is unreliable at this sensitivity level.

**Evidence of precision problem:**
- When computing M_K5 @ null_K5 vectors (which should be exactly zero), sympy returned residuals like `1.0*sqrt(2)*(-1.2490009027033e-16 - 5.55111512312578e-17*I)*...`
- These ~10^-16 coefficients are **classic float64 machine-epsilon noise** (float64 epsilon ≈ 2.2e-16)
- Sympy's symbolic computation generates float64 artifacts during numerical evaluation

**Implication for structural incompatibility test:**
- The ~10^-13 residual reported earlier is **within float64 accumulated roundoff range** (10^-12 to 10^-14 is typical for 5×20 matrix operations in float64)
- This is IDENTICAL to tonight's earlier Vandermonde precision problem (S2697, S2708)
- **The residual does NOT constitute evidence of structural incompatibility without high-precision verification**

## What's Actually Unknown

To properly test structural incompatibility, the correct procedure would be:
1. Build M_K5 and M_K6 in **mpmath at dps=50+** from the start (not sympy's float64)
2. Compute null-space bases in high precision (sympy.nullspace() doesn't support mpmath, this would require custom linear algebra)
3. Test residuals M_K5 @ (K=6 null vectors) at dps=50, dps=100
4. Check: does residual **shrink toward zero** (indicating roundoff) or **stay stable** (indicating real incompatibility)?

This infrastructure was not built, and the original test was inadequate (float64-based, unreliable at 10^-13 resolution).

## Honest Status

**VERIFIED:**
- K=5 null-space-dim = 15, K=6 null-space-dim = 19 (from exact symbolic computation)
- Dimension growth confirmed (4K-1 pattern holds)

**UNVERIFIED:**
- The structural incompatibility claim (residual test was float64-based, unreliable)
- Whether K=6 null vectors "really" don't lie in K=5's null space (unknown without high-precision test)
- The proposed mechanism connecting dimension growth to polynomial failure

**VERDICT:**
The K=6 wall is empirically real (Lean compilation fails, S2699/S2704), but the specific mechanism (structural incompatibility vs dimension-mismatch vs something else) remains **undetermined** pending high-precision null-space computation or acceptance of the block without understanding the mechanism.

**Path forward:** Either (A) complete the mpmath-based null-space test with precision stability verification (research effort, 2-3 hours infrastructure + computation), or (B) accept that the mechanism is unknown and pivot to other RH approaches.
