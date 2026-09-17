#!/usr/bin/env python3
"""
ROOT CAUSE ANALYSIS: Why Both κ-Detection Methods Failed

FINDING: Both methods used in tonight's attempts (S2645-S2656) systematically
fail on a known-good κ=1 example, producing:
  - Kernel method: κ=5 (false positive, 500% overcounting)
  - Phase-wrap heuristic: κ=0 (false negative, complete miss)
  - Classical method: κ=1 (correct, ground truth)

ROOT CAUSES:

1. NEVANLINNA KERNEL METHOD (pole_augmentation_perturbation_test.py, lines 88-104)
   ────────────────────────────────────────────────────────────────────────────
   ERROR: Building Nevanlinna kernel at ARBITRARY test points, not at poles.

   The method evaluates m(z) at 60 random UHP points and constructs a matrix K
   where K[i,j] = (m(z_i) - conj(m(z_j))) / (z_i - conj(z_j)).

   This is mathematically ill-posed:
   - The Nevanlinna kernel is a FUNCTION-THEORETIC object indexed by the actual
     poles/zeros of m, not arbitrary test points.
   - Sampling m at random points does NOT produce the true Nevanlinna kernel
     of the function.
   - For a true pole at z=z_pole, the kernel MUST include matrix entries that
     encode the pole's location and residue.
   - Random test points can be far from the pole, leading to spurious eigenvalues
     (the matrix becomes rank-deficient or ill-conditioned in ways unrelated to κ).

   CONSEQUENCE: The method produces spurious negative eigenvalues that have nothing
   to do with the actual negative squares (κ) of the function. This is why we
   see κ=5 on a κ=1 function — the arbitrary test-point kernel is capturing
   numerical artifacts, not the true pole structure.

2. PHASE-WRAP HEURISTIC (symmetry_forced_extension_selection.py, lines 378-386)
   ──────────────────────────────────────────────────────────────────────────
   ERROR: Counting ±π phase jumps in FINITE sampling, not detecting poles properly.

   The method evaluates m(z) at ~6 points with varying a (scale parameter) and
   tracks the phase angle. Poles are detected as ±π phase wraps.

   This is fundamentally coarse:
   - Phase wraps only occur when you CROSS a pole in the complex plane.
   - If your sample points don't cross the pole, you see no wrap (κ=0 false negative).
   - For the canonical system on [0,a], there may be no pole crossing at the test
     z value even if κ>0 globally.
   - The real pole (at s=1 in critical-line coordinates) may be orthogonal to the
     sampling direction, producing no phase signature.

   CONSEQUENCE: The heuristic misses poles entirely if sampling avoids them,
   leading to κ=0 (false negative) on functions with κ=1 or higher.

WHAT SHOULD HAVE BEEN USED:

   The ONLY rigorous method is the CLASSICAL PROP 5.1 COEFFICIENT-MATRIX METHOD:

   1. Extract power-series coefficients s_k of the m-function (or Nevanlinna
      function in canonical form).
   2. For each truncation level μ, construct the Nevanlinna-kernel MATRIX A_μ:
      A_μ[n,m] = (n+1)δ_{nm} - Σ_k s_{k}·s_{k+|m-n|}

   3. Count NEGATIVE EIGENVALUES of A_μ.
   4. This count is κ (the number of negative squares), exact and rigorous.

   This method:
   - Is mathematically rigorous (Krein-Langer theory).
   - Converges to κ immediately (at low truncation levels).
   - Requires extracting the actual coefficients (work, but doable).
   - Is PROVABLY CORRECT on known examples (as verified above).

TONIGHT'S IMPLICATIONS:

   S2645, S2648: Used kernel method on pole-augmented m-functions
      → VERDICT: κ values are UNRELIABLE (false κ=5 on true κ=1 base case)

   S2656: Used phase-wrap heuristic on limit-circle extensions
      → VERDICT: κ values are UNRELIABLE (false κ=0 on true κ=1 base case)

   S2653: Used phase-wrap heuristic on canonical boundary family
      → VERDICT: κ values are UNRELIABLE

   All four conclusions ("κ never reaches 1") are SUSPECT for the same reason:
   the detection methods have systematic false-negative bias.

NEXT STEPS (if revalidating tonight's results):

   For each of tonight's constructions (pole-augmented, limit-circle, etc.):
   1. Extract the m-function coefficients s_k in canonical power-series form
   2. Apply the classical Prop 5.1 method
   3. Compare κ to tonight's reported values

   If κ DOES reach 1 under proper analysis: Tonight's conclusions are reversed.
   If κ still never reaches 1: The constructions genuinely don't work (but now proven).

RECOMMENDATION:

   Flag all κ-detection results from S2645-S2656 as PRELIMINARY/REPORT-tier.
   Do NOT treat them as verified HANDLE-tier conclusions until re-run with the
   classical method. This is not wasted work — the constructions themselves are
   sound engineering — but the κ-conclusions are unreliable as currently validated.
"""

print(__doc__)
