# Homotopy & Universality Class Investigation — Findings & Direction

## Status
**Date:** 2026-09-16 (session start)  
**Goal:** Find a fresh analytic approach to the c^(-0.44N) exponent that respects full coupling  
**Prior work:** S2642 definitively proved τ is NOT block-diagonal; 49.8-89.6% off-diagonal coupling

## Findings from Eigenvalue Trajectory Analysis

### 1. Homotopy (Eigenvalue Trajectory) Analysis
Tracked λ(t) for τ(t) = (W02-WR) - t·Wp, t ∈ [0,1].

**Observation at float64:** eigenvalues hit machine precision noise near t=1. The bottom of the spectrum clusters at ~10^(-15) to 10^(-18), making detailed analysis impossible at float64.

**Key result:** The trajectories show smooth decay from negative values at t=0 to near-zero at t=1, but float64 cannot resolve the actual minimum margin.

**Status:** Homotopy approach is conceptually sound but requires very high precision tracking of small eigenvalues, which is expensive. Not pursued further at this stage.

### 2. Inertia Tracking (LDL^T Factorization)
Used Sylvester's criterion to count negative/positive eigenvalues as t varies. This is numerically stable even for small eigenvalues.

**Key observations:**
- **c=30, N=4:** 6 negative eigs at t=0 → 5 negative at t≈0.67 → 0 at t=1
- **c=30, N=8:** 6 negative eigs at t=0 → stays 6 until t≈0.95 → 0 at t=1
- **c=100, N=8:** 8 negative eigs at t=0 → stays 8 until t≈0.95 → 0 at t=1

**Interpretation:** The primes don't flip all negative eigenvalues at once. Instead, the cancellation mechanism shows a **staged process** where:
1. For t ∈ [0, t_cross), most negative eigenvalues remain stable
2. At t ≈ t_cross (very close to 1), all remaining negative eigenvalues flip to positive nearly simultaneously
3. The minimum eigenvalue at t=1 is the tiniest positive value resulting from this final flip

**Mechanism insight:** This staged flipping suggests that the smallest eigenvalue arises from an **avoided crossing** between a pair of nearly-degenerate eigenvalues that approach zero and repel each other just at t=1. The repulsion pushes one down nearly to zero.

### 3. Spectral Structure Analysis
Computed eigenvalue spectra of components:

**Background (W02-WR) at c=30, N=8:**
- λ_min = -3.27
- λ_max = +5.62 (spectral norm)
- n_neg = 6
- condition number κ ≈ 1.72

**Primes (Wp):**
- Spectral norm ≈ 5.54 (comparable to background!)
- This is NOT a small perturbation regime

**Ratio:** ||Wp||_2 / |λ_min(BG)| ≈ 1.7 (order-1 perturbation)

**Key insight:** Standard perturbation theory (which assumes ΔA << A) breaks down. The system is in a **strong-coupling regime** where the perturbation is the same order as the background.

## Structural Observations

### Not Block-Diagonal (S2642 Confirmed)
- Off-diagonal coupling between even/odd modes: 49.8-89.6% of max entry
- Minimum eigenvector is 50/50 mix of both parities
- This rules out any analysis based on sector decomposition

### Strong Coupling, Not Weak Perturbation
- Perturbation strength (Wp spectral norm) ~ 1.5-2x the background negativity
- Standard first-order perturbation theory doesn't apply
- Eigenvectors reorganize globally as t changes (overlaps drop to 0.001-0.05)

### Transition is Sharp Near t=1
- Inertia change (6→0 negative eigs) happens in final ~5% of the t range
- This suggests a **near-degeneracy** mechanism where two eigenvalues approach zero and then repel

## Hypothesis: Avoided Crossing at Level Crossing

The structure suggests the following mechanism:

1. **At t=0:** Background has 6 negative eigenvalues in a tight range (e.g., -3 to -0.1)
2. **For t ∈ [0,0.95]:** These eigenvalues move toward zero slowly, following smooth perturbative trajectories
3. **Near t≈0.99-1.00:** Two eigenvalues (one negative, one positive?) nearly collide at zero
4. **At t=1:** They repel due to level repulsion: one pushes down to -4e-16 (numerical noise), the other pushes up to +λ_min (the measured tiny positive value)

This is **a classical avoided-crossing scenario** from quantum/structural mechanics. When two energy levels approach each other, they typically repel rather than cross (in systems with Hermitian matrices).

If this is correct, the 0.44 exponent arises from understanding:
- **How rapidly** the negative eigenvalues approach zero as t→1
- **When** and **how strongly** the avoided crossing occurs
- **The repulsion law:** λ_⊥ ~ (gap_min)^β for some exponent β

## Direction: Toward Universality Class Theory

Rather than detailed perturbative tracking, the key is to identify which **universality class** this matrix belongs to.

### Candidates:
1. **Toeplitz-plus-low-rank:** τ has Toeplitz-like structure (W02, WR from smooth kernels) plus M ~ c/ln(c) rank corrections (primes)
2. **Near-singular Hermitian perturbation:** The background is nearly singular (has negative eigenvalues close to zero), and the perturbation is large enough to perturb the null/near-null space
3. **Level-crossing avoiding in random matrix ensembles:** If τ belongs to a known RMT class, the smallest-eigenvalue statistics might be known

### Classical Results to Consult:
- **Bai-Silverstein:** Extreme eigenvalues of Toeplitz-plus-perturbation
- **Edelman-Guionnet-Péché:** Spike models (low-rank perturbation of random matrices)
- **Horn-Johnson:** Perturbation bounds when perturbation and original matrix are comparable
- **Trefethen-Embree:** Pseudospectra and stability of near-singular matrices

### Concrete Next Steps:

1. **Verify avoided-crossing hypothesis:** Use Arb to compute the two eigenvalues nearest to zero at several t values close to t=1 (e.g., t = 0.98, 0.99, 0.995, 0.999, 1.0). Plot their trajectories and measure the gap at closest approach.

2. **Identify universality class:** Check if the spectrum of τ (at t=1) matches any known ensemble:
   - Compute spectral density (histogram of all eigenvalues)
   - Check for Marchenko-Pastur law (if related to random-matrix perturbations)
   - Compare with Toeplitz matrix spectral distributions

3. **Scaling of repulsion:** If two eigenvalues λ₁(t) < 0 and λ₂(t) > 0 approach each other as t→1, measure:
   - The gap g(t) = λ₂(t) - λ₁(t) as a function of t
   - Fit g(t) ~ (1-t)^γ or exp(-a/(1-t)) to extract the repulsion exponent
   - Relate this to c and N

4. **Connection to 0.44 exponent:** If avoided crossing explains the mechanism, derive how the repulsion strength and gap at t=1 depend on c and N:
   - Does gap ~ c^(-α*N)? Extract α.
   - Or does gap ~ f(c,N) for some other function?

## Current Code Artifacts
- `homotopy_eigenvalue_trajectory.py` — tracks λ(t), limited by float64 precision but conceptually sound
- `inertia_tracking_via_ldlt.py` — MOST ROBUST; reveals staged sign-flip via LDL^T inertia counting
- `universality_class_investigation.py` — spectral properties of components; reveals strong coupling
- `arb_minimal_eigenvalue_scaling.py` — (running in background) extracts λ_min at high Arb precision for multiple (c,N)

## Assessment

**What works:**
- ✓ Confirmed τ is NOT block-diagonal (S2642)
- ✓ Identified strong-coupling regime (perturbation ~ background in magnitude)
- ✓ Inertia tracking shows staged flip mechanism, suggesting avoided crossing
- ✓ Structural understanding clearer than before

**What's open:**
- ? Exact avoided-crossing dynamics (need high-precision tracking near t=1)
- ? Which classical universality class applies
- ? Exact derivation of 0.44 exponent from crossing dynamics + c,N scaling

**Momentum:**
The avoided-crossing hypothesis is concrete and testable. Next pass: (1) verify it numerically with high precision near t=1, (2) identify the universality class, (3) derive the exponent from crossing + perturbation scaling laws.

Not stuck, not dead — just shifted from perturbation-decomposition approach to avoided-crossing + universality-class approach. Both honor the full coupling; the second is more aligned with classical matrix theory.
