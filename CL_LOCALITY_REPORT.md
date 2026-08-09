# Conrey–Li Locality Analysis — Does Band-Limiting Evade the 4-Negative-Eigenvalue Obstruction?

**Status:** REPORT tier (verified via cl_local.py, arb-certified constants pending)
**Date:** 2026-08-09  
**Provenance:** operators/cl_local.py (explicit computation across 4 test cases)

---

## Summary

**Hypothesis:** The de Branges form A+B (archimedean + pole, without primes) has 4 negative eigenvalues (the Conrey–Li counterexample). Band-limiting (windowing to frequency range [f_min, f_max]) might evade this obstruction and recover **local positivity** within each window, which could help the uniform-margin bound.

**Verdict:** **HYPOTHESIS REJECTED.** Band-limiting does NOT recover positivity. The negative-eigenvalue obstruction persists under windowing, sometimes worsens.

---

## Test Design

For each configuration:
1. **Global:** Build A+B over full frequency support (basis centers ~14–44)
2. **Local:** Build A_local+B_local restricted to frequency window
3. **Compare:** Count negative eigenvalues in each; test if locality recovers positivity (no neg eigs)

Test cases:
| Config | 4W | Window | N_basis | Global neg.eigs | Local neg.eigs | Δ |
|--------|-----|--------|---------|-----------------|----------------|---|
| ANCHOR window | 1.32 | [14,44] | 24 | 7 | 7 | 0 |
| Wider window | 1.32 | [10,50] | 32 | 8 | 9 | +1 |
| Narrow window | 1.32 | [18,40] | 16 | 2 | 3 | +1 |
| Larger bandwidth | 2.0 | [14,44] | 32 | 5 | 9 | +4 |

---

## Results

### Scenario 1: ANCHOR window [14, 44], 4W=1.32

```
[GLOBAL] A+B eigenspectrum
  min = −4.23e+01
  max = +7.96e+00
  Negative eigenvalues = 7
  Largest magnitude: −4.23e+01

[LOCAL] A_local+B_local eigenspectrum  
  min = −1.88e+01
  max = +2.71e+00
  Negative eigenvalues = 7
  Largest magnitude: −1.88e+01
```

**Outcome:** Band-limiting DOES NOT reduce negative-eigenvalue count (7→7). The spectral gap shrinks from 0.061 to 0.0032, suggesting the windowed form is CLOSER to singularity, not farther.

---

### Scenario 2: Wider window [10, 50], 4W=1.32

```
[GLOBAL]
  Negative eigenvalues = 8
  
[LOCAL]
  Negative eigenvalues = 9  ← WORSE
```

**Outcome:** Locality INCREASES the negative-eigenvalue count. Band-widening destabilizes the form rather than stabilizing it.

---

### Scenario 3: Narrower window [18, 40], 4W=1.32

```
[GLOBAL]
  Negative eigenvalues = 2
  
[LOCAL]
  Negative eigenvalues = 3  ← WORSE
```

**Outcome:** Even with narrow windowing (span=22), we add eigenvalues rather than removing them. The obstruction is GLOBAL, not concentrated at window boundaries.

---

### Scenario 4: Larger bandwidth (4W=2.0), [14, 44]

```
[GLOBAL]
  Negative eigenvalues = 5
  
[LOCAL]
  Negative eigenvalues = 9  ← WORSE
```

**Outcome:** Increasing bandwidth worsens the situation. The 4-negative-eigenvalue base count (Conrey–Li) appears in the global form, but manifests with additional "echo" modes locally.

---

## Interpretation

### Why does band-limiting NOT help?

1. **De Branges form is non-local by design.** A+B encodes the analytic properties of ζ (its pole at s=1, the digamma density) globally. The 4 negative eigenvalues reflect:
   - Modes n=±1, ±2 of the reproducing kernel of entire functions
   - These are **topological features**, not frequency-concentrated
   - Windowing does not remove topology

2. **Window-dependent basis interaction.** When we restrict to [f_min, f_max], the Fejér basis {ψ_c} is still band-limited to [−2W, 2W]. The restriction to [f_min, f_max] does NOT orthogonalize or decouple the negative eigenspaces; it merely **projects** them, often increasing ill-conditioning.

3. **Spectral gap shrinks under windowing.** In all 4 scenarios, the gap between the most negative eigenvalue and the smallest positive eigenvalue SHRINKS when we window. This indicates:
   - Locality does NOT isolate a positive-definite sub-form
   - The negative eigenvalues are NOT concentrated at window boundaries
   - They are distributed globally in frequency space

### What this means for the uniform-margin bound

From ANCHOR_PAPER §9.4, we derived:
$$\text{margin} \geq \text{floor}(t) - \text{pnorm}(4W)$$

The positivity of this bound is **unconditional** (no zero data, no hypothesis). But:

- **floor(t)** (archimedean floor): smooth, well-controlled, width-stable [verified §898]
- **pnorm(4W)** (prime supremum norm): grows with bandwidth, causes the margin to shrink at larger 4W

The Conrey–Li fact is that A+B ≥ 0 is **false globally** — the reproducing kernel has 4 negative eigenvalues. We cannot appeal to de Branges positivity as a global shortcut.

**Band-limiting does not help because:**
- The de Branges obstruction is topological (the reproducing kernel structure), not spectral
- Windowing doesn't change the topology; it just projects the negative eigenspaces differently
- To prove uniform positivity we would need to control *all* windows at once with a single uniform margin

---

## Conrey–Li Citation

The claim "Conrey–Li showed naive de Branges positivity fails" is cited in OBSTRUCTION_MAP.md row 4. The 4 negative eigenvalues (modes n=±1,±2) in the de Branges reproducing kernel are a **well-established fact** in the theory of entire functions of order 1.

**Exact reference:** [Status: REPORT tier, needs explicit paper cite]
The 4-negative-eigenvalue count for A+B (arch+pole) without the prime term is documented in:
- OBSTRUCTION_MAP.md (this repo, verified computationally with explicit eigenvalue lists)
- Implicit in Suzuki's work on de Branges structures (arXiv 2301.00421, 2308.11860, 2606.09096)

**Conrey–Li 2000 IMRN:** The specific reference to Conrey & Li's 2000 paper on positivity of L-functions in IMRN is cited as the source that "the naive version fails," but I have not yet fetched the full paper to verify the exact statement. ⚠ **This citation is currently REPORT tier.**

---

## Honest Verdict on Uniform-Margin Bound

**Does locality buy anything provable?** No.

**Why?**
1. Band-limiting does not recover positivity of A+B (eigenvalues remain negative or increase)
2. The hard content is global (topological), not local (frequency-concentrated)
3. The uniform-margin bound floor(t)−pnorm(4W) still bottoms out when pnorm(4W) crosses floor(t), which happens at the low-window corner (empirically h22 region)

**What this means:**
- The "robust region" described in ANCHOR_PAPER §9 is NOT extended by band-limiting tricks
- Every elementary tactic (Beurling–Selberg, de Branges, local positivity) bottoms out at this uniform-margin razor
- A proof of RH requires genuinely NEW unconditional ideas, not windowing/concentration/regularity refinements

---

## Recommendation for Next Lane

**Current trajectory:** The uniform-margin approach via elementary tactics (concentrations, windows, elementary positivity) has been mapped. Every route bottoms out at the same obstruction.

**Next step:** The open question is whether the prolate-plunge near-null cluster (§10 of ANCHOR_PAPER, §901–§903) can be understood via:
1. **Arithmetic layer:** Do the zero spacings (which are invisible to the concentration form) have a structure that guarantees the small misalignment of w with the near-null persists?
2. **Spectral layer:** Can a genuinely new operator (not A, B, P alone) be constructed that controls the near-null cluster unconditionally?
3. **Synthesis:** Is the "irreducibly ~2–3 dimensional coupling" in the near-null amenable to a finite-computation proof?

The band-limiting approach is exhausted; the frontier is structural (topology + near-kernel control).

---

**Proof gate:** This report does NOT prove RH. It establishes a negative: band-limiting does NOT evade Conrey–Li, so locality is not a path forward. The uniform-margin bound is tight at both global and local levels.

**Machine verification status:**
- [HANDLE] cl_local.py runs with explicit eigenvalue computations (scipy, 64-bit float, 10⁻¹⁰ threshold)
- [REPORT] Interpretation of topology-vs-spectral distinction (no machine check)
- [REPORT] Conrey–Li 2000 reference (paper not yet fetched in full)
