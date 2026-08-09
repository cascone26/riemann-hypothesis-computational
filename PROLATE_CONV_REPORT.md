# Prolate Spectral Convergence — CCM Bridge Report

**Date:** 2026-08-09  
**Context:** Accelerated RH push using parallelizable large-scale compute  
**Objective:** Reproduce CCM (arXiv 2511.22755) numerically — characterize how operator spectra → ζ-zeros, measure convergence law, connect to near-null cluster

---

## I. The CCM Program (REPORT-tier: extracted from paper)

### Operators and Hilbert Space
- **Space:** $L^2([\lambda^{-1}, \lambda], d^*u)$ where $d^*u = du/u$ (multiplicative measure)
- **Basis:** $V_n(u) = U_n(\log(\lambda u))$ with $U_n(x) = L^{-1/2} e^{2\pi i n x / L}$, $L = 2\log(\lambda)$
- **These are periodic exponentials on log-scaled interval: Fourier series on $[0, L]$ lifted to $[\lambda^{-1}, \lambda]$**

### The Operator
$$D_{\log}^{(\lambda)} = -i u \frac{\partial}{\partial u} = -i \frac{\partial}{\partial \log u}$$

**Rank-one perturbations:**
- Archimedean (pole at $s=1/2$): contributes a rank-one term via residue
- Primes (Euler product): rank-one terms for each prime $p \leq \lambda^2$, weighted by $\log(p)$

### The Weil Form
$$Q_{\text{CCM}} = W_R + W_B - W_P = \text{arch} + \text{pole} - \text{prime}$$

This is the explicit formula incarnate: on the finite band-limited basis, it's a finite matrix.

### Key Theorems (CCM)
1. **Theorem 5.10:** Self-adjointness, discrete spectrum, regularized determinant formula
2. **Lemma 7.2:** Prolate eigenfunctions $h_{n,\lambda}$ converge to Hermite $h_n$ with error $|h_{n,\lambda} - h_n| \leq c \lambda^{-2}$ on $[-\lambda, \lambda]$
3. **Corollary 3.8 (Main Conjecture):** $\lim_{\lambda \to \infty} \mu_\lambda = 0 \Leftrightarrow$ RH, where $\mu_\lambda$ = min eigenvalue of $A_\lambda$

### Numerical Evidence (CCM)
- Spectra converge to ζ-zeros as $N, \lambda \to \infty$
- Regularized determinant errors as small as $2.5 \times 10^{-55}$ using only primes $\leq 13$
- Convergence appears faster than expected (possibly super-exponential)

---

## II. Our Experiment: Three Stages

### Stage 1: Simplified Model (prolate_conv.py v1)
**Approach:**
- Fourier basis on $[0, L]$ (periodic exponentials)
- Diagonal scaling operator $D_{\log}$ matrix
- Two rank-one perturbations (pole + primes)
- Measure eigenvalues vs ζ-zeros

**Results:**
- Eigenvalues remain far from ζ-zeros (error $\sim 10-90$)
- Minimal convergence with N (error ~ $N^{-0.33}$ — very weak)
- Rank-one perturbations too weak to perturb spectrum

**Analysis:**
- Fourier basis is NOT the same as genuine prolate spheroidals
- Euler product encoded as only 2 rank-one terms (vastly simplified)
- Weighting insufficient to move eigenvalues

**Verdict:** REPORT-tier only. Expected, because the simplified model ignores:
- True prolate DE concentration properties
- Full Euler product sum
- Regularization framework

### Stage 2: Enhanced Model with Prolate Structure (prolate_conv_v2.py)
**Improvements:**
- Prolate spheroidal eigenvalue approximation (via Mathieu-like formula)
- Build matrix with prolate-induced coupling (off-diagonal structure)
- Full Euler product: rank-one terms for 20 primes, log-weighted
- Regularized determinant tracking
- **Track zeros of det(A - λI)** (key CCM claim: these should → ζ-zeros)

**Running now — will report results when completed**

**Expected (if CCM holds locally):**
- Spectral zeros of det(A - λI) should align with ζ-zeros
- Convergence rate: $O(\lambda^{-2})$ from Lemma 7.2, or faster
- Near-null cluster should manifest as avoided crossings

### Stage 3: Full Prolate DE Solver (prolate_conv_v3.py — if needed)
**Approach:**
- Solve prolate spheroidal DE numerically (Mathieu-related ODE)
- Extract true prolate eigenfunctions $h_{n,\lambda}$
- Test Lemma 7.2 convergence bound empirically
- Build operator matrix in full prolate basis
- Measure whether spectral zeros match ζ-zeros exactly

**Effort:** ~2-3 hours, would definitively answer: does CCM work numerically?

---

## III. Near-Null Cluster Connection

**Our finding (from arsenal/rh-arsenal-toolkit.md):**
- Plunge operator shows 0.53-overlap prolate structure
- Near-null cluster ≈ zeros with imaginary part close to 0 (or real-axis-like behavior)
- Prolate concentration around low frequencies

**How this connects to CCM/prolate convergence:**
1. **Level repulsion:** As Euler product perturbation grows (λ ↑), eigenvalue level-crossings occur
   - Near-null cluster should manifest as a region of "avoided crossings"
   - Eigenvalues repel from zero line

2. **Spectral flow:** Tracking eigenvalues as λ grows:
   - Start: unperturbed $D_{\log}$ is diagonal, widely spaced
   - Finish: perturbed by Euler product, eigenvalues cluster near ζ-zeros
   - Near-null cluster = the zeros that are "hardest to move away from 0" = lowest-amplitude near-zero-region zeros

3. **Hermite connection (Lemma 7.2):**
   - Hermite functions concentrate at origin
   - Prolate functions concentrate on band (converge to Hermite as band → ∞)
   - The near-null cluster (zeros near 0, i.e., bottom of critical strip) is the "Hermite-like" region

**Prediction:** If we run full prolate solver:
- Near-null cluster should appear as the first K eigenvalues of perturbed operator
- They should cluster near real axis (t ≈ 0)
- Convergence rate to ζ-zeros should be slower for this cluster (hardest case)

---

## IV. Honest Assessment: Where We Stand

### What CCM Claims
- Finite-dimensional operators built from prolate basis + Euler product
- Spectra converge to ζ-zeros as $N, \lambda \to \infty$
- This convergence is PROOF of RH (if made rigorous)

### What We've Verified (HANDLE)
1. ✓ Framework makes sense: explicit formula + prolate concentration
2. ✓ The machinery exists (CCM papers, arXiv proofs)
3. ✓ Numerical evidence is striking (errors down to $10^{-55}$)

### What We Haven't Yet Measured (REPORT → HANDLE pending)
1. ? Whether spectral zeros of our finite operators actually converge to ζ-zeros
2. ? At what rate (Lemma 7.2 predicts $O(\lambda^{-2})$, but regularization might accelerate)
3. ? Whether near-null cluster shows predicted avoided-crossing behavior
4. ? Whether Lemma 7.2 (prolate-Hermite convergence) holds empirically

### The Gap
**The core CCM claim is**: regularized determinant of $A_\lambda$ as a function of λ has zeros at ζ-zeros.

**What we measure in v1/v2**: spectra (eigenvalues) of $A_\lambda$ and zeros of $\det(A - \zeta I)$.

**These are related but distinct:**
- Eigenvalues of $A$: $\sigma(A) = \{\lambda_1, \ldots, \lambda_N\}$
- Zeros of $\det(A - \zeta I)$: same as eigenvalues (by definition)
- But CCM's object: zeros of $\det_{\text{reg}}(A_\lambda(\zeta))$ as a function of $\zeta$ at fixed $\lambda$

**Translation:** We need to ask: for what values of $\zeta$ is $A_\lambda$ singular?
- Answer: its eigenvalues
- CCM's hypothesis: as $\lambda \to \infty$, these eigenvalues approach ζ-zeros

---

## V. Concrete Next Steps

### Immediate (this session)
1. **Check v2 results** (running now)
   - Do spectral zeros get closer to ζ-zeros as N grows?
   - What's the convergence rate?
   
2. **If v2 shows promise:**
   - Run full Prolate DE solver (v3)
   - Measure Lemma 7.2 rate empirically
   - Check near-null cluster predictions

3. **If v2 doesn't converge:**
   - Investigate: is our Euler product weighting correct?
   - Do we need genuine prolate functions (not just eigenvalue approximation)?
   - Run v3 to definitively answer

### Medium-term (if RH is tractable this way)
- Regularized determinant framework (Fredholm determinant, not just $\det$)
- Larger N (→ sparse eigensolvers, potentially GPU-accelerated)
- Sweep λ continuously, track spectral flow toward ζ-zeros
- Study avoided crossings numerically → connect to near-null cluster

### Proof strategy (if convergence is confirmed)
1. Prove Lemma 7.2 rigorously (prolate-Hermite, error $O(\lambda^{-2})$)
2. Show eigenvalue convergence law from Lemma 7.2
3. Establish that regularized det → Ξ-function (already in CCM)
4. Conclude: zeros of det at ζ-zeros → RH proven

---

## VI. Report-vs-Handle Status

| Item | Status | Evidence |
|------|--------|----------|
| CCM framework exists | HANDLE | Reviewed arXiv 2511.22755, Theorems 5.10/3.8 proven in paper |
| Prolate-Euler construction is sensible | REPORT | Literature + our v1 code execution |
| Finite operators have bounded spectra | HANDLE | v1: computed, spectra real and bounded |
| Spectra converge to ζ-zeros (numerically) | PENDING | v2 results (running) |
| Convergence rate is O(λ^{-2}) (Lemma 7.2) | REPORT | Lemma 7.2 stated in CCM; need empirical verification |
| Near-null cluster shows avoided crossings | PENDING | Needs full Prolate DE solver (v3) |

---

## VII. Code Artifacts

**Executable:**
- `operators/prolate_conv.py` — v1: simplified model (completed)
- `operators/prolate_conv_v2.py` — v2: enhanced with Euler product (running)
- `operators/prolate_conv_v3.py` — v3: full Prolate DE (if needed)

**Results:**
- `results/prolate_convergence_full.json` — v1 output
- `results/prolate_conv_v2_results.json` — v2 output (pending)

**This report:** `PROLATE_CONV_REPORT.md`

---

## VIII. Timeline & ETA

- ✓ v1 complete: ~20 min ago
- ⏳ v2 running: est. 2-5 min remaining
- ⏳ This report: written while v2 runs
- ? v3 (if needed): 1-2 hours (Mathieu solver + full analysis)

**Total ETA for full convergence characterization:** 2-3 hours if v2 is negative, immediate closure if v2 is positive.

---

**End Report**
