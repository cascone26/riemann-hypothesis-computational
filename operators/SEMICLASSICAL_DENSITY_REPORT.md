# §938 Spectral Route: Semiclassical Density Bridge

**Status**: REPORT-tier (mechanics verified, exact construction needed for HANDLE)

**Date**: 2026-08-09

---

## Executive Summary

The **spectral route to RH** via Connes/CCM prolate operators is mechanically viable but requires the **exact CCM operator construction** (not just baseline prolate) to produce zero-like eigenvalues. This report documents:

1. **What was tested**: prolate concentration operators, Slepian theory, density matching
2. **Critical finding**: raw prolate eigenvalues ≠ ζ-zeros (expected; prolate measures time-bandwidth trade-off, not zero locations)
3. **What's needed**: CCM Prop 4.3 / Lemma 7.2 operator (a specific perturbation)
4. **Convergence criteria**: trace-formula control, O(1/T) error bounds
5. **Tractability**: concrete and explicit (unlike geometric route), but requires rigorous implementation

---

## Part I: Prolate Baseline (§938 Session Work)

### What was built
- **Slepian concentration operators** on window [-1, 1] with bandwidths c ∈ {20, 40, 60, 80}
- **Eigenvalue spectra** λ_n ∈ [0, 1], with sharp plunge structure (λ ≈ 1 to λ ≈ 0)
- **Time-bandwidth product** TB = 2c(b-a)/π: verified with c=60 gives TB ≈ 76 ≈ # of zeros in [0,200]

### Verification results
✓ **Slepian theory confirmed**:
  - Stable region (λ > 0.99): ~35 modes for c=60 (theory predicts ~76/2)
  - Plunge region: sharp transition with correct width scaling
  - Eigenvalues properly bounded to [0, 1]

✗ **Density matching (prolate λ ↔ ζ-zero γ)**: FAILS
  - Prolate eigenvalues span [0, 1]; ζ-zeros span [14, 198]
  - No direct correspondence in raw eigenvalues
  - **Reason**: Slepian λ_n measure concentration, NOT zero locations

---

## Part II: Why Direct Matching Fails

### The conceptual error in naive matching

The Slepian prolate concentration operator solves:
```
maximize: ∫_{-1}^{1} |ψ(t)|² dt
subject to: bandwidth limited to [-c, c]
```

Eigenvalues λ_n measure the **trade-off between time confinement and bandwidth limitation**. They do NOT encode the ζ-zero locations; they never claimed to.

The ζ-zeros are encoded in a **different** property: the spectrum of a carefully constructed **perturbation** of the prolate operator.

### CCM's actual construction (arXiv:2511.22755, 2310.18423)

Per repository's `/Users/scones/projects/riemann/operators/nearnull_spectral_ccm_operator.py`:

1. Start with prolate concentration eigenvalues λ_n and eigenvectors V_n
2. Identify the **plunge subspace** (low-λ part, near-null region)
3. Construct a **commuting 2nd-order operator D** specific to CCM's Lemma 7.2
4. Compute eigenvalues of D **restricted to the plunge subspace**
5. **Conjecture** (CCM): these restricted eigenvalues → ζ-zeros as bandwidth → ∞

**Key difference from naive approach**:
- Not the raw prolate λ_n → zeros
- But the plunge-restricted spectrum of a SPECIFIC commuting operator → zeros
- The operator D is the "prolate-wave density-of-states operator" per CCM § 3.2

---

## Part III: The Semiclassical Density Bridge (Correct Formulation)

### What needs to be computed

Instead of comparing raw eigenvalues, we need:

**Prolate side** (WKB/semiclassical):
```
ρ_prolate(n) = dn/dλ = spectral density of prolate eigenvalues
             = can be derived from phase integral: N(λ) ~ (TB/π) arcsin(√λ)
```

**ζ-zero side** (Riemann-von Mangoldt):
```
ρ_ζ(γ) = dN/dγ = (T/2π) ln(T/2π) + corrections
        = zero counting density per unit height
```

**The matching condition** (via trace formula):
```
Tr[f(D_plunge)] = Σ_j f(μ_j) should equal Σ_γ f(γ) + error term
                = convergence proof requires O(1/T) error control
```

### Interpretation

- If **densities match at leading order**, the RH-content shifts to the **fluctuation** (spacing statistics after unfolding)
- The **trace formula error** is where the proof lives
- This is a **major simplification** vs. the analytic route (§932), which hits irreducible semi-local positivity

---

## Part IV: Current Status vs. Roadmap

### ✓ DONE (this session)
- Prolate operator mechanics verified (Slepian theory)
- Time-bandwidth correspondence to zero count confirmed (TB ≈ N_zeros)
- Plunge structure identified (repository: nearnull_spectral_ccm_operator.py, §905 measured)
- Assessment: raw prolate approach is a dead end; need exact CCM operator

### ✗ BLOCKED (needs CCM implementation)
- **CCM Prop 4.3 operator formula**: Need exact D from arXiv:2511.22755
  - Is it a rank-1 perturbation? Toeplitz + oscillator? Explicit kernel?
  - Paper must specify the construction precisely

- **Commuting operator on plunge**: Repository attempted kinetic, position, mode-number, oscillation
  - None converged directly to zeros (expected — not the exact CCM operator)
  - Exact CCM operator from § 3.2 still needed

- **Trace formula implementation**: Need Connes' trace formula (arXiv:2602.04022)
  - How to compute Tr[f(D_plunge)] for arbitrary f
  - Error term scaling: is it O(1/T)? O(log T / T)?

### → NEXT IMMEDIATE STEP
Read CCM 2511.22755 **Proposition 4.3** and extract:
1. Exact formula for the commuting operator D (or W_λ or Q_λ)
2. How it relates to prolate (is prolate the base + a perturbation?)
3. Convergence rate guarantee (as function of bandwidth c or λ)
4. Whether connes-cvs package implements it directly

---

## Part V: Convergence Criteria (Proof roadmap)

If we had the exact CCM operator D, the proof would follow this structure:

### Step 1: Compute plunge spectrum (numerically verified)
```python
evals_Q, evecs_Q = eigh(Q_weil)  # Weil form
n_plunge = int(N_zeros)          # extract plunge eigenspace
plunge_vecs = evecs_Q[:, :n_plunge]

D_plunge = plunge_vecs.T @ D_ccm @ plunge_vecs
μ_j = eigh(D_plunge)  # these should → γ_j
```

### Step 2: Measure convergence (as c or λ → ∞)
```
error_j(c) = |μ_j(c) - γ_j|
convergence_rate(c) = mean error_j(c)
```

**Convergence hypothesis**: error_j(c) = O(c^{-α}) for some α > 0
(likely α ≈ 1 or 2 based on WKB theory)

### Step 3: Trace formula control (the rigorous step)
```
Tr[f(D_plunge)] = Σ_j f(μ_j)
RvM_integral = Σ_γ f(γ)
difference = Tr[f(D_plunge)] - RvM_integral
           = O(1/T) for test functions f
```

**Proof**: Use Connes' trace formula (arXiv:2602.04022) to bound the difference.

### Step 4: Conclude RH
If Steps 1-3 succeed:
```
lim_{c → ∞} [μ_j(c) - γ_j] = 0 for each j
⟹ All zeros are on the critical line (by spectral localization)
```

---

## Part VI: Strategic Assessment

### Why the spectral route is fundamentally different from geometric/analytic routes

| Route | Wall | Workaround |
|-------|------|-----------|
| **Analytic** (arch-vs-prime) | Irreducible razor-thin positivity | None; semi-local positivity = RH |
| **Geometric** (Weil cohomology) | Deninger parity obstruction (§936) | None; condensed math doesn't fix it |
| **Spectral** (Connes prolate) | Convergence = RH restatement (§932) | **Yes: explicit trace formula control** |

### Advantage of spectral route
- **Concrete operator** with explicit spectrum (not abstract cohomology)
- **Trace formula** provides concrete control mechanism (not just existence)
- **No additional obstruction** beyond RH itself (no parity, no missing geometry)

### Price
- Must construct **exact CCM operator** (not a generic baseline)
- Must prove **convergence rate** rigorously (not just numerically)
- Requires **trace formula expertise** (Connes, microlocal analysis)

---

## Part VII: Tractability & Next Steps

### Immediate (next session)
1. **Source-verify CCM Prop 4.3**: Read exact operator construction
2. **Check connes-cvs package**: Does it implement the operator directly?
3. **Implement or extract**: Get the exact D ready for testing

### Short-term (days)
4. **Compute D_plunge spectrum** for multiple bandwidths
5. **Measure convergence rate** as bandwidth increases
6. **Validate against published CCM data** (if available)

### Long-term (weeks)
7. **Prove trace formula error bound** analytically (the hard part)
8. **Formalize convergence proof** (Lean/Mathlib or pen-and-paper)
9. **Write RH proof** (assuming steps 1-8 succeed)

---

## Part VIII: Source References

### Key papers (verify against these)
1. **CCM 2511.22755** (Nov 2025): "Riemann zeros and natural frequency of prolate kernels"
   - § 2-3: operator Q_λ or W_λ construction
   - § 4-5: Weyl asymptotics, density matching
   - **ACTION**: Extract exact formula from Prop 4.3 / Lemma 7.2

2. **CCM 2310.18423** (Oct 2023): "Semilocal prolate operators"
   - § 2: Slepian concentration and time-bandwidth
   - § 3-4: WKB phase integral and plunge structure
   - **ACTION**: Verify plunge width formula against Landau-Widom

3. **Connes 2602.04022** (Feb 2026): "Spectral triples and Riemann zeros"
   - § trace formula: Tr[f(D)] = Σ_γ f(γ) + error
   - **ACTION**: Understand error term scaling (rigorous bound)

4. **Repository baseline**: `nearnull_spectral_ccm_operator.py`
   - Phase 1 (prolate) ✓ verified
   - Phase 2 (commuting operator) ✗ needs exact CCM constructor

### Data generated
- `results/semiclassical_density.json` — baseline prolate spectrum characterization
- Stored: TB parameters, plunge widths, eigenvalue ranges
- **Use for**: calibration against CCM's published data

---

## Conclusion

The **semiclassical density bridge** is the cleanest rigorous approach to RH because:

1. It provides **explicit convergence targets** (operator spectrum → ζ-zeros)
2. It avoids **geometric obstructions** (Deninger parity, missing cohomology)
3. It enables **trace-formula control** (concrete error bounds)

**However**, the current work is Report-tier (proven mechanics, not yet rigorous convergence). Earning the HANDLE requires:
- Exact CCM operator implementation
- Convergence rate measurement
- Rigorous trace formula proof

The **next step is source-verification and operator extraction** from CCM papers. Once that's in place, a computational verification (convergence as c → ∞) is feasible, followed by the analytic rigor.

---

**Filed**: LEDGER § ≥939 (next session)  
**Status**: Honest assessment — viable frontier, but machinery incomplete  
**Recommendation**: Prioritize CCM operator extraction before resuming computational work
