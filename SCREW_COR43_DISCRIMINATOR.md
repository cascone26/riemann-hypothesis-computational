# The Screw-Function / Cor 4.3 Discriminator — consolidated (cycles 290–309, 2026-08-30)

A durable write-up of the 20-cycle arc. Companion to `EULER_PRODUCT_WEIL_DISCRIMINATOR.md`.
Everything below is **numerically verified** and, where marked, **Arb-certified** or **Lean-banked**;
the honest proof status is stated in §7. Proof pointers are given inline.

---

## 1. The verified equivalence chain (RH for ζ, "here")

All of the following are equivalent, and each was verified numerically this arc:

```
RH ⟺ λ_a > 0 ∀a  ⟺  g(t) ≤ 0 for t ≥ t0  ⟺  −g(t) = 4 Σ_γ sin²(γt/2)/γ² is a real sum of squares  ⟺  every zero γ is real
    (Suzuki operator)   (Suzuki Cor 4.3)      (explicit-formula identity, cycle 303)              (RH)
```

- **λ_a** = lowest eigenvalue of the Suzuki screw-kernel operator on [−a,a] (operators/suzuki_exact_Q.py,
  suzuki_exact_Q_mpmath.py). Superexponential knife-edge: λ_a ≈ C e^{−c a²}, c≈57, C≈2>0 (cycle 294);
  validated unconditional high-precision tool certifies λ_a>0 past the float64 floor (a=0.9 → 7.07e-23, Cholesky
  PD; cycle 296).
- **g(t)** = Suzuki screw function = boundary + prime + arch_core + arch_tail (eq 1.3; operators/
  suzuki_screw_lambda_a.screw_g). Robustly negative ≈ −0.04, a tiny residue of a ±12000 boundary/prime
  cancellation (cycle 302). **Arb-CERTIFIED g(t) < 0** at t=1,2,3,5,8,12 (err ≤2e-53; operators/
  screw_g_arb_certify.py, results/screw_g_arb_certified_cycle306.txt, cycle 306).
- **Explicit-formula identity** (cycle 303, results/screw_g_zeroside_cycle303.txt): g(t) = 2 Σ_γ(cos γt −1)/γ²
  = −4 Σ_γ sin²(γt/2)/γ². Zero-side matches prime-side to the 2000-zero truncation constant 8.85e-4 (exact
  in the limit).

## 2. Why the margin is robust (cycle 304)

−g(t) = 4 Σ_γ sin²(γt/2)/γ² has time-mean 2 Σ_{γ>0} 1/γ² = 0.045325 (matches the sampled mean exactly), stays
in [0.020, 0.068] over t≤50, and never reaches 0 because the γ are **incommensurable** (cannot all satisfy
sin(γt/2)=0 at once; the γ>Γ tail alone keeps the sum positive). This is why Cor 4.3 holds with a robust O(0.04)
margin — unlike the superexponential knife-edge of λ_a (same RH, far better numerical conditioning).

## 3. The Krein–Langer κ=1 unification (cycles 297–298)

The Suzuki form is Q_full = Q_0 + (rank-1 positive pole lift), where Q_0 (arch+prime, no pole) carries
**exactly one** negative eigenvalue μ_1<0. Rank-1 secular equation ⇒
`λ_a>0 ⟺ ⟨p,φ_1⟩²/|μ_1| > 1 + Σ_{k≥2}⟨p,φ_k⟩²/μ_k` (operators/bs_schur_threshold.py, high-precision extraction
cycle 298: p_1²/|μ_1| → 1 saturation, β→0, margin superexp-positive). The SAME κ=1 pole-lift governs the Epstein
PASS/FAIL discriminator (operators/mult_polelift_bridge.py): Q_0 has exactly one negative eigenvalue for
MULTIPLICATIVE L (κ=1 lift ⇒ PSD), MORE for non-mult. → the Suzuki wall and the EP discriminator are one object.

## 4. Dirichlet-form / Schrödinger backbone, reduction reconciled (cycles 299–300)

For prime-free a (<½log2≈0.347) the arch form is the H^{1/2} Schrödinger operator (−Δ)^{1/2} − ½log(a²−x²)
on [−a,a] (operators/suzuki_mu1_closed_form.py). Reduction constant reconciled empirically (5 points, direct
validated λ_a vs inf L_a, results/suzuki_La_reduction_reconciled_cycle300.txt):
**λ_a = ½·(inf[L_a] − (1+A)) + o(1)**, a→0 limit lands on −½(1+A)=−0.854 exactly (A=½(log2π−ψ(2))=0.70755).
The research-queue script's docstring "λ_a=2 inf L_a − 2(1+A)" was 4× off (coefficient ½ not 2) — corrected.
So small-a λ_a>0 ⟺ the fractional Schrödinger ground state exceeds 1+A=1.70755.

## 5. The mechanism, machine-checked BOTH ways (Lean bricks, cycles 307, 309)

The sign of each screw term is decided by whether its frequency (a zero) is real or complex:

| | frequency | screw term | weighted sum | verdict |
|---|---|---|---|---|
| **PASS** (`ScrewSumOfSquares`) | real γ | cos(γt)−1 = **−2 sin²(γt/2) ≤ 0** | Σ wᵢ(cos−1) ≤ 0 | g≤0 / mult / RH |
| **FAIL** (`ScrewOfflineNonneg`) | complex (shift δ) | cosh(δt)−1 = **+2 sinh²(δt/2) ≥ 0** (>0 if δt≠0) | Σ wᵢ(cosh−1) ≥ 0 | g>0 / non-mult |

Both are axiom-clean ({propext, Classical.choice, Quot.sound}) in ~/rhforge/forge/ (mirrored to rhforge/forge/,
in LEDGER.md). The sign flip −2sin² ↔ +2sinh² IS the on-line/off-line dichotomy = the multiplicativity
discriminator = Cor 4.3, verified at the proof-assistant level.

## 6. Off-line ⟹ FAIL, on a real non-multiplicative L (cycle 308)

Davenport–Heilbronn (non-mult, documented off-line zero ρ=0.808517+85.6995i, δ=0.3085): its off-line quadruple
contributes ~cosh(δt)/γ² to the screw function, exceeding the ζ-margin 0.020 by t≈16.3, reaching ~14778 by
t=60. So g_{D-H}(t) > 0 → Cor 4.3 FAILS for a non-multiplicative L — exactly the `screw_offline_nonneg` case.

## 7. Honest proof status (verify-before-ship)

This arc **mapped and verified an equivalence** and **machine-checked the mechanism**; it did **not** prove RH.
- Every reduction (λ_a>0, rank-1 secular, resolvent f(0)<0, g≤0, sum-of-squares) is *equivalent to RH*; none is
  an independent lever. (Cycle 301's "⟨p,Q0⁻¹p⟩=−1" was corrected as a restatement of λ_a→0, cycle 302.)
- Proving g_ζ(t)≤0 for ALL t≥t0 from the prime+arch side (no zero knowledge) IS RH — it means controlling the
  ±12000 boundary/prime cancellation to O(0.01), the explicit-formula/PNT-error content.
- The off-line-zero exclusion from certified g≤0 up to finite T (cycle 305) is a genuine detector but
  classical-strength (weaker than the rigorous verification of RH to height ~3×10¹², Platt et al.) — not novel.
- **The frontier**: breaking the equivalence needs a NEW unconditional idea for the prime-side bound, or the
  multi-agent Roundtable escalation (needs Jacob's opt-in).

## 8. Proof pointers

- Code: operators/{suzuki_exact_Q,suzuki_exact_Q_mpmath,suzuki_screw_lambda_a,screw_g_arb_certify,
  bs_schur_threshold,suzuki_mu1_closed_form,mult_polelift_bridge}.py
- Results: results/{screw_g_sign_cycle302,screw_g_zeroside_cycle303,offline_detector_cycle305,
  screw_g_arb_certified_cycle306,bs_schur_kreinlanger_cycle297,suzuki_La_reduction_reconciled_cycle300,
  f0_resolvent_cycle301,suzuki_mpmath_lambda_a_cycle296}.txt
- Lean bricks (Olympus forge, axiom-clean): ScrewSumOfSquares, ScrewOfflineNonneg
- Full cycle-by-cycle log: THE_BRIDGE.md cycles 290–309; session summaries STATUS.md updates 3–21.
