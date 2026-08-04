# The Crown Readout, Derived — mechanism note (2026-08-03, §629–632)

Object: G(u) = (1/N)Σ_k (−1)^k e^{iθ(γ_k)} e^{iuγ_k}  (= −⟨e^{−iπS(γ_k)}e^{iuγ_k}⟩, §614).
Claim (numerically established): G(log n) reads μ(n) with full Selberg-class universality.

## The derivation chain (each step machine-verified)
1. **Oscillator model** (§629): S(t) ≈ Σ_p a_p sin(t·m·log p + ·), a_{p,m} = 1/(π m p^{m/2}).
   Bessel expansion of e^{−iπS}: coefficient at frequency u=log n (n squarefree = Πp):
   (−1)^{ω(n)} Π_{p|n} J₁(1/√p) · Π_{q∤n} J₀(1/√q).
   → μ(n) sign EXACT, support {log n} EXACT, per-factor ladder J₁/J₀ ≈ ½ (measured 0.62 with
   sampling correction), background B(X) = ΠJ₀ ~ (log X)^{−1/4}.
2. **Sampling enhancement** (§631, third by-parts identity — family of §610/§621):
   ∫e^{−iπS}e^{iut}dS = (u/π)∫e^{−iπS}e^{iut}dt + bdry  ⟹  E(u) = 1 + u/θ̄′ + θ′-corr.
   T-trend verified; n^{1/4} coincidence rejected.
3. **Compensation** (§632): ⟨e^{−iπS_X(γ)}e^{iuγ}⟩ is X-STABLE → 97% of G_true at X=10⁴,
   while the dt-average decays (log X)^{−1/4}. The zeros' explicit-formula knowledge of each
   oscillator cancels the J₀ decay exactly. The compensation = dN-oscillator cross-correlation
   (Montgomery-type diagonal).
4. **Fine structure** (§615c): prime-power channels = two-order perturbation series of the same
   expansion (c₂=−0.377, c₃=+0.098, rms 9e−4) — the Bessel model's subleading orders.

## Remaining rigor gaps (named, classical)
(i) Selberg off-diagonal moments → the ~3–15% residuals and oscillator-correlation control.
(ii) The cross-correlation lemma: ⟨dN·oscillator⟩ diagonal — Montgomery/explicit-formula computation.
(iii) Error terms / uniformity in u, T.

## Obstruction map (why this wasn't in the literature)
The weight |ζ′|/ζ′ is non-analytic — no contour representation (§627). Phase-only/alternating
Gonek weights unclaimed (checked: 2601.18025, 1811.04530). The comb is LONG-RANGE ORDER
(§628: C(d) plateaus) — an ergodic-mean statement, not a pair-correlation statement.

## CAPSTONE (§633–634): the zero-free double expansion
⟨e^{iγv}⟩ over zeros = Landau's theorem (unconditional; verified 4 decimals). Therefore
G(u) = Σ_{k} [Π_n J_{k_n}(c_n)] × [1 if u+Σk_nL_n=0; else −(T/2πN)Λ(e^{|·|})/e^{|·|/2}].
Evaluated with NO zero data (u=log2, X=200, order ≤2): −0.3138 vs true −0.335 (94%).
The sampling enhancement E(u) = the Landau comb ⊛ Bessel expansion, term-by-term.
All inputs classical: S-explicit-formula oscillators + Bessel + Landau. Theorem path:
multi-index tail control (factorial Bessel decay) + Landau error terms + Selberg off-diagonal.

## RESUMMED FORM (§649 — the final structure)
A_zeros(u,X) = (1/K)[ ∫e^{iΨ}(θ′/π)dt + [Se^{iΨ}] − i∫S·Ψ′e^{iΨ}dt ],  Ψ = −πS_X + ut  (exact).
Verified vs measurement at all X to ~0.001. Three boxes: (i) smooth — stationary phase;
(ii) S_smooth·Ψ′ — convergent prime-pair bilinear (the E(u) enhancement, first-order EF);
(iii) sawtooth (S − S_smooth)·Ψ′ — zero−prime pair correlation, second-order EF/Montgomery
(deficit measured: −0.0006/−0.0054/−0.0270/−0.0476 at X=300/10³/3·10³/10⁴). Items (i),(ii)
classical-unconditional; (iii) is where pair-correlation inputs (and any conditionality) enter.
