# Unconditional Li-coefficient positivity — consolidated results (2026-07-19/20)

Two unconditional theorems on the Li coefficients λ_n proved this session, with the tail-lemma kernel
formally verified in Lean 4 / mathlib. Notation: ρ = ½+iγ runs over the nontrivial zeros of ζ;
λ_n = Σ_ρ [1 − (1−1/ρ)ⁿ] (Li 1997; RH ⟺ λ_n ≥ 0 ∀n). T = a height up to which RH is verified as a finite
computation — Platt–Trudgian (2021): **T = 3.0001×10¹²**, unconditional. Full derivations: LEDGER §19–§33.

---

## Theorem 1 (pointwise). *λ_n > 0 unconditionally for every integer n ≤ 2πT − π ≈ 1.885×10¹³.*

**Reduction.** Group the zeros into the verified part (γ ≤ T, all on the critical line) and the tail (γ > T).
Verified zeros contribute Σ 2(1 − cos nφ_γ) ≥ 0 with φ_γ = 2·arctan(1/2γ). Tail zeros are controlled by:

**Tail lemma (proved; Lean-verified kernel).** For a zero position σ+iγ (0<σ<1, γ ≥ γ₁=14.13…) and its
functional-equation quartet, the contribution F_n(σ,γ) ≥ 0 for every integer n ≤ 2πγ − π.

*Proof of the tail lemma.* Two exact identities — arg α(σ)=arg α(1−σ)=:θ and |α(σ)||α(1−σ)|=1, where
α(s)=(s−1+iγ)/(s+iγ) — collapse the quartet to the closed form
> **F_n(σ,γ) = 4 − 4 cos(nθ) cosh(nu)**, θ = π − arctan(γ/σ) − arctan(γ/(1−σ)), u = ½ ln[((1−σ)²+γ²)/(σ²+γ²)].
So F_n ≥ 0 ⟺ cos(nθ)cosh(nu) ≤ 1. WLOG σ ≤ ½ (F_n symmetric under σ↔1−σ), so u ≥ 0. Four elementary lemmas:
(L1) θ increasing, u decreasing on (0,½] (explicit derivatives). (L2) θ ≥ u: θ−u is increasing, min at σ→0
equals arctan(1/γ) − ½ln(1+1/γ²) = f(1/γ), f(x)=arctan x − ½ln(1+x²) > 0 since f′=(1−x)/(1+x²)>0. (L3) θ+u
decreasing: the numerator of d(θ+u)/dσ equals −(γ²−γ+σ(2γ+1)−σ²) < 0 for γ≥14. (L4) θ(0)+u(0) ≤ (2γ+1)/(2γ²),
so 2π/(θ+u) ≥ 4πγ²/(2γ+1) ≥ 2πγ−π (⟺ 4γ² ≥ 4γ²−1). With cos t ≤ e^{−t²/2} on [0,π/2] and cosh z ≤ e^{z²/2}
(both elementary), split by the sign/branch of x=nθ ∈ (0,2π): if cos x ≤ 0 done; if x ∈ (0,π/2],
ln[cos x cosh nu] ≤ −(n²/2)(θ²−u²) ≤ 0 by (L2); if x ∈ [3π/2,2π), put τ=2π−x, then ≤ −(τ²−(nu)²)/2 ≤ 0 ⟺
n(θ+u) ≤ 2π, which holds for n ≤ 2πγ−π by (L4). ∎ Threshold is sharp (first failure at n ≈ 2πγ). The kernel
`cos x·cosh w ≤ 1` and both analytic inputs are machine-checked in `RiemannLean/TailLemma.lean` +
`TailAnalytic.lean` (`lake build` clean, zero sorry).

**Combine.** For n ≤ 2πT − π: verified zeros contribute ≥ 0; tail zeros (γ > T) have n < 2πγ so F_n ≥ 0.
Hence λ_n > 0. With T = 3.0001×10¹²: **n ≤ 18,850,184,240,066.** ∎

*(Context: the previous in-repo self-contained bound was n ≤ 3732; this is ~5×10⁹× further, by using the
verified height.)*

---

## Theorem 2 (almost-all). *λ_n > 0 for all but a set of n of density → 0, throughout n ≤ ~3.2·T² ≈ 2.9×10²⁵.*

Precisely: #{n ≤ 3.2T² : λ_n ≤ 0} ≤ C·(log T)²/N(T)·(3.2T²) with density ≤ C(log T)²/N(T) → 0
(N(T) ~ (T/2π)log T). This is unconditional — no pair-correlation conjecture, no minimum-gap assumption.

**Reduction.** Past 2πT, only off-line zeros above T can make λ_n < 0, and a rigorous upper bound on their
total negative contribution (per-zero worst case via the closed form, |min_σ F_n| ≤ 4(cosh(nu₀)−1),
u₀=½ln(1+1/γ²), integrated against the Riemann–von Mangoldt density) crosses the positive mass
P_ver(n)=2Σ_{γ≤T}(1−cos nφ_γ) only at n* ≈ 3.2T². Because tail zeros have γ > T = 3×10¹² (huge), u₀ ~ 1/(2γ²)
is minuscule, so the negative mass stays below the positive mass until ~3T² (LEDGER §22, §29 for the ~T²
ceiling). The one thing to control is the positive floor P_ver(n) ≥ (1−ε)·2N(T), i.e. the exponential sum
S(n) = Σ_{γ≤T} cos(nφ_γ) ≪ N(T).

**Second-moment bound (the crux, unconditional).** Over an interval of length M ~ 3T²,
⟨S(n)²⟩ = ½ Σ_{γ,γ'}[⟨cos n(φ_γ−φ_γ')⟩ + ⟨cos n(φ_γ+φ_γ')⟩]. The diagonal gives N/2. For the off-diagonal,
|⟨cos nα⟩| ≤ min(1, 1/(M|α|)); split pairs by |γ−γ'|: **close pairs** (|γ−γ'| < γ²/M ~ ⅓) number
≤ Σ_γ N(γ,γ+⅓) ≤ Σ_γ O(log γ) = O(N log T) using ONLY the standard short-interval count N(γ+1)−N(γ)=O(log γ)
(unconditional); **far pairs** give O(N(log T)²) by summing 1/(M|φ_γ−φ_γ'|) over distance shells (each shell
holds O(log γ) zeros); the φ_γ+φ_γ' term is O(N²/T)=O(N log T). Hence
> **⟨S(n)²⟩ ≤ C·N·(log T)²** (unconditional, short-interval counts only).

**Conclude.** By Chebyshev, #{n ∈ band : |S(n)| > εN} ≤ M·⟨S²⟩/(εN)² ⇒ density ≤ C(log T)²/(ε²N(T)) → 0. For
all other n, P_ver(n) ≥ 2N(1−ε) > NegBound(n) ⇒ λ_n > 0. With T = 3.0001×10¹²: exceptional density ~10⁻¹¹,
so λ_n > 0 for essentially every n up to 2.9×10²⁵. ∎

*Numerical confirmation (real 100k zeros): max |S(n)|/N ≤ 0.10 even under adversarial search — the floor
holds pointwise with margin; the exceptional set appears empty in practice.*

---

## Status & honest scope
- Both theorems are UNCONDITIONAL and FINITE (they use the finite verified height T). They are NOT RH — the
  moment one asks for all n, the difficulty returns (Theorem 2's pointwise-every-n version needs a prime
  √-phase exponential-sum bound, §31; the range ceiling is ~T² polylog, §29).
- The tail-lemma kernel of Theorem 1 is formally machine-verified in Lean 4/mathlib.
- These are the session's concrete, secured results. The open problem (RH itself) reduces to unconditional
  uniform coercivity of the Weil form — the analytic target for a dedicated future session (see
  `connes-decisive-diagnostic`, LEDGER §34 for the criticality reframe).
