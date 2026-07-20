# Proof notes — the Route-6 tail-safety lemma (2026-07-19)

**Goal.** Prove: for a nontrivial-zero position ρ=σ+iγ (0<σ<1, γ≥γ₁=14.13…) and its functional-equation
quartet, the Li-coefficient contribution F_n(σ,γ) ≥ 0 for every integer n < 2πγ (− O(1)). This makes the
unconditional bound λ_n > 0 for n ≤ ⌊2πT⌋ (T = verification height) fully rigorous instead of numerical.

Notation: α(s) = (s−1+iγ)/(s+iγ); g(s) = 2(1 − Re α(s)ⁿ); F_n(σ,γ) = g(σ) + g(1−σ) (the quartet;
on-line σ=½ it is the conjugate pair g(½) = 2(1−cos nφ_γ) ≥ 0).

---

## Step 1 — EXACT closed form (PROVEN; verified numerically to 1e-12)

Two exact identities:
- **(i) arg α(σ) = arg α(1−σ)** =: θ. Proof: arg α(s) = arg(s−1+iγ) − arg(s+iγ) =
  [π − arctan(γ/(1−s))] − arctan(γ/s), symmetric under s↔1−s. ∴ θ = π − arctan(γ/σ) − arctan(γ/(1−σ)).
- **(ii) |α(σ)|·|α(1−σ)| = 1.** Proof: |α(s)|² = ((s−1)²+γ²)/(s²+γ²); the product telescopes to 1.

Write |α(σ)| = e^{u}, |α(1−σ)| = e^{−u}, u = ½ ln[((1−σ)²+γ²)/(σ²+γ²)]. Then Re α(σ)ⁿ = e^{nu}cos nθ,
Re α(1−σ)ⁿ = e^{−nu}cos nθ, so

    **F_n(σ,γ) = 4 − 4 cos(nθ) cosh(nu).**       (verified: |direct − closed| < 1.1e-12)

∴ F_n ≥ 0  ⟺  **cos(nθ) cosh(nu) ≤ 1.**

---

## Step 2 — Worst case is σ→0 (VERIFIED numerically; analytic monotonicity = remaining gap #1)

The maximum of cos(nθ)cosh(nu) over σ∈(0,½] occurs at σ→0 (checked across g∈{25,100,400}, n up to 2πγ−2;
first-negativity always first appears at σ→0). At σ=0: θ₀ = arctan(1/γ) =: φ, and a clean exact relation
holds — |α(0)| = sec φ (⇔ u₀ = ln sec φ), since |α(0)|² = (1+γ²)/γ² = 1 + tan²φ = sec²φ. So the lemma
reduces to the scalar inequality

    **P(n) := cos(nφ) · cosh(n ln sec φ) ≤ 1,  φ = arctan(1/γ).**

*(Remaining gap #1: prove ∂/∂σ of cos(nθ)cosh(nu) ≥ 0 on (0,½] so the sup is at σ→0. Numerically robust.)*

---

## Step 3 — The scalar inequality (elementary; two of three regions PROVEN)

Let t = nφ (phase) and κ = κ(φ) = (ln sec φ)/φ = −(ln cos φ)/φ. Then P(n) = cos(t)·cosh(κt). Two facts:

- **κ < 1 for γ ≥ 14** (PROVEN / verified): κ = −ln cos φ / φ = φ/2 + φ³/12 + … ; for φ = arctan(1/γ) ≤
  arctan(1/14) ≈ 0.0713, κ ≤ 0.0354 ≪ 1. (Monotone, elementary.)
- Since n < 2πγ ⇒ t = nφ < 2πγ·arctan(1/γ) < 2π, so **t ∈ (0, 2π)** on the whole admissible range.

Split by the sign of cos t:
- **Region A, cos t ≤ 0** (t ∈ [π/2, 3π/2]): cosh(κt) > 0 ⇒ P = cos t·cosh κt ≤ 0 < 1.  **PROVEN.**
- **Region B, t ∈ (0, π/2)**: since κ ≤ 1, cosh(κt) ≤ cosh(t), so P ≤ cos t·cosh t ≤ 1, using the classical
  inequality **cos t · cosh t ≤ 1 on (0,π/2)** (⇔ cos t ≤ sech t; equality only at 0; verified: max = 1.000…).
  **PROVEN.**
- **Region C, t ∈ (3π/2, 2π)**: cos t > 0, the tight region. Put τ = 2π − t ∈ (0, π/2). Need
  cos τ · cosh(κ(2π−τ)) ≤ 1. Expanding to second order (φ small): cosh(κ(2π−τ)) ≈ 1 + ½κ²(2π−τ)²,
  cos τ ≈ 1 − τ²/2; product ≤ 1 ⟺ τ² ≥ κ²(2π−τ)² ⟺ τ ≥ κ(2π−τ) ⟺ τ ≥ 2πκ/(1+κ). With τ = 2π − t =
  2π − nφ and κ ≈ φ/2, this is nφ ≤ 2π − πφ, i.e. **n ≤ 2π/φ − π = 2πγ − π + O(1/γ)**. So P(n) ≤ 1 up to
  the threshold, and the first failure is at n ≈ 2πγ − π (matches numerics: first-bad n/γ → 2π⁻, gap ≈ 3 ≈ π).
  *(Remaining gap #2: replace the 2nd-order expansion in Region C with rigorous Taylor remainder bounds to
  get an explicit sub-threshold, e.g. n ≤ 2πγ − 4. Elementary but must be done carefully — Region C is tight.)*

---

## Status
- **PROVEN rigorously:** the exact closed form (Step 1); Regions A and B of the scalar inequality; κ<1.
- **Verified numerically, analytic finish remaining (2 bounded gaps):** worst-case-σ monotonicity (gap #1);
  the tight Region-C boundary estimate with explicit remainder (gap #2).
- **Consequence when closed:** F_n(σ,γ) ≥ 0 for all σ∈(0,1), γ≥14, n ≤ 2πγ − 4 (say), rigorously ⇒ the
  route-6 bound λ_n > 0 for n ≤ 2πT − O(1) ≈ 1.885×10¹³ (T = Platt–Trudgian 3.0001×10¹²) becomes a THEOREM,
  not a numerically-verified claim. Still finite, still not RH — but a complete, rigorous, unconditional result.

The opaque tail-safety lemma is now an exact closed form + an elementary scalar inequality that is 2/3
proven, with two small bounded analytic gaps. This is the concrete, finishable result flagged as "what's next."
