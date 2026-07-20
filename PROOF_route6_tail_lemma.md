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

---

# COMPLETE PROOF (2026-07-19) — both gaps closed, elementary and rigorous

**Theorem.** For every γ ≥ γ₁ = 14.134…, every σ ∈ (0,1), and every integer n with 1 ≤ n ≤ 2πγ − π:
F_n(σ,γ) ≥ 0.

By the σ↔1−σ symmetry of F_n, assume σ ≤ ½; then u = u(σ) ≥ 0. Two standard bounds are used throughout:
**(B1)** cos t ≤ e^{−t²/2} for t ∈ [0, π/2] (Taylor: difference has non-positive leading term), and
**(B2)** cosh x ≤ e^{x²/2} for all x. Also θ = θ(σ) ∈ (0, π/2), u = u(σ) ≥ 0, and (shown below) t := nθ ∈ (0,2π).

**Lemma 1 (θ increasing, u decreasing on (0,½]).**
dθ/dσ = γ/(σ²+γ²) − γ/((1−σ)²+γ²) > 0 (since σ² < (1−σ)²).
du/dσ = −(1−σ)/((1−σ)²+γ²) − σ/(σ²+γ²) < 0. ∎

**Lemma 2 (θ ≥ u on (0,½]).** By Lemma 1, θ−u is increasing, so its min is at σ→0: θ(0)−u(0) =
arctan(1/γ) − ½ln(1+1/γ²) = f(1/γ) with f(x)=arctan x − ½ln(1+x²). Since f(0)=0 and f'(x)=(1−x)/(1+x²) > 0
for x<1 (here x=1/γ<1), f>0. Hence θ ≥ u. ∎

**Lemma 3 (θ+u decreasing on (0,½]; so θ+u ≤ θ(0)+u(0)).**
d(θ+u)/dσ = (γ−σ)/(σ²+γ²) − (γ+1−σ)/((1−σ)²+γ²). Its sign is that of the numerator
(γ+1−σ)(σ²+γ²) − (γ−σ)((1−σ)²+γ²), which expands EXACTLY to −(γ²−γ+σ(2γ+1)−σ²). For γ≥14, σ∈(0,½]:
γ²−γ = γ(γ−1) > 0 and σ(2γ+1)−σ² > 0, so the numerator is < 0, i.e. d(θ+u)/dσ < 0. ∎

**Lemma 4 (threshold).** θ(0)+u(0) = arctan(1/γ) + ½ln(1+1/γ²) ≤ 1/γ + 1/(2γ²) = (2γ+1)/(2γ²). Hence for
all σ, θ(σ)+u(σ) ≤ (2γ+1)/(2γ²), so 2π/(θ+u) ≥ 4πγ²/(2γ+1) ≥ 2πγ − π (⟺ 4γ² ≥ 4γ²−1). ∎
Also t = nθ ≤ nθ(½) ≤ n·(1/γ) ≤ (2πγ−π)/γ < 2π, so t ∈ (0,2π). ∎

**Main argument.** F_n = 4 − 4cos(nθ)cosh(nu) (exact, Step 1). It suffices to show cos(nθ)cosh(nu) ≤ 1.
- If cos(nθ) ≤ 0: the product is ≤ 0 < 1. ✓
- If cos(nθ) > 0 and t=nθ ∈ (0,π/2]: by (B1),(B2), ln[cos(nθ)cosh(nu)] ≤ −(nθ)²/2 + (nu)²/2 =
  −(n²/2)(θ²−u²) ≤ 0 by Lemma 2 (θ≥u≥0). ✓
- If cos(nθ) > 0 and t=nθ ∈ [3π/2,2π): put τ = 2π−nθ ∈ (0,π/2]. Then cos(nθ)=cos τ, and by (B1),(B2)
  ln[cos(nθ)cosh(nu)] ≤ −τ²/2 + (nu)²/2 ≤ 0 ⟺ τ ≥ nu ⟺ n(θ+u) ≤ 2π ⟺ n ≤ 2π/(θ+u). By Lemma 4,
  2π/(θ+u) ≥ 2πγ−π ≥ n. ✓
(The complementary phase interval t∈(π/2,3π/2) has cos(nθ)≤0, already handled.)

Therefore cos(nθ)cosh(nu) ≤ 1, i.e. **F_n(σ,γ) ≥ 0** for all σ∈(0,1), γ≥14, integer n ≤ 2πγ−π. ∎∎

**Corollary (the theorem the bound needed).** Let T be a height up to which RH is verified (Platt–Trudgian
2021: T = 3.0001×10¹², unconditional finite computation). Zeros with γ ≤ T are on the line ⇒ contribute
2(1−cos nφ_γ) ≥ 0; zeros with γ > T contribute F_n(σ,γ) ≥ 0 for n ≤ 2πγ−π (above; γ>T) whenever
n ≤ 2πT−π. Hence **λ_n = Σ_ρ[1−(1−1/ρ)ⁿ] > 0 unconditionally for every integer n ≤ 2πT − π =
18,850,184,240,066 (≈ 1.885×10¹³).** Finite, not RH — but a complete, rigorous, unconditional theorem.

Every algebraic step verified numerically (closed form to 1e-12; all derivative formulas, the exact
numerator identity, f>0, and the N_C bound confirmed). Elementary: only arctan/log calculus, cos t≤e^{−t²/2},
cosh x≤e^{x²/2}. No unproven lemma remains.
