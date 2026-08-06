# The Guinand–Weil form as a universal Selberg-class Gram matrix

**Status:** verified numerical result + machine-checked structural formalization. **No proof of RH or any
GRH.** This document records what was established (an identity, its consequences, and finite rigorous
certificates) and states honestly what remains open (the positivity from arithmetic alone, which *is* the
respective Riemann hypothesis). Session LEDGER §767–785 (2026-08-06). All claims have a code/proof pointer.

> **Relation to prior work — this is not a novelty claim (LEDGER §785).** After completing the below, a
> literature check found that the substantive findings here are **already published**, independently, by
> **Suzuki (arXiv:2606.09096)** — the screw-function self-adjoint operator whose quadratic form *is* the
> Weil form — and **Kim et al (arXiv:2607.24830)**, whose results include: the **archimedean-law
> universality across the Selberg class** with offset `B₀ = log q − 2log2` (their R2, verified for
> Dirichlet characters and modular forms — covering §3's χ₃/χ₄ and Δ cases); **off-line-zero detection**
> with residual blow-up rate `→2δ` (their R6 — §3's detector); and the **Paley–Wiener zero-density
> ceiling** (their R8 — the Nyquist edge). The "Weil form = Gram matrix of zeros" identity below is the
> classical explicit formula written as a bilinear form. So this document is an **independent
> computational reproduction** — via a different apparatus (Gaussian-cosine test functions + a Gram-matrix
> framing, vs their finite-element / operator-spectrum route) — plus a Lean formalization of the abstract
> structure. That two independent methods agree is a real cross-validation; the **discoveries are theirs,
> not ours.** The genuinely local artifacts are the machine-checked Lean development (Bessel functions —
> absent from mathlib — and the abstract Gram theorems) and the rigorous Arb certificates.

---

## 1. The object

For a test function `φ_a(x) = e^{−x²/2} cos(ax)` (frequency `a > 0`), with Fourier transform
`φ̂_a(r) = √(2π)/2 · (e^{−(r−a)²/2} + e^{−(r+a)²/2})`, the **Guinand–Weil form** of an L-function `L(s)` on
a finite basis `{a_1,…,a_n}` is the `n×n` matrix

```
Q_{jk} = arch(a_j,a_k) + pole(a_j,a_k) − prime(a_j,a_k)
```

- `arch(a,b) = (1/2π) ∫ φ̂_a(r) φ̂_b(r) Ω(r) dr`, where `Ω(r)` is the archimedean density of `L`
  (`Ω_ζ(r) = Re ψ(¼+ir/2) − log π`; for a general L-function the conductor/weight/degree enter here).
- `pole(a,b) = 2 Re[φ̂_a(i/2) φ̂_b(i/2)]` (present only when `L` has a pole, i.e. `ζ`; zero for cusp forms
  and non-principal characters).
- `prime(a,b) = 2 Σ_n Λ_L(n)/√n · g(a,b,log n)`, the L-function's von Mangoldt / Satake prime sum against
  the closed-form cross-correlation kernel `g`.

## 2. The identity (LEDGER §771–772)

**`Q` is the Gram matrix of the zero-evaluation vectors.** Writing `v_j = (√2 · φ̂_{a_j}(γ))_{γ>0}` over the
non-trivial zeros `ρ = ½ + iγ` of `L`,

```
Q_{jk} = 2 Σ_{γ>0} φ̂_{a_j}(γ) φ̂_{a_k}(γ) = ⟨v_j, v_k⟩ = (VᵀV)_{jk}.
```

The left side is computed **from the arithmetic** (digamma integral + primes, no zeros); the right side is
the Gram matrix of the zeros. Verified to machine precision (max |Q_arith − Q_gram| ≈ 3×10⁻¹⁵ for ζ; code
`operators/weil_gram_identity.py`). The diagonal `D(a) = Q_{aa} = 2 Σ_{γ>0} φ̂_a(γ)²` is a **sum of squares**
(LEDGER §771).

**Consequence.** `Q` is positive-semidefinite ⟺ the `v_j` are real vectors ⟺ every `γ` is real ⟺ the
Riemann hypothesis for `L`. An off-line zero `ρ = ½ + ε + iγ` injects the *complex* coordinate
`φ̂((ρ−½)/i) = φ̂(γ − iε)`, turning `Q` into a Gram matrix under a bilinear (not positive) form.

## 3. Three consequences (all verified)

1. **Nyquist certification edge (§774, §776).** `Q = VᵀV` is *positive-definite* (rigorously certifiable,
   min-eig > 0) iff the columns of `V` are linearly independent — i.e. iff the basis is **coarser than the
   mean zero spacing** `2π/log(T/2π)`. A denser (sub-Nyquist) basis makes `V` rank-deficient ⇒ only PSD
   (min-eig → 0⁺), and rigorous certification saturates. The saturation edge = the zero-density Nyquist
   limit, confirmed across windows and heights. Formalized in Lean (`weil_gram_posDef`).

2. **Rigorous finite certificates with teeth (§775, §780).** The form is certified PD on finite bases by
   Arb interval arithmetic + exact-rational LDL^T — a ladder 4×4 → 6×6 → **11×11** (the last at the Nyquist
   edge). These are not merely "RH-consistent": the form **detects** an off-line zero straddled by the
   basis — a 2-function basis at a zero is rank-1 PSD on-line and goes **indefinite off-line**
   (min-eig −1.36 at ε=0.5, growing with ε). So a PD certificate is an *in-band* verification that no
   off-line zero is resolved by that basis's frequency band.

3. **Universality across the Selberg class (§781–784).** The identity and all of the above are **not
   special to ζ**. Only `Ω` (conductor/weight/degree) and the prime data (character twist / Satake
   parameters) change:
   - **ζ** — degree 1 (RH). `Q_arith = Q_gram` to 3×10⁻¹⁵.
   - **Dirichlet `L(s,χ₃)`, `L(s,χ₄)`** — degree 1 (GRH). `Ω_χ = log(q/π) + Re ψ(¾+ir/2)`, χ-twisted
     primes, no pole. Diagonal ratios ≈ 1.000; conductor-dependence `log(q/π)` confirmed on two conductors.
   - **Ramanujan Δ (weight 12)** — degree 2, GL(2) automorphic. `τ(n)` from the q-expansion,
     `Ω_Δ = Re ψ(6+ir) − log(2π)` (from `Γ_ℂ(s+11/2)`), Satake power sums `b_{p^k}=α^k+β^k` via the
     Chebyshev recursion, and the `Γ_ℂ` degree-2 normalization. **`Q_arith = Q_gram`, ratio 1.0000 exact.**
   The off-line detector is L-function-agnostic (the Δ detection numbers are *identical* to ζ's — detection
   depends only on `φ̂` and `ε`).

## 4. What is proven vs open

**Proven / verified:**
- The Gram identity `Q_arith = Q_gram` (the explicit formula), machine-precision, for every L-function tested.
- Finite-basis positive-definiteness, rigorously (Arb + exact LDL^T), on the certificate ladder.
- The detector: an off-line zero forces the finite form indefinite.
- The abstract structural theorems in Lean (`RiemannLean/WeilGram.lean`, trusted-kernel axioms only):
  any real `VᵀV` is PSD; PD iff `V.mulVec` injective; hence `Q = VᵀV` real ⟹ PSD.

**Open (this *is* the Riemann hypothesis):**
- Proving `Q ⪰ 0` for **all** test functions from the **arithmetic side alone** (no assumption on zero
  locations). Every reduction attempted on this family collapses to RH-completeness — including a
  midpoint-factorization shortcut that was **tested and refuted** before banking (§779): `h ≥ 0` does not
  make `e^{−(a−b)²/4} h((a+b)/2)` positive-definite; the true diagonal `D` is PD-compatible (slowly
  varying) *because* the zeros are real, which is exactly the content to be proven.
- The single identified analytic frontier is a Gonek-uniform prime-sum resummation (LEDGER §732), on which
  both the certificate side and the sign-theorem side converge.

## 5. Reproduce

- Identity / Gram: `operators/weil_gram_identity.py` (ζ), `weil_dirichlet_gram.py` (χ₃/χ₄),
  `weil_degree2_delta.py` (Δ).
- Nyquist edge: `operators/weil_saturation_scout.py`. Rigorous certificates: `operators/organ2_certify.py`,
  `organ2_certify_6x6.py`, `organ2_certify_N.py` (parametrized, exact-arb basis).
- Detector: inline in the session log (§780/§784).
- Lean: `RiemannLean/WeilGram.lean` + `./audit_modules.sh` (39 theorems, 0 sorries, trusted axioms only).
- Full provenance: `LEDGER.md` §767–784, `PROCESS.md`.

**Honest one-line summary.** The Guinand–Weil positivity form is, for every Selberg-class L-function tested,
exactly the Gram matrix of that L-function's zeros; it is finitely certifiable and genuinely detects
off-line zeros; and proving its positivity from arithmetic alone remains the Riemann hypothesis.
