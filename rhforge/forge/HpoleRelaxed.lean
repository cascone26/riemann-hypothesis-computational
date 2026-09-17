import HpoleRiemannZeta

/-!
# §1427 — pole bound relaxation: `A` a FREE parameter above the pole constant

§1426's `hpole` slot wants `(-logDeriv ζ σ).re ≤ 1/(σ-1) + A` with `A` a chosen constant, but §1405 delivers it
with `A = ‖logDeriv (LFunctionTrivChar₁ 1) σ‖`.  Relaxing to any `A` above that norm decouples `A` (which sits
inside `K = 3A+4B+C`, hence inside `σ = σK`) from the pole term — mirror of §1411's `B` decoupling.  §1360 bounds
`‖logDeriv (LFunctionTrivChar₁ 1) σ‖` by an absolute constant on the relevant compact `σ`-range, so `A` can be
taken `O(1)`.  EXIT0 / std-axioms.
-/

open Complex ArithmeticFunction LSeries DirichletCharacter

theorem hpole_relaxed {σ A : ℝ} (hσ : 1 < σ)
    (hA : ‖logDeriv (LFunctionTrivChar₁ 1) ((σ : ℂ))‖ ≤ A) :
    (-logDeriv riemannZeta ((σ : ℂ))).re ≤ 1 / (σ - 1) + A :=
  le_trans (hpole_riemannZeta hσ) (by linarith)

#print axioms hpole_relaxed
