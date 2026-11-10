import CCMTailExact

/-!
# S3007 — Pick-function form of zero-side positivity (lane 2)

The far tail is Löwner with generator g(ν) = π·Σ_j w_j/(ν − x_j) (S2972: g = π·Cw) and confluent diagonal.  Positivity of the
tail is therefore a statement about the GENERATOR being a positive combination of Cauchy poles — a Pick/Nevanlinna property:
* `lowner_cauchy_psd` — w ≥ 0, (ν_i − x_j)C_ij = 1, ν injective ⇒ for every x:
  xᵀ·lowner(diag(C diag(w) Cᵀ), π·Cw, ν)·x ≥ 0.
Reading (Report): RH above height T ⇔ the far-tail generator −g_tail/π = Σ_γ w_γ/(γ − ν) has nonnegative pole weights
(w_γ = 4sin²(γL/2)/L ≥ 0 iff γ real) ⇔ it is a Pick function on the window.  Löwner's theorem (operator monotonicity ⇔
PSD Löwner matrices) is the classical bridge; the arithmetic side g_A − π·C_K w_K must therefore be a Pick function —
the lane-2 positivity target is now a one-variable function-theoretic statement.
-/

open Matrix Finset

namespace CCMPick

variable {n m : Type*} [Fintype n] [Fintype m] [DecidableEq n] [DecidableEq m]

theorem lowner_cauchy_psd (ν : n → ℝ) (xz : m → ℝ) (C : Matrix n m ℝ) (w : m → ℝ)
    (hC : ∀ i j, (ν i - xz j) * C i j = 1) (hν : Function.Injective ν) (hw : ∀ j, 0 ≤ w j)
    (x : n → ℝ) :
    0 ≤ x ⬝ᵥ (CCMLowner.lowner (fun i => (C * diagonal w * Cᵀ) i i) (Real.pi • (C *ᵥ w)) ν *ᵥ x) := by
  rw [← CCMCauchy.cauchy_gram_eq_lowner ν xz C w hC hν]
  exact CCMTailExact.tail_gram_nonneg C w hw x

end CCMPick

#print axioms CCMPick.lowner_cauchy_psd
