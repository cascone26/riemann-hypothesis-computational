import Mathlib.Analysis.SpecialFunctions.Complex.Circle
import Mathlib.Analysis.Meromorphic.Divisor

/-!
# §1375 — product lower bound `∏‖s-ρ‖^{m ρ} ≥ r^{Σ m ρ}`

If every factor `‖s - ρ‖ ≥ r ≥ 0` over a finite set `S`, then the weighted product is `≥ r^(Σ m ρ)`.  Applied on
the `R`-sphere with `r = R/2` (§1374) and `S` = the half-disk zeros, this is the lower bound `‖P_half‖ ≥ (R/2)^n`
that §1373 turns into the sphere bound `‖g‖ ≤ Mζ/(R/2)^n`.  EXIT0 / std-axioms.
-/

theorem prod_norm_sub_ge {S : Finset ℂ} {s : ℂ} {r : ℝ} (hr : 0 ≤ r) (m : ℂ → ℕ)
    (h : ∀ ρ ∈ S, r ≤ ‖s - ρ‖) :
    r ^ (∑ ρ ∈ S, m ρ) ≤ ∏ ρ ∈ S, ‖s - ρ‖ ^ (m ρ) := by
  rw [← Finset.prod_pow_eq_pow_sum]
  apply Finset.prod_le_prod
  · intro ρ _; exact pow_nonneg hr _
  · intro ρ hρ; exact pow_le_pow_left₀ hr (h ρ hρ) _

#print axioms prod_norm_sub_ge
