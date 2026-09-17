import ZetaFactorization
import Mathlib.NumberTheory.LSeries.Dirichlet

/-!
# §1423 — ζ factorization at a centre with `1 < Re c` (removes the z₀ bookkeeping)

Specializes §1343 `exists_factorization_riemannZeta` to `z₀ := c` for any disk `ball c R` (`R>0`) centred at a
point with `1 < Re c`: the centre lies in the ball (`mem_ball_self`), `c ≠ 1` (since `Re c > 1 ≠ 1`), and
`ζ c ≠ 0` (`riemannZeta_ne_zero_of_one_lt_re`).  This supplies §1410's factorization data `hga/hgn/hfac_cod`
with no auxiliary nonvanishing-point hypothesis — the last mechanical wiring step before the scalar fixed-point.
EXIT0 / std-axioms.
-/

open Complex Metric Filter

theorem exists_factorization_riemannZeta_at_center {c : ℂ} {R : ℝ} (hR : 0 < R) (hc : 1 < c.re) :
    ∃ g : ℂ → ℂ, AnalyticOnNhd ℂ g (ball c R) ∧ (∀ u : ball c R, g u ≠ 0) ∧
      riemannZeta =ᶠ[codiscreteWithin (ball c R)]
        (∏ᶠ u, (· - u) ^ (MeromorphicOn.divisor riemannZeta (ball c R) u)) • g := by
  refine exists_factorization_riemannZeta (mem_ball_self hR) ?_ (riemannZeta_ne_zero_of_one_lt_re hc)
  intro h
  rw [h] at hc
  simp at hc

#print axioms exists_factorization_riemannZeta_at_center
