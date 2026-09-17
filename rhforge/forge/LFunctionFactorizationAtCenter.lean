import LFunctionFactorization
import Mathlib.NumberTheory.LSeries.Nonvanishing

/-!
# §1501 — L factorization on `ball c R'` from `1 < Re c` (center is the nonvanishing point)

L-analogue of ζ §1423.  When `1 < Re c`, `LFunction χ c ≠ 0` (nonvanishing on `Re ≥ 1`,
`LFunction_ne_zero_of_one_le_re`), so the centre `c` itself serves as the nonvanishing point `z₀` of §1476,
yielding the divisor factorization on `ball c R'` with no external nonvanishing-point hypothesis.  This is the
`g₀`-supplier the sphere-internalized region (the L analogue of ζ §1424) consumes.  EXIT0 / std-axioms.
-/

open Complex DirichletCharacter Metric Filter

theorem exists_factorization_LFunction_at_center {N : ℕ} [NeZero N] {χ : DirichletCharacter ℂ N}
    (hχ : χ ≠ 1) {c : ℂ} {R' : ℝ} (hR' : 0 < R') (hcre : 1 < c.re) :
    ∃ g : ℂ → ℂ, AnalyticOnNhd ℂ g (ball c R') ∧ (∀ u : ball c R', g u ≠ 0) ∧
      LFunction χ =ᶠ[codiscreteWithin (ball c R')]
        (∏ᶠ u, (· - u) ^ (MeromorphicOn.divisor (LFunction χ) (ball c R') u)) • g :=
  exists_factorization_LFunction hχ (mem_ball_self hR')
    (LFunction_ne_zero_of_one_le_re χ (Or.inl hχ) hcre.le)

#print axioms exists_factorization_LFunction_at_center
