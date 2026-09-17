import ZetaDivisorZeroReGtOne
import ZetaDivisorFinite

/-!
# §1354 — divisor-support points satisfy `Re ρ ≤ 1` (§1336's `hleft`)

Every point `ρ` of the ζ divisor-support on `ball c R` (a zero or pole of ζ) has `Re ρ ≤ 1`: contrapositive of §1347
(`divisor = 0` whenever `Re > 1`, since ζ is analytic and nonvanishing there).  Feeds §1336's `hleft`
(`∀ ρ ∈ S, ρ.re < 1 + 1/(2K)`) — with `K > 0`, `Re ρ ≤ 1 < 1 + 1/(2K)`.  EXIT0 / std-axioms.
-/

open Complex Metric

theorem divisor_support_re_le_one {c : ℂ} {R : ℝ} {ρ : ℂ}
    (hρmem : ρ ∈ (divisor_riemannZeta_ball_support_finite c R).toFinset) : ρ.re ≤ 1 := by
  rw [Set.Finite.mem_toFinset, Function.mem_support] at hρmem
  by_contra h
  have hρU : ρ ∈ ball c R :=
    (MeromorphicOn.divisor riemannZeta (ball c R)).supportWithinDomain hρmem
  exact hρmem (divisor_riemannZeta_eq_zero_of_one_lt_re hρU (not_le.mp h))

#print axioms divisor_support_re_le_one
