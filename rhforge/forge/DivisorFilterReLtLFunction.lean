import DivisorSupportReLeOneLFunction

/-!
# §1494 — filtered L divisor-support points satisfy `Re ρ < σK` (`χ ≠ 1`; the region `hleft`)

L-analogue of ζ §1408.  Every filtered `LFunction χ` divisor-support point has `Re ρ ≤ 1` (§1489), so `Re ρ < σK`
whenever `1 < σK = 1+1/(2K)`.  This is the `hleft` hypothesis of the region composition §1336.  EXIT0 / std-axioms.
-/

open Classical DirichletCharacter

theorem divisor_filter_support_re_lt_LFunction {N : ℕ} [NeZero N] {χ : DirichletCharacter ℂ N}
    (hχ : χ ≠ 1) {c : ℂ} {R' σK : ℝ} (hσK : 1 < σK) (p : ℂ → Prop) {ρ : ℂ}
    (hρ : ρ ∈ (divisor_LFunction_ball_support_finite hχ c R').toFinset.filter p) : ρ.re < σK :=
  lt_of_le_of_lt (divisor_support_re_le_one_LFunction hχ (Finset.mem_of_mem_filter ρ hρ)) hσK

#print axioms divisor_filter_support_re_lt_LFunction
