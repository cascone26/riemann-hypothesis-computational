import DivisorSupportReLeOneLFunction

/-!
# §1490 — filtered L divisor-support points satisfy `Re ρ ≤ Re z` (for `1 ≤ Re z`)

L-analogue of ζ §1403.  Every point of the (filtered) `LFunction χ` divisor-support has `Re ρ ≤ 1` (§1489), so
`Re ρ ≤ Re z` whenever `1 ≤ Re z` — in particular at the growth point `z = σ+2iγ`.  This discharges the `hzre`
hypothesis of the drop-zeros growth bound §1399.  EXIT0 / std-axioms.
-/

open Classical DirichletCharacter

theorem divisor_filter_support_re_le_LFunction {N : ℕ} [NeZero N] {χ : DirichletCharacter ℂ N}
    (hχ : χ ≠ 1) {c : ℂ} {R : ℝ} {z : ℂ} (hz : 1 ≤ z.re) (p : ℂ → Prop) {ρ : ℂ}
    (hρ : ρ ∈ (divisor_LFunction_ball_support_finite hχ c R).toFinset.filter p) :
    ρ.re ≤ z.re :=
  le_trans (divisor_support_re_le_one_LFunction hχ (Finset.mem_of_mem_filter ρ hρ)) hz

#print axioms divisor_filter_support_re_le_LFunction
