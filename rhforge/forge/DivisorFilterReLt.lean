import DivisorSupportReLeOne

/-!
# §1408 — filtered divisor-support points satisfy `Re ρ < σK` (§1406's `hleft`)

Every (filtered) divisor-support point has `Re ρ ≤ 1` (§1354), so `Re ρ < σK` whenever `1 < σK` — exactly §1406's
`hleft` at `σK = 1+1/(2K) > 1`.  EXIT0 / std-axioms.
-/

open Classical

theorem divisor_filter_support_re_lt {c : ℂ} {R' σK : ℝ} (hσK : 1 < σK) (p : ℂ → Prop) {ρ : ℂ}
    (hρ : ρ ∈ (divisor_riemannZeta_ball_support_finite c R').toFinset.filter p) : ρ.re < σK :=
  lt_of_le_of_lt (divisor_support_re_le_one (Finset.mem_of_mem_filter ρ hρ)) hσK

#print axioms divisor_filter_support_re_lt
