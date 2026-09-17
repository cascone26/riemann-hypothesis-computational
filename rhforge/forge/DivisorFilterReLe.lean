import DivisorSupportReLeOne

/-!
# §1403 — filtered divisor-support points satisfy `Re ρ ≤ Re z` (for `1 ≤ Re z`)

Every point of the (filtered) ζ divisor-support has `Re ρ ≤ 1` (§1354), so `Re ρ ≤ Re z` whenever `1 ≤ Re z` — in
particular at the growth point `z = σ+2iγ` (`Re z = σ > 1`).  This discharges §1399's `hzre` (needed so
`Re(1/(z-ρ)) ≥ 0` via §1396), giving the growth bound `hgrowth`.  EXIT0 / std-axioms.
-/

open Classical

theorem divisor_filter_support_re_le {c : ℂ} {R : ℝ} {z : ℂ} (hz : 1 ≤ z.re) (p : ℂ → Prop)
    {ρ : ℂ} (hρ : ρ ∈ (divisor_riemannZeta_ball_support_finite c R).toFinset.filter p) :
    ρ.re ≤ z.re :=
  le_trans (divisor_support_re_le_one (Finset.mem_of_mem_filter ρ hρ)) hz

#print axioms divisor_filter_support_re_le
