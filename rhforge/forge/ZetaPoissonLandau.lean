import ZetaBlaschkeLandau
import PoissonBlaschkeDrop
import PoissonAvgBound
import ZetaCircleLogBounds
import ZeroGap
import PoleAvoidGeom

/-!
# S3048 — zeta Poisson–Landau growth bound, kappa = 1/4 (lane 1, Poisson–Landau P5b-2, growth form)

`zeta_poisson_growth`: there is A1 >= 0 such that for 1 < sigma <= 9/8, |gamma| >= 3, with L = log(2+|gamma|):
   Re(-zeta'/zeta(sigma + i gamma)) <= L/4 + 2(A1 + log(L+6) + log(3 + 1/(sigma-1))).
Assembly: factorization on ball c 1; ZeroGap radius 3/8 < r < r' < 1/2 avoiding all zero distances; split zeros by
ball c r; translate to 0; poisson_blaschke_drop (S3042) at radius r; neg_avg_le (S3041) with hup/hlo from S3044;
2(A+D)/(pi r) <= 2(A+D) since pi r > 1.
-/

open Complex Metric Filter Classical

namespace ZetaPoisson

theorem zeta_poisson_growth : ∃ A₁ : ℝ, 0 ≤ A₁ ∧ ∀ σ γ : ℝ, 1 < σ → σ ≤ 9 / 8 → 3 ≤ |γ| →
    (-logDeriv riemannZeta ((σ : ℂ) + I * (γ : ℝ))).re
      ≤ Real.log (2 + |γ|) / 4 + 2 * (A₁ + Real.log (Real.log (2 + |γ|) + 6) + Real.log (3 + 1 / (σ - 1))) := by
  obtain ⟨A₁, hA₁, hup⟩ := ZetaCircleLog.log_zeta_upper
  refine ⟨A₁, hA₁, fun σ γ hσ hσ2 hγ => ?_⟩
  set c : ℂ := (σ : ℂ) + I * (γ : ℝ) with hc
  have hcre : 1 < c.re := by simp [hc]; exact hσ
  have hball : (1 : ℂ) ∉ ball c 1 := one_notMem_ball_test (by linarith)
  obtain ⟨g, hga, hgn, hcod⟩ := exists_factorization_riemannZeta_at_center (by norm_num : (0 : ℝ) < 1) hcre
  set T := (divisor_riemannZeta_ball_support_finite c 1).toFinset with hT
  obtain ⟨r, r', h38r, hrr', hr'12, hgap⟩ :=
    ZeroGap.exists_gap (T.image (fun u => ‖u - c‖)) (by norm_num : (3 : ℝ) / 8 < 1 / 2)
  have hr : 0 < r := by linarith
  have hr'1 : r' ≤ 1 := by linarith
  set p : ℂ → Prop := fun u => u ∈ ball c r with hp
  set mm : ℂ → ℕ := fun u => (MeromorphicOn.divisor riemannZeta (ball c 1) u).toNat with hmm
  set Sin := T.filter p with hSin
  set Tout := T.filter (fun u => ¬ p u) with hTout
  have hcod' := hcod
  rw [smul_factor_eq_split hball p] at hcod'
  have hpt : ∀ z ∈ ball c 1, riemannZeta z
      = (∏ u ∈ Sin, (z - u) ^ mm u) * ((∏ u ∈ Tout, (z - u) ^ mm u) * g z) := by
    intro z hz
    refine eq_of_codiscreteWithin_of_continuousAt isOpen_ball hz hcod' ?_ ?_
    · exact (analyticAt_riemannZeta (fun h => hball (h ▸ hz))).continuousAt
    · exact (split_fun_analyticOnNhd mm hga z hz).continuousAt
  have hζc : riemannZeta c ≠ 0 := riemannZeta_ne_zero_of_one_lt_re hcre
  set F : ℂ → ℂ := fun w => riemannZeta (c + w) with hF
  set h₁ : ℂ → ℂ := fun s => (∏ u ∈ Tout, (s - u) ^ mm u) * g s with hh₁
  set g₀ : ℂ → ℂ := fun w => h₁ (c + w) with hg₀
  set S : Finset ℂ := Sin.image (fun u => u - c) with hS
  set m' : ℂ → ℕ := fun a => mm (a + c) with hm'
  have hshift : ∀ {ρ : ℝ} (w : ℂ), w ∈ ball (0 : ℂ) ρ → c + w ∈ ball c ρ := fun w hw => by
    rw [mem_ball, dist_eq_norm, add_sub_cancel_left]; simpa using hw
  have hfac : ∀ w ∈ ball (0 : ℂ) 1, F w = (∏ a ∈ S, (w - a) ^ m' a) * g₀ w := by
    intro w hw
    simp only [hF, hg₀, hh₁]
    rw [hpt _ (hshift w hw), hS, Finset.prod_image (fun u _ v _ h => by simpa using h)]
    congr 1
    refine Finset.prod_congr rfl (fun u _ => ?_)
    simp only [hm', sub_add_cancel]
    ring_nf
  have hg₀an : AnalyticOnNhd ℂ g₀ (ball 0 1) := by
    have h1 : AnalyticOnNhd ℂ h₁ (ball c 1) := prod_rest_mul_analyticOnNhd mm hga
    intro w hw
    exact (h1 (c + w) (hshift w hw)).comp (analyticAt_const.add analyticAt_id)
  have hg₀ne : ∀ w ∈ ball (0 : ℂ) r', g₀ w ≠ 0 := by
    intro w hw
    have hsr' : c + w ∈ ball c r' := hshift w hw
    simp only [hg₀, hh₁]
    refine mul_ne_zero (Finset.prod_ne_zero_iff.mpr (fun u hu => pow_ne_zero _ ?_))
      (hgn ⟨c + w, ball_subset_ball hr'1 hsr'⟩)
    intro h
    have hcu : c + w = u := sub_eq_zero.mp h
    have hu' := Finset.mem_filter.mp hu
    have hnotr : ¬ ‖u - c‖ < r := by
      intro hlt; apply hu'.2; simp only [hp, mem_ball, dist_eq_norm]; exact hlt
    have hlt' : ‖u - c‖ < r' := by
      rw [← hcu, add_sub_cancel_left]; simpa using hw
    rcases hgap _ (Finset.mem_image.mpr ⟨u, hu'.1, rfl⟩) with h1 | h1 <;> linarith [not_lt.mp hnotr]
  have hSmem : ∀ a ∈ S, ∃ u ∈ Sin, a = u - c := fun a ha => by
    obtain ⟨u, hu, rfl⟩ := Finset.mem_image.mp ha; exact ⟨u, hu, rfl⟩
  have hSnorm : ∀ a ∈ S, ‖a‖ < r := by
    intro a ha
    obtain ⟨u, hu, rfl⟩ := hSmem a ha
    have := (Finset.mem_filter.mp hu).2
    simpa [p, mem_ball, dist_eq_norm] using this
  have hSre : ∀ a ∈ S, a.re < 0 := by
    intro a ha
    obtain ⟨u, hu, rfl⟩ := hSmem a ha
    have := divisor_filter_support_re_lt hcre p hu
    simp; linarith
  have h0 : (0 : ℂ) ∉ S := fun h => by have := hSre 0 h; simp at this
  have hF0 : F 0 ≠ 0 := by simpa [hF] using hζc
  have hPB := PoissonBlaschke.poisson_blaschke_drop m' hr hrr' hr'1 hSnorm h0 (fun a ha => (hSre a ha).le)
    hg₀an hg₀ne hfac hF0
  set L := Real.log (2 + |γ|) with hL
  set A := A₁ + Real.log (L + 6) with hA
  set D := Real.log (3 + 1 / (σ - 1)) with hD
  have hL0 : 0 ≤ L := Real.log_nonneg (by linarith [abs_nonneg γ])
  have hA0 : 0 ≤ A := by have := Real.log_nonneg (by linarith : (1 : ℝ) ≤ L + 6); linarith
  have hD0 : 0 ≤ D := Real.log_nonneg (by
    have : 0 < 1 / (σ - 1) := one_div_pos.mpr (by linarith)
    linarith)
  have hsph_ball : ∀ w ∈ sphere (0 : ℂ) r, c + w ∈ ball c 1 := fun w hw =>
    hshift w (sphere_subset_closedBall.trans (closedBall_subset_ball (by linarith)) hw)
  have hFne : ∀ w ∈ sphere (0 : ℂ) r, F w ≠ 0 := by
    intro w hw h
    have hz := hsph_ball w hw
    have hmem := zeta_zero_mem_divisor_support hball hz h (mem_ball_self (by norm_num)) hζc
    have hn : ‖(c + w) - c‖ = r := by rw [add_sub_cancel_left]; simpa using hw
    rcases hgap _ (Finset.mem_image.mpr ⟨c + w, hmem, rfl⟩) with h1 | h1 <;> linarith
  have hcont : ContinuousOn (fun w => Real.log ‖F w‖) (sphere (0 : ℂ) r) := by
    intro w hw
    have hz := hsph_ball w hw
    have hcz : ContinuousAt F w :=
      ((analyticAt_riemannZeta (fun h => hball (h ▸ hz))).continuousAt).comp
        (continuousAt_const.add continuousAt_id)
    exact (hcz.norm.log (norm_ne_zero_iff.mpr (hFne w hw))).continuousWithinAt
  have hupw : ∀ w ∈ sphere (0 : ℂ) r, Real.log ‖F w‖ ≤ A + L / 2 * max (-w.re) 0 := by
    intro w hw
    have hwr : ‖w‖ = r := by simpa using hw
    have hwn : ‖w‖ ≤ 1 / 2 := by linarith
    have := hup σ γ w hσ (by linarith) hγ hwn
    simp only [hF]; rw [hA]; linarith
  have hlow : ∀ w ∈ sphere (0 : ℂ) r, 0 ≤ w.re → -D ≤ Real.log ‖F w‖ := fun w _ hw =>
    ZetaCircleLog.log_zeta_lower hσ hw
  have hAvg := PoissonAvg.neg_avg_le hr hA0 hD0 hL0 hcont hupw hlow
  have hc1 : c ≠ 1 := by intro h; rw [h] at hcre; simp at hcre
  have hdc : DifferentiableAt ℂ riemannZeta c := (analyticAt_riemannZeta hc1).differentiableAt
  have hld : logDeriv F 0 = logDeriv riemannZeta c := by
    have hFc : F = riemannZeta ∘ (fun w => c + w) := rfl
    rw [hFc, logDeriv_comp (by simpa using hdc) (differentiableAt_const _ |>.add differentiableAt_id)]
    simp
  rw [hld] at hPB
  have hπr : 1 ≤ Real.pi * r := by nlinarith [Real.pi_gt_three]
  have hfrac : 2 * (A + D) / (Real.pi * r) ≤ 2 * (A + D) := by
    rw [div_le_iff₀ (by linarith)]; nlinarith
  linarith

end ZetaPoisson

#print axioms ZetaPoisson.zeta_poisson_growth
