import BlaschkeLandau
import ZetaFactorizationAtCenter
import SmulFactorEqSplit
import EqOfCodiscreteContinuous
import SplitFunAnalytic
import GFullAnalytic
import AnalyticAtZeta
import OneLeDivisorZeta
import ZetaZeroMemFilter
import DivisorFilterReLt

/-!
# S3030 — ζ instantiation of the one-radius Blaschke–Landau lemma (lane 1, κ 2 → 1 program, brick 4b)

For 1 < Re c, R < R', 1 ∉ ball c R', ‖ζ‖ ≤ Ms < Mζ on sphere c R, and a zero ρ of ζ in ball c R:

  Re(−ζ'/ζ(c)) ≤ 2·(log Mζ − log‖ζ(c)‖)/R − (Re 1/(c − ρ) + Re(ρ − c)/R²).

Compare the pipeline's two-radii lemma (hzero_Tz_sh): 2·log(…)/(R/2) − Re 1/(c − ρ) — the growth term is HALVED, at the cost
of the tiny correction Re(ρ − c)/R² = −(σ − β)/R².  Built from: the §1343 factorization on ball c R', split (S2xxx
smul_factor_eq_split) into zeros in ball c R (→ Blaschke product, centred) and the rest (→ g₀), translated to centre 0.
-/

open Complex Metric Filter Classical

theorem zeta_landau_blaschke {c ρ : ℂ} {R R' Ms Mζ : ℝ} (hR : 0 < R) (hRR' : R < R') (hcre : 1 < c.re)
    (hball : (1 : ℂ) ∉ ball c R') (hζsphere : ∀ z ∈ sphere c R, ‖riemannZeta z‖ ≤ Ms) (hsl : Ms < Mζ)
    (hρ : ρ ∈ ball c R) (hζρ : riemannZeta ρ = 0) :
    (-logDeriv riemannZeta c).re
      ≤ 2 * (Real.log Mζ - Real.log ‖riemannZeta c‖) / R - ((1 / (c - ρ)).re + (ρ - c).re / R ^ 2) := by
  have hR' : 0 < R' := hR.trans hRR'
  obtain ⟨g, hga, hgn, hcod⟩ := exists_factorization_riemannZeta_at_center hR' hcre
  set p : ℂ → Prop := fun u => u ∈ ball c R with hp
  set T := (divisor_riemannZeta_ball_support_finite c R').toFinset with hT
  set mm : ℂ → ℕ := fun u => (MeromorphicOn.divisor riemannZeta (ball c R') u).toNat with hmm
  set Sin := T.filter p with hSin
  set Tout := T.filter (fun u => ¬ p u) with hTout
  have hcod' := hcod
  rw [smul_factor_eq_split hball p] at hcod'
  have hpt : ∀ z ∈ ball c R', riemannZeta z
      = (∏ u ∈ Sin, (z - u) ^ mm u) * ((∏ u ∈ Tout, (z - u) ^ mm u) * g z) := by
    intro z hz
    refine eq_of_codiscreteWithin_of_continuousAt isOpen_ball hz hcod' ?_ ?_
    · have hz1 : z ≠ 1 := fun h => hball (h ▸ hz)
      exact (analyticAt_riemannZeta hz1).continuousAt
    · exact (split_fun_analyticOnNhd mm hga z hz).continuousAt
  have hζc : riemannZeta c ≠ 0 := riemannZeta_ne_zero_of_one_lt_re hcre
  -- translated data
  set F : ℂ → ℂ := fun w => riemannZeta (c + w) with hF
  set h₁ : ℂ → ℂ := fun s => (∏ u ∈ Tout, (s - u) ^ mm u) * g s with hh₁
  set g₀ : ℂ → ℂ := fun w => h₁ (c + w) with hg₀
  set S : Finset ℂ := Sin.image (fun u => u - c) with hS
  set m' : ℂ → ℕ := fun a => mm (a + c) with hm'
  have hshift : ∀ w ∈ ball (0 : ℂ) R', c + w ∈ ball c R' := fun w hw => by
    rw [mem_ball, dist_eq_norm, add_sub_cancel_left]; simpa using hw
  have hfac : ∀ w ∈ ball (0 : ℂ) R', F w = (∏ a ∈ S, (w - a) ^ m' a) * g₀ w := by
    intro w hw
    simp only [hF, hg₀, hh₁]
    rw [hpt _ (hshift w hw), hS, Finset.prod_image (fun u _ v _ h => by simpa using h)]
    congr 1
    refine Finset.prod_congr rfl (fun u _ => ?_)
    simp only [hm', sub_add_cancel]
    ring_nf
  have hg₀an : AnalyticOnNhd ℂ g₀ (ball 0 R') := by
    have h1 : AnalyticOnNhd ℂ h₁ (ball c R') := prod_rest_mul_analyticOnNhd mm hga
    intro w hw
    exact (h1 (c + w) (hshift w hw)).comp (analyticAt_const.add analyticAt_id)
  have hg₀ne : ∀ w ∈ ball (0 : ℂ) R, g₀ w ≠ 0 := by
    intro w hw
    have hsR : c + w ∈ ball c R := by
      rw [mem_ball, dist_eq_norm, add_sub_cancel_left]; simpa using hw
    simp only [hg₀, hh₁]
    refine mul_ne_zero (Finset.prod_ne_zero_iff.mpr (fun u hu => pow_ne_zero _ ?_))
      (hgn ⟨c + w, ball_subset_ball hRR'.le hsR⟩)
    intro h
    have : c + w = u := sub_eq_zero.mp h
    exact (Finset.mem_filter.mp hu).2 (this ▸ hsR)
  have hSmem : ∀ a ∈ S, ∃ u ∈ Sin, a = u - c := fun a ha => by
    obtain ⟨u, hu, rfl⟩ := Finset.mem_image.mp ha; exact ⟨u, hu, rfl⟩
  have hSnorm : ∀ a ∈ S, ‖a‖ < R := by
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
  have hρSin : ρ ∈ Sin :=
    zeta_zero_mem_filter hball (ball_subset_ball hRR'.le hρ) hζρ (mem_ball_self hR') hζc p hρ
  have ha₀ : ρ - c ∈ S := Finset.mem_image.mpr ⟨ρ, hρSin, rfl⟩
  have hm1 : 1 ≤ m' (ρ - c) := by
    have := one_le_divisor_riemannZeta_at_zero hball (ball_subset_ball hRR'.le hρ) hζρ (mem_ball_self hR') hζc
    simp only [hm', sub_add_cancel, hmm]
    omega
  have hM : ∀ w ∈ sphere (0 : ℂ) R, ‖F w‖ ≤ Ms := by
    intro w hw
    apply hζsphere
    rw [mem_sphere, dist_eq_norm, add_sub_cancel_left]; simpa using hw
  have hF0 : F 0 ≠ 0 := by simpa [hF] using hζc
  have key := Blaschke.landau_blaschke m' hR hRR' hSnorm h0 (fun a ha => (hSre a ha).le) ha₀ hg₀an hg₀ne
    hfac hM hsl hF0
  -- logDeriv of the translate
  have hc1 : c ≠ 1 := by intro h; rw [h] at hcre; simp at hcre
  have hdc : DifferentiableAt ℂ riemannZeta c := (analyticAt_riemannZeta hc1).differentiableAt
  have hld : logDeriv F 0 = logDeriv riemannZeta c := by
    have hFc : F = riemannZeta ∘ (fun w => c + w) := rfl
    rw [hFc, logDeriv_comp (by simpa using hdc) (differentiableAt_const _ |>.add differentiableAt_id)]
    simp
  have hF0' : F 0 = riemannZeta c := by simp [hF]
  rw [hld, hF0'] at key
  -- drop multiplicity ≥ 1 and rewrite the zero term
  have hX : 0 ≤ (1 / (-(ρ - c))).re + (ρ - c).re / R ^ 2 :=
    Blaschke.blaschke_zero_term_nonneg hR (hSnorm _ ha₀).le (fun h => h0 (h ▸ ha₀)) (hSre _ ha₀).le
  have hneg : -(ρ - c) = c - ρ := by ring
  rw [hneg] at hX key
  have hm1' : (1 : ℝ) ≤ (m' (ρ - c) : ℝ) := by exact_mod_cast hm1
  nlinarith

#print axioms zeta_landau_blaschke
