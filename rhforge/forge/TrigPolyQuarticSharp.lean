import TrigPoly5841SharpC0
import TrigPolyRegionPole

/-!
# S2971 — degree-4 trigonometric polynomial on the sharp pipeline: c₀ ≈ 3.47e-3 (lane 1)

`P(θ) = 14379 + 24480 cos θ + 14900 cos 2θ + 6000 cos 3θ + 1250 cos 4θ = (27 + 120x + 100x²)²`, x = cos θ
(i.e. 10⁴·(9/10 + x)²(3/10 + x)², a Stechkin-type quartic; grid-optimal shape ≈ (0.913+x)²(0.277+x)²).
Score ε(c₁−c₀−ε)/((c₀+ε)Σ_{k≥1}c_k) = 0.028646 vs 0.026991 for 5-8-4-1 (+6.1%).
* `pQ_nonneg` — P ≥ 0.
* `riemannZeta_zero_free_quartic_sharp` — with R = 1/2, R' = 1, ε = 4384, B = C₂ = C₃ = C₄ = (33/4)·log|γ|:
  β ≤ 1 − (25063328/18763)/(14379A₀ + (769395/2)·log|γ|).
* `riemannZeta_zero_free_quartic_clean` — β ≤ 1 − (25063328/7218126100)/log|γ| (≈ 3.4723e-3) for |γ| ≥ T₀.
Ladder: 221/68040 ≈ 3.248e-3 (S2969) → ≈ 3.472e-3 (+6.9%).
-/

open Complex Metric TrigPoly5841Sharp

namespace TrigPolyQuartic

theorem pQ_nonneg (θ : ℝ) :
    0 ≤ 14379 + 24480 * Real.cos θ + 14900 * Real.cos (2 * θ) + 6000 * Real.cos (3 * θ)
      + 1250 * Real.cos (4 * θ) := by
  have h4 : Real.cos (4 * θ) = 2 * Real.cos (2 * θ) ^ 2 - 1 := by
    rw [show 4 * θ = 2 * (2 * θ) by ring, Real.cos_two_mul]
  rw [h4, Real.cos_three_mul, Real.cos_two_mul]
  have : 14379 + 24480 * Real.cos θ + 14900 * (2 * Real.cos θ ^ 2 - 1)
      + 6000 * (4 * Real.cos θ ^ 3 - 3 * Real.cos θ)
      + 1250 * (2 * (2 * Real.cos θ ^ 2 - 1) ^ 2 - 1)
      = (27 + 120 * Real.cos θ + 100 * Real.cos θ ^ 2) ^ 2 := by ring
  rw [this]; positivity

set_option maxHeartbeats 2000000 in
theorem riemannZeta_zero_free_quartic_sharp :
    ∃ A₀ : ℝ, 0 ≤ A₀ ∧ ∀ γ β : ℝ,
      Real.exp 65536 ≤ |γ| →
      8 * (8 * Real.log 4 + 16 * Real.log 15 + 16 * Real.log (3 + 4 * A₀ + 88)) ≤ Real.log |γ| →
      riemannZeta ((β : ℂ) + I * (γ : ℝ)) = 0 →
      β ≤ 1 - (25063328 / 18763) / (14379 * A₀ + 769395 / 2 * Real.log |γ|) := by
  obtain ⟨A₀, hA0, H⟩ := TrigPolyRegion.zeta_region_P_pole
  refine ⟨A₀, hA0, ?_⟩
  intro γ β hxe hthr hζρ₀
  set L := Real.log |γ| with hLdef
  have hL : 65536 ≤ L := by
    calc (65536 : ℝ) = Real.log (Real.exp 65536) := (Real.log_exp 65536).symm
      _ ≤ L := Real.log_le_log (Real.exp_pos 65536) hxe
  have hγ3 : 3 ≤ |γ| := le_trans (by
      have := Real.add_one_le_exp (65536 : ℝ); linarith) hxe
  have hγpos : 0 < |γ| := by linarith
  set κL : ℝ := 33 / 4 * L with hκL
  have hκL0 : 0 ≤ κL := by rw [hκL]; positivity
  set 𝒦 := 14379 * A₀ + 769395 / 2 * L with h𝒦
  have hK : (![14379, 24480, 14900, 6000, 1250] : Fin 5 → ℝ) 0 * A₀
      + (![14379, 24480, 14900, 6000, 1250] : Fin 5 → ℝ) 1 * κL
      + ∑ k : Fin 3, (![14379, 24480, 14900, 6000, 1250] : Fin 5 → ℝ) k.succ.succ
          * (![κL, κL, κL] : Fin 3 → ℝ) k = 𝒦 := by
    simp [Fin.sum_univ_three, h𝒦, hκL]; ring
  have hKpos : 0 < 𝒦 := by rw [h𝒦]; nlinarith
  set σ := 1 + 4384 / 𝒦 with hσdef
  have hq : 0 < 4384 / 𝒦 := by positivity
  have hq4 : 4384 / 𝒦 ≤ 1 / 4 := by rw [div_le_iff₀ hKpos, h𝒦]; nlinarith
  have hσ1 : 1 < σ := by rw [hσdef]; linarith
  have hσ54 : σ ≤ 5 / 4 := by rw [hσdef]; linarith
  have hβσ : β < σ := by
    by_contra hcon
    push Not at hcon
    have hre : 1 < ((β : ℂ) + I * (γ : ℝ)).re := by simp; linarith
    exact riemannZeta_ne_zero_of_one_lt_re hre hζρ₀
  by_cases hfar : β < σ - 1 / 4
  · have h1 : 4384 / 𝒦 + (25063328 / 18763) / 𝒦 ≤ 1 / 4 := by
      rw [← add_div, div_le_iff₀ hKpos, h𝒦]; nlinarith
    linarith
  push Not at hfar
  have hinv : 1 / (σ - 1) = 𝒦 / 4384 := by
    rw [hσdef, show 1 + 4384 / 𝒦 - 1 = 4384 / 𝒦 by ring, one_div_div]
  have hll : 16 * Real.log L ≤ L / 8 := by
    have := mul_log_le_self (C := 128) (x := L) (by norm_num) (by norm_num; linarith)
    linarith
  have hLpos : 0 < L := by linarith
  have hCA : 0 < 3 + 4 * A₀ + 88 := by linarith
  have hD : Real.log (3 + 1 / (σ - 1)) ≤ Real.log (3 + 4 * A₀ + 88) + Real.log L := by
    rw [hinv, ← Real.log_mul hCA.ne' hLpos.ne']
    have hpos : 0 < 3 + 𝒦 / 4384 := by positivity
    refine Real.log_le_log hpos ?_
    rw [h𝒦]
    have : (14379 * A₀ + 769395 / 2 * L) / 4384 ≤ 4 * A₀ + 88 * L := by
      rw [div_le_iff₀ (by norm_num)]; nlinarith
    nlinarith
  have hK15 : Real.log (1 / (1 - (σ - 1 / 2)) + 3 + 4 / (σ - 1 / 2)) ≤ Real.log 15 := by
    have ha : 0 < 1 - (σ - 1 / 2) := by linarith
    have hb : 0 < σ - 1 / 2 := by linarith
    have h1 : 1 / (1 - (σ - 1 / 2)) ≤ 4 := by rw [div_le_iff₀ ha]; linarith
    have h2 : 4 / (σ - 1 / 2) ≤ 8 := by rw [div_le_iff₀ hb]; linarith
    refine Real.log_le_log (by positivity) ?_
    linarith
  have habs : ∀ t : ℝ, 3 ≤ |t| → |t| ≤ 4 * |γ| → sharpRHS σ t (1 / 2) ≤ κL := by
    intro t ht htγ
    unfold sharpRHS
    have htpos : 0 < |t| := by linarith
    have hlt : Real.log |t| ≤ Real.log 4 + L := by
      rw [hLdef, ← Real.log_mul (by norm_num) hγpos.ne']
      exact Real.log_le_log htpos htγ
    rw [show (8 : ℝ) / (1 / 2) = 16 by norm_num, hκL]
    nlinarith [hlt, hK15, hD, hll, hthr]
  have hkγ : ∀ k : ℕ, 1 ≤ k → k ≤ 4 → 3 ≤ |(k : ℝ) * γ| ∧ |(k : ℝ) * γ| ≤ 4 * |γ| := by
    intro k hk1 hk4
    have hk1' : (1 : ℝ) ≤ k := by exact_mod_cast hk1
    have hk4' : (k : ℝ) ≤ 4 := by exact_mod_cast hk4
    rw [abs_mul, abs_of_nonneg (by linarith : (0 : ℝ) ≤ k)]
    constructor <;> nlinarith
  have hT : ∀ k : ℕ, 2 ≤ k → k ≤ 4 → TrigPolyRegion.Tz σ γ k ≤ κL := by
    intro k hk2 hk4
    obtain ⟨h1, h2⟩ := hkγ k (by omega) hk4
    exact TrigPoly5841Sharp.Tz_le_sharpRHS hσ1 (by linarith) (by norm_num) (by norm_num)
      (by linarith) h1 (habs _ h1 h2)
  -- B-term on the sharp sphere bound
  have hN : 1 ≤ ⌊|γ|⌋₊ := Nat.le_floor (by norm_num; linarith)
  set Ms : ℝ := GrowthSharp.sharpMs σ γ (1 / 2) ⌊|γ|⌋₊ with hMs
  have hMs0 : 0 ≤ Ms := TrigPoly5841AllSharp.sharpMs_nonneg hN hσ1 (by norm_num) (by norm_num)
    (by linarith) (by norm_num : (1 / 2 : ℝ) < 1) (by linarith)
  have hMpos : 0 < Ms + 1 := by linarith
  have hBexp : 4 * Real.log ((Ms + 1) / ‖riemannZeta ((σ : ℂ) + I * (γ : ℝ))‖) / (1 / 2 / 2) ≤ κL :=
    le_trans (TrigPoly5841AllSharp.hB_sharp hσ1 (by linarith) (by norm_num) (by norm_num)
      (by linarith) hγ3 hMpos) (habs γ hγ3 (by linarith))
  have hp : ((β : ℂ) + I * (γ : ℝ)) ∈ closedBall ((σ : ℂ) + I * (γ : ℝ)) (1 / 2 / 2) := by
    rw [mem_closedBall, dist_eq_norm, TrigPoly5841SharpC0.norm_sub_center', abs_sub_comm,
      abs_of_pos (by linarith)]
    linarith
  have hρball : ((β : ℂ) + I * (γ : ℝ)) ∈ ball ((σ : ℂ) + I * (γ : ℝ)) 1 := by
    rw [mem_ball, dist_eq_norm, TrigPoly5841SharpC0.norm_sub_center', abs_sub_comm,
      abs_of_pos (by linarith)]
    linarith
  have hσ' : σ = 1 + 4384 / ((![14379, 24480, 14900, 6000, 1250] : Fin 5 → ℝ) 0 * A₀
      + (![14379, 24480, 14900, 6000, 1250] : Fin 5 → ℝ) 1 * κL
      + ∑ k : Fin 3, (![14379, 24480, 14900, 6000, 1250] : Fin 5 → ℝ) k.succ.succ
          * (![κL, κL, κL] : Fin 3 → ℝ) k) := by
    rw [hK]
  have hgap := H (γ := γ) (β := β) (ε := 4384) (B := κL) (K := 3)
    ![14379, 24480, 14900, 6000, 1250] ![κL, κL, κL]
    (fun θ => by
      have := pQ_nonneg θ
      simp [Fin.sum_univ_five]; linarith)
    (by norm_num) (by norm_num) (fun k => by fin_cases k <;> norm_num) (by norm_num) hκL0
    (fun k => by fin_cases k <;> simpa using hκL0)
    (by rw [hK]; exact hKpos) (by rw [hK, h𝒦]; nlinarith) (by rw [← hσ']; exact hβσ)
    (by
      rw [← hσ']
      exact TrigPolyRegion.hzero_Tz (by norm_num : (0 : ℝ) < 1 / 2) (by norm_num) hσ1 hβσ
        (one_notMem_ball_test (by linarith)) hρball hζρ₀ hp
        (fun z hz => GrowthSharp.sphere_le_sharpMs hN hσ1 (by norm_num) (by norm_num)
          (by linarith) (by norm_num : (1 / 2 : ℝ) < 1) (by linarith) hz)
        (by linarith) hMpos hBexp)
    (by
      rw [← hσ']
      intro k
      fin_cases k
      · exact le_trans (by simpa using hT 2 (by norm_num) (by norm_num)) le_rfl
      · exact le_trans (by simpa using hT 3 (by norm_num) (by norm_num)) le_rfl
      · exact le_trans (by simpa using hT 4 (by norm_num) (by norm_num)) le_rfl)
  rw [hK] at hgap
  have e : (4384 : ℝ) * ((![14379, 24480, 14900, 6000, 1250] : Fin 5 → ℝ) 1
      - (![14379, 24480, 14900, 6000, 1250] : Fin 5 → ℝ) 0 - 4384)
      / (((![14379, 24480, 14900, 6000, 1250] : Fin 5 → ℝ) 0 + 4384) * 𝒦)
      = (25063328 / 18763) / 𝒦 := by
    simp; field_simp; ring
  rw [e] at hgap
  exact hgap

theorem riemannZeta_zero_free_quartic_clean :
    ∃ T₀ : ℝ, ∀ γ β : ℝ, T₀ ≤ |γ| →
      riemannZeta ((β : ℂ) + I * (γ : ℝ)) = 0 → β ≤ 1 - (25063328 / 7218126100) / Real.log |γ| := by
  obtain ⟨A₀, hA₀nn, hmain⟩ := riemannZeta_zero_free_quartic_sharp
  set Cst := 8 * (8 * Real.log 4 + 16 * Real.log 15 + 16 * Real.log (3 + 4 * A₀ + 88)) with hCst
  refine ⟨max (max (Real.exp 65536) (Real.exp Cst)) (Real.exp (14379 * 2 * A₀ / 5)), ?_⟩
  intro γ β hγ hζ
  have h1 : Real.exp 65536 ≤ |γ| := le_trans (le_trans (le_max_left _ _) (le_max_left _ _)) hγ
  have hexp2 : Real.exp Cst ≤ |γ| := le_trans (le_trans (le_max_right _ _) (le_max_left _ _)) hγ
  have hexp3 : Real.exp (14379 * 2 * A₀ / 5) ≤ |γ| := le_trans (le_max_right _ _) hγ
  have h2 : Cst ≤ Real.log |γ| := by
    calc Cst = Real.log (Real.exp Cst) := (Real.log_exp _).symm
      _ ≤ Real.log |γ| := Real.log_le_log (Real.exp_pos _) hexp2
  have hL : 65536 ≤ Real.log |γ| := by
    calc (65536 : ℝ) = Real.log (Real.exp 65536) := (Real.log_exp 65536).symm
      _ ≤ Real.log |γ| := Real.log_le_log (Real.exp_pos 65536) h1
  have hA : 14379 * 2 * A₀ / 5 ≤ Real.log |γ| := by
    calc 14379 * 2 * A₀ / 5 = Real.log (Real.exp (14379 * 2 * A₀ / 5)) := (Real.log_exp _).symm
      _ ≤ Real.log |γ| := Real.log_le_log (Real.exp_pos _) hexp3
  have hmain' := hmain γ β h1 h2 hζ
  have hKpos : 0 < 14379 * A₀ + 769395 / 2 * Real.log |γ| := by nlinarith
  have hcmp : 14379 * A₀ + 769395 / 2 * Real.log |γ| ≤ 384700 * Real.log |γ| := by nlinarith
  have hdiv : (25063328 / 18763) / (384700 * Real.log |γ|)
      ≤ (25063328 / 18763) / (14379 * A₀ + 769395 / 2 * Real.log |γ|) :=
    div_le_div_of_nonneg_left (by norm_num) hKpos hcmp
  have he : (25063328 / 7218126100 : ℝ) / Real.log |γ|
      = (25063328 / 18763) / (384700 * Real.log |γ|) := by
    field_simp; ring
  rw [he]; linarith

end TrigPolyQuartic

#print axioms TrigPolyQuartic.pQ_nonneg
#print axioms TrigPolyQuartic.riemannZeta_zero_free_quartic_sharp
#print axioms TrigPolyQuartic.riemannZeta_zero_free_quartic_clean
