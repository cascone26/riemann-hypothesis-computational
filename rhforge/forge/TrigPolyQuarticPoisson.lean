import TrigPolyQuarticSharp
import HzeroAtPointSh
import ZetaPoissonLandau
import MulLogLeSelf

/-!
# S3051 — quartic zero-free region on the POISSON–LANDAU bound: κ = 1/4 + 2⁻¹², c₀ = 25063328/218964210 ≈ 0.11446 (lane 1)

Inputs (S3048/S3050, no sphere sup any more): Re(−ζ'/ζ(σ+iγ)) ≤ L/4 + 2(A₁ + log(L+6) + log(3+1/(σ−1))) and the zero form
with + 2 − 1/(σ−β).  Leading coefficient 1/4 (S3032: 1, S3020: 2, S3014: 4).  κ = 1025/4096, log|γ| ≥ 2³²:
β ≤ 1 − (25063328/18763)/(14379A₀ + (23897875/2048)·log|γ|);  clean: β ≤ 1 − (25063328/218964210)/log|γ| ≈ 1 − 0.11446/log|γ|.
-/

open Complex Metric

namespace TrigPolyQuarticPoisson

set_option maxHeartbeats 2000000 in
theorem riemannZeta_zero_free_quartic_poisson :
    ∃ A₀ C : ℝ, 0 ≤ A₀ ∧ 0 ≤ C ∧ ∀ γ β : ℝ,
      Real.exp 4294967296 ≤ |γ| →
      8192 * (Real.log 6 / 4 + 2 * C + 2 * Real.log 2 + 2 * Real.log (3 + 4 * A₀ + 44) + 2) ≤ Real.log |γ| →
      riemannZeta ((β : ℂ) + I * (γ : ℝ)) = 0 →
      β ≤ 1 - (25063328 / 18763) / (14379 * A₀ + 23897875 / 2048 * Real.log |γ|) := by
  obtain ⟨A₀, hA0, H⟩ := TrigPolyRegion.zeta_region_P_pole
  obtain ⟨Ag, hAg, Hg⟩ := ZetaPoisson.zeta_poisson_growth
  obtain ⟨Az, hAz, Hz⟩ := ZetaPoisson.zeta_poisson_zero
  set C := Ag + Az with hCdef
  have hC0 : 0 ≤ C := by rw [hCdef]; linarith
  refine ⟨A₀, C, hA0, hC0, ?_⟩
  intro γ β hxe hthr hζρ₀
  set L := Real.log |γ| with hLdef
  have hL : 4294967296 ≤ L := by
    calc (4294967296 : ℝ) = Real.log (Real.exp 4294967296) := (Real.log_exp 4294967296).symm
      _ ≤ L := Real.log_le_log (Real.exp_pos 4294967296) hxe
  have hγ3 : 3 ≤ |γ| := le_trans (by
      have := Real.add_one_le_exp (4294967296 : ℝ); linarith) hxe
  have hγpos : 0 < |γ| := by linarith
  set κL : ℝ := 1025 / 4096 * L with hκL
  have hκL0 : 0 ≤ κL := by rw [hκL]; positivity
  set 𝒦 := 14379 * A₀ + 23897875 / 2048 * L with h𝒦
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
  have hq8 : 4384 / 𝒦 ≤ 1 / 8 := by rw [div_le_iff₀ hKpos, h𝒦]; nlinarith
  have hσ98 : σ ≤ 9 / 8 := by rw [hσdef]; linarith
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
  have hll : 4 * Real.log L ≤ L / 8192 := by
    have := mul_log_le_self (C := 32768) (x := L) (by norm_num) (by norm_num; linarith)
    linarith
  have hLpos : 0 < L := by linarith
  have hCA : 0 < 3 + 4 * A₀ + 44 := by linarith
  have hD : Real.log (3 + 1 / (σ - 1)) ≤ Real.log (3 + 4 * A₀ + 44) + Real.log L := by
    rw [hinv, ← Real.log_mul hCA.ne' hLpos.ne']
    have hpos : 0 < 3 + 𝒦 / 4384 := by positivity
    refine Real.log_le_log hpos ?_
    rw [h𝒦]
    have : (14379 * A₀ + 23897875 / 2048 * L) / 4384 ≤ 4 * A₀ + 44 * L := by
      rw [div_le_iff₀ (by norm_num)]; nlinarith
    nlinarith
  have hll2 : Real.log (2 * L) = Real.log 2 + Real.log L := Real.log_mul (by norm_num) hLpos.ne'
  have habs : ∀ t : ℝ, 3 ≤ |t| → |t| ≤ 4 * |γ| →
      Real.log (2 + |t|) / 4 + 2 * (C + Real.log (Real.log (2 + |t|) + 6) + Real.log (3 + 1 / (σ - 1))) + 2 ≤ κL := by
    intro t ht htγ
    have htpos : 0 < 2 + |t| := by positivity
    have h6 : 2 + |t| ≤ 6 * |γ| := by linarith
    have hlt : Real.log (2 + |t|) ≤ Real.log 6 + L := by
      rw [hLdef, ← Real.log_mul (by norm_num) hγpos.ne']
      exact Real.log_le_log htpos h6
    have hl6 : Real.log 6 ≤ 6 := by
      have := Real.log_le_sub_one_of_pos (show (0:ℝ) < 6 by norm_num); linarith
    have hLL : Real.log (Real.log (2 + |t|) + 6) ≤ Real.log 2 + Real.log L := by
      rw [← hll2]
      have hpos : 0 < Real.log (2 + |t|) + 6 := by
        have := Real.log_nonneg (show (1:ℝ) ≤ 2 + |t| by linarith [abs_nonneg t]); linarith
      exact Real.log_le_log hpos (by linarith)
    rw [hκL]
    nlinarith [hlt, hLL, hD, hll, hthr]
  have hkγ : ∀ k : ℕ, 1 ≤ k → k ≤ 4 → 3 ≤ |(k : ℝ) * γ| ∧ |(k : ℝ) * γ| ≤ 4 * |γ| := by
    intro k hk1 hk4
    have hk1' : (1 : ℝ) ≤ k := by exact_mod_cast hk1
    have hk4' : (k : ℝ) ≤ 4 := by exact_mod_cast hk4
    rw [abs_mul, abs_of_nonneg (by linarith : (0 : ℝ) ≤ k)]
    constructor <;> nlinarith
  have hT : ∀ k : ℕ, 2 ≤ k → k ≤ 4 → TrigPolyRegion.Tz σ γ k ≤ κL := by
    intro k hk2 hk4
    obtain ⟨h1, h2⟩ := hkγ k (by omega) hk4
    have hg := Hg σ ((k : ℝ) * γ) hσ1 hσ98 h1
    have hb := habs _ h1 h2
    unfold TrigPolyRegion.Tz
    linarith
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
      have := TrigPolyQuartic.pQ_nonneg θ
      simp [Fin.sum_univ_five]; linarith)
    (by norm_num) (by norm_num) (fun k => by fin_cases k <;> norm_num) (by norm_num) hκL0
    (fun k => by fin_cases k <;> simpa using hκL0)
    (by rw [hK]; exact hKpos) (by rw [hK, h𝒦]; nlinarith) (by rw [← hσ']; exact hβσ)
    (by
      rw [← hσ']
      have hz := Hz σ β γ hσ1 hσ98 hγ3 hβσ hfar hζρ₀
      have hab := habs γ hγ3 (by linarith)
      rw [TrigPolyRegion.Tz_one_sh]
      linarith)
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

theorem riemannZeta_zero_free_quartic_poisson_clean :
    ∃ T₀ : ℝ, ∀ γ β : ℝ, T₀ ≤ |γ| →
      riemannZeta ((β : ℂ) + I * (γ : ℝ)) = 0 → β ≤ 1 - (25063328 / 218964210) / Real.log |γ| := by
  obtain ⟨A₀, C, hA₀nn, hC0, hmain⟩ := riemannZeta_zero_free_quartic_poisson
  set Cst := 8192 * (Real.log 6 / 4 + 2 * C + 2 * Real.log 2 + 2 * Real.log (3 + 4 * A₀ + 44) + 2) with hCst
  refine ⟨max (max (Real.exp 4294967296) (Real.exp Cst)) (Real.exp (14379 * A₀)), ?_⟩
  intro γ β hγ hζ
  have h1 : Real.exp 4294967296 ≤ |γ| := le_trans (le_trans (le_max_left _ _) (le_max_left _ _)) hγ
  have hexp2 : Real.exp Cst ≤ |γ| := le_trans (le_trans (le_max_right _ _) (le_max_left _ _)) hγ
  have hexp3 : Real.exp (14379 * A₀) ≤ |γ| := le_trans (le_max_right _ _) hγ
  have h2 : Cst ≤ Real.log |γ| := by
    calc Cst = Real.log (Real.exp Cst) := (Real.log_exp _).symm
      _ ≤ Real.log |γ| := Real.log_le_log (Real.exp_pos _) hexp2
  have hL : 4294967296 ≤ Real.log |γ| := by
    calc (4294967296 : ℝ) = Real.log (Real.exp 4294967296) := (Real.log_exp 4294967296).symm
      _ ≤ Real.log |γ| := Real.log_le_log (Real.exp_pos 4294967296) h1
  have hA : 14379 * A₀ ≤ Real.log |γ| := by
    calc 14379 * A₀ = Real.log (Real.exp (14379 * A₀)) := (Real.log_exp _).symm
      _ ≤ Real.log |γ| := Real.log_le_log (Real.exp_pos _) hexp3
  have hmain' := hmain γ β h1 h2 hζ
  have hKpos : 0 < 14379 * A₀ + 23897875 / 2048 * Real.log |γ| := by nlinarith
  have hcmp : 14379 * A₀ + 23897875 / 2048 * Real.log |γ| ≤ 11670 * Real.log |γ| := by nlinarith
  have hdiv : (25063328 / 18763) / (11670 * Real.log |γ|)
      ≤ (25063328 / 18763) / (14379 * A₀ + 23897875 / 2048 * Real.log |γ|) :=
    div_le_div_of_nonneg_left (by norm_num) hKpos hcmp
  have he : (25063328 / 218964210 : ℝ) / Real.log |γ|
      = (25063328 / 18763) / (11670 * Real.log |γ|) := by
    field_simp; ring
  rw [he]; linarith

end TrigPolyQuarticPoisson

#print axioms TrigPolyQuarticPoisson.riemannZeta_zero_free_quartic_poisson
#print axioms TrigPolyQuarticPoisson.riemannZeta_zero_free_quartic_poisson_clean
