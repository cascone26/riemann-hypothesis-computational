import TrigPolyQuarticSharp
import GrowthConvexClean
import HBexplicitOfGeSh
import HgrowthRiemannZetaScalarSh
import HzeroAtPointSh

/-!
# S3020 — quartic zero-free region on the SHARP Borel–Carathéodory chain: κ = 2 + 2⁻¹², c₀ ≈ 1.432e-2 (lane 1)

S3018 proved the sharp BC derivative bound (constant 2, via Schwarz); S3020 regenerates the ten-file growth chain
(§1378 … §1444, `hgrowth_riemannZeta_scalar_sh`, `TrigPolyRegion.hzero_Tz_sh`) with it, so every growth term is
2·log(…)/(R/2) instead of 4·log(…)/(R/2).  On the convexity sphere bound (S3008) the growth scalar becomes
convRHS2 = 2·log(2+|t|) + (4/R)(log(C+1) + log(log(2+|t|)+6) + log(3+1/(σ−1))), leading coefficient 2 (was 4).
With R = 1/2, κ = 8193/4096, log|γ| ≥ 2³⁶:
β ≤ 1 − (25063328/18763)/(14379A₀ + (191019795/2048)·log|γ|);  clean: β ≤ 1 − (25063328/1750081299)/log|γ| ≈ 1 − 0.014321/log|γ|.
(S3016: 7.161e-3 — DOUBLED.)
-/

open Complex Metric

namespace TrigPolyQuarticSharpBC

noncomputable def convRHS (C σ t R : ℝ) : ℝ :=
  2 * Real.log (2 + |t|) + 4 / R * (Real.log (C + 1) + Real.log (Real.log (2 + |t|) + 6)
    + Real.log (3 + 1 / (σ - 1)))

theorem scalar_le_convRHS {C σ t R : ℝ} (hC : 0 ≤ C) (hσ : 1 < σ) (hR0 : σ - 1 ≤ R) (hR : 0 < R) :
    2 * Real.log ((C * (2 + |t|) ^ ((1 - σ + R) / 2) * (Real.log (2 + |t|) + 6) + 1)
      * (3 + 1 / (σ - 1))) / (R / 2) ≤ convRHS C σ t R := by
  unfold convRHS
  set X := 2 + |t| with hX
  set p := (1 - σ + R) / 2 with hp
  set Lg := Real.log X + 6 with hLg
  set D := 3 + 1 / (σ - 1) with hD
  have hX1 : 1 ≤ X := by rw [hX]; linarith [abs_nonneg t]
  have hX0 : 0 < X := by linarith
  have hp0 : 0 ≤ p := by rw [hp]; linarith
  have hpR : p ≤ R / 2 := by rw [hp]; linarith
  have hlogX : 0 ≤ Real.log X := Real.log_nonneg hX1
  have hLg6 : 6 ≤ Lg := by rw [hLg]; linarith
  have hD0 : 0 < D := by
    rw [hD]; have : 0 < 1 / (σ - 1) := by apply div_pos one_pos; linarith
    linarith
  have hXp : 1 ≤ X ^ p := Real.one_le_rpow hX1 hp0
  have hXp0 : 0 < X ^ p := by linarith
  have hM1 : C * X ^ p * Lg + 1 ≤ (C + 1) * (X ^ p * Lg) := by nlinarith
  have hMpos : 0 < C * X ^ p * Lg + 1 := by positivity
  have hlog : Real.log ((C * X ^ p * Lg + 1) * D)
      ≤ Real.log (C + 1) + p * Real.log X + Real.log Lg + Real.log D := by
    calc Real.log ((C * X ^ p * Lg + 1) * D)
        ≤ Real.log (((C + 1) * (X ^ p * Lg)) * D) :=
          Real.log_le_log (mul_pos hMpos hD0) (mul_le_mul_of_nonneg_right hM1 hD0.le)
      _ = Real.log (C + 1) + p * Real.log X + Real.log Lg + Real.log D := by
          rw [Real.log_mul (by positivity) hD0.ne', Real.log_mul (by positivity) (by positivity),
            Real.log_mul hXp0.ne' (by positivity), Real.log_rpow hX0]; ring
  have hsplit : 2 * Real.log ((C * X ^ p * Lg + 1) * D) / (R / 2)
      = 4 / R * Real.log ((C * X ^ p * Lg + 1) * D) := by
    rw [div_div_eq_mul_div]; ring
  rw [hsplit]
  have h8R : 0 ≤ 4 / R := by positivity
  have step1 := mul_le_mul_of_nonneg_left hlog h8R
  have step2 : 4 / R * (p * Real.log X) ≤ 2 * Real.log X := by
    have : 4 / R * (p * Real.log X) = Real.log X * (4 * p / R) := by ring
    rw [this]
    have hq : 4 * p / R ≤ 2 := by rw [div_le_iff₀ hR]; linarith
    nlinarith
  have hexp : 4 / R * (Real.log (C + 1) + p * Real.log X + Real.log Lg + Real.log D)
      = 4 / R * (p * Real.log X) + 4 / R * (Real.log (C + 1) + Real.log Lg + Real.log D) := by ring
  linarith

theorem Tz_le_conv {C : ℝ} (hC0 : 0 ≤ C)
    (hC : ∀ (σ₀ γ R : ℝ) (z : ℂ), 1 < σ₀ → σ₀ - 1 ≤ R → R ≤ 1 → σ₀ + R ≤ 2 → 3 ≤ |γ| →
      z ∈ sphere ((σ₀ : ℂ) + (γ : ℂ) * I) R →
      ‖riemannZeta z‖ ≤ C * (2 + |γ|) ^ ((1 - σ₀ + R) / 2) * (Real.log (2 + |γ|) + 6))
    {σ t : ℝ} (hσ : 1 < σ) (hσR : σ - 1 ≤ 1 / 2) (hσ2 : σ + 1 / 2 ≤ 2) (ht3 : 3 ≤ |t|) :
    (-logDeriv riemannZeta ((σ : ℂ) + I * (t : ℂ))).re ≤ convRHS C σ t (1 / 2) := by
  set c : ℂ := (σ : ℂ) + I * (t : ℂ) with hc
  have hcre_eq : c.re = σ := by rw [hc]; simp [Complex.add_re, Complex.mul_re]
  have hcre : 1 < c.re := by rw [hcre_eq]; exact hσ
  have hball : (1 : ℂ) ∉ ball c 1 := by rw [hc]; exact one_notMem_ball_test (by linarith)
  have hc' : c = (σ : ℂ) + (t : ℂ) * I := by rw [hc]; ring
  have key := hgrowth_riemannZeta_scalar_sh
    (Ms := C * (2 + |t|) ^ ((1 - σ + 1 / 2) / 2) * (Real.log (2 + |t|) + 6))
    (B := 2 * Real.log ((C * (2 + |t|) ^ ((1 - σ + 1 / 2) / 2) * (Real.log (2 + |t|) + 6) + 1)
      * (3 + 1 / (c.re - 1))) / (1 / 2 / 2)) (by norm_num : (0 : ℝ) < 1 / 2) (by norm_num) hcre hball
    (fun z hz => by
      rw [hc'] at hz
      exact hC σ t (1 / 2) z hσ hσR (by norm_num) hσ2 ht3 hz) le_rfl
  rw [hcre_eq] at key
  exact le_trans key (scalar_le_convRHS hC0 hσ hσR (by norm_num))

set_option maxHeartbeats 2000000 in
theorem riemannZeta_zero_free_quartic_sharpbc :
    ∃ A₀ C : ℝ, 0 ≤ A₀ ∧ 0 ≤ C ∧ ∀ γ β : ℝ,
      Real.exp 68719476736 ≤ |γ| →
      8192 * (2 * Real.log 6 + 8 * Real.log (C + 1) + 8 * Real.log 2 + 8 * Real.log (3 + 4 * A₀ + 44)) ≤ Real.log |γ| →
      riemannZeta ((β : ℂ) + I * (γ : ℝ)) = 0 →
      β ≤ 1 - (25063328 / 18763) / (14379 * A₀ + 191019795 / 2048 * Real.log |γ|) := by
  obtain ⟨A₀, hA0, H⟩ := TrigPolyRegion.zeta_region_P_pole
  obtain ⟨C, hC1, hC⟩ := ConvexSphereMs.convex_sphere_Ms
  have hC0 : 0 ≤ C := by linarith
  refine ⟨A₀, C, hA0, hC0, ?_⟩
  intro γ β hxe hthr hζρ₀
  set L := Real.log |γ| with hLdef
  have hL : 68719476736 ≤ L := by
    calc (68719476736 : ℝ) = Real.log (Real.exp 68719476736) := (Real.log_exp 68719476736).symm
      _ ≤ L := Real.log_le_log (Real.exp_pos 68719476736) hxe
  have hγ3 : 3 ≤ |γ| := le_trans (by
      have := Real.add_one_le_exp (68719476736 : ℝ); linarith) hxe
  have hγpos : 0 < |γ| := by linarith
  set κL : ℝ := 8193 / 4096 * L with hκL
  have hκL0 : 0 ≤ κL := by rw [hκL]; positivity
  set 𝒦 := 14379 * A₀ + 191019795 / 2048 * L with h𝒦
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
  have hll : 16 * Real.log L ≤ L / 8192 := by
    have := mul_log_le_self (C := 131072) (x := L) (by norm_num) (by norm_num; linarith)
    linarith
  have hLpos : 0 < L := by linarith
  have hCA : 0 < 3 + 4 * A₀ + 44 := by linarith
  have hD : Real.log (3 + 1 / (σ - 1)) ≤ Real.log (3 + 4 * A₀ + 44) + Real.log L := by
    rw [hinv, ← Real.log_mul hCA.ne' hLpos.ne']
    have hpos : 0 < 3 + 𝒦 / 4384 := by positivity
    refine Real.log_le_log hpos ?_
    rw [h𝒦]
    have : (14379 * A₀ + 191019795 / 2048 * L) / 4384 ≤ 4 * A₀ + 44 * L := by
      rw [div_le_iff₀ (by norm_num)]; nlinarith
    nlinarith
  have hll2 : Real.log (2 * L) = Real.log 2 + Real.log L := Real.log_mul (by norm_num) hLpos.ne'
  have habs : ∀ t : ℝ, 3 ≤ |t| → |t| ≤ 4 * |γ| → convRHS C σ t (1 / 2) ≤ κL := by
    intro t ht htγ
    unfold convRHS
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
    rw [show (4 : ℝ) / (1 / 2) = 8 by norm_num, hκL]
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
    exact le_trans (Tz_le_conv hC0 hC hσ1 (by linarith) (by linarith) h1) (habs _ h1 h2)
  -- B-term on the convexity sphere bound
  set Ms : ℝ := C * (2 + |γ|) ^ ((1 - σ + 1 / 2) / 2) * (Real.log (2 + |γ|) + 6) with hMs
  have hMs0 : 0 ≤ Ms := by
    rw [hMs]
    have := Real.log_nonneg (show (1:ℝ) ≤ 2 + |γ| by linarith [abs_nonneg γ])
    positivity
  have hMpos : 0 < Ms + 1 := by linarith
  have hcre : ((σ : ℂ) + I * (γ : ℝ)).re = σ := by simp
  have hBexp : 2 * Real.log ((Ms + 1) / ‖riemannZeta ((σ : ℂ) + I * (γ : ℝ))‖) / (1 / 2 / 2) ≤ κL := by
    refine hBexplicit_of_ge_sh (by rw [hcre]; exact hσ1) hMpos (by norm_num) ?_
    rw [hcre]
    exact le_trans (scalar_le_convRHS hC0 hσ1 (by linarith) (by norm_num)) (habs γ hγ3 (by linarith))
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
      exact TrigPolyRegion.hzero_Tz_sh (by norm_num : (0 : ℝ) < 1 / 2) (by norm_num) hσ1 hβσ
        (one_notMem_ball_test (by linarith)) hρball hζρ₀ hp
        (fun z hz => by
          have hz' : z ∈ sphere ((σ : ℂ) + (γ : ℂ) * I) (1 / 2) := by
            rw [show ((σ : ℂ) + (γ : ℂ) * I) = (σ : ℂ) + I * (γ : ℝ) by ring]; exact hz
          exact hC σ γ (1 / 2) z hσ1 (by linarith) (by norm_num) (by linarith) hγ3 hz')
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

theorem riemannZeta_zero_free_quartic_sharpbc_clean :
    ∃ T₀ : ℝ, ∀ γ β : ℝ, T₀ ≤ |γ| →
      riemannZeta ((β : ℂ) + I * (γ : ℝ)) = 0 → β ≤ 1 - (25063328 / 1750081299) / Real.log |γ| := by
  obtain ⟨A₀, C, hA₀nn, hC0, hmain⟩ := riemannZeta_zero_free_quartic_sharpbc
  set Cst := 8192 * (2 * Real.log 6 + 8 * Real.log (C + 1) + 8 * Real.log 2 + 8 * Real.log (3 + 4 * A₀ + 44)) with hCst
  refine ⟨max (max (Real.exp 68719476736) (Real.exp Cst)) (Real.exp (14379 * A₀)), ?_⟩
  intro γ β hγ hζ
  have h1 : Real.exp 68719476736 ≤ |γ| := le_trans (le_trans (le_max_left _ _) (le_max_left _ _)) hγ
  have hexp2 : Real.exp Cst ≤ |γ| := le_trans (le_trans (le_max_right _ _) (le_max_left _ _)) hγ
  have hexp3 : Real.exp (14379 * A₀) ≤ |γ| := le_trans (le_max_right _ _) hγ
  have h2 : Cst ≤ Real.log |γ| := by
    calc Cst = Real.log (Real.exp Cst) := (Real.log_exp _).symm
      _ ≤ Real.log |γ| := Real.log_le_log (Real.exp_pos _) hexp2
  have hL : 68719476736 ≤ Real.log |γ| := by
    calc (68719476736 : ℝ) = Real.log (Real.exp 68719476736) := (Real.log_exp 68719476736).symm
      _ ≤ Real.log |γ| := Real.log_le_log (Real.exp_pos 68719476736) h1
  have hA : 14379 * A₀ ≤ Real.log |γ| := by
    calc 14379 * A₀ = Real.log (Real.exp (14379 * A₀)) := (Real.log_exp _).symm
      _ ≤ Real.log |γ| := Real.log_le_log (Real.exp_pos _) hexp3
  have hmain' := hmain γ β h1 h2 hζ
  have hKpos : 0 < 14379 * A₀ + 191019795 / 2048 * Real.log |γ| := by nlinarith
  have hcmp : 14379 * A₀ + 191019795 / 2048 * Real.log |γ| ≤ 93273 * Real.log |γ| := by nlinarith
  have hdiv : (25063328 / 18763) / (93273 * Real.log |γ|)
      ≤ (25063328 / 18763) / (14379 * A₀ + 191019795 / 2048 * Real.log |γ|) :=
    div_le_div_of_nonneg_left (by norm_num) hKpos hcmp
  have he : (25063328 / 1750081299 : ℝ) / Real.log |γ|
      = (25063328 / 18763) / (93273 * Real.log |γ|) := by
    field_simp; ring
  rw [he]; linarith

end TrigPolyQuarticSharpBC

#print axioms TrigPolyQuarticSharpBC.scalar_le_convRHS
#print axioms TrigPolyQuarticSharpBC.Tz_le_conv
#print axioms TrigPolyQuarticSharpBC.riemannZeta_zero_free_quartic_sharpbc
#print axioms TrigPolyQuarticSharpBC.riemannZeta_zero_free_quartic_sharpbc_clean
