import TrigPoly5841Sharp
import HBexplicitOfGe

/-!
# S2967 — 5-8-4-1 region with ALL THREE growth terms (B, C₂, C₃) on the sharp N-term bound (lane 1)

S2965 moved C₂, C₃ to the sharp growth `sharpRHS` (leading 8·log|t|).  Here the iγ term B — fed by
`4·log(Mζ/|ζ(σ+iγ)|)/(R/2)` with the sphere bound Ms — is moved too, using `Ms := sharpMs σ γ R ⌊|γ|⌋₊`:
* `scalar_le_sharpRHS` — `4·log((sharpMs+1)(3+1/(σ-1)))/(R/2) ≤ sharpRHS σ t R`  (scalar core of S2963).
* `hB_sharp` — `4·log((sharpMs+1)/‖ζ(σ+iγ)‖)/(R/2) ≤ sharpRHS σ γ R`  (via `hBexplicit_of_ge`).
* `zeta_region_5841_allsharp` — region with `sharpRHS σ (kγ) R ≤ B, C₂, C₃` (k = 1, 2, 3) and only real
  side conditions (σ ≤ 2, 0 < R < R' ≤ |γ|, R < 1, σ − R < 1, |γ| ≥ 3).
All three growth inputs now have leading coefficient 8 (was 8/R for each).
-/

open Complex Metric

namespace TrigPoly5841AllSharp

open TrigPoly5841Sharp

theorem scalar_le_sharpRHS {σ t R : ℝ} (hσ : 1 < σ) (hσ2 : σ ≤ 2) (hR : 0 < R)
    (hR1 : R < 1) (haR : σ - R < 1) (ht : 3 ≤ |t|) :
    4 * Real.log ((GrowthSharp.sharpMs σ t R ⌊|t|⌋₊ + 1) * (3 + 1 / (σ - 1))) / (R / 2)
      ≤ sharpRHS σ t R := by
  have hN : 1 ≤ ⌊|t|⌋₊ := Nat.le_floor (by norm_num; linarith)
  unfold sharpRHS
  set x : ℝ := |t| with hxdef
  set a : ℝ := σ - R with hadef
  set p : ℝ := 1 - a with hpdef
  set M : ℝ := GrowthSharp.sharpMs σ t R ⌊|t|⌋₊ with hMdef
  set K : ℝ := 1 / p + 3 + 4 / a with hKdef
  set D : ℝ := 3 + 1 / (σ - 1) with hDdef
  have hapos : 0 < a := by rw [hadef]; linarith
  have hppos : 0 < p := by rw [hpdef]; linarith
  have hpR : p ≤ R := by rw [hpdef, hadef]; linarith
  have hx0 : 0 < x := by linarith
  have hlogx : 0 ≤ Real.log x := Real.log_nonneg (by linarith)
  have hK : 0 < K := by rw [hKdef]; positivity
  have hD : 0 < D := by
    rw [hDdef]; have : 0 < 1 / (σ - 1) := by apply div_pos one_pos; linarith
    linarith
  have hX : 0 < x ^ p := Real.rpow_pos_of_pos hx0 p
  have hM1 : M + 1 ≤ x ^ p * K := by
    rw [hMdef, hKdef, hpdef, hadef, hxdef]
    exact GrowthSharpClean.sharpMs_floor_add_one_le hσ hσ2 hR hR1 haR ht
  have hMpos : 0 < M + 1 := by
    have hN0 : (1 : ℝ) ≤ (⌊|t|⌋₊ : ℝ) := by exact_mod_cast hN
    have h1 : 0 ≤ ((⌊|t|⌋₊ : ℝ) ^ (1 - (σ - R)) - 1) / (1 - (σ - R)) := by
      apply div_nonneg _ (by linarith)
      have := Real.one_le_rpow hN0 (show 0 ≤ 1 - (σ - R) by linarith)
      linarith
    have h3 : 0 ≤ (⌊|t|⌋₊ : ℝ) ^ (1 - (σ - R)) / (|t| - R) :=
      div_nonneg (Real.rpow_nonneg (by linarith) _) (by linarith)
    have h4 : 0 ≤ (σ + |t| + R) * ((⌊|t|⌋₊ : ℝ) ^ (-(σ - R)) / (σ - R)) :=
      mul_nonneg (by linarith [abs_nonneg t])
        (div_nonneg (Real.rpow_nonneg (by linarith) _) (by linarith))
    rw [hMdef]; unfold GrowthSharp.sharpMs; linarith
  have hlog : Real.log ((M + 1) * D) ≤ p * Real.log x + Real.log K + Real.log D := by
    calc Real.log ((M + 1) * D)
        ≤ Real.log ((x ^ p * K) * D) :=
          Real.log_le_log (mul_pos hMpos hD) (mul_le_mul_of_nonneg_right hM1 hD.le)
      _ = p * Real.log x + Real.log K + Real.log D := by
          rw [Real.log_mul (by positivity) hD.ne', Real.log_mul hX.ne' hK.ne',
            Real.log_rpow hx0]
  have hsplit : 4 * Real.log ((M + 1) * D) / (R / 2) = 8 / R * Real.log ((M + 1) * D) := by
    rw [div_div_eq_mul_div]; ring
  rw [hsplit]
  have h8R : 0 ≤ 8 / R := by positivity
  have step1 : 8 / R * Real.log ((M + 1) * D)
      ≤ 8 / R * (p * Real.log x + Real.log K + Real.log D) :=
    mul_le_mul_of_nonneg_left hlog h8R
  have step2 : 8 / R * (p * Real.log x) ≤ 8 * Real.log x := by
    have hpR' : p / R ≤ 1 := (div_le_one hR).mpr hpR
    have : 8 / R * (p * Real.log x) = 8 * Real.log x * (p / R) := by ring
    rw [this]
    nlinarith
  have hexp : 8 / R * (p * Real.log x + Real.log K + Real.log D)
      = 8 / R * (p * Real.log x) + 8 / R * (Real.log K + Real.log D) := by ring
  linarith

theorem sharpMs_nonneg {σ t R R' : ℝ} {N : ℕ} (hN : 1 ≤ N) (hσ : 1 < σ) (hR : 0 < R)
    (hR1 : R < 1) (haR : σ - R < 1) (hRR' : R < R') (hR't : R' ≤ |t|) :
    0 ≤ GrowthSharp.sharpMs σ t R N := by
  have hz : (σ : ℂ) + I * (t : ℂ) + (R : ℂ) ∈ sphere ((σ : ℂ) + I * (t : ℂ)) R := by
    rw [mem_sphere, dist_eq]
    simp only [add_sub_cancel_left, Complex.norm_real, Real.norm_eq_abs, abs_of_pos hR]
  exact le_trans (norm_nonneg _) (GrowthSharp.sphere_le_sharpMs hN hσ hR hR1 haR hRR' hR't hz)

theorem hB_sharp {σ γ R : ℝ} (hσ : 1 < σ) (hσ2 : σ ≤ 2) (hR : 0 < R) (hR1 : R < 1)
    (haR : σ - R < 1) (hγ : 3 ≤ |γ|) (hMpos : 0 < GrowthSharp.sharpMs σ γ R ⌊|γ|⌋₊ + 1) :
    4 * Real.log ((GrowthSharp.sharpMs σ γ R ⌊|γ|⌋₊ + 1)
        / ‖riemannZeta ((σ : ℂ) + I * (γ : ℂ))‖) / (R / 2) ≤ sharpRHS σ γ R := by
  have hcre : ((σ : ℂ) + I * (γ : ℂ)).re = σ := by simp
  refine hBexplicit_of_ge (by rw [hcre]; exact hσ) hMpos hR ?_
  rw [hcre]
  exact scalar_le_sharpRHS hσ hσ2 hR hR1 haR hγ

theorem zeta_region_5841_allsharp :
    ∃ A₀ : ℝ, 0 ≤ A₀ ∧ ∀ {γ β ε B C₂ C₃ σ R R' : ℝ},
      0 < ε → 0 ≤ B → 0 ≤ C₂ → 0 ≤ C₃ →
      ε ≤ 5 * A₀ + 8 * B + 4 * C₂ + C₃ →
      σ = 1 + ε / (5 * A₀ + 8 * B + 4 * C₂ + C₃) →
      β < σ →
      0 < R → R < R' → R' ≤ |γ| → R < 1 → σ ≤ 2 → σ - R < 1 → 3 ≤ |γ| →
      ((β : ℂ) + I * (γ : ℝ)) ∈ ball ((σ : ℂ) + I * (γ : ℝ)) R' →
      riemannZeta ((β : ℂ) + I * (γ : ℝ)) = 0 →
      ((β : ℂ) + I * (γ : ℝ)) ∈ closedBall ((σ : ℂ) + I * (γ : ℝ)) (R / 2) →
      sharpRHS σ γ R ≤ B → sharpRHS σ (2 * γ) R ≤ C₂ → sharpRHS σ (3 * γ) R ≤ C₃ →
      β ≤ 1 - ε * (3 - ε) / ((5 + ε) * (5 * A₀ + 8 * B + 4 * C₂ + C₃)) := by
  obtain ⟨A₀, hA0, H⟩ := zeta_region_5841_sharp
  refine ⟨A₀, hA0, ?_⟩
  intro γ β ε B C₂ C₃ σ R R' hε hB hC2 hC3 hεK hσdef hβσ hR hRR' hR'γ hR1 hσ2 haR hγ3
    hρball hζρ hp hGB hG2 hG3
  have hKpos : 0 < 5 * A₀ + 8 * B + 4 * C₂ + C₃ := lt_of_lt_of_le hε hεK
  have hσ : 1 < σ := by
    rw [hσdef]; have : 0 < ε / (5 * A₀ + 8 * B + 4 * C₂ + C₃) := div_pos hε hKpos; linarith
  have hN : 1 ≤ ⌊|γ|⌋₊ := Nat.le_floor (by norm_num; linarith)
  set Ms : ℝ := GrowthSharp.sharpMs σ γ R ⌊|γ|⌋₊ with hMs
  have hMs0 : 0 ≤ Ms := sharpMs_nonneg hN hσ hR hR1 haR hRR' hR'γ
  have hMpos : 0 < Ms + 1 := by linarith
  have h2γ : 3 ≤ |2 * γ| := by rw [abs_mul]; norm_num; linarith
  exact H (Ms := Ms) (Mζ := Ms + 1) (Rg := R) hε hB hC2 hC3 hεK hσdef hβσ hR hRR'
    (one_notMem_ball_test hR'γ) hρball hζρ hp
    (fun z hz => GrowthSharp.sphere_le_sharpMs hN hσ hR hR1 haR hRR' hR'γ hz)
    (by linarith) hMpos
    (le_trans (hB_sharp hσ hσ2 hR hR1 haR hγ3 hMpos) hGB)
    hσ2 hR hR1 haR h2γ hG2 hG3

end TrigPoly5841AllSharp

#print axioms TrigPoly5841AllSharp.scalar_le_sharpRHS
#print axioms TrigPoly5841AllSharp.hB_sharp
#print axioms TrigPoly5841AllSharp.zeta_region_5841_allsharp
