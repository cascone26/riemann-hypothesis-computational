import HdecompHgBRelaxedOfSphereSh
import HzeroZpow
import DivisorNonnegAll
import OneLeDivisorZeta
import DivisorFilterReLt
import ZetaZeroMemFilter
import TrigPolyRegionPole

/-! [S3020 SHARP COPY: Borel–Carathéodory constant 4 → 2 via SharpBC (S3018); names suffixed _sh]

# S2934 — point-generic ζ ZERO input at σ+iγ (lane 1)

The 3-4-1 tower derived the zero input `Re(−ζ'/ζ)(σK+iγ) ≤ B − 1/(σK−β)` only inside the monolithic
§1425/§1409 chain at the fixed abscissa σK = 1+1/(2K).  Here it is extracted at ANY σ > 1:
given a zero ρ₀ = β+iγ with β < σ in closedBall(c, R/2), c = σ+iγ, the BC sphere data
(‖ζ‖ ≤ Ms < Mζ on the R-sphere, 4·log(Mζ/‖ζ c‖)/(R/2) ≤ B) and 1 ∉ ball(c, R'), we get
  Re(−ζ'/ζ)(σ+iγ) ≤ B − 1/(σ−β).
`hzero_Tz_sh` restates it as the `hzero` hypothesis of S2926/S2932 (`Tz σ γ 1`), so the P-generic region
(5-8-4-1 included) now has all three analytic inputs available as compiled lemmas (pole S2932, zero here,
growth S2928 + S2930 log absorption).
-/

open Complex Metric

theorem hzero_at_point_sh {σ β γ R R' Mζ Ms B : ℝ} (hR : 0 < R) (hRR' : R < R') (hσ : 1 < σ)
    (hβσ : β < σ)
    (hball : (1 : ℂ) ∉ ball ((σ : ℂ) + I * (γ : ℝ)) R')
    (hρ₀ball : ((β : ℂ) + I * (γ : ℝ)) ∈ ball ((σ : ℂ) + I * (γ : ℝ)) R')
    (hζρ₀ : riemannZeta ((β : ℂ) + I * (γ : ℝ)) = 0)
    (hp : ((β : ℂ) + I * (γ : ℝ)) ∈ closedBall ((σ : ℂ) + I * (γ : ℝ)) (R / 2))
    (hζsphere : ∀ z ∈ sphere ((σ : ℂ) + I * (γ : ℝ)) R, ‖riemannZeta z‖ ≤ Ms)
    (hsl : Ms < Mζ) (hMζ : 0 < Mζ)
    (hBexplicit : 2 * Real.log (Mζ / ‖riemannZeta ((σ : ℂ) + I * (γ : ℝ))‖) / (R / 2) ≤ B) :
    (-logDeriv riemannZeta ((σ : ℂ) + I * (γ : ℝ))).re ≤ B - 1 / (σ - β) := by
  have hcre : 1 < ((σ : ℂ) + I * (γ : ℝ)).re := by simpa using hσ
  obtain ⟨g, hdec, hgB⟩ :=
    hdecomp_and_hgB_relaxed_of_sphere_sh hR hRR' hcre hball hζsphere hsl hMζ hBexplicit
  have hz₀ : ((σ : ℂ) + I * (γ : ℝ)) ∈ ball ((σ : ℂ) + I * (γ : ℝ)) R' :=
    mem_ball_self (lt_trans hR hRR')
  have hz₀ne := riemannZeta_ne_zero_of_one_lt_re hcre
  exact re_neg_le_of_zpow_decomp
    (zeta_zero_mem_filter hball hρ₀ball hζρ₀ hz₀ hz₀ne _ hp)
    (fun ρ hρ => divisor_filter_support_re_lt hσ _ hρ)
    (fun ρ _ => divisor_riemannZeta_nonneg_all hball ρ)
    (one_le_divisor_riemannZeta_at_zero hball hρ₀ball hζρ₀ hz₀ hz₀ne)
    hβσ hgB hdec

namespace TrigPolyRegion

theorem Tz_one_sh (σ γ : ℝ) : Tz σ γ 1 = (-logDeriv riemannZeta ((σ : ℂ) + I * (γ : ℝ))).re := by
  unfold Tz; simp

theorem hzero_Tz_sh {σ β γ R R' Mζ Ms B : ℝ} (hR : 0 < R) (hRR' : R < R') (hσ : 1 < σ)
    (hβσ : β < σ)
    (hball : (1 : ℂ) ∉ ball ((σ : ℂ) + I * (γ : ℝ)) R')
    (hρ₀ball : ((β : ℂ) + I * (γ : ℝ)) ∈ ball ((σ : ℂ) + I * (γ : ℝ)) R')
    (hζρ₀ : riemannZeta ((β : ℂ) + I * (γ : ℝ)) = 0)
    (hp : ((β : ℂ) + I * (γ : ℝ)) ∈ closedBall ((σ : ℂ) + I * (γ : ℝ)) (R / 2))
    (hζsphere : ∀ z ∈ sphere ((σ : ℂ) + I * (γ : ℝ)) R, ‖riemannZeta z‖ ≤ Ms)
    (hsl : Ms < Mζ) (hMζ : 0 < Mζ)
    (hBexplicit : 2 * Real.log (Mζ / ‖riemannZeta ((σ : ℂ) + I * (γ : ℝ))‖) / (R / 2) ≤ B) :
    Tz σ γ 1 ≤ B - 1 / (σ - β) := by
  rw [Tz_one_sh]
  exact hzero_at_point_sh hR hRR' hσ hβσ hball hρ₀ball hζρ₀ hp hζsphere hsl hMζ hBexplicit

end TrigPolyRegion

#print axioms hzero_at_point_sh
#print axioms TrigPolyRegion.hzero_Tz_sh
