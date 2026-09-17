import PerronFVlineIntegrable
import Mathlib.MeasureTheory.Integral.IntegralEqImproper

set_option maxHeartbeats 1000000

open Complex MeasureTheory Filter Topology

/-- §1644 — the truncated vertical integral converges to the full line integral: for `σ > 1`,
`∫_{-T}^{T} F(σ+it) dt → ∫_ℝ F(σ+it) dt` as `T → ∞`.  Since `F` is integrable on the whole line (§1625), the
symmetric truncation exhausts it — the truncated-rectangle vertical sides approach the untruncated Perron line. -/
theorem perron_F_vline_truncated_tendsto {x : ℝ} (hx : 0 < x) {σ : ℝ} (hσ : 1 < σ) :
    Tendsto (fun T : ℝ => ∫ t : ℝ in (-T)..T,
        -logDeriv riemannZeta ((σ : ℂ) + (t : ℝ) * I) *
            (x : ℂ) ^ (((σ : ℂ) + (t : ℝ) * I) + 1) /
          (((σ : ℂ) + (t : ℝ) * I) * (((σ : ℂ) + (t : ℝ) * I) + 1)))
      atTop (𝓝 (∫ t : ℝ,
        -logDeriv riemannZeta ((σ : ℂ) + (t : ℝ) * I) *
            (x : ℂ) ^ (((σ : ℂ) + (t : ℝ) * I) + 1) /
          (((σ : ℂ) + (t : ℝ) * I) * (((σ : ℂ) + (t : ℝ) * I) + 1)))) :=
  intervalIntegral_tendsto_integral (integrable_perron_F_vline hx hσ)
    tendsto_neg_atTop_atBot tendsto_id

#print axioms perron_F_vline_truncated_tendsto
