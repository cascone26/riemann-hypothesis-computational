import VlineResidueLimit

set_option maxHeartbeats 1000000

open Complex intervalIntegral Filter Topology

/-!
# §1662 ★ — THE PERRON POLE-TERM SHIFT: crossing `s = 1` deposits the PNT main term `2π·x²/2`

Specialization of §1661 to the Perron pole term `(x²/2)/(s-1)` (§1614/§1623/§1626 spelling),
`p = 1`: as `T → ∞`, moving the truncated Perron pole integral from any line `σ₀ > 1` to any line
`σ₁ < 1` (into the zero-free region) deposits exactly `2π·(x²/2)` — the residue of the smoothed
main term. This is the analytic mechanism turning the §1624 line split `∫F = ∫pole + ∫R` into
`Ψ(x) = x²/2 + (shifted remainder)`: the pole side is now fully machine-checked.
-/

theorem perron_pole_shift_tendsto {x σ₁ σ₀ : ℝ} (h1 : σ₁ < 1) (h2 : 1 < σ₀) :
    Tendsto (fun T : ℝ =>
        (∫ y in (-T)..T, (x : ℂ) ^ (2 : ℂ) / 2 / (((σ₀ : ℂ) + (y : ℝ) * Complex.I) - 1)) -
          (∫ y in (-T)..T, (x : ℂ) ^ (2 : ℂ) / 2 / (((σ₁ : ℂ) + (y : ℝ) * Complex.I) - 1)))
      atTop (𝓝 (2 * (Real.pi : ℂ) * ((x : ℂ) ^ (2 : ℂ) / 2))) := by
  have h := vline_diff_residue_tendsto ((x : ℂ) ^ (2 : ℂ) / 2) 1 (σ₁ := σ₁) (σ₀ := σ₀)
    (by simpa using h1) (by simpa using h2)
  simpa only [div_eq_mul_inv] using h

#print axioms perron_pole_shift_tendsto
