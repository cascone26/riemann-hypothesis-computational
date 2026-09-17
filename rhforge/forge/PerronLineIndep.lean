import PerronRectZero
import PerronHorizontalEdge
import PerronHorizEdgeVanish
import PerronHorizEdgeBottom
import PerronFVlineTruncated

set_option maxHeartbeats 1000000

open Complex Filter Topology

/-- §1648 — CONTOUR INDEPENDENCE: the full Perron line integral `∫_ℝ F(σ+it) dt` is the same on every
vertical line `σ > 1` (`F s = -ζ'/ζ(s)·x^{s+1}/(s(s+1))`).  Cauchy–Goursat §1611 on `[σ₁,σ₀]×[-T,T]` gives a
4-term identity `≡ 0`; as `T→∞` the horizontal edges vanish (§1637/§1646) and the vertical sides converge to the
full lines (§1644), forcing the two line integrals to agree.  EXIT0 / std-axioms. -/
theorem perron_line_integral_indep {x σ₁ σ₀ : ℝ} (hx : 0 < x) (hσ₁ : 1 < σ₁) (hσ₀ : σ₁ ≤ σ₀) :
    (∫ t : ℝ, -logDeriv riemannZeta ((σ₀ : ℂ) + (t : ℝ) * I) *
        (x : ℂ) ^ (((σ₀ : ℂ) + (t : ℝ) * I) + 1) /
      (((σ₀ : ℂ) + (t : ℝ) * I) * (((σ₀ : ℂ) + (t : ℝ) * I) + 1)))
      = ∫ t : ℝ, -logDeriv riemannZeta ((σ₁ : ℂ) + (t : ℝ) * I) *
          (x : ℂ) ^ (((σ₁ : ℂ) + (t : ℝ) * I) + 1) /
        (((σ₁ : ℂ) + (t : ℝ) * I) * (((σ₁ : ℂ) + (t : ℝ) * I) + 1)) := by
  have hσ₀gt : 1 < σ₀ := lt_of_lt_of_le hσ₁ hσ₀
  set A : ℝ → ℂ := fun T => ∫ σ' : ℝ in σ₁..σ₀,
    -logDeriv riemannZeta ((σ' : ℂ) + (-T : ℝ) * I) *
        (x : ℂ) ^ (((σ' : ℂ) + (-T : ℝ) * I) + 1) /
      (((σ' : ℂ) + (-T : ℝ) * I) * (((σ' : ℂ) + (-T : ℝ) * I) + 1)) with hAdef
  set B : ℝ → ℂ := fun T => ∫ σ' : ℝ in σ₁..σ₀,
    -logDeriv riemannZeta ((σ' : ℂ) + (T : ℝ) * I) *
        (x : ℂ) ^ (((σ' : ℂ) + (T : ℝ) * I) + 1) /
      (((σ' : ℂ) + (T : ℝ) * I) * (((σ' : ℂ) + (T : ℝ) * I) + 1)) with hBdef
  set C : ℝ → ℂ := fun T => I • ∫ y : ℝ in (-T)..T,
    -logDeriv riemannZeta ((σ₀ : ℂ) + (y : ℝ) * I) *
        (x : ℂ) ^ (((σ₀ : ℂ) + (y : ℝ) * I) + 1) /
      (((σ₀ : ℂ) + (y : ℝ) * I) * (((σ₀ : ℂ) + (y : ℝ) * I) + 1)) with hCdef
  set D : ℝ → ℂ := fun T => I • ∫ y : ℝ in (-T)..T,
    -logDeriv riemannZeta ((σ₁ : ℂ) + (y : ℝ) * I) *
        (x : ℂ) ^ (((σ₁ : ℂ) + (y : ℝ) * I) + 1) /
      (((σ₁ : ℂ) + (y : ℝ) * I) * (((σ₁ : ℂ) + (y : ℝ) * I) + 1)) with hDdef
  have hid : ∀ T : ℝ, A T - B T + C T - D T = 0 := fun T =>
    perron_rect_boundary_zero (x := x) hx (z := ⟨σ₁, -T⟩) (w := ⟨σ₀, T⟩) hσ₁ hσ₀gt
  have hA : Tendsto A atTop (𝓝 0) := perron_horizontal_edge_bottom_tendsto_zero hx hσ₁ hσ₀
  have hB : Tendsto B atTop (𝓝 0) := perron_horizontal_edge_tendsto_zero hx hσ₁ hσ₀
  have hC : Tendsto C atTop (𝓝 (I • ∫ y : ℝ,
      -logDeriv riemannZeta ((σ₀ : ℂ) + (y : ℝ) * I) *
          (x : ℂ) ^ (((σ₀ : ℂ) + (y : ℝ) * I) + 1) /
        (((σ₀ : ℂ) + (y : ℝ) * I) * (((σ₀ : ℂ) + (y : ℝ) * I) + 1)))) :=
    (perron_F_vline_truncated_tendsto hx hσ₀gt).const_smul I
  have hD : Tendsto D atTop (𝓝 (I • ∫ y : ℝ,
      -logDeriv riemannZeta ((σ₁ : ℂ) + (y : ℝ) * I) *
          (x : ℂ) ^ (((σ₁ : ℂ) + (y : ℝ) * I) + 1) /
        (((σ₁ : ℂ) + (y : ℝ) * I) * (((σ₁ : ℂ) + (y : ℝ) * I) + 1)))) :=
    (perron_F_vline_truncated_tendsto hx hσ₁).const_smul I
  have hlim := ((hA.sub hB).add hC).sub hD
  have hlim0 : Tendsto (fun T => A T - B T + C T - D T) atTop (𝓝 0) := by
    simp only [hid]; exact tendsto_const_nhds
  have heq := tendsto_nhds_unique hlim0 hlim
  have h0 : I • (∫ y : ℝ,
        -logDeriv riemannZeta ((σ₀ : ℂ) + (y : ℝ) * I) *
            (x : ℂ) ^ (((σ₀ : ℂ) + (y : ℝ) * I) + 1) /
          (((σ₀ : ℂ) + (y : ℝ) * I) * (((σ₀ : ℂ) + (y : ℝ) * I) + 1)))
      - I • (∫ y : ℝ,
        -logDeriv riemannZeta ((σ₁ : ℂ) + (y : ℝ) * I) *
            (x : ℂ) ^ (((σ₁ : ℂ) + (y : ℝ) * I) + 1) /
          (((σ₁ : ℂ) + (y : ℝ) * I) * (((σ₁ : ℂ) + (y : ℝ) * I) + 1))) = 0 := by
    simpa using heq.symm
  rw [← smul_sub] at h0
  rcases smul_eq_zero.mp h0 with hI | hval
  · exact absurd hI Complex.I_ne_zero
  · exact sub_eq_zero.mp hval

#print axioms perron_line_integral_indep
