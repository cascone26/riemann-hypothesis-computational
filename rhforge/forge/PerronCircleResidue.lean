import PerronResidueDeposit
import ZeroFreeBall

set_option maxHeartbeats 1000000

/-!
# §1618 — the Perron integrand's circle residue is unconditional: `∮_{C(1,R)} F = 2πi·(x²/2)`

Combining the residue deposit §1616 with the zero-free disc §1617 (which discharges its hypotheses) and rewriting
`Ĥ(z)/(z-1) = F(z)` on the sphere (where `z ≠ 1`, so `(z-1)` cancels), the Perron integrand `F` itself has circle
residue `2πi·(x²/2)` at `s=1` — for a concrete `R > 0`.  This is the PNT main term realized as an unconditional
contour integral of the actual integrand.  EXIT0 / std-axioms.
-/

open Complex

theorem perron_circle_residue {x : ℝ} (hx : 0 < x) :
    ∃ R : ℝ, 0 < R ∧
      (∮ z in C((1 : ℂ), R), -logDeriv riemannZeta z * (x : ℂ) ^ (z + 1) / (z * (z + 1)))
        = 2 * ↑Real.pi * I * ((x : ℂ) ^ (2 : ℂ) / 2) := by
  obtain ⟨R, hR, hζ, h0, hm1⟩ := exists_zero_free_ball
  refine ⟨R, hR, ?_⟩
  rw [← perron_residue_deposit hx hR hζ h0 hm1]
  apply circleIntegral.integral_congr hR.le
  intro z hz
  rw [Metric.mem_sphere, Complex.dist_eq] at hz
  have hz1 : z ≠ 1 := by
    intro h; rw [h, sub_self, norm_zero] at hz; linarith
  have hsub : z - 1 ≠ 0 := sub_ne_zero.mpr hz1
  simp only [Function.update_of_ne hz1]
  generalize -logDeriv riemannZeta z * (x : ℂ) ^ (z + 1) / (z * (z + 1)) = A
  field_simp

#print axioms perron_circle_residue
