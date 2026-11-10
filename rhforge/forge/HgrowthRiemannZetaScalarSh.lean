import ReNegLogDerivLeBOfSphereSh
import HBexplicitOfGeSh

set_option maxHeartbeats 1000000

/-! [S3020 SHARP COPY: Borel–Carathéodory constant 4 → 2 via SharpBC (S3018); names suffixed _sh]

# §1555 — scalar ζ growth `Re(-logDeriv ζ c) ≤ B` from sphere data + `hball`, BC discharged

ζ-analogue of §1540 (the entire-`LFunction` growth), but pole-aware: uses §1429 `re_neg_logDeriv_le_B_of_sphere_sh`
(which needs `hball : 1 ∉ ball c R'`, the pole avoided) and §1444 `hBexplicit_of_ge_sh` (BC scalar from the ζ lower
bound §1443).  Setting `Mζ = Ms+1`, the growth value bound reduces to a pure scalar `hByscalar`.  Applied at
`c = σK+2iγ` (with `hball` from `R' ≤ |2γ|` via §1371) this supplies the ζ part of the QUADRATIC (`χ²=1`) growth
term (the other part = the finite-product correction §1554).  `Ms ≥ 0` is extracted from a sphere point.
EXIT0 / std-axioms.
-/

open Complex Metric

theorem hgrowth_riemannZeta_scalar_sh {c : ℂ} {R R' Ms B : ℝ} (hR : 0 < R) (hRR' : R < R')
    (hcre : 1 < c.re) (hball : (1 : ℂ) ∉ ball c R')
    (hζsphere : ∀ z ∈ sphere c R, ‖riemannZeta z‖ ≤ Ms)
    (hByscalar : 2 * Real.log ((Ms + 1) * (3 + 1 / (c.re - 1))) / (R / 2) ≤ B) :
    (-logDeriv riemannZeta c).re ≤ B := by
  have hz0 : c + (R : ℂ) ∈ sphere c R := by
    rw [mem_sphere, dist_eq]
    simp only [add_sub_cancel_left, Complex.norm_real, Real.norm_eq_abs, abs_of_pos hR]
  have hMs : 0 ≤ Ms := le_trans (norm_nonneg _) (hζsphere _ hz0)
  have hMspos : 0 < Ms + 1 := by linarith
  refine re_neg_logDeriv_le_B_of_sphere_sh hR hRR' hcre hball hζsphere (by linarith) hMspos ?_
  exact hBexplicit_of_ge_sh hcre hMspos hR hByscalar

#print axioms hgrowth_riemannZeta_scalar_sh
