import HgrowthLFunctionFromSphere
import LFunctionSphereN1
import HBexplicitOfGeLFunction

set_option maxHeartbeats 1000000

/-!
# §1540 — `Re(-logDeriv (LFunction ψ) c) ≤ B` from a PURE SCALAR bound on `B`

Combines §1510 (growth from sphere data), §1535 (scalar N=1 sphere bound `Ms = 1 + N·2^{-(c.re-R)} +
(‖c‖+R)·N/(c.re-R)`) and §1539 (BC scalar via the L lower bound) at a nontrivial character `ψ`.  Setting
`Mζ = Ms+1`, `R' = 2R`, every complex-analytic hypothesis of §1510 is discharged, leaving a single scalar
inequality `hByscalar` (no norms, no `L ψ c`).  Applied at `ψ = χ²`, `c = σK+2iγ` this discharges §1509's
`hgrowth`.  EXIT0 / std-axioms.
-/

open Complex DirichletCharacter

theorem hgrowth_LFunction_scalar {N : ℕ} [NeZero N] {ψ : DirichletCharacter ℂ N} (hψ : ψ ≠ 1)
    {c : ℂ} {R B : ℝ} (hR : 0 < R) (hcre : 1 < c.re) (hcR : 0 < c.re - R)
    (hByscalar : 4 * Real.log ((1 + ((N : ℝ) * (2 : ℝ) ^ (-(c.re - R))
          + (‖c‖ + R) * (N : ℝ) / (c.re - R)) + 1)
        * (3 + 1 / (c.re - 1))) / (R / 2) ≤ B) :
    (-logDeriv (LFunction ψ) c).re ≤ B := by
  set Ms : ℝ := 1 + ((N : ℝ) * (2 : ℝ) ^ (-(c.re - R)) + (‖c‖ + R) * (N : ℝ) / (c.re - R)) with hMs
  have hMspos : 0 < Ms + 1 := by
    have ha : (0 : ℝ) ≤ ‖c‖ + R := add_nonneg (norm_nonneg c) hR.le
    have hb : (0 : ℝ) ≤ (N : ℝ) := Nat.cast_nonneg N
    have hnn : 0 ≤ (‖c‖ + R) * (N : ℝ) / (c.re - R) := div_nonneg (mul_nonneg ha hb) hcR.le
    have hnn2 : 0 ≤ (N : ℝ) * (2 : ℝ) ^ (-(c.re - R)) := mul_nonneg hb (Real.rpow_nonneg (by norm_num) _)
    rw [hMs]; linarith
  refine hgrowth_LFunction_from_sphere (Mζ := Ms + 1) (Ms := Ms) hψ hR (by linarith : R < 2 * R)
    hcre ?_ (by linarith) hMspos ?_
  · intro z hz
    exact norm_LFunction_le_sphere_N1 hψ hR hcR hz
  · exact hBexplicit_of_ge_LFunction hψ hcre hMspos hR hByscalar

#print axioms hgrowth_LFunction_scalar
