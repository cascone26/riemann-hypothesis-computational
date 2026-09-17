import Mathlib.NumberTheory.LSeries.Dirichlet
import Mathlib.NumberTheory.LSeries.DirichletContinuation

/-!
# §1536 — L lower bound from a twisted-Möbius upper bound: `‖L(↗χ·↗μ) s‖ ≤ M ⟹ M⁻¹ ≤ ‖LFunction χ s‖`

L-analogue of ζ §1420 (`norm_riemannZeta_ge_inv`).  Mathlib's `LSeries.mul_mu_eq_one` gives
`L ↗χ s · L (↗χ·↗μ) s = 1` for `1 < Re s`, and `LFunction_eq_LSeries` identifies `LFunction χ s = L ↗χ s` there,
so `‖LFunction χ s‖ · ‖L (↗χ·↗μ) s‖ = 1`.  Any upper bound `M` on the twisted-Möbius series norm therefore yields
the lower bound `‖LFunction χ s‖ ≥ M⁻¹` — the reusable core of the classical `|L(σ+it,χ)| ≥ 1/ζ(σ)` estimate (next
brick supplies `M = ζ(σ)` via `|χ(n)μ(n)| ≤ 1`).  This is the last input the BC scalar `hBexplicit` needs.
EXIT0 / std-axioms.
-/

open ArithmeticFunction DirichletCharacter
open scoped LSeries.notation ArithmeticFunction.Moebius

theorem norm_LFunction_ge_inv {N : ℕ} [NeZero N] (χ : DirichletCharacter ℂ N) {s : ℂ} {M : ℝ}
    (hs : 1 < s.re) (hM : ‖L (↗χ * ↗μ) s‖ ≤ M) :
    M⁻¹ ≤ ‖LFunction χ s‖ := by
  have hprod : LFunction χ s * L (↗χ * ↗μ) s = 1 := by
    rw [LFunction_eq_LSeries χ hs]; exact LSeries.mul_mu_eq_one χ hs
  have hnorm : ‖LFunction χ s‖ * ‖L (↗χ * ↗μ) s‖ = 1 := by
    rw [← norm_mul, hprod, norm_one]
  have hmu_pos : 0 < ‖L (↗χ * ↗μ) s‖ := by
    rcases (norm_nonneg (L (↗χ * ↗μ) s)).lt_or_eq with h | h
    · exact h
    · rw [← h, mul_zero] at hnorm; exact absurd hnorm (by norm_num)
  have hMpos : 0 < M := lt_of_lt_of_le hmu_pos hM
  have key : 1 ≤ M * ‖LFunction χ s‖ := by
    nlinarith [hnorm, mul_nonneg (norm_nonneg (LFunction χ s))
      (by linarith : (0 : ℝ) ≤ M - ‖L (↗χ * ↗μ) s‖)]
  calc M⁻¹ = 1 * M⁻¹ := (one_mul _).symm
    _ ≤ (M * ‖LFunction χ s‖) * M⁻¹ :=
        mul_le_mul_of_nonneg_right key (le_of_lt (inv_pos.mpr hMpos))
    _ = ‖LFunction χ s‖ := by field_simp

#print axioms norm_LFunction_ge_inv
