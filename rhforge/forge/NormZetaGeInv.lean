import Mathlib.NumberTheory.LSeries.Dirichlet

/-!
# §1420 — ζ lower bound from a Möbius upper bound: `‖L↗μ s‖ ≤ M ⟹ M⁻¹ ≤ ‖ζ s‖`

`ζ(s)·(1/ζ(s)) = 1` with `1/ζ(s) = L↗μ s` (Mathlib `LSeries_zeta_mul_Lseries_moebius` +
`LSeries_zeta_eq_riemannZeta`), so `‖ζ s‖·‖L↗μ s‖ = 1`.  Any upper bound `M` on the Möbius series norm therefore
gives the lower bound `‖ζ s‖ ≥ M⁻¹`.  This is the reusable core of the classical `|ζ(σ+it)| ≥ 1/ζ(σ)` estimate
(next brick supplies `M = ζ(σ)`), feeding §1419's `‖ζ c‖`.  EXIT0 / std-axioms.
-/

open ArithmeticFunction
open scoped LSeries.notation Moebius

theorem norm_riemannZeta_ge_inv {s : ℂ} {M : ℝ} (hs : 1 < s.re) (hM : ‖L ↗μ s‖ ≤ M) :
    M⁻¹ ≤ ‖riemannZeta s‖ := by
  have hprod : riemannZeta s * L ↗μ s = 1 := by
    rw [← LSeries_zeta_eq_riemannZeta hs]; exact LSeries_zeta_mul_Lseries_moebius hs
  have hnorm : ‖riemannZeta s‖ * ‖L ↗μ s‖ = 1 := by
    rw [← norm_mul, hprod, norm_one]
  have hmu_pos : 0 < ‖L ↗μ s‖ := by
    rcases (norm_nonneg (L ↗μ s)).lt_or_eq with h | h
    · exact h
    · rw [← h, mul_zero] at hnorm; exact absurd hnorm (by norm_num)
  have hMpos : 0 < M := lt_of_lt_of_le hmu_pos hM
  have key : 1 ≤ M * ‖riemannZeta s‖ := by
    nlinarith [hnorm, mul_nonneg (norm_nonneg (riemannZeta s)) (by linarith : (0:ℝ) ≤ M - ‖L ↗μ s‖)]
  calc M⁻¹ = 1 * M⁻¹ := (one_mul _).symm
    _ ≤ (M * ‖riemannZeta s‖) * M⁻¹ :=
        mul_le_mul_of_nonneg_right key (le_of_lt (inv_pos.mpr hMpos))
    _ = ‖riemannZeta s‖ := by field_simp

#print axioms norm_riemannZeta_ge_inv
