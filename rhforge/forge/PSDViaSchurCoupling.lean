import Mathlib

open Matrix

/-- **PATH D coupling mechanism: the defect is paired away, not eliminated.**

Certifying the Weil form's positivity in the strong-block reduction (RH ⟺ Q_SS ⪰ 0): the strong block
`A` (where the archimedean and prime forms coincide, `A ≻ 0` with a razor-thin margin) couples to the
residual block `D` through `B`. Even when the residual block `D` is NOT positive semidefinite on its
own — the negative-symbol defect — the FULL block form is PSD provided the coupling-corrected **Schur
complement** `D - Bᴴ A⁻¹ B` is PSD. The coupling `Bᴴ A⁻¹ B` absorbs (pairs away) the defect. This is
the algebraic heart of the numerical corner certification (PATH D, [14,44] certified Q_SS>0): what a
naive concentration/uncertainty bound cannot see (§888, margin stays negative), the Schur coupling
closes. -/
theorem psd_via_schur_coupling {m n : Type*} [Fintype m] [Fintype n] [DecidableEq m]
    (A : Matrix m m ℝ) (B : Matrix m n ℝ) (D : Matrix n n ℝ)
    (hA : A.PosDef) (hSchur : (D - Bᴴ * A⁻¹ * B).PosSemidef) :
    (fromBlocks A B Bᴴ D).PosSemidef := by
  haveI : Invertible A := hA.isUnit.invertible
  exact (Matrix.PosDef.fromBlocks₁₁ B D hA).mpr hSchur
