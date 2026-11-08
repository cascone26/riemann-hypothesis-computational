/-
CauchyDetPos.lean (HP RH center, LEDGER S2914, 2026-11-08) — Cauchy determinant POSITIVITY.
x, y strictly increasing, all x_i + y_j > 0 ⟹ det [1/(x_i+y_j)] > 0. Every square submatrix of such a Cauchy matrix is again
of this form (increasing sub-sequences), so Cauchy matrices with these orderings are TOTALLY POSITIVE — the core of S2859's
TP-band theorem (parity blocks of τ are Cauchy-Gram in x = ν², poles γ²).
-/
import Mathlib
import CauchyDetStep

open Matrix

namespace CauchyDet

theorem det_cauchy_pos : ∀ {n : ℕ} (x y : Fin n → ℝ), StrictMono x → StrictMono y →
    (∀ i j, 0 < x i + y j) → 0 < (cauchy x y).det := by
  intro n
  induction n with
  | zero => intro x y _ _ _; simp [cauchy, Matrix.det_isEmpty]
  | succ n ih =>
    intro x y hx hy hpos
    rw [det_cauchy_succ x y (fun i j => (hpos i j).ne')]
    have hlast : ∀ i : Fin n, i.castSucc < Fin.last n := fun i => Fin.castSucc_lt_last i
    refine mul_pos (inv_pos.mpr (hpos _ _)) (mul_pos (mul_pos ?_ ?_) ?_)
    · exact Finset.prod_pos fun i _ => div_pos (sub_pos.mpr (hx (hlast i))) (hpos _ _)
    · exact ih (x ∘ Fin.castSucc) (y ∘ Fin.castSucc) (hx.comp Fin.strictMono_castSucc)
        (hy.comp Fin.strictMono_castSucc) (fun i j => hpos _ _)
    · exact Finset.prod_pos fun j _ => div_pos (sub_pos.mpr (hy (hlast j))) (hpos _ _)

end CauchyDet

#print axioms CauchyDet.det_cauchy_pos
