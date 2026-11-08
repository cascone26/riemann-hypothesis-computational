/-
CauchyTP.lean (HP RH center, LEDGER S2915, 2026-11-08) — TOTAL POSITIVITY of Cauchy matrices.
x : Fin n → ℝ, y : Fin m → ℝ strictly increasing, all x_i + y_j > 0. Then every square minor with strictly increasing row
selection r and column selection c is POSITIVE:  det [1/(x_{r a} + y_{c b})]_{a,b} > 0.
Also the "pole" form used by S2859: for nodes u strictly increasing, poles v strictly DECREASING... handled by the substitution
x = u, y = −v ∘ rev if needed; here the + form is the canonical statement.
-/
import Mathlib
import CauchyDetPos

open Matrix

namespace CauchyTP

noncomputable def cauchyRect {n m : ℕ} (x : Fin n → ℝ) (y : Fin m → ℝ) : Matrix (Fin n) (Fin m) ℝ :=
  Matrix.of fun i j => (x i + y j)⁻¹

theorem minor_pos {n m k : ℕ} (x : Fin n → ℝ) (y : Fin m → ℝ) (hx : StrictMono x) (hy : StrictMono y)
    (hpos : ∀ i j, 0 < x i + y j) (r : Fin k → Fin n) (c : Fin k → Fin m)
    (hr : StrictMono r) (hc : StrictMono c) :
    0 < ((cauchyRect x y).submatrix r c).det := by
  have : (cauchyRect x y).submatrix r c = CauchyDet.cauchy (x ∘ r) (y ∘ c) := by
    ext a b; simp [cauchyRect, CauchyDet.cauchy]
  rw [this]
  exact CauchyDet.det_cauchy_pos (x ∘ r) (y ∘ c) (hx.comp hr) (hy.comp hc) (fun a b => hpos _ _)

/-- All entries positive (1×1 minors) — the k = 1 case, stated directly. -/
theorem entry_pos {n m : ℕ} (x : Fin n → ℝ) (y : Fin m → ℝ) (hpos : ∀ i j, 0 < x i + y j) (i : Fin n) (j : Fin m) :
    0 < cauchyRect x y i j := by
  simp [cauchyRect]; exact hpos i j

end CauchyTP

#print axioms CauchyTP.minor_pos
