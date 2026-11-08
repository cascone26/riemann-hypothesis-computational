/-
SGramSchmidt.lean (HP RH center, LEDGER S2907, 2026-11-07)
Gram–Schmidt in the S-inner product (S ≻ 0), supplying S2884's family hypothesis from positivity alone:
for generators f : Fin r → (n → ℝ) (e.g. od(k) in frequency order), O := Gram–Schmidt of f w.r.t. ⟪x,y⟫_S satisfies
  O j ∈ span (f '' Iic j)    and    ∀ i < j, ∀ w ∈ span (f '' Iic i),  w ⬝ᵥ (S *ᵥ O j) = 0 .
-/
import Mathlib

open Matrix

namespace SGS

variable {n : Type*} [Fintype n] [DecidableEq n] {r : ℕ}

/-- Gram–Schmidt with respect to the inner product induced by a positive-definite S. -/
noncomputable def gsS (S : Matrix n n ℝ) (hS : S.PosDef) (f : Fin r → n → ℝ) : Fin r → n → ℝ :=
  letI := S.toNormedAddCommGroup hS
  letI := S.toInnerProductSpace hS.posSemidef
  InnerProductSpace.gramSchmidt ℝ f

theorem gsS_mem_span (S : Matrix n n ℝ) (hS : S.PosDef) (f : Fin r → n → ℝ) (j : Fin r) :
    gsS S hS f j ∈ Submodule.span ℝ (f '' Set.Iic j) := by
  letI := S.toNormedAddCommGroup hS
  letI := S.toInnerProductSpace hS.posSemidef
  exact InnerProductSpace.gramSchmidt_mem_span ℝ f le_rfl

theorem gsS_orth (S : Matrix n n ℝ) (hS : S.PosDef) (f : Fin r → n → ℝ) {i j : Fin r} (hij : i < j)
    (w : n → ℝ) (hw : w ∈ Submodule.span ℝ (f '' Set.Iic i)) :
    w ⬝ᵥ (S *ᵥ gsS S hS f j) = 0 := by
  letI := S.toNormedAddCommGroup hS
  letI := S.toInnerProductSpace hS.posSemidef
  -- span of f up to i = span of the GS vectors up to i
  rw [← InnerProductSpace.span_gramSchmidt_Iic ℝ f i] at hw
  have key : ∀ v ∈ Submodule.span ℝ (InnerProductSpace.gramSchmidt ℝ f '' Set.Iic i),
      (S *ᵥ InnerProductSpace.gramSchmidt ℝ f j) ⬝ᵥ v = 0 := by
    intro v hv
    induction hv using Submodule.span_induction with
    | mem x hx =>
      obtain ⟨k, hk, rfl⟩ := hx
      have ho := InnerProductSpace.gramSchmidt_orthogonal ℝ f (ne_of_lt (lt_of_le_of_lt hk hij))
      have e : @inner ℝ (n → ℝ) _ (InnerProductSpace.gramSchmidt ℝ f k) (InnerProductSpace.gramSchmidt ℝ f j)
          = (S *ᵥ InnerProductSpace.gramSchmidt ℝ f j) ⬝ᵥ star (InnerProductSpace.gramSchmidt ℝ f k) := rfl
      rw [e] at ho
      simpa using ho
    | zero => exact dotProduct_zero _
    | add x y _ _ hx hy => rw [dotProduct_add, hx, hy, add_zero]
    | smul a x _ hx => rw [dotProduct_smul, hx, smul_zero]
  have h := key w hw
  simpa [gsS, dotProduct_comm] using h

end SGS

#print axioms SGS.gsS_mem_span
#print axioms SGS.gsS_orth
