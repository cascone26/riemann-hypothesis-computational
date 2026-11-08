/-
CauchyInterlacing.lean (HP RH center, LEDGER S2899, 2026-11-07)
Cauchy interlacing for a codimension-1 compression (absent from Mathlib), from Courant–Fischer halves
(RayleighEigenTail S2897, CourantFischerMinMax S2898). The last missing ingredient for a full Lean S2884.
V ⊆ E, finrank V = m, finrank E = m+1, S on V with ⟪S v, w⟫ = ⟪T v, w⟫ (compression of T). Then
  λ_{k+1}(T) ≤ μ_k(S) ≤ λ_k(T)   (eigenvalues sorted decreasingly).
-/
import Mathlib
import RayleighEigenTail
import CourantFischerMinMax

open Module

namespace CauchyInterlace

variable {E : Type*} [NormedAddCommGroup E] [InnerProductSpace ℝ E] [FiniteDimensional ℝ E]
variable {m : ℕ} {T : E →ₗ[ℝ] E} {V : Submodule ℝ E} {S : V →ₗ[ℝ] V}

/-- the span of eigenvectors indexed by a finset is orthogonal to the others -/
lemma inner_eq_zero_of_mem_span {F : Type*} [NormedAddCommGroup F] [InnerProductSpace ℝ F]
    {ι : Type*} [Fintype ι] [DecidableEq ι] (c : OrthonormalBasis ι ℝ F) (A : Finset ι) (i : ι)
    (hi : i ∉ A) (v : F) (hv : v ∈ Submodule.span ℝ (Set.range fun j : A => c j)) :
    @inner ℝ F _ (c i) v = 0 := by
  have hle : Submodule.span ℝ (Set.range fun j : A => c j) ≤ LinearMap.ker (innerₛₗ ℝ (c i)) := by
    rw [Submodule.span_le]
    rintro _ ⟨j, rfl⟩
    have hij : i ≠ (j : ι) := fun h => hi (h ▸ j.2)
    simp [c.orthonormal.inner_eq_zero hij]
  simpa using hle hv

lemma finrank_span_sub {F : Type*} [NormedAddCommGroup F] [InnerProductSpace ℝ F]
    {ι : Type*} [Fintype ι] (c : OrthonormalBasis ι ℝ F) (A : Finset ι) :
    finrank ℝ (Submodule.span ℝ (Set.range fun j : A => c j)) = A.card := by
  rw [finrank_span_eq_card]
  · simp
  · exact c.orthonormal.linearIndependent.comp _ Subtype.val_injective

lemma finrank_map_subtype (W : Submodule ℝ V) : finrank ℝ (W.map V.subtype) = finrank ℝ W :=
  (Submodule.equivMapOfInjective V.subtype Subtype.val_injective W).finrank_eq.symm

theorem interlace_upper (hT : T.IsSymmetric) (hS : S.IsSymmetric)
    (hST : ∀ v w : V, @inner ℝ V _ (S v) w = @inner ℝ E _ (T v) w)
    (hE : finrank ℝ E = m + 1) (hV : finrank ℝ V = m) (k : Fin m) :
    hS.eigenvalues hV k ≤ hT.eigenvalues hE k.castSucc := by
  classical
  set c := hS.eigenvectorBasis hV
  set H : Submodule ℝ V := Submodule.span ℝ (Set.range fun j : (Finset.Iic k) => c j)
  have hHdim : finrank ℝ H = k + 1 := by rw [finrank_span_sub, Fin.card_Iic]
  obtain ⟨x, hxW, hx0, hx⟩ := CourantFischer.exists_le_eigenvalue hT hE k.castSucc (H.map V.subtype)
    (by rw [finrank_map_subtype, hHdim]; simp)
  obtain ⟨v, hvH, rfl⟩ := Submodule.mem_map.mp hxW
  have hhead := RayleighTail.rayleigh_head_ge hS hV k v (fun i hi =>
    inner_eq_zero_of_mem_span c (Finset.Iic k) i (by simp; exact hi) v hvH)
  rw [hST] at hhead
  have hn : (0 : ℝ) < ‖(v : E)‖ ^ 2 := by
    have : (v : E) ≠ 0 := hx0
    positivity
  have hnv : ‖v‖ = ‖(v : E)‖ := rfl
  rw [hnv] at hhead
  simp only [Submodule.subtype_apply] at hx
  nlinarith

theorem interlace_lower (hT : T.IsSymmetric) (hS : S.IsSymmetric)
    (hST : ∀ v w : V, @inner ℝ V _ (S v) w = @inner ℝ E _ (T v) w)
    (hE : finrank ℝ E = m + 1) (hV : finrank ℝ V = m) (k : Fin m) :
    hT.eigenvalues hE k.succ ≤ hS.eigenvalues hV k := by
  classical
  set c := hS.eigenvectorBasis hV
  set H : Submodule ℝ V := Submodule.span ℝ (Set.range fun j : (Finset.Ici k) => c j)
  have hHdim : finrank ℝ H = m - k := by rw [finrank_span_sub, Fin.card_Ici]
  obtain ⟨x, hxW, hx0, hx⟩ := CourantFischer.exists_ge_eigenvalue hT hE k.succ (H.map V.subtype)
    (by rw [finrank_map_subtype, hHdim]; simp)
  obtain ⟨v, hvH, rfl⟩ := Submodule.mem_map.mp hxW
  have htail := RayleighTail.rayleigh_tail_le hS hV k v (fun i hi =>
    inner_eq_zero_of_mem_span c (Finset.Ici k) i (by simp; exact hi) v hvH)
  rw [hST] at htail
  have hn : (0 : ℝ) < ‖(v : E)‖ ^ 2 := by
    have : (v : E) ≠ 0 := hx0
    positivity
  have hnv : ‖v‖ = ‖(v : E)‖ := rfl
  rw [hnv] at htail
  simp only [Submodule.subtype_apply] at hx
  nlinarith

end CauchyInterlace

#print axioms CauchyInterlace.interlace_upper
#print axioms CauchyInterlace.interlace_lower
