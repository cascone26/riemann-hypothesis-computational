/-
CauchyInterlacingIso.lean (HP RH center, LEDGER S2901, 2026-11-07)
Cauchy interlacing through a linear isometry ι : F →ₗᵢ E (generalizing S2899's subspace version), so it applies
directly to matrices and their principal submatrices (ι = extension by zero) with no eigenvalue transport.
finrank F = m, finrank E = m+1, S on F with ⟪S v, w⟫ = ⟪T (ι v), ι w⟫:  λ_{k+1}(T) ≤ μ_k(S) ≤ λ_k(T).
-/
import Mathlib
import RayleighEigenTail
import CourantFischerMinMax
import CauchyInterlacing

open Module

namespace CauchyInterlaceIso

variable {E F : Type*} [NormedAddCommGroup E] [InnerProductSpace ℝ E] [FiniteDimensional ℝ E]
  [NormedAddCommGroup F] [InnerProductSpace ℝ F] [FiniteDimensional ℝ F]
variable {m : ℕ} {T : E →ₗ[ℝ] E} {S : F →ₗ[ℝ] F} (ι : F →ₗᵢ[ℝ] E)

lemma finrank_map_iso (W : Submodule ℝ F) : finrank ℝ (W.map ι.toLinearMap) = finrank ℝ W :=
  (Submodule.equivMapOfInjective ι.toLinearMap ι.injective W).finrank_eq.symm

theorem interlace_upper (hT : T.IsSymmetric) (hS : S.IsSymmetric)
    (hST : ∀ v w : F, @inner ℝ F _ (S v) w = @inner ℝ E _ (T (ι v)) (ι w))
    (hE : finrank ℝ E = m + 1) (hF : finrank ℝ F = m) (k : Fin m) :
    hS.eigenvalues hF k ≤ hT.eigenvalues hE k.castSucc := by
  classical
  set c := hS.eigenvectorBasis hF
  set H : Submodule ℝ F := Submodule.span ℝ (Set.range fun j : (Finset.Iic k) => c j)
  have hHdim : finrank ℝ H = k + 1 := by rw [CauchyInterlace.finrank_span_sub, Fin.card_Iic]
  obtain ⟨x, hxW, hx0, hx⟩ := CourantFischer.exists_le_eigenvalue hT hE k.castSucc
    (H.map ι.toLinearMap) (by rw [finrank_map_iso, hHdim]; simp)
  obtain ⟨v, hvH, rfl⟩ := Submodule.mem_map.mp hxW
  have hhead := RayleighTail.rayleigh_head_ge hS hF k v (fun i hi =>
    CauchyInterlace.inner_eq_zero_of_mem_span c (Finset.Iic k) i (by simp; exact hi) v hvH)
  rw [hST] at hhead
  have hv0 : v ≠ 0 := fun h => hx0 (by simp [h])
  have hn : (0 : ℝ) < ‖v‖ ^ 2 := by positivity
  have hnorm : ‖ι v‖ = ‖v‖ := ι.norm_map v
  simp only [LinearIsometry.coe_toLinearMap] at hx
  rw [hnorm] at hx
  nlinarith

theorem interlace_lower (hT : T.IsSymmetric) (hS : S.IsSymmetric)
    (hST : ∀ v w : F, @inner ℝ F _ (S v) w = @inner ℝ E _ (T (ι v)) (ι w))
    (hE : finrank ℝ E = m + 1) (hF : finrank ℝ F = m) (k : Fin m) :
    hT.eigenvalues hE k.succ ≤ hS.eigenvalues hF k := by
  classical
  set c := hS.eigenvectorBasis hF
  set H : Submodule ℝ F := Submodule.span ℝ (Set.range fun j : (Finset.Ici k) => c j)
  have hHdim : finrank ℝ H = m - k := by rw [CauchyInterlace.finrank_span_sub, Fin.card_Ici]
  obtain ⟨x, hxW, hx0, hx⟩ := CourantFischer.exists_ge_eigenvalue hT hE k.succ
    (H.map ι.toLinearMap) (by rw [finrank_map_iso, hHdim]; simp)
  obtain ⟨v, hvH, rfl⟩ := Submodule.mem_map.mp hxW
  have htail := RayleighTail.rayleigh_tail_le hS hF k v (fun i hi =>
    CauchyInterlace.inner_eq_zero_of_mem_span c (Finset.Ici k) i (by simp; exact hi) v hvH)
  rw [hST] at htail
  have hv0 : v ≠ 0 := fun h => hx0 (by simp [h])
  have hn : (0 : ℝ) < ‖v‖ ^ 2 := by positivity
  have hnorm : ‖ι v‖ = ‖v‖ := ι.norm_map v
  simp only [LinearIsometry.coe_toLinearMap] at hx
  rw [hnorm] at hx
  nlinarith

end CauchyInterlaceIso

#print axioms CauchyInterlaceIso.interlace_upper
#print axioms CauchyInterlaceIso.interlace_lower
