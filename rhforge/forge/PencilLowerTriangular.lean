/-
PencilLowerTriangular.lean (HP RH center, LEDGER S2905, 2026-11-07)
ALGEBRAIC HALF OF S2884: Shi's unshifted contrast pencil is lower-triangular in nested parity bases.
For ANY families E_i ∈ span{evc(k) : k ∈ K_i}, O_j ∈ span{od(k)} that are S-orthogonal to span{od(k) : k ∈ K_i} (i<j),
with displacement DS − SD = b1ᵀ − 1bᵀ (Shi 3.31), D = diag ν (ν odd, ν(center)=0):
    B_ij := E_iᵀ (S D) O_j = 0   for i < j.
(Gram–Schmidt in the S-inner product, frequency-ordered, is one instance.) Combined with PencilThetaInterlace (S2903)
this completes S2884 in Lean.
-/
import Mathlib
import PencilDisplacementSwap
import ParityDisplacementBasis

open Matrix

namespace PencilTri

variable {N : ℕ}

/-- the coordinate-sum functional -/
def sumL : (Fin (2 * N + 1) → ℝ) →ₗ[ℝ] ℝ := ∑ i, LinearMap.proj i

lemma sumL_apply (v : Fin (2 * N + 1) → ℝ) : sumL v = ∑ i, v i := by simp [sumL]

lemma sum_zero_of_mem_span_evc (c : Fin (2 * N + 1)) (K : Set (Fin (2 * N + 1))) (v : Fin (2 * N + 1) → ℝ)
    (hv : v ∈ Submodule.span ℝ (ParityDisp.evc c '' K)) : ∑ i, v i = 0 := by
  have hle : Submodule.span ℝ (ParityDisp.evc c '' K) ≤ LinearMap.ker (sumL (N := N)) := by
    rw [Submodule.span_le]; rintro _ ⟨k, -, rfl⟩
    simp [sumL_apply, ParityDisp.sum_evc]
  simpa [sumL_apply] using hle hv

lemma sum_zero_of_mem_span_od (K : Set (Fin (2 * N + 1))) (v : Fin (2 * N + 1) → ℝ)
    (hv : v ∈ Submodule.span ℝ (ParityDisp.od '' K)) : ∑ i, v i = 0 := by
  have hle : Submodule.span ℝ (ParityDisp.od '' K) ≤ LinearMap.ker (sumL (N := N)) := by
    rw [Submodule.span_le]; rintro _ ⟨k, -, rfl⟩
    simp [sumL_apply, ParityDisp.sum_od]
  simpa [sumL_apply] using hle hv

theorem lower_triangular {r : ℕ} (S : Matrix (Fin (2 * N + 1)) (Fin (2 * N + 1)) ℝ)
    (ν : Fin (2 * N + 1) → ℝ) (hodd : ∀ i, ν i.rev = -ν i) (c : Fin (2 * N + 1)) (hc : ν c = 0)
    (b : Fin (2 * N + 1) → ℝ)
    (hdisp : diagonal ν * S - S * diagonal ν = vecMulVec b (fun _ => 1) - vecMulVec (fun _ => 1) b)
    (K : Fin r → Set (Fin (2 * N + 1)))
    (E O : Fin r → Fin (2 * N + 1) → ℝ)
    (hE : ∀ i, E i ∈ Submodule.span ℝ (ParityDisp.evc c '' K i))
    (hOmem : ∀ j, ∃ Kj : Set (Fin (2 * N + 1)), O j ∈ Submodule.span ℝ (ParityDisp.od '' Kj))
    (hO : ∀ i j, i < j → ∀ w ∈ Submodule.span ℝ (ParityDisp.od '' K i), w ⬝ᵥ (S *ᵥ O j) = 0)
    {i j : Fin r} (hij : i < j) :
    E i ⬝ᵥ ((S * diagonal ν) *ᵥ O j) = 0 := by
  obtain ⟨Kj, hKj⟩ := hOmem j
  refine PencilSwap.form_vanishes S (diagonal ν) b hdisp (Submodule.span ℝ (ParityDisp.od '' K i))
    (E i) (O j) (sum_zero_of_mem_span_evc c (K i) (E i) (hE i)) (sum_zero_of_mem_span_od Kj (O j) hKj)
    ?_ (hO i j hij)
  rw [diagonal_transpose]
  have hmap := ParityDisp.D_span_le ν hodd c hc (K i)
  exact hmap (Submodule.mem_map_of_mem (f := Matrix.mulVecLin (diagonal ν)) (hE i))

end PencilTri

#print axioms PencilTri.lower_triangular
