/-
PencilTriangularPosDef.lean (HP RH center, LEDGER S2908, 2026-11-07)
S2884 triangularity from POSITIVITY ALONE: if S ≻ 0 (certifiable per (λ,N), S2880/81/88), then for any frequency
ordering freq : Fin r → Fin (2N+1), the canonical S-Gram–Schmidt families E := gsS(evc ∘ freq), O := gsS(od ∘ freq)
give a lower-triangular pencil block  E_iᵀ (S D) O_j = 0  (i < j). No family hypotheses remain.
-/
import Mathlib
import PencilLowerTriangular
import SGramSchmidt

open Matrix

namespace PencilPD

variable {N r : ℕ}

theorem lower_triangular_of_posDef (S : Matrix (Fin (2 * N + 1)) (Fin (2 * N + 1)) ℝ) (hS : S.PosDef)
    (ν : Fin (2 * N + 1) → ℝ) (hodd : ∀ i, ν i.rev = -ν i) (c : Fin (2 * N + 1)) (hc : ν c = 0)
    (b : Fin (2 * N + 1) → ℝ)
    (hdisp : diagonal ν * S - S * diagonal ν = vecMulVec b (fun _ => 1) - vecMulVec (fun _ => 1) b)
    (freq : Fin r → Fin (2 * N + 1)) {i j : Fin r} (hij : i < j) :
    SGS.gsS S hS (ParityDisp.evc c ∘ freq) i ⬝ᵥ
      ((S * diagonal ν) *ᵥ SGS.gsS S hS (ParityDisp.od ∘ freq) j) = 0 := by
  set K : Fin r → Set (Fin (2 * N + 1)) := fun i => freq '' Set.Iic i
  have himg : ∀ (g : Fin (2 * N + 1) → Fin (2 * N + 1) → ℝ) (i : Fin r),
      (g ∘ freq) '' Set.Iic i = g '' K i := by
    intro g i; exact (Set.image_image g freq _).symm
  refine PencilTri.lower_triangular S ν hodd c hc b hdisp K
    (SGS.gsS S hS (ParityDisp.evc c ∘ freq)) (SGS.gsS S hS (ParityDisp.od ∘ freq)) ?_ ?_ ?_ hij
  · intro i
    have := SGS.gsS_mem_span S hS (ParityDisp.evc c ∘ freq) i
    rwa [himg] at this
  · intro j
    exact ⟨K j, by have := SGS.gsS_mem_span S hS (ParityDisp.od ∘ freq) j; rwa [himg] at this⟩
  · intro i j hij w hw
    rw [← himg] at hw
    exact SGS.gsS_orth S hS (ParityDisp.od ∘ freq) hij w hw

end PencilPD

#print axioms PencilPD.lower_triangular_of_posDef
