/-
ParityDisplacementBasis.lean (HP RH center, LEDGER S2904, 2026-11-07)
Algebraic half of S2884, concrete parity step: with an odd frequency map ν on Fin (2N+1) (ν(rev i) = −ν i, ν(center)=0),
D = diag ν maps even parity vectors to odd ones FREQUENCY-WISE:
  D·ev(k) = ν_k·od(k),  D·e_c = 0,  D·(ev(k) − 2 e_c) = ν_k·od(k)  (zero-mean correction is invisible to D),
hence D(span{evc(k) : k ∈ K}) ⊆ span{od(k) : k ∈ K} — the step that makes the pencil block lower-triangular.
-/
import Mathlib

open Matrix

namespace ParityDisp

variable {N : ℕ}

def ev (k : Fin (2 * N + 1)) : Fin (2 * N + 1) → ℝ := Pi.single k 1 + Pi.single k.rev 1
def od (k : Fin (2 * N + 1)) : Fin (2 * N + 1) → ℝ := Pi.single k 1 - Pi.single k.rev 1
def evc (c k : Fin (2 * N + 1)) : Fin (2 * N + 1) → ℝ := ev k - (2 : ℝ) • Pi.single c 1

variable (ν : Fin (2 * N + 1) → ℝ) (hodd : ∀ i, ν i.rev = -ν i)

lemma diag_single (j : Fin (2 * N + 1)) :
    diagonal ν *ᵥ Pi.single j 1 = ν j • Pi.single j 1 := by
  ext i; by_cases h : i = j <;> simp [mulVec, dotProduct, diagonal, Pi.single_apply, h]

include hodd in
theorem D_ev (k : Fin (2 * N + 1)) : diagonal ν *ᵥ ev k = ν k • od k := by
  simp only [ev, od, mulVec_add, diag_single, hodd, neg_smul, smul_sub]
  abel

theorem D_center (c : Fin (2 * N + 1)) (hc : ν c = 0) : diagonal ν *ᵥ Pi.single c 1 = 0 := by
  rw [diag_single, hc, zero_smul]

include hodd in
theorem D_evc (c : Fin (2 * N + 1)) (hc : ν c = 0) (k : Fin (2 * N + 1)) :
    diagonal ν *ᵥ evc c k = ν k • od k := by
  rw [evc, mulVec_sub, mulVec_smul, D_center ν c hc, smul_zero, sub_zero, D_ev ν hodd]

include hodd in
/-- D maps the span of zero-mean even contrasts at frequencies K into the span of odd vectors at the same frequencies. -/
theorem D_span_le (c : Fin (2 * N + 1)) (hc : ν c = 0) (K : Set (Fin (2 * N + 1))) :
    (Submodule.span ℝ (evc c '' K)).map (Matrix.mulVecLin (diagonal ν))
      ≤ Submodule.span ℝ (od '' K) := by
  rw [Submodule.map_span_le]
  rintro _ ⟨k, hk, rfl⟩
  rw [mulVecLin_apply, D_evc ν hodd c hc]
  exact Submodule.smul_mem _ _ (Submodule.subset_span ⟨k, hk, rfl⟩)

/-- zero-mean: every odd vector and every corrected even vector has zero total sum (when k ≠ center). -/
theorem sum_od (k : Fin (2 * N + 1)) : ∑ i, od k i = 0 := by
  simp [od, Finset.sum_sub_distrib]

theorem sum_evc (c k : Fin (2 * N + 1)) : ∑ i, evc c k i = 0 := by
  simp [evc, ev, Finset.sum_add_distrib, Finset.sum_sub_distrib, ← Finset.mul_sum]; norm_num

end ParityDisp

#print axioms ParityDisp.D_span_le
#print axioms ParityDisp.sum_evc
