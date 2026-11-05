/-
SignChangeCount.lean  (HP RH center, LEDGER S2862, 2026-11-04)

Summed form of the parity lemma (PartialFractionSignChange, S2861): for strictly increasing poles
x₀<…<x_m and nonzero residues u, with R(y)=∑ₖ uₖ/(y−xₖ):
  #{j : u_j·u_{j+1} < 0}  ≥  m − #{j : R has a zero in (x_j, x_{j+1})}.
This is the combinatorial half of the S2860 zero-counting sign-change law for the CCM/Shi ground state.
-/
import Mathlib
import PartialFractionSignChange

open Set

namespace PFSignChange

lemma consecutive_out {m : ℕ} (x : Fin (m+1) → ℝ) (hx : StrictMono x) (j : Fin m) :
    ∀ k, k ≠ j.castSucc → k ≠ j.succ → x k ∉ Icc (x j.castSucc) (x j.succ) := by
  intro k h1 h2 hk
  have hv1 : k.val ≠ j.val := fun h => h1 (Fin.ext (by simpa using h))
  have hv2 : k.val ≠ j.val + 1 := fun h => h2 (Fin.ext (by simpa using h))
  rcases lt_or_gt_of_ne hv1 with h | h
  · have : k < j.castSucc := by rw [Fin.lt_def]; simpa using h
    exact absurd hk.1 (not_le.mpr (hx this))
  · have h' : j.val + 1 < k.val := by omega
    have : j.succ < k := by rw [Fin.lt_def]; simpa using h'
    exact absurd hk.2 (not_le.mpr (hx this))

open Classical in
/-- **Sign-change count ≥ (#intervals) − (#intervals carrying a zero of R).** -/
theorem signChanges_ge {m : ℕ} (x u : Fin (m+1) → ℝ) (hx : StrictMono x) (hu : ∀ k, u k ≠ 0) :
    m - (Finset.univ.filter
          (fun j : Fin m => ∃ y ∈ Ioo (x j.castSucc) (x j.succ), R x u y = 0)).card
      ≤ (Finset.univ.filter (fun j : Fin m => u j.castSucc * u j.succ < 0)).card := by
  set B := Finset.univ.filter
    (fun j : Fin m => ∃ y ∈ Ioo (x j.castSucc) (x j.succ), R x u y = 0)
  set S := Finset.univ.filter (fun j : Fin m => u j.castSucc * u j.succ < 0)
  have hsub : Finset.univ.filter
      (fun j : Fin m => ¬ ∃ y ∈ Ioo (x j.castSucc) (x j.succ), R x u y = 0) ⊆ S := by
    intro j hj
    simp only [Finset.mem_filter, Finset.mem_univ, true_and] at hj
    simp only [S, Finset.mem_filter, Finset.mem_univ, true_and]
    have hlt : x j.castSucc < x j.succ := hx j.castSucc_lt_succ
    exact mul_neg_of_no_zero_between x u j.castSucc j.succ hlt (consecutive_out x hx j)
      (hu _) (hu _) (fun y hy h0 => hj ⟨y, hy, h0⟩)
  have hcard : B.card + (Finset.univ.filter
      (fun j : Fin m => ¬ ∃ y ∈ Ioo (x j.castSucc) (x j.succ), R x u y = 0)).card = m := by
    have := Finset.card_filter_add_card_filter_not (s := (Finset.univ : Finset (Fin m)))
      (fun j : Fin m => ∃ y ∈ Ioo (x j.castSucc) (x j.succ), R x u y = 0)
    simpa [B] using this
  have hle := Finset.card_le_card hsub
  omega

end PFSignChange

#print axioms PFSignChange.signChanges_ge
