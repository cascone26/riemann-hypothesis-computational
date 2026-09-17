import Mathlib.Tactic

/-!
# §1807 — Interval converse core: a real root outside [−2,2] breaks the weighted form

Companion to §1791 (roots in [−2,2] ⟹ `(4−r²)`-weighted Hankel form PSD).  Here the CONVERSE seed: if a real
root `r₀` lies OUTSIDE `[−2,2]` (`r₀² > 4`, so the weight `4−r₀² < 0`) and a real vector isolates it
(`P_c(r₀)=1`, `P_c=0` on the rest), the weighted form equals `4−r₀² < 0` — not PSD.  With §1791 this brackets
the interval both ways.  Unlike the complex case (§1803) no conjugate pairing is needed — the weight sign does
the work.  EXIT0 / std-axioms.
-/

open Finset

variable {d : ℕ}

/-- **Interval converse core.**  An isolated real root outside `[−2,2]` drives the `(4−r²)`-weighted Hankel
form strictly negative. -/
theorem interval_isolated_neg (rest : Multiset ℝ) (c : Fin (d + 1) → ℝ) (r₀ : ℝ)
    (hr₀ : (∑ i : Fin (d + 1), c i * r₀ ^ (i : ℕ)) = 1)
    (hrest : ∀ s ∈ rest, (∑ i : Fin (d + 1), c i * s ^ (i : ℕ)) = 0)
    (hout : 4 < r₀ ^ 2) :
    ((r₀ ::ₘ rest).map (fun r => (4 - r ^ 2) * (∑ i : Fin (d + 1), c i * r ^ (i : ℕ)) ^ 2)).sum < 0 := by
  have hval : ((r₀ ::ₘ rest).map
      (fun r => (4 - r ^ 2) * (∑ i : Fin (d + 1), c i * r ^ (i : ℕ)) ^ 2)).sum = 4 - r₀ ^ 2 := by
    simp only [Multiset.map_cons, Multiset.sum_cons]
    rw [hr₀]
    have hz : (rest.map (fun r => (4 - r ^ 2) * (∑ i : Fin (d + 1), c i * r ^ (i : ℕ)) ^ 2)).sum = 0 := by
      have hzero : rest.map (fun r => (4 - r ^ 2) * (∑ i : Fin (d + 1), c i * r ^ (i : ℕ)) ^ 2)
          = rest.map (fun _ => (0 : ℝ)) := by
        refine Multiset.map_congr rfl (fun s hs => ?_)
        rw [hrest s hs]; ring
      rw [hzero]; simp
    rw [hz]; ring
  rw [hval]; linarith
