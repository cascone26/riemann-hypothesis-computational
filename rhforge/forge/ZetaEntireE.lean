import Mathlib

/-!
# S2979 — E(s) = (s−1)ζ(s) is ENTIRE (lane 1: Phragmén–Lindelöf input)

Mathlib's `PhragmenLindelof.vertical_strip` needs `DiffContOnCl ℂ f (re ⁻¹' Ioo a b)`; ζ itself has the pole at 1.
* `zetaE := update (fun s ↦ (s − 1)·ζ(s)) 1 1`.
* `zetaE_eq` — s ≠ 1 ⇒ zetaE s = (s − 1)·ζ(s);  `zetaE_one` — zetaE 1 = 1.
* `zetaE_differentiable` — `Differentiable ℂ zetaE`  (removable singularity: `riemannZeta_residue_one` +
  `differentiableOn_compl_singleton_and_continuousAt_iff`, the same pattern as Mathlib's Nonvanishing).
* `zetaE_diffContOnCl` — DiffContOnCl on every set, in particular on the critical strip.
Next: the growth hypothesis exp(B·exp(c|t|)) on 0 < Re s < 1 (forge strip bounds are polynomial) and the
normaliser for the convexity exponent.
-/

open Complex Filter Topology Set

namespace ZetaEntire

noncomputable def zetaE : ℂ → ℂ := Function.update (fun s => (s - 1) * riemannZeta s) 1 1

theorem zetaE_eq {s : ℂ} (hs : s ≠ 1) : zetaE s = (s - 1) * riemannZeta s := by
  simp [zetaE, Function.update_of_ne hs]

theorem zetaE_one : zetaE 1 = 1 := by simp [zetaE]

theorem zetaE_differentiable : Differentiable ℂ zetaE := by
  have hon : DifferentiableOn ℂ zetaE (univ : Set ℂ) := by
    refine (differentiableOn_compl_singleton_and_continuousAt_iff (c := 1) Filter.univ_mem).mp
      ⟨?_, ?_⟩
    · intro x hx
      have hx1 : x ≠ 1 := hx.2
      have hd : DifferentiableAt ℂ (fun s => (s - 1) * riemannZeta s) x :=
        (differentiableAt_id.sub_const 1).mul (differentiableAt_riemannZeta hx1)
      have heq : (fun s => (s - 1) * riemannZeta s) =ᶠ[𝓝 x] zetaE := by
        filter_upwards [isOpen_ne.mem_nhds hx1] with y hy
        exact (zetaE_eq hy).symm
      exact (hd.congr_of_eventuallyEq heq.symm).differentiableWithinAt
    · simpa only [zetaE, continuousAt_update_same] using riemannZeta_residue_one
  exact differentiableOn_univ.mp hon

theorem zetaE_diffContOnCl (S : Set ℂ) : DiffContOnCl ℂ zetaE S :=
  zetaE_differentiable.diffContOnCl

end ZetaEntire

#print axioms ZetaEntire.zetaE_eq
#print axioms ZetaEntire.zetaE_one
#print axioms ZetaEntire.zetaE_differentiable
#print axioms ZetaEntire.zetaE_diffContOnCl
