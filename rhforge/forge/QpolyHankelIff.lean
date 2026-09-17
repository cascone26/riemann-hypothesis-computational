import Mathlib.Analysis.Complex.Polynomial.Basic
import HankelIff
import QpolyMapConj

/-!
# §1817 — Qpoly-level Hankel criterion (on-circle ⊗ real-rootedness, CAPSTONE)

For real palindromic data `v`, the complex Hankel form over the Chebyshev-image polynomial's roots
`(Qpoly v m).roots` is positive semidefinite ⟺ all those roots are real.  Composes §1816 (`Qpoly` real-
coefficiented ⟹ its roots are conjugation-closed, via `roots_map` over ℂ) with §1806 (`hankel_psd_iff_all_real`).
Via §1814 (P-root ⟺ Qpoly-root at `ρ+ρ⁻¹`) and §1812/§1781 (`|ρ|=1 ⟺ ρ+ρ⁻¹` real in `[−2,2]`), this ties the
on-circle module to the Hankel-PSD iff: the CS/Chebyshev mechanism realized through a checkable positivity
condition.  EXIT0 / std-axioms.
-/

open Polynomial Finset

private theorem conj_inj : Function.Injective (starRingEnd ℂ) :=
  Function.LeftInverse.injective (fun z => Complex.conj_conj z)

/-- **Qpoly-level Hankel criterion.**  Complex Hankel form over `Qpoly`'s roots PSD ⟺ all roots real. -/
theorem qpoly_hankel_psd_iff_roots_real (v : ℕ → ℂ) (m : ℕ) (hv : ∀ i, (v i).im = 0)
    (hsimple : ∀ r ∈ (Qpoly v m).roots, (Qpoly v m).roots.count r = 1) :
    (∀ (n : ℕ) (c : Fin (n + 1) → ℝ),
        0 ≤ (((Qpoly v m).roots.map
          (fun r => (∑ i : Fin (n + 1), (c i : ℂ) * r ^ (i : ℕ)) ^ 2)).sum).re)
      ↔ (∀ r ∈ (Qpoly v m).roots, r.im = 0) := by
  have hcard : (Qpoly v m).roots.card = (Qpoly v m).natDegree :=
    IsAlgClosed.card_roots_eq_natDegree
  have hcc : (Qpoly v m).roots.map (starRingEnd ℂ) = (Qpoly v m).roots := by
    rw [roots_map_of_injective_of_card_eq_natDegree conj_inj hcard, Qpoly_map_conj v m hv]
  exact hankel_psd_iff_all_real (Qpoly v m).roots hcc hsimple
