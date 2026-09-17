import DicksonSum

/-!
# §1815 — Dickson polynomials have real (conjugation-fixed) coefficients

`Polynomial.map conj (dickson n) = dickson n` over ℂ: the Dickson recurrence `D₀=2, D₁=X, D_{n+2}=X·D_{n+1}−D_n`
has integer coefficients, so complex conjugation fixes them.  Feeds §1816 (`Qpoly` real-coefficiented when the
palindromic data `v` is real ⟹ `Qpoly.roots` conjugation-closed ⟹ §1806 applies).  EXIT0 / std-axioms.
-/

open Polynomial

/-- **Dickson coefficients are conjugation-fixed.** -/
theorem dickson_map_conj (n : ℕ) :
    Polynomial.map (starRingEnd ℂ) (dickson n) = dickson n := by
  induction n using Nat.twoStepInduction with
  | zero => simp [dickson, map_ofNat]
  | one => simp [dickson]
  | more k ih0 ih1 =>
    rw [dickson, Polynomial.map_sub, Polynomial.map_mul, Polynomial.map_X, ih0, ih1]
