import Psi1RBridge

set_option maxHeartbeats 1000000

/-!
# §1609 — the de-smoothing squeeze: `(Ψℝ(x)-Ψℝ(x-h))/h ≤ ψ(x) ≤ (Ψℝ(x+h)-Ψℝ(x))/h`

`ψ` is trapped between the backward and forward difference quotients of the smoothed `Ψℝ = ∫₁·ψ`.  This is the exact
two-sided inequality that converts a `Ψ`-asymptotic (`Ψℝ(x) = x²/2 + O(error)`, from the contour shift depositing
the §1608 residue) into a `ψ`-asymptotic (`ψ(x) = x + O(error/h + h)` after optimizing `h`) — the last analytic
step of the prime number theorem with error term.  Both halves come from the monotone sandwich §1604/§1605 applied
at `x` (forward) and `x-h` (backward).  Contour-free.  EXIT0 / std-axioms.
-/

open Chebyshev

theorem psi_desmooth_squeeze {x h : ℝ} (hh : 0 < h) :
    (psi1R x - psi1R (x - h)) / h ≤ ψ x ∧ ψ x ≤ (psi1R (x + h) - psi1R x) / h := by
  constructor
  · -- backward quotient ≤ ψ x
    rw [div_le_iff₀ hh]
    have hback := (psi1R_sandwich (x := x - h) hh).2
    rw [sub_add_cancel] at hback
    rw [mul_comm] at hback
    exact hback
  · -- ψ x ≤ forward quotient
    rw [le_div_iff₀ hh]
    have hfwd := (psi1R_sandwich (x := x) hh).1
    rw [mul_comm] at hfwd
    exact hfwd

#print axioms psi_desmooth_squeeze
