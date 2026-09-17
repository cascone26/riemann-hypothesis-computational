import PsiDesmoothSqueeze

set_option maxHeartbeats 1000000

open Chebyshev

/-- §1642 — backward de-smoothing error: the backward difference quotient of `Ψℝ` approximates `ψ(x)` from below
with error at most `ψ(x) - ψ(x-h)`.  `0 ≤ ψ(x) - (Ψℝ(x)-Ψℝ(x-h))/h ≤ ψ(x) - ψ(x-h)`.  Companion of §1640;
together they two-side the de-smoothing of `ψ`. -/
theorem psi_desmooth_backward_error {x h : ℝ} (hh : 0 < h) :
    0 ≤ ψ x - (psi1R x - psi1R (x - h)) / h
      ∧ ψ x - (psi1R x - psi1R (x - h)) / h ≤ ψ x - ψ (x - h) := by
  have hhi : (psi1R x - psi1R (x - h)) / h ≤ ψ x := (psi_desmooth_squeeze hh).1
  have hlo : ψ (x - h) ≤ (psi1R x - psi1R (x - h)) / h := by
    rw [le_div_iff₀ hh, mul_comm]
    have := (psi1R_sandwich (x := x - h) hh).1
    rwa [sub_add_cancel] at this
  exact ⟨by linarith, by linarith⟩

#print axioms psi_desmooth_backward_error
