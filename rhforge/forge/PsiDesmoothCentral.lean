import PsiDesmoothError
import PsiDesmoothBackward

set_option maxHeartbeats 1000000

open Chebyshev

/-- §1651 — central-difference de-smoothing error: the *symmetric* difference quotient of the
smoothed `Ψℝ` approximates `ψ(x)` with a two-sided error controlled by half the forward/backward
jumps: `-(ψ x - ψ(x-h))/2 ≤ (Ψℝ(x+h)-Ψℝ(x-h))/(2h) - ψ x ≤ (ψ(x+h)-ψ x)/2`. Combines the forward
(§1640) and backward (§1642) error bounds via `(Q_f + Q_b)/2 = Q_central` — the central quotient
converges to `ψ(x)` at least twice as tightly (both jump terms halved) as either one-sided quotient,
the standard reason central differences are preferred for de-smoothing/numerical differentiation. -/
theorem psi_desmooth_central_error {x h : ℝ} (hh : 0 < h) :
    -(ψ x - ψ (x - h)) / 2 ≤ (psi1R (x + h) - psi1R (x - h)) / (2 * h) - ψ x
      ∧ (psi1R (x + h) - psi1R (x - h)) / (2 * h) - ψ x ≤ (ψ (x + h) - ψ x) / 2 := by
  have hf := psi_desmooth_forward_error (x := x) hh
  have hb := psi_desmooth_backward_error (x := x) hh
  have hcentral : (psi1R (x + h) - psi1R (x - h)) / (2 * h)
      = ((psi1R (x + h) - psi1R x) / h + (psi1R x - psi1R (x - h)) / h) / 2 := by
    field_simp
    ring
  rw [hcentral]
  exact ⟨by linarith [hf.1, hf.2, hb.1, hb.2], by linarith [hf.1, hf.2, hb.1, hb.2]⟩

#print axioms psi_desmooth_central_error
