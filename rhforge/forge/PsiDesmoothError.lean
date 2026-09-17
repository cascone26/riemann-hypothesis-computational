import PsiDesmoothSqueeze

set_option maxHeartbeats 1000000

open Chebyshev

/-- §1640 — the de-smoothing error bound: the forward difference quotient of the smoothed `Ψℝ` approximates `ψ(x)`
from above with error at most `ψ(x+h) - ψ(x)`.  Concretely `0 ≤ (Ψℝ(x+h)-Ψℝ(x))/h - ψ(x) ≤ ψ(x+h) - ψ(x)`.
Combined with a smoothed asymptotic `Ψℝ(x)=x²/2+O(err)` and optimizing `h`, this converts to the `ψ`-asymptotic
(PNT with error). From §1609 (ψ(x) ≤ forward quotient) + §1605 (forward quotient ≤ ψ(x+h)). -/
theorem psi_desmooth_forward_error {x h : ℝ} (hh : 0 < h) :
    0 ≤ (psi1R (x + h) - psi1R x) / h - ψ x
      ∧ (psi1R (x + h) - psi1R x) / h - ψ x ≤ ψ (x + h) - ψ x := by
  have hlo : ψ x ≤ (psi1R (x + h) - psi1R x) / h := (psi_desmooth_squeeze hh).2
  have hhi : (psi1R (x + h) - psi1R x) / h ≤ ψ (x + h) := by
    rw [div_le_iff₀ hh, mul_comm]; exact (psi1R_sandwich hh).2
  exact ⟨by linarith, by linarith⟩

#print axioms psi_desmooth_forward_error
