import PsiDesmoothCentral

set_option maxHeartbeats 1000000

open Chebyshev

/-- §1653 — the central de-smoothing error, in clean absolute-value form:
`|Q_central - ψ x| ≤ (ψ(x+h) - ψ(x-h)) / 2`. Sharpens §1651's two-sided asymmetric bound into a
single usable bound by combining it with `ψ`'s monotonicity (`Chebyshev.psi_mono`, Mathlib): since
`ψ(x-h) ≤ ψ x ≤ ψ(x+h)`, both of §1651's half-widths `(ψx-ψ(x-h))/2` and `(ψ(x+h)-ψx)/2` are each
`≤ (ψ(x+h)-ψ(x-h))/2`, so the sharp asymmetric interval sits inside the symmetric one. -/
theorem psi_desmooth_central_error_abs {x h : ℝ} (hh : 0 < h) :
    |(psi1R (x + h) - psi1R (x - h)) / (2 * h) - ψ x| ≤ (ψ (x + h) - ψ (x - h)) / 2 := by
  have hmono1 : ψ (x - h) ≤ ψ x := psi_mono (by linarith)
  have hmono2 : ψ x ≤ ψ (x + h) := psi_mono (by linarith)
  have hc := psi_desmooth_central_error (x := x) hh
  rw [abs_le]
  exact ⟨by linarith [hc.1], by linarith [hc.2]⟩

#print axioms psi_desmooth_central_error_abs
