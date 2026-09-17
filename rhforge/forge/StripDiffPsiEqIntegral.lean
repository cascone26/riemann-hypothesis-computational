import CpowDiffEqIntegral
import PsiStripIntegral

set_option maxHeartbeats 1000000

/-!
# §1581 — per-strip combine: `((n+1)^{-s} - n^{-s})·ψ(n) = -s·∫_n^{n+1} ψ(t)t^{-s-1} dt`

Fourth Mellin brick: §1579 (the difference is `-s` times the strip integral of `t^{-s-1}`) composed with §1580
(`ψ` factors out of the strip integral).  Each Abel-difference term of §1578 is now a genuine integral of
`ψ(t)t^{-s-1}` over its strip — ready to be glued into `∫_1^N` by adjacent-interval additivity.
EXIT0 / std-axioms.
-/

open Complex Chebyshev

theorem strip_diff_psi_eq_integral {s : ℂ} (hs : s ≠ 0) {i : ℕ} (hi : 0 < i) :
    ((((i + 1 : ℕ) : ℂ)) ^ (-s) - (i : ℂ) ^ (-s)) * ((ψ (i : ℝ) : ℝ) : ℂ)
      = -s * ∫ t : ℝ in (i : ℝ)..((i : ℝ) + 1), ((ψ t : ℝ) : ℂ) * (t : ℂ) ^ (-s - 1) := by
  rw [psi_strip_integral hi]
  have h79 := cpow_diff_eq_integral hs hi
  have hbnd : (((i + 1 : ℕ) : ℝ)) = (i : ℝ) + 1 := by push_cast; ring
  rw [hbnd] at h79
  rw [h79]
  ring

#print axioms strip_diff_psi_eq_integral
