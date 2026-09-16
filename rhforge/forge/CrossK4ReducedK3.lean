import CrossK4LeadingCoeff
import CrossK4ThirteenthCoeff

set_option maxHeartbeats 4000000
set_option maxRecDepth 4000

/-!
# K=4 decoupling: eliminate D3_4 from the degree-3 equation via the degree-15 equation

Direct 4-block extension of the K=3 ReducedK{...} cascade, continuing S2536-2538. Sympy-derived AND
verified BEFORE writing this Lean (derive_and_verify_k4_reductions.py): combo = crossK4ThirteenthCoeff - (octic-style
coefficient, b=0) * crossK4LeadingCoeff eliminates D3_4 completely (verified via the 2x2-per-block
linear-system solve showing D3_4's multiplier is EXACTLY 0), and the resulting closed form (grouped by
each surviving block's coeff1/coeff3, with pure-z multiplier polynomials machine-translated from the
verified sympy expressions via sympy_to_lean.py -- not hand-transcribed, avoiding the S2537 bug class)
was independently reconstructed and confirmed to exactly match the raw combination (diff=0).

Regenerated 2026-09-16 (S2561) after the original was lost to a scratchpad rotation -- reconstructed
sympy_to_lean.py from its functional spec (LEDGER S2539) and re-derived via this same generic script;
output diffed byte-for-byte against the one surviving original (CrossK4ReducedK3.lean) to validate the
reconstruction before trusting K13/K11/K9/K7/K5's fresh output. -/

theorem crossK4ReducedK3 (N : ℕ) (hN : 8 ≤ N) (L : ℝ) (hL : L ≠ 0) (z1 z2 z3 z4 : ℂ)
    (hz1_im : z1.im ≠ 0) (hz2_im : z2.im ≠ 0) (hz3_im : z3.im ≠ 0) (hz4_im : z4.im ≠ 0)
    (p1 q1 r1 s1 p2 q2 r2 s2 p3 q3 r3 s3 p4 q4 r4 s4 : ℂ)
    (hsum : p1 • vPlusLC N L z1 + q1 • vMinusLC N L z1
        + r1 • vPlusLC N L (starRingEnd ℂ z1) + s1 • vMinusLC N L (starRingEnd ℂ z1)
        + p2 • vPlusLC N L z2 + q2 • vMinusLC N L z2
        + r2 • vPlusLC N L (starRingEnd ℂ z2) + s2 • vMinusLC N L (starRingEnd ℂ z2)
        + p3 • vPlusLC N L z3 + q3 • vMinusLC N L z3
        + r3 • vPlusLC N L (starRingEnd ℂ z3) + s3 • vMinusLC N L (starRingEnd ℂ z3)
        + p4 • vPlusLC N L z4 + q4 • vMinusLC N L z4
        + r4 • vPlusLC N L (starRingEnd ℂ z4) + s4 • vMinusLC N L (starRingEnd ℂ z4) = 0) :
    ((q1 - p1) * (((Real.sqrt 2 : ℝ) : ℂ) * Complex.sin (L * z1 / 2)) * (starRingEnd ℂ z1) ^ 2 + (s1 - r1) * (starRingEnd ℂ (((Real.sqrt 2 : ℝ) : ℂ) * Complex.sin (L * z1 / 2))) * z1 ^ 2)
        * (-z2^2*(starRingEnd ℂ z2)^2*z3^2*(starRingEnd ℂ z3)^2*z4^2 - z2^2*(starRingEnd ℂ z2)^2*z3^2*(starRingEnd ℂ z3)^2*(starRingEnd ℂ z4)^2 - z2^2*(starRingEnd ℂ z2)^2*z3^2*z4^2*(starRingEnd ℂ z4)^2 - z2^2*(starRingEnd ℂ z2)^2*(starRingEnd ℂ z3)^2*z4^2*(starRingEnd ℂ z4)^2 - z2^2*z3^2*(starRingEnd ℂ z3)^2*z4^2*(starRingEnd ℂ z4)^2 - (starRingEnd ℂ z2)^2*z3^2*(starRingEnd ℂ z3)^2*z4^2*(starRingEnd ℂ z4)^2)
      + ((q2 - p2) * (((Real.sqrt 2 : ℝ) : ℂ) * Complex.sin (L * z2 / 2)) * (starRingEnd ℂ z2) ^ 2 + (s2 - r2) * (starRingEnd ℂ (((Real.sqrt 2 : ℝ) : ℂ) * Complex.sin (L * z2 / 2))) * z2 ^ 2)
        * (-z1^2*(starRingEnd ℂ z1)^2*z3^2*(starRingEnd ℂ z3)^2*z4^2 - z1^2*(starRingEnd ℂ z1)^2*z3^2*(starRingEnd ℂ z3)^2*(starRingEnd ℂ z4)^2 - z1^2*(starRingEnd ℂ z1)^2*z3^2*z4^2*(starRingEnd ℂ z4)^2 - z1^2*(starRingEnd ℂ z1)^2*(starRingEnd ℂ z3)^2*z4^2*(starRingEnd ℂ z4)^2 - z1^2*z3^2*(starRingEnd ℂ z3)^2*z4^2*(starRingEnd ℂ z4)^2 - (starRingEnd ℂ z1)^2*z3^2*(starRingEnd ℂ z3)^2*z4^2*(starRingEnd ℂ z4)^2)
      + ((q3 - p3) * (((Real.sqrt 2 : ℝ) : ℂ) * Complex.sin (L * z3 / 2)) * (starRingEnd ℂ z3) ^ 2 + (s3 - r3) * (starRingEnd ℂ (((Real.sqrt 2 : ℝ) : ℂ) * Complex.sin (L * z3 / 2))) * z3 ^ 2)
        * (-z1^2*(starRingEnd ℂ z1)^2*z2^2*(starRingEnd ℂ z2)^2*z4^2 - z1^2*(starRingEnd ℂ z1)^2*z2^2*(starRingEnd ℂ z2)^2*(starRingEnd ℂ z4)^2 - z1^2*(starRingEnd ℂ z1)^2*z2^2*z4^2*(starRingEnd ℂ z4)^2 - z1^2*(starRingEnd ℂ z1)^2*(starRingEnd ℂ z2)^2*z4^2*(starRingEnd ℂ z4)^2 - z1^2*z2^2*(starRingEnd ℂ z2)^2*z4^2*(starRingEnd ℂ z4)^2 - (starRingEnd ℂ z1)^2*z2^2*(starRingEnd ℂ z2)^2*z4^2*(starRingEnd ℂ z4)^2)
      + ((q4 - p4) * (((Real.sqrt 2 : ℝ) : ℂ) * Complex.sin (L * z4 / 2)) * (starRingEnd ℂ z4) ^ 2 + (s4 - r4) * (starRingEnd ℂ (((Real.sqrt 2 : ℝ) : ℂ) * Complex.sin (L * z4 / 2))) * z4 ^ 2)
        * (-z1^2*(starRingEnd ℂ z1)^2*z2^2*(starRingEnd ℂ z2)^2*z3^2 - z1^2*(starRingEnd ℂ z1)^2*z2^2*(starRingEnd ℂ z2)^2*(starRingEnd ℂ z3)^2 - z1^2*(starRingEnd ℂ z1)^2*z2^2*z3^2*(starRingEnd ℂ z3)^2 - z1^2*(starRingEnd ℂ z1)^2*(starRingEnd ℂ z2)^2*z3^2*(starRingEnd ℂ z3)^2 - z1^2*z2^2*(starRingEnd ℂ z2)^2*z3^2*(starRingEnd ℂ z3)^2 - (starRingEnd ℂ z1)^2*z2^2*(starRingEnd ℂ z2)^2*z3^2*(starRingEnd ℂ z3)^2)
      + ((p1 - q1) * (((Real.sqrt 2 : ℝ) : ℂ) * Complex.sin (L * z1 / 2)) + (r1 - s1) * (starRingEnd ℂ (((Real.sqrt 2 : ℝ) : ℂ) * Complex.sin (L * z1 / 2))))
        * (-z1^2*(starRingEnd ℂ z1)^2*z2^2*(starRingEnd ℂ z2)^2*z3^2*(starRingEnd ℂ z3)^2 + z2^2*(starRingEnd ℂ z2)^2*z3^2*(starRingEnd ℂ z3)^2*z4^2*(starRingEnd ℂ z4)^2)
      + ((p2 - q2) * (((Real.sqrt 2 : ℝ) : ℂ) * Complex.sin (L * z2 / 2)) + (r2 - s2) * (starRingEnd ℂ (((Real.sqrt 2 : ℝ) : ℂ) * Complex.sin (L * z2 / 2))))
        * (-z1^2*(starRingEnd ℂ z1)^2*z2^2*(starRingEnd ℂ z2)^2*z3^2*(starRingEnd ℂ z3)^2 + z1^2*(starRingEnd ℂ z1)^2*z3^2*(starRingEnd ℂ z3)^2*z4^2*(starRingEnd ℂ z4)^2)
      + ((p3 - q3) * (((Real.sqrt 2 : ℝ) : ℂ) * Complex.sin (L * z3 / 2)) + (r3 - s3) * (starRingEnd ℂ (((Real.sqrt 2 : ℝ) : ℂ) * Complex.sin (L * z3 / 2))))
        * (-z1^2*(starRingEnd ℂ z1)^2*z2^2*(starRingEnd ℂ z2)^2*z3^2*(starRingEnd ℂ z3)^2 + z1^2*(starRingEnd ℂ z1)^2*z2^2*(starRingEnd ℂ z2)^2*z4^2*(starRingEnd ℂ z4)^2) = 0 := by
  have hE15 := crossK4LeadingCoeff N hN L hL z1 z2 z3 z4 hz1_im hz2_im hz3_im hz4_im p1 q1 r1 s1 p2 q2 r2 s2 p3 q3 r3 s3 p4 q4 r4 s4 hsum
  have hE3 := crossK4ThirteenthCoeff N hN L hL z1 z2 z3 z4 hz1_im hz2_im hz3_im hz4_im p1 q1 r1 s1 p2 q2 r2 s2 p3 q3 r3 s3 p4 q4 r4 s4 hsum
  linear_combination hE3 - (z1^2*(starRingEnd ℂ z1)^2*z2^2*(starRingEnd ℂ z2)^2*z3^2*(starRingEnd ℂ z3)^2) * hE15

#print axioms crossK4ReducedK3
