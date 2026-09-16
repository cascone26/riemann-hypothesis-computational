import CrossK5Key
import NuInjectivePolyVanish

open Polynomial

set_option maxRecDepth 4000
set_option maxHeartbeats 4000000

/-!
# K=5 step 2: the degree-19 numerator is the ZERO polynomial

S2592: direct 5-block extension of CrossK4PolyVanish, same proven pattern (generalize the repeated
quartic-in-ν products before the final linear_combination, avoiding the OOM class of bug entirely). -/

theorem crossK5PolyVanish (N : ℕ) (hN : 10 ≤ N) (L : ℝ) (hL : L ≠ 0) (z1 z2 z3 z4 z5 : ℂ)
    (hz1_im : z1.im ≠ 0) (hz2_im : z2.im ≠ 0) (hz3_im : z3.im ≠ 0) (hz4_im : z4.im ≠ 0) (hz5_im : z5.im ≠ 0)
    (p1 q1 r1 s1 p2 q2 r2 s2 p3 q3 r3 s3 p4 q4 r4 s4 p5 q5 r5 s5 : ℂ)
    (hsum : p1 • vPlusLC N L z1 + q1 • vMinusLC N L z1
        + r1 • vPlusLC N L (starRingEnd ℂ z1) + s1 • vMinusLC N L (starRingEnd ℂ z1)
        + p2 • vPlusLC N L z2 + q2 • vMinusLC N L z2
        + r2 • vPlusLC N L (starRingEnd ℂ z2) + s2 • vMinusLC N L (starRingEnd ℂ z2)
        + p3 • vPlusLC N L z3 + q3 • vMinusLC N L z3
        + r3 • vPlusLC N L (starRingEnd ℂ z3) + s3 • vMinusLC N L (starRingEnd ℂ z3)
        + p4 • vPlusLC N L z4 + q4 • vMinusLC N L z4
        + r4 • vPlusLC N L (starRingEnd ℂ z4) + s4 • vMinusLC N L (starRingEnd ℂ z4)
        + p5 • vPlusLC N L z5 + q5 • vMinusLC N L z5
        + r5 • vPlusLC N L (starRingEnd ℂ z5) + s5 • vMinusLC N L (starRingEnd ℂ z5) = 0) :
    ((C ((p1 + q1) * (((Real.sqrt 2 : ℝ) : ℂ) * Complex.sin (L * z1 / 2)) * (z1 * (starRingEnd ℂ z1) ^ 2)
            + (r1 + s1) * (starRingEnd ℂ (((Real.sqrt 2 : ℝ) : ℂ) * Complex.sin (L * z1 / 2))) * (z1 ^ 2 * (starRingEnd ℂ z1)))
        + C ((q1 - p1) * (((Real.sqrt 2 : ℝ) : ℂ) * Complex.sin (L * z1 / 2)) * (starRingEnd ℂ z1) ^ 2
          + (s1 - r1) * (starRingEnd ℂ (((Real.sqrt 2 : ℝ) : ℂ) * Complex.sin (L * z1 / 2))) * z1 ^ 2) * X
        + C (-(p1 + q1) * (((Real.sqrt 2 : ℝ) : ℂ) * Complex.sin (L * z1 / 2)) * z1
          - (r1 + s1) * (starRingEnd ℂ (((Real.sqrt 2 : ℝ) : ℂ) * Complex.sin (L * z1 / 2))) * (starRingEnd ℂ z1)) * X ^ 2
        + C ((p1 - q1) * (((Real.sqrt 2 : ℝ) : ℂ) * Complex.sin (L * z1 / 2))
          + (r1 - s1) * (starRingEnd ℂ (((Real.sqrt 2 : ℝ) : ℂ) * Complex.sin (L * z1 / 2)))) * X ^ 3)
        * ((C z2 + X) * (C z2 - X) * (C (starRingEnd ℂ z2) + X) * (C (starRingEnd ℂ z2) - X))
        * ((C z3 + X) * (C z3 - X) * (C (starRingEnd ℂ z3) + X) * (C (starRingEnd ℂ z3) - X))
        * ((C z4 + X) * (C z4 - X) * (C (starRingEnd ℂ z4) + X) * (C (starRingEnd ℂ z4) - X))
        * ((C z5 + X) * (C z5 - X) * (C (starRingEnd ℂ z5) + X) * (C (starRingEnd ℂ z5) - X))
      + (C ((p2 + q2) * (((Real.sqrt 2 : ℝ) : ℂ) * Complex.sin (L * z2 / 2)) * (z2 * (starRingEnd ℂ z2) ^ 2)
            + (r2 + s2) * (starRingEnd ℂ (((Real.sqrt 2 : ℝ) : ℂ) * Complex.sin (L * z2 / 2))) * (z2 ^ 2 * (starRingEnd ℂ z2)))
        + C ((q2 - p2) * (((Real.sqrt 2 : ℝ) : ℂ) * Complex.sin (L * z2 / 2)) * (starRingEnd ℂ z2) ^ 2
          + (s2 - r2) * (starRingEnd ℂ (((Real.sqrt 2 : ℝ) : ℂ) * Complex.sin (L * z2 / 2))) * z2 ^ 2) * X
        + C (-(p2 + q2) * (((Real.sqrt 2 : ℝ) : ℂ) * Complex.sin (L * z2 / 2)) * z2
          - (r2 + s2) * (starRingEnd ℂ (((Real.sqrt 2 : ℝ) : ℂ) * Complex.sin (L * z2 / 2))) * (starRingEnd ℂ z2)) * X ^ 2
        + C ((p2 - q2) * (((Real.sqrt 2 : ℝ) : ℂ) * Complex.sin (L * z2 / 2))
          + (r2 - s2) * (starRingEnd ℂ (((Real.sqrt 2 : ℝ) : ℂ) * Complex.sin (L * z2 / 2)))) * X ^ 3)
        * ((C z1 + X) * (C z1 - X) * (C (starRingEnd ℂ z1) + X) * (C (starRingEnd ℂ z1) - X))
        * ((C z3 + X) * (C z3 - X) * (C (starRingEnd ℂ z3) + X) * (C (starRingEnd ℂ z3) - X))
        * ((C z4 + X) * (C z4 - X) * (C (starRingEnd ℂ z4) + X) * (C (starRingEnd ℂ z4) - X))
        * ((C z5 + X) * (C z5 - X) * (C (starRingEnd ℂ z5) + X) * (C (starRingEnd ℂ z5) - X))
      + (C ((p3 + q3) * (((Real.sqrt 2 : ℝ) : ℂ) * Complex.sin (L * z3 / 2)) * (z3 * (starRingEnd ℂ z3) ^ 2)
            + (r3 + s3) * (starRingEnd ℂ (((Real.sqrt 2 : ℝ) : ℂ) * Complex.sin (L * z3 / 2))) * (z3 ^ 2 * (starRingEnd ℂ z3)))
        + C ((q3 - p3) * (((Real.sqrt 2 : ℝ) : ℂ) * Complex.sin (L * z3 / 2)) * (starRingEnd ℂ z3) ^ 2
          + (s3 - r3) * (starRingEnd ℂ (((Real.sqrt 2 : ℝ) : ℂ) * Complex.sin (L * z3 / 2))) * z3 ^ 2) * X
        + C (-(p3 + q3) * (((Real.sqrt 2 : ℝ) : ℂ) * Complex.sin (L * z3 / 2)) * z3
          - (r3 + s3) * (starRingEnd ℂ (((Real.sqrt 2 : ℝ) : ℂ) * Complex.sin (L * z3 / 2))) * (starRingEnd ℂ z3)) * X ^ 2
        + C ((p3 - q3) * (((Real.sqrt 2 : ℝ) : ℂ) * Complex.sin (L * z3 / 2))
          + (r3 - s3) * (starRingEnd ℂ (((Real.sqrt 2 : ℝ) : ℂ) * Complex.sin (L * z3 / 2)))) * X ^ 3)
        * ((C z1 + X) * (C z1 - X) * (C (starRingEnd ℂ z1) + X) * (C (starRingEnd ℂ z1) - X))
        * ((C z2 + X) * (C z2 - X) * (C (starRingEnd ℂ z2) + X) * (C (starRingEnd ℂ z2) - X))
        * ((C z4 + X) * (C z4 - X) * (C (starRingEnd ℂ z4) + X) * (C (starRingEnd ℂ z4) - X))
        * ((C z5 + X) * (C z5 - X) * (C (starRingEnd ℂ z5) + X) * (C (starRingEnd ℂ z5) - X))
      + (C ((p4 + q4) * (((Real.sqrt 2 : ℝ) : ℂ) * Complex.sin (L * z4 / 2)) * (z4 * (starRingEnd ℂ z4) ^ 2)
            + (r4 + s4) * (starRingEnd ℂ (((Real.sqrt 2 : ℝ) : ℂ) * Complex.sin (L * z4 / 2))) * (z4 ^ 2 * (starRingEnd ℂ z4)))
        + C ((q4 - p4) * (((Real.sqrt 2 : ℝ) : ℂ) * Complex.sin (L * z4 / 2)) * (starRingEnd ℂ z4) ^ 2
          + (s4 - r4) * (starRingEnd ℂ (((Real.sqrt 2 : ℝ) : ℂ) * Complex.sin (L * z4 / 2))) * z4 ^ 2) * X
        + C (-(p4 + q4) * (((Real.sqrt 2 : ℝ) : ℂ) * Complex.sin (L * z4 / 2)) * z4
          - (r4 + s4) * (starRingEnd ℂ (((Real.sqrt 2 : ℝ) : ℂ) * Complex.sin (L * z4 / 2))) * (starRingEnd ℂ z4)) * X ^ 2
        + C ((p4 - q4) * (((Real.sqrt 2 : ℝ) : ℂ) * Complex.sin (L * z4 / 2))
          + (r4 - s4) * (starRingEnd ℂ (((Real.sqrt 2 : ℝ) : ℂ) * Complex.sin (L * z4 / 2)))) * X ^ 3)
        * ((C z1 + X) * (C z1 - X) * (C (starRingEnd ℂ z1) + X) * (C (starRingEnd ℂ z1) - X))
        * ((C z2 + X) * (C z2 - X) * (C (starRingEnd ℂ z2) + X) * (C (starRingEnd ℂ z2) - X))
        * ((C z3 + X) * (C z3 - X) * (C (starRingEnd ℂ z3) + X) * (C (starRingEnd ℂ z3) - X))
        * ((C z5 + X) * (C z5 - X) * (C (starRingEnd ℂ z5) + X) * (C (starRingEnd ℂ z5) - X))
      + (C ((p5 + q5) * (((Real.sqrt 2 : ℝ) : ℂ) * Complex.sin (L * z5 / 2)) * (z5 * (starRingEnd ℂ z5) ^ 2)
            + (r5 + s5) * (starRingEnd ℂ (((Real.sqrt 2 : ℝ) : ℂ) * Complex.sin (L * z5 / 2))) * (z5 ^ 2 * (starRingEnd ℂ z5)))
        + C ((q5 - p5) * (((Real.sqrt 2 : ℝ) : ℂ) * Complex.sin (L * z5 / 2)) * (starRingEnd ℂ z5) ^ 2
          + (s5 - r5) * (starRingEnd ℂ (((Real.sqrt 2 : ℝ) : ℂ) * Complex.sin (L * z5 / 2))) * z5 ^ 2) * X
        + C (-(p5 + q5) * (((Real.sqrt 2 : ℝ) : ℂ) * Complex.sin (L * z5 / 2)) * z5
          - (r5 + s5) * (starRingEnd ℂ (((Real.sqrt 2 : ℝ) : ℂ) * Complex.sin (L * z5 / 2))) * (starRingEnd ℂ z5)) * X ^ 2
        + C ((p5 - q5) * (((Real.sqrt 2 : ℝ) : ℂ) * Complex.sin (L * z5 / 2))
          + (r5 - s5) * (starRingEnd ℂ (((Real.sqrt 2 : ℝ) : ℂ) * Complex.sin (L * z5 / 2)))) * X ^ 3)
        * ((C z1 + X) * (C z1 - X) * (C (starRingEnd ℂ z1) + X) * (C (starRingEnd ℂ z1) - X))
        * ((C z2 + X) * (C z2 - X) * (C (starRingEnd ℂ z2) + X) * (C (starRingEnd ℂ z2) - X))
        * ((C z3 + X) * (C z3 - X) * (C (starRingEnd ℂ z3) + X) * (C (starRingEnd ℂ z3) - X))
        * ((C z4 + X) * (C z4 - X) * (C (starRingEnd ℂ z4) + X) * (C (starRingEnd ℂ z4) - X))
      : Polynomial ℂ) = 0 := by
  have hkey := crossK5Key N L hL z1 z2 z3 z4 z5 hz1_im hz2_im hz3_im hz4_im hz5_im p1 q1 r1 s1 p2 q2 r2 s2 p3 q3 r3 s3 p4 q4 r4 s4 p5 q5 r5 s5 hsum
  apply poly_vanish_at_nu_nodes N L hL _ ?_ ?_
  · refine lt_of_le_of_lt (?_ : _ ≤ 19) (by omega)
    compute_degree
  · intro n
    have hkn := hkey n
    dsimp only at hkn
    simp only [eval_add, eval_mul, eval_sub, eval_C, eval_X, eval_pow]
    push_cast at hkn ⊢
    revert hkn
    generalize hQ1 : (z1 + (2 * Real.pi * (((n : ℕ) : ℝ) - (N : ℝ)) / L : ℂ)) * (z1 - (2 * Real.pi * (((n : ℕ) : ℝ) - (N : ℝ)) / L : ℂ)) * (starRingEnd ℂ z1 + (2 * Real.pi * (((n : ℕ) : ℝ) - (N : ℝ)) / L : ℂ))
        * (starRingEnd ℂ z1 - (2 * Real.pi * (((n : ℕ) : ℝ) - (N : ℝ)) / L : ℂ)) = Q1
    generalize hQ2 : (z2 + (2 * Real.pi * (((n : ℕ) : ℝ) - (N : ℝ)) / L : ℂ)) * (z2 - (2 * Real.pi * (((n : ℕ) : ℝ) - (N : ℝ)) / L : ℂ)) * (starRingEnd ℂ z2 + (2 * Real.pi * (((n : ℕ) : ℝ) - (N : ℝ)) / L : ℂ))
        * (starRingEnd ℂ z2 - (2 * Real.pi * (((n : ℕ) : ℝ) - (N : ℝ)) / L : ℂ)) = Q2
    generalize hQ3 : (z3 + (2 * Real.pi * (((n : ℕ) : ℝ) - (N : ℝ)) / L : ℂ)) * (z3 - (2 * Real.pi * (((n : ℕ) : ℝ) - (N : ℝ)) / L : ℂ)) * (starRingEnd ℂ z3 + (2 * Real.pi * (((n : ℕ) : ℝ) - (N : ℝ)) / L : ℂ))
        * (starRingEnd ℂ z3 - (2 * Real.pi * (((n : ℕ) : ℝ) - (N : ℝ)) / L : ℂ)) = Q3
    generalize hQ4 : (z4 + (2 * Real.pi * (((n : ℕ) : ℝ) - (N : ℝ)) / L : ℂ)) * (z4 - (2 * Real.pi * (((n : ℕ) : ℝ) - (N : ℝ)) / L : ℂ)) * (starRingEnd ℂ z4 + (2 * Real.pi * (((n : ℕ) : ℝ) - (N : ℝ)) / L : ℂ))
        * (starRingEnd ℂ z4 - (2 * Real.pi * (((n : ℕ) : ℝ) - (N : ℝ)) / L : ℂ)) = Q4
    generalize hQ5 : (z5 + (2 * Real.pi * (((n : ℕ) : ℝ) - (N : ℝ)) / L : ℂ)) * (z5 - (2 * Real.pi * (((n : ℕ) : ℝ) - (N : ℝ)) / L : ℂ)) * (starRingEnd ℂ z5 + (2 * Real.pi * (((n : ℕ) : ℝ) - (N : ℝ)) / L : ℂ))
        * (starRingEnd ℂ z5 - (2 * Real.pi * (((n : ℕ) : ℝ) - (N : ℝ)) / L : ℂ)) = Q5
    intro hkn
    linear_combination hkn

#print axioms crossK5PolyVanish
