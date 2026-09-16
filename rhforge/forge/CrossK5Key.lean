import GroupClearIdentity
import VPlusVMinusFullKey

set_option maxRecDepth 4000
set_option maxHeartbeats 8000000

/-!
# K=5 cross-block combined numerator identity -- extending the concrete arc to a third instance

S2591: direct 5-block extension of CrossK4Key, using the SAME proven toolkit (generalize for true
ring1 opacity + congrArg for zero-transcription per-block scaling + push_cast for cast normalization)
that closed the entire K=4 odd-group arc tonight. K=5's D-cancellation was verified via sympy first
(S2582) before writing this Lean, matching the project's verify-before-build discipline. -/

theorem crossK5Key (N : ℕ) (L : ℝ) (hL : L ≠ 0) (z1 z2 z3 z4 z5 : ℂ)
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
    ∀ n : Fin (2 * N + 1),
      let A1 := ((Real.sqrt 2 : ℝ) : ℂ) * Complex.sin (L * z1 / 2)
      let A2 := ((Real.sqrt 2 : ℝ) : ℂ) * Complex.sin (L * z2 / 2)
      let A3 := ((Real.sqrt 2 : ℝ) : ℂ) * Complex.sin (L * z3 / 2)
      let A4 := ((Real.sqrt 2 : ℝ) : ℂ) * Complex.sin (L * z4 / 2)
      let A5 := ((Real.sqrt 2 : ℝ) : ℂ) * Complex.sin (L * z5 / 2)
      let cubic1 := (((p1 + q1) * A1 * (z1 * (starRingEnd ℂ z1) ^ 2)
          + (r1 + s1) * (starRingEnd ℂ A1) * (z1 ^ 2 * (starRingEnd ℂ z1))
          + ((q1 - p1) * A1 * (starRingEnd ℂ z1) ^ 2 + (s1 - r1) * (starRingEnd ℂ A1) * z1 ^ 2) * (2 * Real.pi * (((n : ℕ) : ℝ) - (N : ℝ)) / L : ℂ)
          + (-(p1 + q1) * A1 * z1 - (r1 + s1) * (starRingEnd ℂ A1) * (starRingEnd ℂ z1)) * (2 * Real.pi * (((n : ℕ) : ℝ) - (N : ℝ)) / L : ℂ) ^ 2
          + ((p1 - q1) * A1 + (r1 - s1) * (starRingEnd ℂ A1)) * (2 * Real.pi * (((n : ℕ) : ℝ) - (N : ℝ)) / L : ℂ) ^ 3))
      let cubic2 := (((p2 + q2) * A2 * (z2 * (starRingEnd ℂ z2) ^ 2)
          + (r2 + s2) * (starRingEnd ℂ A2) * (z2 ^ 2 * (starRingEnd ℂ z2))
          + ((q2 - p2) * A2 * (starRingEnd ℂ z2) ^ 2 + (s2 - r2) * (starRingEnd ℂ A2) * z2 ^ 2) * (2 * Real.pi * (((n : ℕ) : ℝ) - (N : ℝ)) / L : ℂ)
          + (-(p2 + q2) * A2 * z2 - (r2 + s2) * (starRingEnd ℂ A2) * (starRingEnd ℂ z2)) * (2 * Real.pi * (((n : ℕ) : ℝ) - (N : ℝ)) / L : ℂ) ^ 2
          + ((p2 - q2) * A2 + (r2 - s2) * (starRingEnd ℂ A2)) * (2 * Real.pi * (((n : ℕ) : ℝ) - (N : ℝ)) / L : ℂ) ^ 3))
      let cubic3 := (((p3 + q3) * A3 * (z3 * (starRingEnd ℂ z3) ^ 2)
          + (r3 + s3) * (starRingEnd ℂ A3) * (z3 ^ 2 * (starRingEnd ℂ z3))
          + ((q3 - p3) * A3 * (starRingEnd ℂ z3) ^ 2 + (s3 - r3) * (starRingEnd ℂ A3) * z3 ^ 2) * (2 * Real.pi * (((n : ℕ) : ℝ) - (N : ℝ)) / L : ℂ)
          + (-(p3 + q3) * A3 * z3 - (r3 + s3) * (starRingEnd ℂ A3) * (starRingEnd ℂ z3)) * (2 * Real.pi * (((n : ℕ) : ℝ) - (N : ℝ)) / L : ℂ) ^ 2
          + ((p3 - q3) * A3 + (r3 - s3) * (starRingEnd ℂ A3)) * (2 * Real.pi * (((n : ℕ) : ℝ) - (N : ℝ)) / L : ℂ) ^ 3))
      let cubic4 := (((p4 + q4) * A4 * (z4 * (starRingEnd ℂ z4) ^ 2)
          + (r4 + s4) * (starRingEnd ℂ A4) * (z4 ^ 2 * (starRingEnd ℂ z4))
          + ((q4 - p4) * A4 * (starRingEnd ℂ z4) ^ 2 + (s4 - r4) * (starRingEnd ℂ A4) * z4 ^ 2) * (2 * Real.pi * (((n : ℕ) : ℝ) - (N : ℝ)) / L : ℂ)
          + (-(p4 + q4) * A4 * z4 - (r4 + s4) * (starRingEnd ℂ A4) * (starRingEnd ℂ z4)) * (2 * Real.pi * (((n : ℕ) : ℝ) - (N : ℝ)) / L : ℂ) ^ 2
          + ((p4 - q4) * A4 + (r4 - s4) * (starRingEnd ℂ A4)) * (2 * Real.pi * (((n : ℕ) : ℝ) - (N : ℝ)) / L : ℂ) ^ 3))
      let cubic5 := (((p5 + q5) * A5 * (z5 * (starRingEnd ℂ z5) ^ 2)
          + (r5 + s5) * (starRingEnd ℂ A5) * (z5 ^ 2 * (starRingEnd ℂ z5))
          + ((q5 - p5) * A5 * (starRingEnd ℂ z5) ^ 2 + (s5 - r5) * (starRingEnd ℂ A5) * z5 ^ 2) * (2 * Real.pi * (((n : ℕ) : ℝ) - (N : ℝ)) / L : ℂ)
          + (-(p5 + q5) * A5 * z5 - (r5 + s5) * (starRingEnd ℂ A5) * (starRingEnd ℂ z5)) * (2 * Real.pi * (((n : ℕ) : ℝ) - (N : ℝ)) / L : ℂ) ^ 2
          + ((p5 - q5) * A5 + (r5 - s5) * (starRingEnd ℂ A5)) * (2 * Real.pi * (((n : ℕ) : ℝ) - (N : ℝ)) / L : ℂ) ^ 3))
      let quartic1 := ((z1 + (2 * Real.pi * (((n : ℕ) : ℝ) - (N : ℝ)) / L : ℂ)) * (z1 - (2 * Real.pi * (((n : ℕ) : ℝ) - (N : ℝ)) / L : ℂ)) * (starRingEnd ℂ z1 + (2 * Real.pi * (((n : ℕ) : ℝ) - (N : ℝ)) / L : ℂ)) * (starRingEnd ℂ z1 - (2 * Real.pi * (((n : ℕ) : ℝ) - (N : ℝ)) / L : ℂ)))
      let quartic2 := ((z2 + (2 * Real.pi * (((n : ℕ) : ℝ) - (N : ℝ)) / L : ℂ)) * (z2 - (2 * Real.pi * (((n : ℕ) : ℝ) - (N : ℝ)) / L : ℂ)) * (starRingEnd ℂ z2 + (2 * Real.pi * (((n : ℕ) : ℝ) - (N : ℝ)) / L : ℂ)) * (starRingEnd ℂ z2 - (2 * Real.pi * (((n : ℕ) : ℝ) - (N : ℝ)) / L : ℂ)))
      let quartic3 := ((z3 + (2 * Real.pi * (((n : ℕ) : ℝ) - (N : ℝ)) / L : ℂ)) * (z3 - (2 * Real.pi * (((n : ℕ) : ℝ) - (N : ℝ)) / L : ℂ)) * (starRingEnd ℂ z3 + (2 * Real.pi * (((n : ℕ) : ℝ) - (N : ℝ)) / L : ℂ)) * (starRingEnd ℂ z3 - (2 * Real.pi * (((n : ℕ) : ℝ) - (N : ℝ)) / L : ℂ)))
      let quartic4 := ((z4 + (2 * Real.pi * (((n : ℕ) : ℝ) - (N : ℝ)) / L : ℂ)) * (z4 - (2 * Real.pi * (((n : ℕ) : ℝ) - (N : ℝ)) / L : ℂ)) * (starRingEnd ℂ z4 + (2 * Real.pi * (((n : ℕ) : ℝ) - (N : ℝ)) / L : ℂ)) * (starRingEnd ℂ z4 - (2 * Real.pi * (((n : ℕ) : ℝ) - (N : ℝ)) / L : ℂ)))
      let quartic5 := ((z5 + (2 * Real.pi * (((n : ℕ) : ℝ) - (N : ℝ)) / L : ℂ)) * (z5 - (2 * Real.pi * (((n : ℕ) : ℝ) - (N : ℝ)) / L : ℂ)) * (starRingEnd ℂ z5 + (2 * Real.pi * (((n : ℕ) : ℝ) - (N : ℝ)) / L : ℂ)) * (starRingEnd ℂ z5 - (2 * Real.pi * (((n : ℕ) : ℝ) - (N : ℝ)) / L : ℂ)))
      cubic1 * quartic2 * quartic3 * quartic4 * quartic5
        + cubic2 * quartic1 * quartic3 * quartic4 * quartic5
        + cubic3 * quartic1 * quartic2 * quartic4 * quartic5
        + cubic4 * quartic1 * quartic2 * quartic3 * quartic5
        + cubic5 * quartic1 * quartic2 * quartic3 * quartic4
        = 0 := by
  intro n
  have hLC : (L : ℂ) ≠ 0 := by exact_mod_cast hL
  have hpt := congrFun hsum n
  simp only [Pi.add_apply, Pi.smul_apply, smul_eq_mul, Pi.zero_apply] at hpt
  simp only [vPlusLC_conj, vMinusLC_conj] at hpt
  unfold vPlusLC vMinusLC at hpt
  set A1 := ((Real.sqrt 2 : ℝ) : ℂ) * Complex.sin (L * z1 / 2) with hA1def
  set A2 := ((Real.sqrt 2 : ℝ) : ℂ) * Complex.sin (L * z2 / 2) with hA2def
  set A3 := ((Real.sqrt 2 : ℝ) : ℂ) * Complex.sin (L * z3 / 2) with hA3def
  set A4 := ((Real.sqrt 2 : ℝ) : ℂ) * Complex.sin (L * z4 / 2) with hA4def
  set A5 := ((Real.sqrt 2 : ℝ) : ℂ) * Complex.sin (L * z5 / 2) with hA5def
  simp only [map_add, map_sub, map_mul, map_div₀, Complex.conj_ofReal, map_ofNat] at hpt
  push_cast at hpt
  have h1 := groupClearIdentity L hL z1 hz1_im
    (2 * Real.pi * (((n : ℕ) : ℝ) - (N : ℝ)) / L) A1 (starRingEnd ℂ A1) p1 q1 r1 s1
  have h2 := groupClearIdentity L hL z2 hz2_im
    (2 * Real.pi * (((n : ℕ) : ℝ) - (N : ℝ)) / L) A2 (starRingEnd ℂ A2) p2 q2 r2 s2
  have h3 := groupClearIdentity L hL z3 hz3_im
    (2 * Real.pi * (((n : ℕ) : ℝ) - (N : ℝ)) / L) A3 (starRingEnd ℂ A3) p3 q3 r3 s3
  have h4 := groupClearIdentity L hL z4 hz4_im
    (2 * Real.pi * (((n : ℕ) : ℝ) - (N : ℝ)) / L) A4 (starRingEnd ℂ A4) p4 q4 r4 s4
  have h5 := groupClearIdentity L hL z5 hz5_im
    (2 * Real.pi * (((n : ℕ) : ℝ) - (N : ℝ)) / L) A5 (starRingEnd ℂ A5) p5 q5 r5 s5
  push_cast at h1 h2 h3 h4 h5
  have hL3ne : (L : ℂ) ^ 3 ≠ 0 := pow_ne_zero 3 hLC
  dsimp only
  revert h1
  generalize hQ1 : (z1 + (2 * Real.pi * (((n : ℕ) : ℝ) - (N : ℝ)) / L : ℂ)) * (z1 - (2 * Real.pi * (((n : ℕ) : ℝ) - (N : ℝ)) / L : ℂ)) * (starRingEnd ℂ z1 + (2 * Real.pi * (((n : ℕ) : ℝ) - (N : ℝ)) / L : ℂ))
      * (starRingEnd ℂ z1 - (2 * Real.pi * (((n : ℕ) : ℝ) - (N : ℝ)) / L : ℂ)) = Q1
  intro h1
  revert h2
  generalize hQ2 : (z2 + (2 * Real.pi * (((n : ℕ) : ℝ) - (N : ℝ)) / L : ℂ)) * (z2 - (2 * Real.pi * (((n : ℕ) : ℝ) - (N : ℝ)) / L : ℂ)) * (starRingEnd ℂ z2 + (2 * Real.pi * (((n : ℕ) : ℝ) - (N : ℝ)) / L : ℂ))
      * (starRingEnd ℂ z2 - (2 * Real.pi * (((n : ℕ) : ℝ) - (N : ℝ)) / L : ℂ)) = Q2
  intro h2
  revert h3
  generalize hQ3 : (z3 + (2 * Real.pi * (((n : ℕ) : ℝ) - (N : ℝ)) / L : ℂ)) * (z3 - (2 * Real.pi * (((n : ℕ) : ℝ) - (N : ℝ)) / L : ℂ)) * (starRingEnd ℂ z3 + (2 * Real.pi * (((n : ℕ) : ℝ) - (N : ℝ)) / L : ℂ))
      * (starRingEnd ℂ z3 - (2 * Real.pi * (((n : ℕ) : ℝ) - (N : ℝ)) / L : ℂ)) = Q3
  intro h3
  revert h4
  generalize hQ4 : (z4 + (2 * Real.pi * (((n : ℕ) : ℝ) - (N : ℝ)) / L : ℂ)) * (z4 - (2 * Real.pi * (((n : ℕ) : ℝ) - (N : ℝ)) / L : ℂ)) * (starRingEnd ℂ z4 + (2 * Real.pi * (((n : ℕ) : ℝ) - (N : ℝ)) / L : ℂ))
      * (starRingEnd ℂ z4 - (2 * Real.pi * (((n : ℕ) : ℝ) - (N : ℝ)) / L : ℂ)) = Q4
  intro h4
  revert h5
  generalize hQ5 : (z5 + (2 * Real.pi * (((n : ℕ) : ℝ) - (N : ℝ)) / L : ℂ)) * (z5 - (2 * Real.pi * (((n : ℕ) : ℝ) - (N : ℝ)) / L : ℂ)) * (starRingEnd ℂ z5 + (2 * Real.pi * (((n : ℕ) : ℝ) - (N : ℝ)) / L : ℂ))
      * (starRingEnd ℂ z5 - (2 * Real.pi * (((n : ℕ) : ℝ) - (N : ℝ)) / L : ℂ)) = Q5
  intro h5
  have e1 := congrArg (fun x => (-(Q2*Q3*Q4*Q5*(L:ℂ)^16)) * x) h1
  have e2 := congrArg (fun x => (-(Q1*Q3*Q4*Q5*(L:ℂ)^16)) * x) h2
  have e3 := congrArg (fun x => (-(Q1*Q2*Q4*Q5*(L:ℂ)^16)) * x) h3
  have e4 := congrArg (fun x => (-(Q1*Q2*Q3*Q5*(L:ℂ)^16)) * x) h4
  have e5 := congrArg (fun x => (-(Q1*Q2*Q3*Q4*(L:ℂ)^16)) * x) h5
  have hfinal : (L:ℂ)^19 * (
      (((p1 + q1) * A1 * (z1 * (starRingEnd ℂ z1) ^ 2)
          + (r1 + s1) * (starRingEnd ℂ A1) * (z1 ^ 2 * (starRingEnd ℂ z1))
          + ((q1 - p1) * A1 * (starRingEnd ℂ z1) ^ 2 + (s1 - r1) * (starRingEnd ℂ A1) * z1 ^ 2) * (2 * Real.pi * (((n : ℕ) : ℝ) - (N : ℝ)) / L : ℂ)
          + (-(p1 + q1) * A1 * z1 - (r1 + s1) * (starRingEnd ℂ A1) * (starRingEnd ℂ z1)) * (2 * Real.pi * (((n : ℕ) : ℝ) - (N : ℝ)) / L : ℂ) ^ 2
          + ((p1 - q1) * A1 + (r1 - s1) * (starRingEnd ℂ A1)) * (2 * Real.pi * (((n : ℕ) : ℝ) - (N : ℝ)) / L : ℂ) ^ 3)) * Q2*Q3*Q4*Q5
      + (((p2 + q2) * A2 * (z2 * (starRingEnd ℂ z2) ^ 2)
          + (r2 + s2) * (starRingEnd ℂ A2) * (z2 ^ 2 * (starRingEnd ℂ z2))
          + ((q2 - p2) * A2 * (starRingEnd ℂ z2) ^ 2 + (s2 - r2) * (starRingEnd ℂ A2) * z2 ^ 2) * (2 * Real.pi * (((n : ℕ) : ℝ) - (N : ℝ)) / L : ℂ)
          + (-(p2 + q2) * A2 * z2 - (r2 + s2) * (starRingEnd ℂ A2) * (starRingEnd ℂ z2)) * (2 * Real.pi * (((n : ℕ) : ℝ) - (N : ℝ)) / L : ℂ) ^ 2
          + ((p2 - q2) * A2 + (r2 - s2) * (starRingEnd ℂ A2)) * (2 * Real.pi * (((n : ℕ) : ℝ) - (N : ℝ)) / L : ℂ) ^ 3)) * Q1*Q3*Q4*Q5
      + (((p3 + q3) * A3 * (z3 * (starRingEnd ℂ z3) ^ 2)
          + (r3 + s3) * (starRingEnd ℂ A3) * (z3 ^ 2 * (starRingEnd ℂ z3))
          + ((q3 - p3) * A3 * (starRingEnd ℂ z3) ^ 2 + (s3 - r3) * (starRingEnd ℂ A3) * z3 ^ 2) * (2 * Real.pi * (((n : ℕ) : ℝ) - (N : ℝ)) / L : ℂ)
          + (-(p3 + q3) * A3 * z3 - (r3 + s3) * (starRingEnd ℂ A3) * (starRingEnd ℂ z3)) * (2 * Real.pi * (((n : ℕ) : ℝ) - (N : ℝ)) / L : ℂ) ^ 2
          + ((p3 - q3) * A3 + (r3 - s3) * (starRingEnd ℂ A3)) * (2 * Real.pi * (((n : ℕ) : ℝ) - (N : ℝ)) / L : ℂ) ^ 3)) * Q1*Q2*Q4*Q5
      + (((p4 + q4) * A4 * (z4 * (starRingEnd ℂ z4) ^ 2)
          + (r4 + s4) * (starRingEnd ℂ A4) * (z4 ^ 2 * (starRingEnd ℂ z4))
          + ((q4 - p4) * A4 * (starRingEnd ℂ z4) ^ 2 + (s4 - r4) * (starRingEnd ℂ A4) * z4 ^ 2) * (2 * Real.pi * (((n : ℕ) : ℝ) - (N : ℝ)) / L : ℂ)
          + (-(p4 + q4) * A4 * z4 - (r4 + s4) * (starRingEnd ℂ A4) * (starRingEnd ℂ z4)) * (2 * Real.pi * (((n : ℕ) : ℝ) - (N : ℝ)) / L : ℂ) ^ 2
          + ((p4 - q4) * A4 + (r4 - s4) * (starRingEnd ℂ A4)) * (2 * Real.pi * (((n : ℕ) : ℝ) - (N : ℝ)) / L : ℂ) ^ 3)) * Q1*Q2*Q3*Q5
      + (((p5 + q5) * A5 * (z5 * (starRingEnd ℂ z5) ^ 2)
          + (r5 + s5) * (starRingEnd ℂ A5) * (z5 ^ 2 * (starRingEnd ℂ z5))
          + ((q5 - p5) * A5 * (starRingEnd ℂ z5) ^ 2 + (s5 - r5) * (starRingEnd ℂ A5) * z5 ^ 2) * (2 * Real.pi * (((n : ℕ) : ℝ) - (N : ℝ)) / L : ℂ)
          + (-(p5 + q5) * A5 * z5 - (r5 + s5) * (starRingEnd ℂ A5) * (starRingEnd ℂ z5)) * (2 * Real.pi * (((n : ℕ) : ℝ) - (N : ℝ)) / L : ℂ) ^ 2
          + ((p5 - q5) * A5 + (r5 - s5) * (starRingEnd ℂ A5)) * (2 * Real.pi * (((n : ℕ) : ℝ) - (N : ℝ)) / L : ℂ) ^ 3)) * Q1*Q2*Q3*Q4
      ) = 0 := by
    linear_combination (norm := (push_cast; ring1)) e1 + e2 + e3 + e4 + e5 +
      ((L:ℂ)^20 * (Q1*Q2*Q3*Q4*Q5)) * hpt
  have hgoal := (mul_eq_zero.mp hfinal).resolve_left (pow_ne_zero 19 hLC)
  linear_combination hgoal

#print axioms crossK5Key
