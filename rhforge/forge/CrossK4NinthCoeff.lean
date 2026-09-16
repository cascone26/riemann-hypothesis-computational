import CrossK4PolyVanish
import QuarticTripleProductExpand

open Polynomial

set_option maxRecDepth 4000
set_option maxHeartbeats 8000000

/-!
# K=4 coefficient extraction: the degree-7 coefficient vanishes (NinthCoeff)

Direct 4-block extension of the K=3 `crossK3NinthCoeff` pattern (S2467-2504 series), using
`quarticTripleProductExpand` (S2517) in place of `quarticProductExpand`. Goal formula sympy-verified
BEFORE writing this Lean (verify_k4_coeffs.py, deg7 claim, diff=0). Same set-cubic/set-quartic/
triple-product-expand/targeted-coeff_mul skeleton as `crossK4LeadingCoeff` (S2519), varying only the
target degree and the antidiagonal-pair combination extracted from it.

Search-before-build: grepped `crossK4.*NinthCoeff`/`crossK4.*Coeff7` -- fresh; NOT YET COMPILED (drafted
while Olympus connectivity was flapping, S2520) -- REPORT-tier until #print axioms confirms EXIT0. -/

theorem crossK4NinthCoeff (N : ℕ) (hN : 8 ≤ N) (L : ℝ) (hL : L ≠ 0) (z1 z2 z3 z4 : ℂ)
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
    (-((q1 - p1) * (((Real.sqrt 2 : ℝ) : ℂ) * Complex.sin (L * z1 / 2)) * (starRingEnd ℂ z1) ^ 2
          + (s1 - r1) * (starRingEnd ℂ (((Real.sqrt 2 : ℝ) : ℂ) * Complex.sin (L * z1 / 2))) * z1 ^ 2) * ((z2 ^ 2 * (starRingEnd ℂ z2) ^ 2) * (z3 ^ 2 + (starRingEnd ℂ z3) ^ 2) + (z2 ^ 2 * (starRingEnd ℂ z2) ^ 2) * (z4 ^ 2 + (starRingEnd ℂ z4) ^ 2) + (z3 ^ 2 * (starRingEnd ℂ z3) ^ 2) * (z2 ^ 2 + (starRingEnd ℂ z2) ^ 2)
          + (z3 ^ 2 * (starRingEnd ℂ z3) ^ 2) * (z4 ^ 2 + (starRingEnd ℂ z4) ^ 2) + (z4 ^ 2 * (starRingEnd ℂ z4) ^ 2) * (z2 ^ 2 + (starRingEnd ℂ z2) ^ 2) + (z4 ^ 2 * (starRingEnd ℂ z4) ^ 2) * (z3 ^ 2 + (starRingEnd ℂ z3) ^ 2)
          + (z2 ^ 2 + (starRingEnd ℂ z2) ^ 2) * (z3 ^ 2 + (starRingEnd ℂ z3) ^ 2) * (z4 ^ 2 + (starRingEnd ℂ z4) ^ 2))
        + ((p1 - q1) * (((Real.sqrt 2 : ℝ) : ℂ) * Complex.sin (L * z1 / 2)) + (r1 - s1) * (starRingEnd ℂ (((Real.sqrt 2 : ℝ) : ℂ) * Complex.sin (L * z1 / 2)))) * ((z2 ^ 2 * (starRingEnd ℂ z2) ^ 2) * (z3 ^ 2 * (starRingEnd ℂ z3) ^ 2) + (z2 ^ 2 * (starRingEnd ℂ z2) ^ 2) * (z4 ^ 2 * (starRingEnd ℂ z4) ^ 2) + (z3 ^ 2 * (starRingEnd ℂ z3) ^ 2) * (z4 ^ 2 * (starRingEnd ℂ z4) ^ 2)
          + (z2 ^ 2 * (starRingEnd ℂ z2) ^ 2) * (z3 ^ 2 + (starRingEnd ℂ z3) ^ 2) * (z4 ^ 2 + (starRingEnd ℂ z4) ^ 2) + (z3 ^ 2 * (starRingEnd ℂ z3) ^ 2) * (z2 ^ 2 + (starRingEnd ℂ z2) ^ 2) * (z4 ^ 2 + (starRingEnd ℂ z4) ^ 2)
          + (z4 ^ 2 * (starRingEnd ℂ z4) ^ 2) * (z2 ^ 2 + (starRingEnd ℂ z2) ^ 2) * (z3 ^ 2 + (starRingEnd ℂ z3) ^ 2)))
      + (-((q2 - p2) * (((Real.sqrt 2 : ℝ) : ℂ) * Complex.sin (L * z2 / 2)) * (starRingEnd ℂ z2) ^ 2
          + (s2 - r2) * (starRingEnd ℂ (((Real.sqrt 2 : ℝ) : ℂ) * Complex.sin (L * z2 / 2))) * z2 ^ 2) * ((z1 ^ 2 * (starRingEnd ℂ z1) ^ 2) * (z3 ^ 2 + (starRingEnd ℂ z3) ^ 2) + (z1 ^ 2 * (starRingEnd ℂ z1) ^ 2) * (z4 ^ 2 + (starRingEnd ℂ z4) ^ 2) + (z3 ^ 2 * (starRingEnd ℂ z3) ^ 2) * (z1 ^ 2 + (starRingEnd ℂ z1) ^ 2)
          + (z3 ^ 2 * (starRingEnd ℂ z3) ^ 2) * (z4 ^ 2 + (starRingEnd ℂ z4) ^ 2) + (z4 ^ 2 * (starRingEnd ℂ z4) ^ 2) * (z1 ^ 2 + (starRingEnd ℂ z1) ^ 2) + (z4 ^ 2 * (starRingEnd ℂ z4) ^ 2) * (z3 ^ 2 + (starRingEnd ℂ z3) ^ 2)
          + (z1 ^ 2 + (starRingEnd ℂ z1) ^ 2) * (z3 ^ 2 + (starRingEnd ℂ z3) ^ 2) * (z4 ^ 2 + (starRingEnd ℂ z4) ^ 2))
        + ((p2 - q2) * (((Real.sqrt 2 : ℝ) : ℂ) * Complex.sin (L * z2 / 2)) + (r2 - s2) * (starRingEnd ℂ (((Real.sqrt 2 : ℝ) : ℂ) * Complex.sin (L * z2 / 2)))) * ((z1 ^ 2 * (starRingEnd ℂ z1) ^ 2) * (z3 ^ 2 * (starRingEnd ℂ z3) ^ 2) + (z1 ^ 2 * (starRingEnd ℂ z1) ^ 2) * (z4 ^ 2 * (starRingEnd ℂ z4) ^ 2) + (z3 ^ 2 * (starRingEnd ℂ z3) ^ 2) * (z4 ^ 2 * (starRingEnd ℂ z4) ^ 2)
          + (z1 ^ 2 * (starRingEnd ℂ z1) ^ 2) * (z3 ^ 2 + (starRingEnd ℂ z3) ^ 2) * (z4 ^ 2 + (starRingEnd ℂ z4) ^ 2) + (z3 ^ 2 * (starRingEnd ℂ z3) ^ 2) * (z1 ^ 2 + (starRingEnd ℂ z1) ^ 2) * (z4 ^ 2 + (starRingEnd ℂ z4) ^ 2)
          + (z4 ^ 2 * (starRingEnd ℂ z4) ^ 2) * (z1 ^ 2 + (starRingEnd ℂ z1) ^ 2) * (z3 ^ 2 + (starRingEnd ℂ z3) ^ 2)))
      + (-((q3 - p3) * (((Real.sqrt 2 : ℝ) : ℂ) * Complex.sin (L * z3 / 2)) * (starRingEnd ℂ z3) ^ 2
          + (s3 - r3) * (starRingEnd ℂ (((Real.sqrt 2 : ℝ) : ℂ) * Complex.sin (L * z3 / 2))) * z3 ^ 2) * ((z1 ^ 2 * (starRingEnd ℂ z1) ^ 2) * (z2 ^ 2 + (starRingEnd ℂ z2) ^ 2) + (z1 ^ 2 * (starRingEnd ℂ z1) ^ 2) * (z4 ^ 2 + (starRingEnd ℂ z4) ^ 2) + (z2 ^ 2 * (starRingEnd ℂ z2) ^ 2) * (z1 ^ 2 + (starRingEnd ℂ z1) ^ 2)
          + (z2 ^ 2 * (starRingEnd ℂ z2) ^ 2) * (z4 ^ 2 + (starRingEnd ℂ z4) ^ 2) + (z4 ^ 2 * (starRingEnd ℂ z4) ^ 2) * (z1 ^ 2 + (starRingEnd ℂ z1) ^ 2) + (z4 ^ 2 * (starRingEnd ℂ z4) ^ 2) * (z2 ^ 2 + (starRingEnd ℂ z2) ^ 2)
          + (z1 ^ 2 + (starRingEnd ℂ z1) ^ 2) * (z2 ^ 2 + (starRingEnd ℂ z2) ^ 2) * (z4 ^ 2 + (starRingEnd ℂ z4) ^ 2))
        + ((p3 - q3) * (((Real.sqrt 2 : ℝ) : ℂ) * Complex.sin (L * z3 / 2)) + (r3 - s3) * (starRingEnd ℂ (((Real.sqrt 2 : ℝ) : ℂ) * Complex.sin (L * z3 / 2)))) * ((z1 ^ 2 * (starRingEnd ℂ z1) ^ 2) * (z2 ^ 2 * (starRingEnd ℂ z2) ^ 2) + (z1 ^ 2 * (starRingEnd ℂ z1) ^ 2) * (z4 ^ 2 * (starRingEnd ℂ z4) ^ 2) + (z2 ^ 2 * (starRingEnd ℂ z2) ^ 2) * (z4 ^ 2 * (starRingEnd ℂ z4) ^ 2)
          + (z1 ^ 2 * (starRingEnd ℂ z1) ^ 2) * (z2 ^ 2 + (starRingEnd ℂ z2) ^ 2) * (z4 ^ 2 + (starRingEnd ℂ z4) ^ 2) + (z2 ^ 2 * (starRingEnd ℂ z2) ^ 2) * (z1 ^ 2 + (starRingEnd ℂ z1) ^ 2) * (z4 ^ 2 + (starRingEnd ℂ z4) ^ 2)
          + (z4 ^ 2 * (starRingEnd ℂ z4) ^ 2) * (z1 ^ 2 + (starRingEnd ℂ z1) ^ 2) * (z2 ^ 2 + (starRingEnd ℂ z2) ^ 2)))
      + (-((q4 - p4) * (((Real.sqrt 2 : ℝ) : ℂ) * Complex.sin (L * z4 / 2)) * (starRingEnd ℂ z4) ^ 2
          + (s4 - r4) * (starRingEnd ℂ (((Real.sqrt 2 : ℝ) : ℂ) * Complex.sin (L * z4 / 2))) * z4 ^ 2) * ((z1 ^ 2 * (starRingEnd ℂ z1) ^ 2) * (z2 ^ 2 + (starRingEnd ℂ z2) ^ 2) + (z1 ^ 2 * (starRingEnd ℂ z1) ^ 2) * (z3 ^ 2 + (starRingEnd ℂ z3) ^ 2) + (z2 ^ 2 * (starRingEnd ℂ z2) ^ 2) * (z1 ^ 2 + (starRingEnd ℂ z1) ^ 2)
          + (z2 ^ 2 * (starRingEnd ℂ z2) ^ 2) * (z3 ^ 2 + (starRingEnd ℂ z3) ^ 2) + (z3 ^ 2 * (starRingEnd ℂ z3) ^ 2) * (z1 ^ 2 + (starRingEnd ℂ z1) ^ 2) + (z3 ^ 2 * (starRingEnd ℂ z3) ^ 2) * (z2 ^ 2 + (starRingEnd ℂ z2) ^ 2)
          + (z1 ^ 2 + (starRingEnd ℂ z1) ^ 2) * (z2 ^ 2 + (starRingEnd ℂ z2) ^ 2) * (z3 ^ 2 + (starRingEnd ℂ z3) ^ 2))
        + ((p4 - q4) * (((Real.sqrt 2 : ℝ) : ℂ) * Complex.sin (L * z4 / 2)) + (r4 - s4) * (starRingEnd ℂ (((Real.sqrt 2 : ℝ) : ℂ) * Complex.sin (L * z4 / 2)))) * ((z1 ^ 2 * (starRingEnd ℂ z1) ^ 2) * (z2 ^ 2 * (starRingEnd ℂ z2) ^ 2) + (z1 ^ 2 * (starRingEnd ℂ z1) ^ 2) * (z3 ^ 2 * (starRingEnd ℂ z3) ^ 2) + (z2 ^ 2 * (starRingEnd ℂ z2) ^ 2) * (z3 ^ 2 * (starRingEnd ℂ z3) ^ 2)
          + (z1 ^ 2 * (starRingEnd ℂ z1) ^ 2) * (z2 ^ 2 + (starRingEnd ℂ z2) ^ 2) * (z3 ^ 2 + (starRingEnd ℂ z3) ^ 2) + (z2 ^ 2 * (starRingEnd ℂ z2) ^ 2) * (z1 ^ 2 + (starRingEnd ℂ z1) ^ 2) * (z3 ^ 2 + (starRingEnd ℂ z3) ^ 2)
          + (z3 ^ 2 * (starRingEnd ℂ z3) ^ 2) * (z1 ^ 2 + (starRingEnd ℂ z1) ^ 2) * (z2 ^ 2 + (starRingEnd ℂ z2) ^ 2))) = 0 := by
  have hpoly := crossK4PolyVanish N hN L hL z1 z2 z3 z4 hz1_im hz2_im hz3_im hz4_im p1 q1 r1 s1 p2 q2 r2 s2 p3 q3 r3 s3 p4 q4 r4 s4 hsum
  set cubic1 := (C ((p1 + q1) * (((Real.sqrt 2 : ℝ) : ℂ) * Complex.sin (L * z1 / 2))
            * (z1 * (starRingEnd ℂ z1) ^ 2)
          + (r1 + s1) * (starRingEnd ℂ (((Real.sqrt 2 : ℝ) : ℂ) * Complex.sin (L * z1 / 2)))
            * (z1 ^ 2 * (starRingEnd ℂ z1)))
        + C ((q1 - p1) * (((Real.sqrt 2 : ℝ) : ℂ) * Complex.sin (L * z1 / 2)) * (starRingEnd ℂ z1) ^ 2
          + (s1 - r1) * (starRingEnd ℂ (((Real.sqrt 2 : ℝ) : ℂ) * Complex.sin (L * z1 / 2))) * z1 ^ 2) * X
        + C (-(p1 + q1) * (((Real.sqrt 2 : ℝ) : ℂ) * Complex.sin (L * z1 / 2)) * z1
          - (r1 + s1) * (starRingEnd ℂ (((Real.sqrt 2 : ℝ) : ℂ) * Complex.sin (L * z1 / 2)))
            * (starRingEnd ℂ z1)) * X ^ 2
        + C ((p1 - q1) * (((Real.sqrt 2 : ℝ) : ℂ) * Complex.sin (L * z1 / 2))
          + (r1 - s1) * (starRingEnd ℂ (((Real.sqrt 2 : ℝ) : ℂ) * Complex.sin (L * z1 / 2)))) * X ^ 3) with hcubic1def
  set cubic2 := (C ((p2 + q2) * (((Real.sqrt 2 : ℝ) : ℂ) * Complex.sin (L * z2 / 2))
            * (z2 * (starRingEnd ℂ z2) ^ 2)
          + (r2 + s2) * (starRingEnd ℂ (((Real.sqrt 2 : ℝ) : ℂ) * Complex.sin (L * z2 / 2)))
            * (z2 ^ 2 * (starRingEnd ℂ z2)))
        + C ((q2 - p2) * (((Real.sqrt 2 : ℝ) : ℂ) * Complex.sin (L * z2 / 2)) * (starRingEnd ℂ z2) ^ 2
          + (s2 - r2) * (starRingEnd ℂ (((Real.sqrt 2 : ℝ) : ℂ) * Complex.sin (L * z2 / 2))) * z2 ^ 2) * X
        + C (-(p2 + q2) * (((Real.sqrt 2 : ℝ) : ℂ) * Complex.sin (L * z2 / 2)) * z2
          - (r2 + s2) * (starRingEnd ℂ (((Real.sqrt 2 : ℝ) : ℂ) * Complex.sin (L * z2 / 2)))
            * (starRingEnd ℂ z2)) * X ^ 2
        + C ((p2 - q2) * (((Real.sqrt 2 : ℝ) : ℂ) * Complex.sin (L * z2 / 2))
          + (r2 - s2) * (starRingEnd ℂ (((Real.sqrt 2 : ℝ) : ℂ) * Complex.sin (L * z2 / 2)))) * X ^ 3) with hcubic2def
  set cubic3 := (C ((p3 + q3) * (((Real.sqrt 2 : ℝ) : ℂ) * Complex.sin (L * z3 / 2))
            * (z3 * (starRingEnd ℂ z3) ^ 2)
          + (r3 + s3) * (starRingEnd ℂ (((Real.sqrt 2 : ℝ) : ℂ) * Complex.sin (L * z3 / 2)))
            * (z3 ^ 2 * (starRingEnd ℂ z3)))
        + C ((q3 - p3) * (((Real.sqrt 2 : ℝ) : ℂ) * Complex.sin (L * z3 / 2)) * (starRingEnd ℂ z3) ^ 2
          + (s3 - r3) * (starRingEnd ℂ (((Real.sqrt 2 : ℝ) : ℂ) * Complex.sin (L * z3 / 2))) * z3 ^ 2) * X
        + C (-(p3 + q3) * (((Real.sqrt 2 : ℝ) : ℂ) * Complex.sin (L * z3 / 2)) * z3
          - (r3 + s3) * (starRingEnd ℂ (((Real.sqrt 2 : ℝ) : ℂ) * Complex.sin (L * z3 / 2)))
            * (starRingEnd ℂ z3)) * X ^ 2
        + C ((p3 - q3) * (((Real.sqrt 2 : ℝ) : ℂ) * Complex.sin (L * z3 / 2))
          + (r3 - s3) * (starRingEnd ℂ (((Real.sqrt 2 : ℝ) : ℂ) * Complex.sin (L * z3 / 2)))) * X ^ 3) with hcubic3def
  set cubic4 := (C ((p4 + q4) * (((Real.sqrt 2 : ℝ) : ℂ) * Complex.sin (L * z4 / 2))
            * (z4 * (starRingEnd ℂ z4) ^ 2)
          + (r4 + s4) * (starRingEnd ℂ (((Real.sqrt 2 : ℝ) : ℂ) * Complex.sin (L * z4 / 2)))
            * (z4 ^ 2 * (starRingEnd ℂ z4)))
        + C ((q4 - p4) * (((Real.sqrt 2 : ℝ) : ℂ) * Complex.sin (L * z4 / 2)) * (starRingEnd ℂ z4) ^ 2
          + (s4 - r4) * (starRingEnd ℂ (((Real.sqrt 2 : ℝ) : ℂ) * Complex.sin (L * z4 / 2))) * z4 ^ 2) * X
        + C (-(p4 + q4) * (((Real.sqrt 2 : ℝ) : ℂ) * Complex.sin (L * z4 / 2)) * z4
          - (r4 + s4) * (starRingEnd ℂ (((Real.sqrt 2 : ℝ) : ℂ) * Complex.sin (L * z4 / 2)))
            * (starRingEnd ℂ z4)) * X ^ 2
        + C ((p4 - q4) * (((Real.sqrt 2 : ℝ) : ℂ) * Complex.sin (L * z4 / 2))
          + (r4 - s4) * (starRingEnd ℂ (((Real.sqrt 2 : ℝ) : ℂ) * Complex.sin (L * z4 / 2)))) * X ^ 3) with hcubic4def
  set quartic1 := (C z1 + X) * (C z1 - X) * (C (starRingEnd ℂ z1) + X) * (C (starRingEnd ℂ z1) - X) with hquartic1def
  set quartic2 := (C z2 + X) * (C z2 - X) * (C (starRingEnd ℂ z2) + X) * (C (starRingEnd ℂ z2) - X) with hquartic2def
  set quartic3 := (C z3 + X) * (C z3 - X) * (C (starRingEnd ℂ z3) + X) * (C (starRingEnd ℂ z3) - X) with hquartic3def
  set quartic4 := (C z4 + X) * (C z4 - X) * (C (starRingEnd ℂ z4) + X) * (C (starRingEnd ℂ z4) - X) with hquartic4def
  rw [quarticPolyExpand] at hquartic1def hquartic2def hquartic3def hquartic4def
  have htpe1 := quarticTripleProductExpand z2 z3 z4
  have htpe2 := quarticTripleProductExpand z1 z3 z4
  have htpe3 := quarticTripleProductExpand z1 z2 z4
  have htpe4 := quarticTripleProductExpand z1 z2 z3
  rw [← hquartic2def, ← hquartic3def, ← hquartic4def] at htpe1
  rw [← hquartic1def, ← hquartic3def, ← hquartic4def] at htpe2
  rw [← hquartic1def, ← hquartic2def, ← hquartic4def] at htpe3
  rw [← hquartic1def, ← hquartic2def, ← hquartic3def] at htpe4
  set triple1 := C ((z2 ^ 2 * (starRingEnd ℂ z2) ^ 2) * (z3 ^ 2 * (starRingEnd ℂ z3) ^ 2) * (z4 ^ 2 * (starRingEnd ℂ z4) ^ 2))
      - C ((z2 ^ 2 * (starRingEnd ℂ z2) ^ 2) * (z3 ^ 2 * (starRingEnd ℂ z3) ^ 2) * (z4 ^ 2 + (starRingEnd ℂ z4) ^ 2) + (z2 ^ 2 * (starRingEnd ℂ z2) ^ 2) * (z4 ^ 2 * (starRingEnd ℂ z4) ^ 2) * (z3 ^ 2 + (starRingEnd ℂ z3) ^ 2)
          + (z3 ^ 2 * (starRingEnd ℂ z3) ^ 2) * (z4 ^ 2 * (starRingEnd ℂ z4) ^ 2) * (z2 ^ 2 + (starRingEnd ℂ z2) ^ 2)) * X ^ 2
      + C ((z2 ^ 2 * (starRingEnd ℂ z2) ^ 2) * (z3 ^ 2 * (starRingEnd ℂ z3) ^ 2) + (z2 ^ 2 * (starRingEnd ℂ z2) ^ 2) * (z4 ^ 2 * (starRingEnd ℂ z4) ^ 2) + (z3 ^ 2 * (starRingEnd ℂ z3) ^ 2) * (z4 ^ 2 * (starRingEnd ℂ z4) ^ 2)
          + (z2 ^ 2 * (starRingEnd ℂ z2) ^ 2) * (z3 ^ 2 + (starRingEnd ℂ z3) ^ 2) * (z4 ^ 2 + (starRingEnd ℂ z4) ^ 2) + (z3 ^ 2 * (starRingEnd ℂ z3) ^ 2) * (z2 ^ 2 + (starRingEnd ℂ z2) ^ 2) * (z4 ^ 2 + (starRingEnd ℂ z4) ^ 2)
          + (z4 ^ 2 * (starRingEnd ℂ z4) ^ 2) * (z2 ^ 2 + (starRingEnd ℂ z2) ^ 2) * (z3 ^ 2 + (starRingEnd ℂ z3) ^ 2)) * X ^ 4
      - C ((z2 ^ 2 * (starRingEnd ℂ z2) ^ 2) * (z3 ^ 2 + (starRingEnd ℂ z3) ^ 2) + (z2 ^ 2 * (starRingEnd ℂ z2) ^ 2) * (z4 ^ 2 + (starRingEnd ℂ z4) ^ 2) + (z3 ^ 2 * (starRingEnd ℂ z3) ^ 2) * (z2 ^ 2 + (starRingEnd ℂ z2) ^ 2)
          + (z3 ^ 2 * (starRingEnd ℂ z3) ^ 2) * (z4 ^ 2 + (starRingEnd ℂ z4) ^ 2) + (z4 ^ 2 * (starRingEnd ℂ z4) ^ 2) * (z2 ^ 2 + (starRingEnd ℂ z2) ^ 2) + (z4 ^ 2 * (starRingEnd ℂ z4) ^ 2) * (z3 ^ 2 + (starRingEnd ℂ z3) ^ 2)
          + (z2 ^ 2 + (starRingEnd ℂ z2) ^ 2) * (z3 ^ 2 + (starRingEnd ℂ z3) ^ 2) * (z4 ^ 2 + (starRingEnd ℂ z4) ^ 2)) * X ^ 6
      + C ((z2 ^ 2 + (starRingEnd ℂ z2) ^ 2) * (z3 ^ 2 + (starRingEnd ℂ z3) ^ 2) + (z2 ^ 2 + (starRingEnd ℂ z2) ^ 2) * (z4 ^ 2 + (starRingEnd ℂ z4) ^ 2) + (z3 ^ 2 + (starRingEnd ℂ z3) ^ 2) * (z4 ^ 2 + (starRingEnd ℂ z4) ^ 2)
          + (z2 ^ 2 * (starRingEnd ℂ z2) ^ 2) + (z3 ^ 2 * (starRingEnd ℂ z3) ^ 2) + (z4 ^ 2 * (starRingEnd ℂ z4) ^ 2)) * X ^ 8
      - C ((z2 ^ 2 + (starRingEnd ℂ z2) ^ 2) + (z3 ^ 2 + (starRingEnd ℂ z3) ^ 2) + (z4 ^ 2 + (starRingEnd ℂ z4) ^ 2)) * X ^ 10
      + X ^ 12 with htriple1def
  set triple2 := C ((z1 ^ 2 * (starRingEnd ℂ z1) ^ 2) * (z3 ^ 2 * (starRingEnd ℂ z3) ^ 2) * (z4 ^ 2 * (starRingEnd ℂ z4) ^ 2))
      - C ((z1 ^ 2 * (starRingEnd ℂ z1) ^ 2) * (z3 ^ 2 * (starRingEnd ℂ z3) ^ 2) * (z4 ^ 2 + (starRingEnd ℂ z4) ^ 2) + (z1 ^ 2 * (starRingEnd ℂ z1) ^ 2) * (z4 ^ 2 * (starRingEnd ℂ z4) ^ 2) * (z3 ^ 2 + (starRingEnd ℂ z3) ^ 2)
          + (z3 ^ 2 * (starRingEnd ℂ z3) ^ 2) * (z4 ^ 2 * (starRingEnd ℂ z4) ^ 2) * (z1 ^ 2 + (starRingEnd ℂ z1) ^ 2)) * X ^ 2
      + C ((z1 ^ 2 * (starRingEnd ℂ z1) ^ 2) * (z3 ^ 2 * (starRingEnd ℂ z3) ^ 2) + (z1 ^ 2 * (starRingEnd ℂ z1) ^ 2) * (z4 ^ 2 * (starRingEnd ℂ z4) ^ 2) + (z3 ^ 2 * (starRingEnd ℂ z3) ^ 2) * (z4 ^ 2 * (starRingEnd ℂ z4) ^ 2)
          + (z1 ^ 2 * (starRingEnd ℂ z1) ^ 2) * (z3 ^ 2 + (starRingEnd ℂ z3) ^ 2) * (z4 ^ 2 + (starRingEnd ℂ z4) ^ 2) + (z3 ^ 2 * (starRingEnd ℂ z3) ^ 2) * (z1 ^ 2 + (starRingEnd ℂ z1) ^ 2) * (z4 ^ 2 + (starRingEnd ℂ z4) ^ 2)
          + (z4 ^ 2 * (starRingEnd ℂ z4) ^ 2) * (z1 ^ 2 + (starRingEnd ℂ z1) ^ 2) * (z3 ^ 2 + (starRingEnd ℂ z3) ^ 2)) * X ^ 4
      - C ((z1 ^ 2 * (starRingEnd ℂ z1) ^ 2) * (z3 ^ 2 + (starRingEnd ℂ z3) ^ 2) + (z1 ^ 2 * (starRingEnd ℂ z1) ^ 2) * (z4 ^ 2 + (starRingEnd ℂ z4) ^ 2) + (z3 ^ 2 * (starRingEnd ℂ z3) ^ 2) * (z1 ^ 2 + (starRingEnd ℂ z1) ^ 2)
          + (z3 ^ 2 * (starRingEnd ℂ z3) ^ 2) * (z4 ^ 2 + (starRingEnd ℂ z4) ^ 2) + (z4 ^ 2 * (starRingEnd ℂ z4) ^ 2) * (z1 ^ 2 + (starRingEnd ℂ z1) ^ 2) + (z4 ^ 2 * (starRingEnd ℂ z4) ^ 2) * (z3 ^ 2 + (starRingEnd ℂ z3) ^ 2)
          + (z1 ^ 2 + (starRingEnd ℂ z1) ^ 2) * (z3 ^ 2 + (starRingEnd ℂ z3) ^ 2) * (z4 ^ 2 + (starRingEnd ℂ z4) ^ 2)) * X ^ 6
      + C ((z1 ^ 2 + (starRingEnd ℂ z1) ^ 2) * (z3 ^ 2 + (starRingEnd ℂ z3) ^ 2) + (z1 ^ 2 + (starRingEnd ℂ z1) ^ 2) * (z4 ^ 2 + (starRingEnd ℂ z4) ^ 2) + (z3 ^ 2 + (starRingEnd ℂ z3) ^ 2) * (z4 ^ 2 + (starRingEnd ℂ z4) ^ 2)
          + (z1 ^ 2 * (starRingEnd ℂ z1) ^ 2) + (z3 ^ 2 * (starRingEnd ℂ z3) ^ 2) + (z4 ^ 2 * (starRingEnd ℂ z4) ^ 2)) * X ^ 8
      - C ((z1 ^ 2 + (starRingEnd ℂ z1) ^ 2) + (z3 ^ 2 + (starRingEnd ℂ z3) ^ 2) + (z4 ^ 2 + (starRingEnd ℂ z4) ^ 2)) * X ^ 10
      + X ^ 12 with htriple2def
  set triple3 := C ((z1 ^ 2 * (starRingEnd ℂ z1) ^ 2) * (z2 ^ 2 * (starRingEnd ℂ z2) ^ 2) * (z4 ^ 2 * (starRingEnd ℂ z4) ^ 2))
      - C ((z1 ^ 2 * (starRingEnd ℂ z1) ^ 2) * (z2 ^ 2 * (starRingEnd ℂ z2) ^ 2) * (z4 ^ 2 + (starRingEnd ℂ z4) ^ 2) + (z1 ^ 2 * (starRingEnd ℂ z1) ^ 2) * (z4 ^ 2 * (starRingEnd ℂ z4) ^ 2) * (z2 ^ 2 + (starRingEnd ℂ z2) ^ 2)
          + (z2 ^ 2 * (starRingEnd ℂ z2) ^ 2) * (z4 ^ 2 * (starRingEnd ℂ z4) ^ 2) * (z1 ^ 2 + (starRingEnd ℂ z1) ^ 2)) * X ^ 2
      + C ((z1 ^ 2 * (starRingEnd ℂ z1) ^ 2) * (z2 ^ 2 * (starRingEnd ℂ z2) ^ 2) + (z1 ^ 2 * (starRingEnd ℂ z1) ^ 2) * (z4 ^ 2 * (starRingEnd ℂ z4) ^ 2) + (z2 ^ 2 * (starRingEnd ℂ z2) ^ 2) * (z4 ^ 2 * (starRingEnd ℂ z4) ^ 2)
          + (z1 ^ 2 * (starRingEnd ℂ z1) ^ 2) * (z2 ^ 2 + (starRingEnd ℂ z2) ^ 2) * (z4 ^ 2 + (starRingEnd ℂ z4) ^ 2) + (z2 ^ 2 * (starRingEnd ℂ z2) ^ 2) * (z1 ^ 2 + (starRingEnd ℂ z1) ^ 2) * (z4 ^ 2 + (starRingEnd ℂ z4) ^ 2)
          + (z4 ^ 2 * (starRingEnd ℂ z4) ^ 2) * (z1 ^ 2 + (starRingEnd ℂ z1) ^ 2) * (z2 ^ 2 + (starRingEnd ℂ z2) ^ 2)) * X ^ 4
      - C ((z1 ^ 2 * (starRingEnd ℂ z1) ^ 2) * (z2 ^ 2 + (starRingEnd ℂ z2) ^ 2) + (z1 ^ 2 * (starRingEnd ℂ z1) ^ 2) * (z4 ^ 2 + (starRingEnd ℂ z4) ^ 2) + (z2 ^ 2 * (starRingEnd ℂ z2) ^ 2) * (z1 ^ 2 + (starRingEnd ℂ z1) ^ 2)
          + (z2 ^ 2 * (starRingEnd ℂ z2) ^ 2) * (z4 ^ 2 + (starRingEnd ℂ z4) ^ 2) + (z4 ^ 2 * (starRingEnd ℂ z4) ^ 2) * (z1 ^ 2 + (starRingEnd ℂ z1) ^ 2) + (z4 ^ 2 * (starRingEnd ℂ z4) ^ 2) * (z2 ^ 2 + (starRingEnd ℂ z2) ^ 2)
          + (z1 ^ 2 + (starRingEnd ℂ z1) ^ 2) * (z2 ^ 2 + (starRingEnd ℂ z2) ^ 2) * (z4 ^ 2 + (starRingEnd ℂ z4) ^ 2)) * X ^ 6
      + C ((z1 ^ 2 + (starRingEnd ℂ z1) ^ 2) * (z2 ^ 2 + (starRingEnd ℂ z2) ^ 2) + (z1 ^ 2 + (starRingEnd ℂ z1) ^ 2) * (z4 ^ 2 + (starRingEnd ℂ z4) ^ 2) + (z2 ^ 2 + (starRingEnd ℂ z2) ^ 2) * (z4 ^ 2 + (starRingEnd ℂ z4) ^ 2)
          + (z1 ^ 2 * (starRingEnd ℂ z1) ^ 2) + (z2 ^ 2 * (starRingEnd ℂ z2) ^ 2) + (z4 ^ 2 * (starRingEnd ℂ z4) ^ 2)) * X ^ 8
      - C ((z1 ^ 2 + (starRingEnd ℂ z1) ^ 2) + (z2 ^ 2 + (starRingEnd ℂ z2) ^ 2) + (z4 ^ 2 + (starRingEnd ℂ z4) ^ 2)) * X ^ 10
      + X ^ 12 with htriple3def
  set triple4 := C ((z1 ^ 2 * (starRingEnd ℂ z1) ^ 2) * (z2 ^ 2 * (starRingEnd ℂ z2) ^ 2) * (z3 ^ 2 * (starRingEnd ℂ z3) ^ 2))
      - C ((z1 ^ 2 * (starRingEnd ℂ z1) ^ 2) * (z2 ^ 2 * (starRingEnd ℂ z2) ^ 2) * (z3 ^ 2 + (starRingEnd ℂ z3) ^ 2) + (z1 ^ 2 * (starRingEnd ℂ z1) ^ 2) * (z3 ^ 2 * (starRingEnd ℂ z3) ^ 2) * (z2 ^ 2 + (starRingEnd ℂ z2) ^ 2)
          + (z2 ^ 2 * (starRingEnd ℂ z2) ^ 2) * (z3 ^ 2 * (starRingEnd ℂ z3) ^ 2) * (z1 ^ 2 + (starRingEnd ℂ z1) ^ 2)) * X ^ 2
      + C ((z1 ^ 2 * (starRingEnd ℂ z1) ^ 2) * (z2 ^ 2 * (starRingEnd ℂ z2) ^ 2) + (z1 ^ 2 * (starRingEnd ℂ z1) ^ 2) * (z3 ^ 2 * (starRingEnd ℂ z3) ^ 2) + (z2 ^ 2 * (starRingEnd ℂ z2) ^ 2) * (z3 ^ 2 * (starRingEnd ℂ z3) ^ 2)
          + (z1 ^ 2 * (starRingEnd ℂ z1) ^ 2) * (z2 ^ 2 + (starRingEnd ℂ z2) ^ 2) * (z3 ^ 2 + (starRingEnd ℂ z3) ^ 2) + (z2 ^ 2 * (starRingEnd ℂ z2) ^ 2) * (z1 ^ 2 + (starRingEnd ℂ z1) ^ 2) * (z3 ^ 2 + (starRingEnd ℂ z3) ^ 2)
          + (z3 ^ 2 * (starRingEnd ℂ z3) ^ 2) * (z1 ^ 2 + (starRingEnd ℂ z1) ^ 2) * (z2 ^ 2 + (starRingEnd ℂ z2) ^ 2)) * X ^ 4
      - C ((z1 ^ 2 * (starRingEnd ℂ z1) ^ 2) * (z2 ^ 2 + (starRingEnd ℂ z2) ^ 2) + (z1 ^ 2 * (starRingEnd ℂ z1) ^ 2) * (z3 ^ 2 + (starRingEnd ℂ z3) ^ 2) + (z2 ^ 2 * (starRingEnd ℂ z2) ^ 2) * (z1 ^ 2 + (starRingEnd ℂ z1) ^ 2)
          + (z2 ^ 2 * (starRingEnd ℂ z2) ^ 2) * (z3 ^ 2 + (starRingEnd ℂ z3) ^ 2) + (z3 ^ 2 * (starRingEnd ℂ z3) ^ 2) * (z1 ^ 2 + (starRingEnd ℂ z1) ^ 2) + (z3 ^ 2 * (starRingEnd ℂ z3) ^ 2) * (z2 ^ 2 + (starRingEnd ℂ z2) ^ 2)
          + (z1 ^ 2 + (starRingEnd ℂ z1) ^ 2) * (z2 ^ 2 + (starRingEnd ℂ z2) ^ 2) * (z3 ^ 2 + (starRingEnd ℂ z3) ^ 2)) * X ^ 6
      + C ((z1 ^ 2 + (starRingEnd ℂ z1) ^ 2) * (z2 ^ 2 + (starRingEnd ℂ z2) ^ 2) + (z1 ^ 2 + (starRingEnd ℂ z1) ^ 2) * (z3 ^ 2 + (starRingEnd ℂ z3) ^ 2) + (z2 ^ 2 + (starRingEnd ℂ z2) ^ 2) * (z3 ^ 2 + (starRingEnd ℂ z3) ^ 2)
          + (z1 ^ 2 * (starRingEnd ℂ z1) ^ 2) + (z2 ^ 2 * (starRingEnd ℂ z2) ^ 2) + (z3 ^ 2 * (starRingEnd ℂ z3) ^ 2)) * X ^ 8
      - C ((z1 ^ 2 + (starRingEnd ℂ z1) ^ 2) + (z2 ^ 2 + (starRingEnd ℂ z2) ^ 2) + (z3 ^ 2 + (starRingEnd ℂ z3) ^ 2)) * X ^ 10
      + X ^ 12 with htriple4def
  have hpoly' : cubic1 * triple1 + cubic2 * triple2 + cubic3 * triple3 + cubic4 * triple4 = 0 := by
    linear_combination hpoly - cubic1 * htpe1 - cubic2 * htpe2 - cubic3 * htpe3 - cubic4 * htpe4
  have hc7 := congrArg (fun p => Polynomial.coeff p 7) hpoly'
  simp only [coeff_add, coeff_zero] at hc7
  rw [hcubic1def, htriple1def, hcubic2def, htriple2def, hcubic3def, htriple3def, hcubic4def, htriple4def] at hc7
  rw [coeff_mul, coeff_mul, coeff_mul, coeff_mul] at hc7
  simp only [Finset.Nat.sum_antidiagonal_eq_sum_range_succ_mk] at hc7
  simp only [Finset.sum_range_succ, Finset.sum_range_zero, coeff_add, coeff_sub, coeff_C_mul,
    coeff_C, coeff_X_pow, coeff_X, zero_add] at hc7
  norm_num at hc7
  simp only [map_mul, Complex.conj_ofReal] at hc7 ⊢
  linear_combination hc7

#print axioms crossK4NinthCoeff
