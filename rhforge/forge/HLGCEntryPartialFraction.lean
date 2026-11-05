import HLGCPartialFraction

open scoped MatrixOrder Matrix

/-!
# Off-critical-line lane: the FULL row-entry closed form factors as prefactor × partial-fraction split

Step 4-a of §2391's 4-step plan (assembling toward the matrix-level statement). §2391's
`quadratic_partial_fraction_split` handles the bare algebraic factor `(z²+a·u)/(z²-u²)`; this brick
pulls out §2386's `H_L_G_C_entry_closed_form`'s full common prefactor `4sin²(Lz/2)/(L²(z²-ν_m²))`
(which does NOT depend on the column frequency `ν_n`) and shows the WHOLE row-entry closed form
factors as `prefactor × [(z+ν_m)/(2(z-ν_n)) + (z-ν_m)/(2(z+ν_n))]` — i.e. multiplying
`quadratic_partial_fraction_split` through by the prefactor commutes with the division, which
`field_simp`/`ring` confirm directly given the same four pole-avoidance facts. This is the exact
per-term identity the row-sum `∑_n H^G_{mn,L}(z_k)·e^{off}_n` needs before its terms can be split
into the two `c_off_null_identity`/`c_off_null_identity_reflected` sums (the remaining honest work:
sum both sides over `n`, distribute, factor the now-n-independent prefactor and `(z±ν_m)/2` factors
out of their respective sums, invoke the two null identities). Search-before-build: grepped
`H_L_G_C_entry_partial_fraction`/`entry_partial_fraction` — fresh. Fully general, no zeta/RH
dependence. -/
theorem H_L_G_C_entry_partial_fraction (L : ℝ) (z u a : ℂ)
    (hzu : z ≠ u) (hzu' : z ≠ -u) :
    (4 * (Complex.sin (L * z / 2)) ^ 2 * (z ^ 2 + a * u))
        / (L ^ 2 * (z ^ 2 - a ^ 2) * (z ^ 2 - u ^ 2))
      = (4 * (Complex.sin (L * z / 2)) ^ 2 / (L ^ 2 * (z ^ 2 - a ^ 2)))
          * ((z + a) / (2 * (z - u)) + (z - a) / (2 * (z + u))) := by
  rw [← quadratic_partial_fraction_split z u a hzu hzu']
  by_cases ha : z ^ 2 - a ^ 2 = 0
  · simp [ha]
  · field_simp

#print axioms H_L_G_C_entry_partial_fraction
