import P2KLNodalBridge

open scoped MatrixOrder Matrix

/-!
# Off-critical-line lane: the two-term partial-fraction split of `H^G_L`'s node-dependence

Hand-derivation (worked out before writing any Lean): fixing the evaluation point `z` and the row
index `ν_m`, and viewing §2386's `H^G_L`-entry closed form as a function of the COLUMN frequency
`u:=ν_n`, the factor `(z²+ν_m·u)/(z²-u²)` (stripped of the common `4sin²(Lz/2)/L²` prefactor) admits
the clean two-term partial-fraction split `(z²+ν_m·u)/(z²-u²) = (z+ν_m)/(2(z-u)) + (z-ν_m)/(2(z+u))`
— checked by clearing denominators: `(z+ν_m)(z+u)+(z-ν_m)(z-u) = 2z²+2ν_m u`, matching twice the
numerator. This is the KEY missing algebraic step for the deferred matrix-level null-vector
connection: it turns the ROW SUM `∑_n H^G_{mn,L}(z_k)·e^{off}_n` (a sum over a two-variable rational
function) into TWO copies of the SAME single-variable sum `∑_n e^{off}_n/(z_k∓ν_n)` that
§2380-2383's `Q_off_K_null_identity` already shows vanishes — modulo the STILL-OPEN question of
whether the "+ν_n" copy also vanishes (needs a reflection-symmetry argument using `Q_off_K_even`,
correctly NOT attempted this cycle). Search-before-build: grepped for `partial_fraction_two_term`/
`quadratic_partial_fraction` — fresh. Fully general, no zeta/RH dependence. -/
theorem quadratic_partial_fraction_split (z u a : ℂ) (hzu : z ≠ u) (hzu' : z ≠ -u) :
    (z ^ 2 + a * u) / (z ^ 2 - u ^ 2) = (z + a) / (2 * (z - u)) + (z - a) / (2 * (z + u)) := by
  have hzmu : z - u ≠ 0 := sub_ne_zero_of_ne hzu
  have hzpu : z + u ≠ 0 := by
    intro h
    apply hzu'
    linear_combination h
  have hz2u2 : z ^ 2 - u ^ 2 ≠ 0 := by
    have heq : z ^ 2 - u ^ 2 = (z - u) * (z + u) := by ring
    rw [heq]; exact mul_ne_zero hzmu hzpu
  field_simp
  ring

#print axioms quadratic_partial_fraction_split
