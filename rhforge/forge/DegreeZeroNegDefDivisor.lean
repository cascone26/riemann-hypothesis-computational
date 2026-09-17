import HodgeIndexInequalityDivisor
import ZtoRModelEmbed
import HodgeFormModel

set_option maxHeartbeats 4000000

open IsDedekindDomain NumberField

noncomputable section

variable {K : Type*} [Field K] [NumberField K]

/-!
# §1723 ★★★★★ — THE HODGE FORM IS NEGATIVE-DEFINITE ON DEGREE-ZERO ARAKELOV DIVISORS

The arithmetic analog of the engine of the function-field RH proof: on a curve `C/𝔽_q`, the Weil
bound comes from the intersection form being **negative-definite on the degree-zero (numerically
primitive) part** of the Néron–Severi space — the negative side of the Hodge-index signature
`(1, ρ−1)`. §1720 gave that signature on the ℝ-model; §1713 gave negdef on `ker degLin` there.
This brick lands the negative-definiteness on **genuine ℤ-model Arakelov divisors**: any Arakelov
divisor `D` of Arakelov degree zero with nonzero ℝ-image has `hodgeForm D D < 0`.

Transport is via the pairing/degree-compatible embedding §1717 (`toVecR`) and the §1722 bridge
`hodgeForm_eq_hodgeBilin`, then §1713's `hodgeBilin_self_neg_of_deg_zero`. Together with §1722
(Lorentzian inequality, positive on a hodge-positive direction) this makes the arithmetic
Hodge-index dichotomy — positive on one direction, negative on the degree-zero locus — hold on
actual Arakelov divisors: the function-field-RH object, fully self-contained on the ℤ-model.
-/

/-- **Negative-definiteness on the degree-zero locus (Arakelov divisors).**
For a genuine Arakelov divisor `D` of Arakelov degree zero whose ℝ-image is nonzero, the arithmetic
Hodge form is strictly negative: `hodgeForm D D < 0`. This is the arithmetic analog of "numerically
trivial ⟹ negative self-intersection" — the negative side of the Hodge index that powers the
function-field Riemann Hypothesis. -/
theorem hodgeForm_neg_of_degree_zero (D : ArakelovDivisor K)
    (hdeg : arakDeg D = 0) (hne : toVecR D ≠ 0) :
    hodgeForm D D < 0 := by
  rw [hodgeForm_eq_hodgeBilin]
  apply ArakelovVecR.hodgeBilin_self_neg_of_deg_zero
  · rw [ArakelovVecR.degLin_apply, degFun_toVecR, hdeg]
  · exact hne

/-- Every principal Arakelov divisor with nonzero ℝ-image is Hodge-negative: principal divisors
are degree zero (§1711 `arakDeg_principalArakelov`), so this is the degree-zero negdef specialized
to the principal (numerically-trivial-by-construction) case. -/
theorem hodgeForm_neg_principal (x : Kˣ)
    (hne : toVecR (principalArakelov (K := K) x) ≠ 0) :
    hodgeForm (principalArakelov (K := K) x) (principalArakelov (K := K) x) < 0 :=
  hodgeForm_neg_of_degree_zero _ (arakDeg_principalArakelov x) hne

#print axioms hodgeForm_neg_of_degree_zero
#print axioms hodgeForm_neg_principal
