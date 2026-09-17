import HodgeIndexInequalityOnModel
import ZtoRModelEmbed
import HodgeFormModel

set_option maxHeartbeats 4000000

open IsDedekindDomain NumberField QuadraticMap Submodule

noncomputable section

variable {K : Type*} [Field K] [NumberField K]

/-!
# §1722 ★★★★★ — THE ARITHMETIC HODGE INDEX INEQUALITY on actual ARAKELOV DIVISORS
Transports §1721 (ℝ-model) to genuine Arakelov divisors via the pairing/degree-compatible embedding
§1717. For Arakelov divisors `D, E` whose ℝ-images lie in a subspace with a hodge-positive vector,
with `hodgeForm D D > 0`:  `hodgeForm D D · hodgeForm E E ≤ hodgeForm D E ²` — the arithmetic
intersection form's Lorentzian positivity on real divisors, the function-field-RH analog.
-/

/-- The ℤ-model Hodge form IS the ℝ-model Hodge bilinear form under the embedding (via §1717). -/
theorem hodgeForm_eq_hodgeBilin (D E : ArakelovDivisor K) :
    hodgeForm D E = ArakelovVecR.hodgeBilin K (toVecR D) (toVecR E) := by
  rw [ArakelovVecR.hodgeBilin_apply, degFun_toVecR, degFun_toVecR, pairingFun_toVecR]
  rfl

theorem hodge_index_inequality_divisor
    (W : Submodule ℝ (ArakelovVecR K)) [FiniteDimensional ℝ W]
    (wp : ArakelovVecR K) (hwpW : wp ∈ W) (hwppos : 0 < ArakelovVecR.hodgeBilin K wp wp)
    (D E : ArakelovDivisor K) (hDW : toVecR D ∈ W) (hEW : toVecR E ∈ W)
    (hDpos : 0 < hodgeForm D D) :
    hodgeForm D D * hodgeForm E E ≤ hodgeForm D E ^ 2 := by
  rw [hodgeForm_eq_hodgeBilin, hodgeForm_eq_hodgeBilin, hodgeForm_eq_hodgeBilin]
  rw [hodgeForm_eq_hodgeBilin] at hDpos
  exact ArakelovVecR.hodge_index_inequality_model K W wp hwpW hwppos
    (toVecR D) (toVecR E) hDW hEW hDpos

#print axioms hodge_index_inequality_divisor
