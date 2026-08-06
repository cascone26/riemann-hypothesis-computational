import Mathlib.LinearAlgebra.Matrix.PosDef
import Mathlib.Data.Real.StarOrdered
open Matrix
/-- **Weil-certificate reduction** (isolates the one transcendental step): if the rationalized
lower-bound matrix `M` is PosSemidef (proven by the exact-rational LDL^T certificate,
`organ2_weil_cert_posSemidef`) and the true Weil matrix `Q` dominates it — `Q - M` PosSemidef,
the rigorous Arb interval enclosure (§755) — then the true Weil form `Q` is PosSemidef.
So Weil positivity for the basis reduces to the single labeled interval-bound `hQM`. -/
theorem weil_posSemidef_of_dominates {n : Type*} [Fintype n] [DecidableEq n]
    {Q M : Matrix n n ℝ} (hM : M.PosSemidef) (hQM : (Q - M).PosSemidef) :
    Q.PosSemidef := by
  have : Q = M + (Q - M) := by abel
  rw [this]
  exact hM.add hQM

#print axioms weil_posSemidef_of_dominates
