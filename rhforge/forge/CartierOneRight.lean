import IntersectionMultRight

set_option maxHeartbeats 1000000

/-- One (a unit) on the right meets nothing: `(f)·(1) = 0` — a direct specialization of §1629. -/
theorem cartierInterMult_one_right {R : Type*} [CommRing R] (f : R) :
    cartierInterMult f 1 = 0 :=
  cartierInterMult_unit_right isUnit_one f

#print axioms cartierInterMult_one_right
