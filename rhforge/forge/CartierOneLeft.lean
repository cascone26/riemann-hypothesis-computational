import IntersectionMult

set_option maxHeartbeats 1000000

/-- One (a unit) on the left meets nothing: `(1)·(g) = 0` (symmetric partner of §1631). -/
theorem cartierInterMult_one_left {R : Type*} [CommRing R] (g : R) :
    cartierInterMult 1 g = 0 :=
  cartierInterMult_unit_left isUnit_one g

#print axioms cartierInterMult_one_left
