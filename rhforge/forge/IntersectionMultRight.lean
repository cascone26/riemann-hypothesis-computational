import IntersectionMult

set_option maxHeartbeats 1000000

/-!
# §1629 — symmetric completion of the local intersection-multiplicity API (F₁/arith-Hodge lane)

Right-hand-side companions of the proven `interMult`/`cartierInterMult` lemmas, obtained by symmetry
(`interMult_comm`): the pairing vanishes against the whole ring on the left, is antitone in its *second* ideal,
and a unit on the *right* meets nothing.  These complete the symmetric API of the arithmetic intersection pairing
`D·E` on the affine core of an arithmetic surface — genuine commutative algebra, no RH-crux content.
EXIT0 / std-axioms.
-/

variable {R : Type*} [CommRing R]

/-- The whole ring (empty divisor) meets nothing from the left: `(⊤)·J = 0`. -/
theorem interMult_top_left (J : Ideal R) : interMult (⊤ : Ideal R) J = 0 := by
  rw [interMult_comm]; exact interMult_top J

/-- Intersection multiplicity is antitone in the *second* ideal: enlarging the second divisor ideal
    can only cut a shorter intersection. -/
theorem interMult_antitone_right (I : Ideal R) {J J2 : Ideal R} (h : J2 ≤ J) :
    interMult I J ≤ interMult I J2 := by
  rw [interMult_comm I J, interMult_comm I J2]
  exact interMult_antitone_left I h

/-- A unit (empty divisor) on the right meets nothing: `(f)·(unit) = 0`. -/
theorem cartierInterMult_unit_right {g : R} (hg : IsUnit g) (f : R) :
    cartierInterMult f g = 0 := by
  rw [cartierInterMult_comm]; exact cartierInterMult_unit_left hg f

#print axioms interMult_top_left
#print axioms interMult_antitone_right
#print axioms cartierInterMult_unit_right
