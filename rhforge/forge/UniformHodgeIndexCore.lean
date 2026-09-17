import HodgeIndexCore

set_option maxHeartbeats 1000000

open QuadraticMap Submodule

variable {F M : Type*} [Field F] [LinearOrder F] [IsStrictOrderedRing F]
  [AddCommGroup M] [Module F M] [FiniteDimensional F M]

/-- **Uniform arithmetic Hodge index (abstract core).** If a quadratic form `Q` is positive on one
vector `w` off the kernel of a functional `φ`, and negative-definite on `ker φ`, then it has exactly
ONE positive square: `sigPos Q = 1`. (P = span{w} posdef line; N = ker φ negdef hyperplane.) Stated
over a general ordered field, so instantiating at ℝ carries no instance-diamond. -/
theorem sigPos_eq_one_of_pos_and_negdef_ker
    (Q : QuadraticForm F M) (φ : M →ₗ[F] F)
    (w : M) (hw : 0 < Q w) (hwφ : φ w ≠ 0)
    (hneg : ∀ z ∈ LinearMap.ker φ, z ≠ 0 → Q z < 0) :
    sigPos Q = 1 := by
  have hw0 : w ≠ 0 := by
    intro h; rw [h, QuadraticMap.map_zero] at hw; exact lt_irrefl 0 hw
  set P : Submodule F M := span F {w} with hPdef
  -- P is a positive-definite line
  have hP1 : Module.finrank F P = 1 := finrank_span_singleton hw0
  have hPpos : (Q.restrict P).PosDef := by
    intro z hz
    obtain ⟨t, ht⟩ := mem_span_singleton.mp z.2
    have ht0 : t ≠ 0 := by
      intro h0; apply hz; apply Subtype.ext
      rw [ZeroMemClass.coe_zero, ← ht, h0, zero_smul]
    have : (Q.restrict P) z = t * t * Q w := by
      have h1 : (Q.restrict P) z = Q (z : M) := restrict_apply _ _ _
      rw [h1, ← ht, QuadraticMap.map_smul]
      ring
    rw [this]; exact mul_pos (mul_self_pos.mpr ht0) hw
  -- N = ker φ is a negative-definite hyperplane
  set N : Submodule F M := LinearMap.ker φ with hNdef
  have hNneg : ((-Q).restrict N).PosDef := by
    refine fun z hz => ?_
    have hz0 : (z : M) ≠ 0 := by
      intro h; apply hz; exact Subtype.ext h
    rw [restrict_apply, neg_apply]
    have := hneg (z : M) z.2 hz0
    linarith
  have hφ0 : φ ≠ 0 := by intro h; rw [h] at hwφ; exact hwφ rfl
  have hrange : Module.finrank F (LinearMap.range φ) = 1 := by
    have hsurj : LinearMap.range φ = ⊤ := by
      rw [LinearMap.range_eq_top]
      intro y
      refine ⟨(y / φ w) • w, ?_⟩
      rw [map_smul, smul_eq_mul, div_mul_cancel₀ _ hwφ]
    rw [hsurj, finrank_top, Module.finrank_self]
  have hNrank : Module.finrank F N = Module.finrank F M - 1 := by
    have h := LinearMap.finrank_range_add_finrank_ker φ
    rw [hrange, ← hNdef] at h; omega
  have hM1 : 1 ≤ Module.finrank F M := by
    rw [← hP1]; exact Submodule.finrank_le P
  exact (hodgeIndex_signature Q P N hP1 hPpos hNneg hNrank hM1).1

#print axioms sigPos_eq_one_of_pos_and_negdef_ker

/-- **Uniform arithmetic Hodge index — full signature.** Same hypotheses as `sigPos_eq_one_of_
pos_and_negdef_ker`, additional conclusion `sigNeg Q = dim M - 1`: the form has signature `(1, dim−1)`. -/
theorem sig_full_of_pos_and_negdef_ker
    (Q : QuadraticForm F M) (φ : M →ₗ[F] F)
    (w : M) (hw : 0 < Q w) (hwφ : φ w ≠ 0)
    (hneg : ∀ z ∈ LinearMap.ker φ, z ≠ 0 → Q z < 0) :
    sigPos Q = 1 ∧ sigNeg Q = Module.finrank F M - 1 := by
  have hw0 : w ≠ 0 := by
    intro h; rw [h, QuadraticMap.map_zero] at hw; exact lt_irrefl 0 hw
  set P : Submodule F M := span F {w} with hPdef
  have hP1 : Module.finrank F P = 1 := finrank_span_singleton hw0
  have hPpos : (Q.restrict P).PosDef := by
    intro z hz
    obtain ⟨t, ht⟩ := mem_span_singleton.mp z.2
    have ht0 : t ≠ 0 := by
      intro h0; apply hz; apply Subtype.ext
      rw [ZeroMemClass.coe_zero, ← ht, h0, zero_smul]
    have : (Q.restrict P) z = t * t * Q w := by
      have h1 : (Q.restrict P) z = Q (z : M) := restrict_apply _ _ _
      rw [h1, ← ht, QuadraticMap.map_smul]; ring
    rw [this]; exact mul_pos (mul_self_pos.mpr ht0) hw
  set N : Submodule F M := LinearMap.ker φ with hNdef
  have hNneg : ((-Q).restrict N).PosDef := by
    refine fun z hz => ?_
    have hz0 : (z : M) ≠ 0 := by intro h; apply hz; exact Subtype.ext h
    rw [restrict_apply, neg_apply]
    have := hneg (z : M) z.2 hz0; linarith
  have hφ0 : φ ≠ 0 := by intro h; rw [h] at hwφ; exact hwφ rfl
  have hrange : Module.finrank F (LinearMap.range φ) = 1 := by
    have hsurj : LinearMap.range φ = ⊤ := by
      rw [LinearMap.range_eq_top]; intro y
      exact ⟨(y / φ w) • w, by rw [map_smul, smul_eq_mul, div_mul_cancel₀ _ hwφ]⟩
    rw [hsurj, finrank_top, Module.finrank_self]
  have hNrank : Module.finrank F N = Module.finrank F M - 1 := by
    have h := LinearMap.finrank_range_add_finrank_ker φ
    rw [hrange, ← hNdef] at h; omega
  have hM1 : 1 ≤ Module.finrank F M := by rw [← hP1]; exact Submodule.finrank_le P
  exact hodgeIndex_signature Q P N hP1 hPpos hNneg hNrank hM1

#print axioms sig_full_of_pos_and_negdef_ker
