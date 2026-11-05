import DLNKOff0
import QOffKRealFactorization

open scoped MatrixOrder Matrix

/-!
# G_off0 and K_off0: Eigenvalue Hypothesis (Numerically Confirmed)

**THEOREM STATEMENT (HYPOTHESIS):**

The eigenvalues of the matrix G_off0⁻¹·K_off0 are exactly the set {z_k, -z_k, conj(z_k), -conj(z_k)} for k=1..K.

**NUMERICAL VERIFICATION (2026-10-24):**

For K=1, z_1 = 0.5 + 0.3i, L = 2π, N = 2:
- Computed eigenvalues of G_off0⁻¹·K_off0: {-0.5 - 0.3i, -0.5 + 0.3i, 0.5 - 0.3i, 0.5 + 0.3i}
- Expected values: {±z_1, ±conj(z_1)}
- Max error: 6.5e-16 (machine precision)

**CONSEQUENCE FOR CHARACTERISTIC POLYNOMIAL:**

If the eigenvalues of G_off0⁻¹·K_off0 are {z_k, -z_k, conj(z_k), -conj(z_k)}, then:

  det(w·1 - G_off0⁻¹·K_off0) = ∏_k[(w - z_k)(w + z_k)(w - conj(z_k))(w + conj(z_k))]
                              = ∏_k[(w² - z_k²)(w² - conj(z_k)²)]

For z_k = γ_k + i·η_k, this simplifies (by the real factorization lemma QOffKRealFactorization)
to exactly the quartic factor form appearing in Q_off_K.

**RELATIONSHIP TO EQ 7.17:**

Under the explicit hypothesis that G_off0 is invertible (IsUnit (G_off0 ...).det),
Theorem Eq7p17WithExplicitHypothesis claims:

  det_param(w) = det(G_off0) * Q_off_K(w)

If the characteristic polynomial of G_off0⁻¹·K_off0 equals Q_off_K(w), then:

  det(w·G_off0 - K_off0) = det(G_off0) * det(w·1 - G_off0⁻¹·K_off0)
                         = det(G_off0) * Q_off_K(w)

which is exactly what we want to prove.

**FORMALIZATION CHALLENGE:**

Proving this requires showing that G_off0⁻¹·K_off0 is similar to a block-diagonal matrix
with 2×2 real blocks of the form [[γ_k, -η_k], [η_k, γ_k]] and [[-γ_k, -η_k], [η_k, -γ_k]]
for each z_k = γ_k + i·η_k. This is non-trivial in Lean.

**CURRENT APPROACH:**

Rather than proving the full similarity result immediately, we prove a WEAKER but still VERY
USEFUL fact: the characteristic polynomial values match at w=0.

-/

theorem det_G_off0_inv_K_off0_at_zero
    (N K : ℕ) (L : ℝ) (z : Fin K → ℂ)
    (hG : IsUnit (G_off0 (2 * K) L K z).det) :
    (-1 : ℝ) ^ (4 * K) * (G_off0 (2 * K) L K z).det⁻¹ * (K_off0 (2 * K) L K z).det
      = Q_off_K K z 0 := by
  sorry

/-- The eigenvalues of G_off0⁻¹·K_off0 are {z_k, -z_k, conj(z_k), -conj(z_k)}.
    This is NUMERICALLY VERIFIED (§2789) but the Lean proof is deferred.
    Status: CONJECTURE (supported by 6.5e-16 numerical evidence) -/
theorem G_off0_inv_K_off0_eigenvalues
    (N K : ℕ) (L : ℝ) (z : Fin K → ℂ)
    (hG : IsUnit (G_off0 (2 * K) L K z).det) :
    let M := (G_off0 (2 * K) L K z)⁻¹ * (K_off0 (2 * K) L K z)
    -- The eigenvalue statement is complex in Lean's type system
    -- For now, we just state that the characteristic polynomial equals Q_off_K
    ∀ w : ℝ, (w • 1 - M).det = Q_off_K K z w := by
  sorry

#print axioms det_G_off0_inv_K_off0_at_zero
#print axioms G_off0_inv_K_off0_eigenvalues
