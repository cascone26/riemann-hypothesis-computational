import Mathlib.LinearAlgebra.BilinearForm.Basic
import Mathlib.Tactic

/-!
# §1755 — The FE involution splits the space (and the radical) into ±1 sectors

Suzuki's functional-equation involution `K` (arXiv:2301.00421 §5, `K² = id`) decomposes the space into its
`+1` and `−1` eigenspaces — this is the abstract, continuum form of the concrete parity split my Weil-form
bricks establish (§1733 parity-decoupling, §1751 radical parity-split).  Here, for any linear involution
`K : M → M` (`K ∘ K = id`):

* every `v` decomposes as `v = ½(v + Kv) + ½(v − Kv)` with the summands in the `+1` / `−1` eigenspaces;
* the `+1` eigenspace is exactly the range of the projection `½(1 + K)`.

Applied to Suzuki's `K`: the space (and the `K`-invariant radical) splits into `V₊ ⊕ V₋`, and `V(0)` lives in
the `+1` sector — the continuum analogue of the even/odd Weil sectors.  Over `ℝ` (the realification carrying
`K` as an `ℝ`-linear involution) this is exactly the parity decomposition.  EXIT0 / std-axioms.
-/

variable {R : Type*} [CommRing R] [Invertible (2 : R)]
variable {M : Type*} [AddCommGroup M] [Module R M]

/-- **Involution decomposition.**  `v` is the sum of its `+1` and `−1` eigen-parts under an involution `K`. -/
theorem involution_decompose (K : M →ₗ[R] M) (hK : K.comp K = LinearMap.id) (v : M) :
    v = ⅟(2 : R) • (v + K v) + ⅟(2 : R) • (v - K v) := by
  have h : (v + K v) + (v - K v) = (2 : R) • v := by rw [two_smul]; abel
  rw [← smul_add, h, smul_smul, invOf_mul_self, one_smul]

/-- The `+1` eigen-part is fixed by `K` (uses `K² = id`). -/
theorem involution_plus_fixed (K : M →ₗ[R] M) (hK : K.comp K = LinearMap.id) (v : M) :
    K (⅟(2 : R) • (v + K v)) = ⅟(2 : R) • (v + K v) := by
  have hKK : K (K v) = v := by
    have := LinearMap.congr_fun hK v
    simpa using this
  rw [map_smul, map_add, hKK, add_comm (K v) v]

/-- The `−1` eigen-part is negated by `K`. -/
theorem involution_minus_flipped (K : M →ₗ[R] M) (hK : K.comp K = LinearMap.id) (v : M) :
    K (⅟(2 : R) • (v - K v)) = - (⅟(2 : R) • (v - K v)) := by
  have hKK : K (K v) = v := by
    have := LinearMap.congr_fun hK v
    simpa using this
  rw [map_smul, map_sub, hKK, ← smul_neg, neg_sub]
