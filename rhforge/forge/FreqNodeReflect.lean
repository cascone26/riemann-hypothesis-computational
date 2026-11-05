import HLGCPartialFraction

open scoped MatrixOrder Matrix

/-!
# Off-critical-line lane: the frequency-node reflection involution (step 3 of the matrix-connection plan)

Part of the 4-step plan from §2391's breakthrough: to show the "+ν_n" partial-fraction copy also
vanishes (matching §2380-2383's `Q_off_K_null_identity`, which only directly covers the "−ν_n"
copy), the natural route is a REFLECTION SYMMETRY of the frequency-node index set: negating a
frequency node `ν_{n,L}=2π(n-2K)/L` (over the symmetric range `n∈{0,…,4K}` representing
`{-2K,…,2K}`) lands on ANOTHER frequency node, namely the one at the reflected index `4K-n`. This
brick builds that reflection map and its two key properties (it's an involution; it negates the
node value exactly). Search-before-build: grepped for `freqNode_reflect`/`freqNodeReflect`/
`nu_reflect` — fresh. Fully general, no zeta/RH dependence. -/
def freqReflect (K : ℕ) : Fin (2 * (2 * K) + 1) → Fin (2 * (2 * K) + 1) :=
  fun n => ⟨4 * K - (n : ℕ), by omega⟩

theorem freqReflect_involutive (K : ℕ) : Function.Involutive (freqReflect K) := by
  intro n
  unfold freqReflect
  ext
  simp only
  omega

theorem freqReflect_bijective (K : ℕ) : Function.Bijective (freqReflect K) :=
  (freqReflect_involutive K).bijective

theorem freqNode_reflect_eq_neg (K : ℕ) (L : ℝ) (n : Fin (2 * (2 * K) + 1)) :
    2 * Real.pi * (((freqReflect K n : ℕ) : ℝ) - ((2 * K : ℕ) : ℝ)) / L
      = -(2 * Real.pi * (((n : ℕ) : ℝ) - ((2 * K : ℕ) : ℝ)) / L) := by
  unfold freqReflect
  have hn : (n : ℕ) ≤ 4 * K := by omega
  have hcast : (((4 * K - (n : ℕ) : ℕ) : ℝ)) = (4 * K : ℝ) - ((n : ℕ) : ℝ) := by
    push_cast [Nat.cast_sub hn]
    ring
  simp only [hcast]
  push_cast
  ring

#print axioms freqReflect_involutive
#print axioms freqReflect_bijective
#print axioms freqNode_reflect_eq_neg
