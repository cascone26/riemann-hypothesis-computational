import P2KLDerivEven
import EOffInterpVector
import FreqNodeReflect

open scoped MatrixOrder Matrix

/-!
# Off-critical-line lane: `c_off`/`e_off` are invariant under the frequency-node reflection

Step 3-d of §2391's 4-step plan: assembles §2392-2394's `freqReflect`/`freqNode_reflect_eq_neg`
with §2396-2398's `P_2KL_deriv_even` and `QOffKPoly.lean`'s `Q_off_K_even` to show the interpolation
coefficient `c^{off}_{n,2K,L}` (Eq 7.13) takes the SAME value at a frequency-node index and at its
reflection: both the numerator (`Q_off_K`, even) and denominator (`deriv (P_2KL K L)`, even) are
unchanged by negating their argument, and `freqNode_reflect_eq_neg` says reflecting the INDEX
negates the node VALUE — so the whole ratio is unchanged. Search-before-build: grepped
`c_off_freqReflect`/`c_off_reflect`/`e_off_freqReflect` — fresh. Fully general, no zeta/RH
dependence. -/
theorem c_off_freqReflect_invariant (K : ℕ) (L : ℝ) (z : Fin K → ℂ)
    (n : Fin (2 * (2 * K) + 1)) :
    c_off K L z (freqReflect K n) = c_off K L z n := by
  unfold c_off
  rw [freqNode_reflect_eq_neg K L n, Q_off_K_even, P_2KL_deriv_even]

theorem e_off_freqReflect_invariant (K : ℕ) (L : ℝ) (z : Fin K → ℂ)
    (n : Fin (2 * (2 * K) + 1)) :
    e_off K L z (freqReflect K n) = e_off K L z n := by
  unfold e_off
  rw [c_off_freqReflect_invariant]

#print axioms c_off_freqReflect_invariant
#print axioms e_off_freqReflect_invariant
