import ShiftOverlapParity

/-!
# §1739 — The prime piece is parity-block-diagonal (any weights)

The prime contribution to the Weil form is `Q^prime_kl = −∑ₙ Λ(n)/√n · S_kl(log n)` (verified on the exact
Suzuki form, 2026-08-24: coefficient exactly `Λ(pᵏ)/p^{k/2}`).  Each shift-overlap `S_kl(log n)` vanishes for
`k+l` odd by §1733.  Therefore the ENTIRE prime-weighted sum vanishes for `k+l` odd — **regardless of the
weights** `Λ(n)/√n` and regardless of which prime powers appear:

`∑ₙ (weight n)·S_kl(shift n) = 0`   whenever `k+l` is odd.

So the prime piece carries the SAME even/odd parity block-decoupling as the archimedean kernel: the full Weil
form `arch − prime` is parity-block-diagonal, and its Hodge signature splits over the two parity sectors.  This
is the structural reason the RH-relevant positivity can be analysed sector by sector.  EXIT0 / std-axioms.
-/

open Real intervalIntegral

/-- **Prime-piece parity decoupling.**  For `k+l` odd, any weighted sum of symmetrised shift-overlaps
(the off-diagonal prime piece) vanishes. -/
theorem prime_piece_parity_decouple (a b L : ℝ) (k l : ℤ)
    (hd : a - b ≠ 0) (hs : a + b ≠ 0)
    (hLa : a * L = (k : ℝ) * π) (hLb : b * L = (l : ℝ) * π) (hodd : Odd (k + l))
    (N : Finset ℕ) (weight shift : ℕ → ℝ) :
    (∑ n ∈ N, weight n *
      ((∫ u in (0:ℝ)..(L - shift n), Real.sin (a * (u + shift n)) * Real.sin (b * u))
        + (∫ u in (0:ℝ)..(L - shift n), Real.sin (b * (u + shift n)) * Real.sin (a * u)))) = 0 := by
  apply Finset.sum_eq_zero
  intro n _
  rw [shift_overlap_parity_decouple a b (shift n) L k l hd hs hLa hLb hodd]
  ring
