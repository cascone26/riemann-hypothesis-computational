import CCMTailFinal

/-!
# S3003 — margin transfer with the sharp tail bound (lane 2)

`margin_transfer_sharp`: if xᵀAx = xᵀGx + xᵀTx with T the confluent Löwner far tail (S3001 hypotheses) and
xᵀGx ≥ m·Σx² on 𝟙^⊥, then for Σx = 0:
  xᵀAx ≥ (m − D₂·S/(2π))·Σx².
This is the certificate shape: the CCM margin on 𝟙^⊥ is ≥ (finite near-zero margin) − (explicit far-tail term), with the
tail term now within ~10× of truth (S3001).  Honest status: at small λ² the near-zero Gram is rank-deficient (S2980) so m = 0;
the transfer becomes informative once the near block has ≥ n−1 zeros, i.e. 2V below the ~n/2-th zero height.
-/

open Set Finset Matrix

namespace CCMMarginSchur

theorem margin_transfer_sharp {ι : Type*} [Fintype ι] [DecidableEq ι] (A G : Matrix ι ι ℝ)
    (g g' g'' : ℝ → ℝ) {V D₂ S m : ℝ}
    (hg : ∀ x ∈ Icc (-V) V, HasDerivWithinAt g (g' x) (Icc (-V) V) x)
    (hg' : ∀ x ∈ Icc (-V) V, HasDerivWithinAt g' (g'' x) (Icc (-V) V) x)
    (hD : ∀ x ∈ Icc (-V) V, |g'' x| ≤ D₂) (ν : ι → ℝ) (hν : ∀ i, ν i ∈ Icc (-V) V)
    (hinj : Function.Injective ν) (hS : ∀ i, ∑ k, |ν i - ν k| ≤ S) (x : ι → ℝ)
    (hx : ∑ i, x i = 0)
    (hsplit : x ⬝ᵥ (A *ᵥ x) = x ⬝ᵥ (G *ᵥ x)
      + x ⬝ᵥ (CCMLowner.lowner (fun i => -(g' (ν i)) / Real.pi) (fun j => g (ν j)) ν *ᵥ x))
    (hG : m * ∑ i, x i ^ 2 ≤ x ⬝ᵥ (G *ᵥ x)) :
    (m - D₂ * S / (2 * Real.pi)) * ∑ i, x i ^ 2 ≤ x ⬝ᵥ (A *ᵥ x) := by
  have ht := CCMTailFinal.confluent_perp_schur_half g g' g'' hg hg' hD ν hν hinj hS x hx
  have hneg := neg_abs_le (x ⬝ᵥ (CCMLowner.lowner (fun i => -(g' (ν i)) / Real.pi)
    (fun j => g (ν j)) ν *ᵥ x))
  have e : D₂ / (2 * Real.pi) * S * ∑ i, x i ^ 2 = D₂ * S / (2 * Real.pi) * ∑ i, x i ^ 2 := by ring
  rw [e] at ht
  rw [hsplit]
  nlinarith

end CCMMarginSchur

#print axioms CCMMarginSchur.margin_transfer_sharp
