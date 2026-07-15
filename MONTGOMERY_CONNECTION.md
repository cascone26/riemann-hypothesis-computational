# Montgomery Pair Correlation and sc ~ L^{-2}

## The Rate L^{-1.9}

From our computations:
  sc ~ 3.073e-3 / L^{1.9151}

The exponent 1.9151 ≈ 2 to within 4%. This is remarkably close to L^{-2}.

## Montgomery's Pair Correlation Hypothesis

Montgomery (1973) conjectured:
  Σ_{γ,γ'} f(γ-γ') ~ (T/(2π)) ∫ f(u) |F(u)|² du

where F(u) = 1 - sin²(πu)/(πu)² is the GUE pair correlation.

**Key insight**: The sub-cancel metric sc measures ||V_null^T (A+B-P) V_null||_F / ||V_null^T (A+B) V_null||_F.

The numerator is approximately ||V_null^T P V_null - V_null^T (A+B) V_null||_F,
where by the Weil explicit formula this equals:
  ||V_null^T (Σ_ρ |ĝ(ρ)|²) V_null||_F ≈ Σ_{γ,γ'} |V_null^T M_{γ,γ'} V_null|

where M_{γ,γ'} is the outer product of mode-space projections onto zeros γ and γ'.

If the zeros are uncorrelated (Poisson): Σ_{γ≠γ'} terms add incoherently → O(1/L).
If the zeros follow GUE pair correlation (Montgomery): diagonal terms dominate,
off-diagonal terms partially cancel → net O(1/L²).

**The extra factor of L^{-1}** in sc ~ L^{-2} vs L^{-1} may come from:
  1. Montgomery's pair correlation giving destructive interference among zero-pair terms
  2. The "1 - sin²(πu)/(πu)²" factor creating exact cancellation at u=0

## Test Prediction

Under the Montgomery hypothesis:
  sc(L) ~ C × L^{-2} × (correction factors from diagonal terms)

Under Poisson (zeros uncorrelated):
  sc(L) ~ C × L^{-1} × (no pair-correlation enhancement)

If we saw sc ~ L^{-2.0} exactly, this would be strong evidence that:
1. Zeros follow GUE pair correlation (Montgomery's hypothesis)
2. RH holds (implied by Connes)

Our measurement of L^{-1.9151} is closer to L^{-2} than to L^{-1}, suggesting
pair correlation is at play.

## The Critical Calculation

The sub-cancel sc in terms of the Weil explicit formula:

sc² ≈ (1/||(A+B)||²) × Σ_{k,l} |<ψ_k, Σ_ρ|ĝ(ρ)|² ψ_l>|²

where ψ_k are the near-null eigenvectors of τ.

If zeros are on the critical line (RH), the Riemann-Siegel formula gives:
  Z(t) = Σ_n (π² n² t²)^{-1/4} cos(θ(t) - ... )

and the modes ψ_k ≈ sinc functions near n=0 give:
  <ψ_k, M_{γ,γ'} ψ_l> ~ sinc(n_k - γL/2π) × sinc(n_l - γ'L/2π)

The cross terms Σ_{γ≠γ'} include Montgomery's pair correlation function:
  Σ_{0 < γ < γ', γ-γ'=t} 1 ~ (T/2π) (1 - sin²(πt)/(πt)²)

The "hole" at t=0 (GUE repulsion) reduces the near-diagonal contributions by a factor
proportional to t², and integrating gives an extra 1/L factor relative to uncorrelated zeros.

## Status: Speculative

This is a theoretical conjecture, not proven. The connection to Montgomery's hypothesis
would be a major theoretical result. Verification would require:
1. Computing sc at much larger L (L=20-50) to see L^{-2.0} vs L^{-1.9} precisely
2. Computing the pair-correlation function directly from our near-null eigenvectors
3. Comparing to the Montgomery-Odlyzko law numerically

## What This Means if True

If sc ~ L^{-2} follows from Montgomery's pair correlation:
  - The sub-cancel convergence would be a consequence of zero statistics (GUE)
  - NOT specifically RH (since Montgomery only CONDITIONALLY implies GUE under RH)
  - It would prove: the zeros seen in our bandwidth follow GUE statistics
  - It would NOT prove: all zeros are on the critical line (only those in bandwidth)

**This is a significant finding regardless**: the rate L^{-2} vs L^{-1} encodes
the pair correlation structure of the Riemann zeros.
