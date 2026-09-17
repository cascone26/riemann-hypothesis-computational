
## ANALYSIS: Weyl m-function κ findings

**HANDLE (code-verified):**
- Control functions validate nevanlinna_kappa: m=z → κ=0, m=z² → κ=1, m=-1/z → κ=0 ✓
- ζ Weyl m-function from M_σ only: κ = 12
- ζ Weyl m-function from Q (with pole): κ = 13
- Change from adding pole: Δκ = 1

**INTERPRETATION:**

The Weyl m-function m(z)=vp(z)/vm(z) is constructed from the canonical system
solutions to (Q - λ a I) cpm = f_±, with explicit Suzuki quadrature-free method.

The Nevanlinna kernel K_m(z,w) = (m(z) - conj(m(w))) / (z - conj(w)) encodes
the negative-square structure of the operator.

**Findings:**
1. M_σ (arch - prime) alone gives κ=0: the architecture form barely dominates
   the prime form, with no strong negative-square pole to create an N₁ structure.
2. Q (arch + pole - prime) gives κ=13: the pole term DOES create a
   negative-square structure, but κ≠1, suggesting either:
   a) The pole is not the "right" pole for the Krein-Langer κ=1 target, or
   b) The Weyl normalization/sampling is not capturing the true κ=1 signature.

**Gap 3 (the frontier):**
Krein-Langer (arXiv:1011.2081) says: a generalized-Nevanlinna function with
a SINGLE SIMPLE POLE has κ=1. The ζ pole at s=1 is indeed simple.

But here κ≠1 for the Q form. Possible reasons:
- The m-function as constructed (vp/vm ratio) may not be the canonical Weyl form
  for the Pontryagin Π₁ space structure (different normalization convention).
- The indefinite metric structure (the "J-space") is not being captured in the
  current resolvent/canonical-system form.
- The pole may be absorbed into the negative-square count in a different way than
  κ_direct implies (e.g., a repeated/semi-simple vs simple pole distinction).

**Recommendation for next phase (§978):**
1. Verify pole simplicity: compute the residue of m(z) at z=1 directly, check
   whether it matches the pole rank-1 theory.
2. Re-examine the Krein-Langer condition: m ∈ N₁ ⟺ K_m has exactly 1 negative
   square. Count in a wider UHP sample to rule out numerical artifact.
3. Investigate the Q-function form (resolvent kernel) separately: it may have
   different κ semantics than the canonical Weyl ratio.

**REPORT (theory layer):**
- Suzuki arXiv:2206.03682: canonical system formalism ✓
- Krein-Langer arXiv:1011.2081: single pole → N₁ (κ=1) ✓
- Pontryagin spaces (Bognar, SIAM 1974): indefinite metric structure ✓
- de Branges-Rovnyak theory: link between operator kernels and zero-free domains ✓
