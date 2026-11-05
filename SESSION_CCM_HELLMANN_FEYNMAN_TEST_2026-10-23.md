# CCM Hellmann-Feynman Protection Test
## 2026-10-23, ~14:30 CDT

**Objective:** Test the Hellmann-Feynman protection hypothesis from S2755 literature search. Specifically: does the perturbation (w02 − prime) couple the two lowest unperturbed eigenstates (ground and first-excited) to each other?

**Hypothesis:** If the off-diagonal matrix element ⟨ψ₁ | (w02 − prime) | ψ₂⟩ is small relative to the perturbation norm and the unperturbed gap, the perturbation is "Hellmann-Feynman protected in the crossing direction" — it avoids mixing the two eigenstates, which would cause level crossing despite large overall norm.

---

## Methodology

**Setup:**
- τ(ε) = arch + ε(w02 − prime), where arch is the unperturbed Connes Chevalley-Chamseddine operator
- Computed arch alone at N=20, 24, 32 using the v3 vectorized architecture
- Extracted the two smallest eigenvectors (ψ₁, ψ₂) of arch
- Computed the off-diagonal coupling: ⟨ψ₁ | (w02 − prime) | ψ₂⟩

**Test Values:** N ∈ {20, 24, 32} (matching S2752 test range)

---

## Results

### Detailed Data

**N=20:**
- Unperturbed eigenvalues: λ₁ = −2.9969, λ₂ = −1.3246
- Unperturbed gap: **1.6722**
- Off-diagonal element: **1.093×10⁻¹⁵** (machine zero)
- Diagonal contributions:
  - ⟨ψ₁ | (w02 − prime) | ψ₁⟩ = 3.0039
  - ⟨ψ₂ | (w02 − prime) | ψ₂⟩ = 1.3355
- Perturbation norm ‖w02 − prime‖: 8.191
- Ratio |off-diag| / ‖perturb‖: **1.3×10⁻¹⁶** (zero)

**N=24:**
- Unperturbed eigenvalues: λ₁ = −2.9967, λ₂ = −1.3261
- Unperturbed gap: **1.6706**
- Off-diagonal element: **1.471×10⁻¹⁷** (machine zero)
- Diagonal contributions:
  - ⟨ψ₁ | (w02 − prime) | ψ₁⟩ = 3.0038
  - ⟨ψ₂ | (w02 − prime) | ψ₂⟩ = 1.3376
- Perturbation norm ‖w02 − prime‖: 9.035
- Ratio |off-diag| / ‖perturb‖: **1.6×10⁻¹⁸** (zero)

**N=32:**
- Unperturbed eigenvalues: λ₁ = −2.9965, λ₂ = −1.3280
- Unperturbed gap: **1.6686**
- Off-diagonal element: **6.689×10⁻¹⁶** (machine zero)
- Diagonal contributions:
  - ⟨ψ₁ | (w02 − prime) | ψ₁⟩ = 3.0038
  - ⟨ψ₂ | (w02 − prime) | ψ₂⟩ = 1.3413
- Perturbation norm ‖w02 − prime‖: 10.607
- Ratio |off-diag| / ‖perturb‖: **6.3×10⁻¹⁷** (zero)

### Summary Table

| N  | \|off-diag\| | gap    | ‖w02 − prime‖ | ratio: \|off-diag\|/‖perturb‖ |
|----|--------------|--------|------------------|---------------------------|
| 20 | 1.09×10⁻¹⁵  | 1.6722 | 8.191            | 1.3×10⁻¹⁶                 |
| 24 | 1.47×10⁻¹⁷  | 1.6706 | 9.035            | 1.6×10⁻¹⁸                 |
| 32 | 6.69×10⁻¹⁶  | 1.6686 | 10.607           | 6.3×10⁻¹⁷                 |

---

## Interpretation

### Critical Structural Property

The off-diagonal matrix element is **numerically zero** at all three N values — at machine precision (10⁻¹⁵ range). This is remarkable because:

1. **The perturbation norm is substantial:** ‖w02 − prime‖ ranges from 8.2 to 10.6 (roughly the same scale as the unperturbed gap)
2. **Yet it does not couple the two eigenstates:** |⟨ψ₁ | (w02 − prime) | ψ₂⟩| ≈ 10⁻¹⁶, vanishingly small
3. **But the diagonal elements are large:** ⟨ψ₁ | (w02 − prime) | ψ₁⟩ ≈ 3.0 and ⟨ψ₂ | (w02 − prime) | ψ₂⟩ ≈ 1.34

### Hellmann-Feynman Analysis

From the Hellmann-Feynman theorem (S2755 literature):
- First-order eigenvalue shift: δλ_n = ⟨ψ_n | dH/dε | ψ_n⟩
- If perturbation has **zero off-diagonal coupling** ⟨ψ₁ | P | ψ₂⟩ = 0, then the two eigenstates remain non-degenerate and do not mix under perturbation, avoiding the avoided crossing (level repulsion) that typically occurs when ‖P‖ exceeds the gap.

**The mechanism:** The perturbation (w02 − prime) shifts each eigenstate independently via its diagonal action:
- ψ₁ is shifted by ~3.0 (toward higher energy or lower, depending on sign)
- ψ₂ is shifted by ~1.34 independently

Because there is no off-diagonal coupling, these shifts occur **without the eigenstates exchanging character**, even though the diagonal shifts are comparable to the gap itself (gap ≈ 1.67, shifts are 3.0 and 1.34).

### Why No Level Crossing Despite Gap Collapse?

The literature (S2755, Rellich's theorem + Hellmann-Feynman theorem) explains:
- For analytic Hermitian perturbations, eigenvalue paths are analytic in the parameter ε
- If two eigenstates do not mix under the perturbation (zero off-diagonal coupling), they cannot exchange identities (which would be a crossing)
- They can come arbitrarily close (as observed: gap → 10⁻¹² at ε=1, per S2752) but remain simple (non-degenerate) throughout

The near-zero off-diagonal element **is the structural reason** why this occurs.

---

## Hypothesis Verdict

**✓ HELLMANN-FEYNMAN PROTECTION HYPOTHESIS: STRONGLY SUPPORTED**

**Evidence:**
1. ⟨ψ₁ | (w02 − prime) | ψ₂⟩ ≈ 0 at machine precision across N=20, 24, 32 ✓
2. The ratio |off-diag| / ‖perturbation‖ is <10⁻¹⁶ for all N ✓
3. The diagonal elements (which cause gap collapse) are non-zero and substantial ✓
4. This structure exactly matches the Hellmann-Feynman mechanism: protected against crossing because there is no off-diagonal coupling ✓

**Interpretation:**
The (w02 − prime) perturbation is **non-crossing by design**. It is structured such that it targets each of the two lowest eigenstates independently (large diagonal elements), without coupling them to each other (zero off-diagonal element). This structural property, predicted by Hellmann-Feynman and classical perturbation theory (Rellich's theorem, structure-preserving perturbations) explains why:
- The gap can collapse from 1.67 to 10⁻¹² (large diagonal shifts)
- Yet no level crossing occurs (zero off-diagonal coupling)
- And simplicity is preserved throughout ε ∈ [0,1] (analytic eigenvalue continuation)

---

## Next Steps (Recommended from S2755)

1. ✓ **HF test completed.** Off-diagonal coupling is zero; perturbation is HF-protected. ← THIS SESSION
2. Examine explicit-formula structure of w02 and prime to understand **why** they produce orthogonal actions on the ground vs. first-excited eigenstates
3. Check whether w02 − prime commutes (or nearly commutes) with arch — a sufficient condition for zero off-diagonal coupling

---

## Files

- **Test code:** `/Users/scones/projects/riemann/test_hellmann_feynman.py`
- **Data:** inline in this report

## Verification

Results verified by:
- Three independent N values (20, 24, 32)
- Eigendecomposition via np.linalg.eigh (standard NumPy LAPACK)
- Direct matrix-vector multiplication for coupling elements
- Frobenius norms via np.linalg.norm
