import QOffKPoly
import Mathlib

open scoped MatrixOrder Matrix

/-!
# Analysis: What are the roots of Q_off_K?

MATHEMATICAL ANALYSIS (PART 2 — FINDINGS)
==========================================

Q_off_K DEFINITION:
  Q_off_K(w) = ∏_k (w⁴ - 2(γ_k² - η_k²)w² + (γ_k² + η_k²)²)

where z_k = γ_k + iη_k (complex parameters with nonzero imaginary part).

Each quartic factor is a BIQUADRATIC:
  f_k(w) = w⁴ - aw² + b = 0  (u = w² substitution)
  
where:
  a = 2(γ_k² - η_k²) = 2·Re(z_k²)
  b = (γ_k² + η_k²)² = |z_k|⁴

Solving u² - au + b = 0:
  u = (a ± √(a² - 4b)) / 2

CLOSED-FORM EVALUATION:
  Δ = a² - 4b = 4(γ_k² - η_k²)² - 4(γ_k² + η_k²)²
              = 4[(γ_k² - η_k²)² - (γ_k² + η_k²)²]
              = 4[((γ_k² - η_k²) - (γ_k² + η_k²))((γ_k² - η_k²) + (γ_k² + η_k²))]
              = 4[(-2η_k²)(2γ_k²)]
              = -16γ_k² η_k²  (negative, so √Δ is imaginary)

  √Δ = ±4iγ_k η_k

Therefore:
  u₁ = (a + 4iγ_k η_k)/2 = γ_k² - η_k² + 2iγ_k η_k = (γ_k + iη_k)² = z_k²
  u₂ = (a - 4iγ_k η_k)/2 = γ_k² - η_k² - 2iγ_k η_k = (γ_k - iη_k)² = (conj z_k)²

ROOTS OF THE QUARTIC:
  If u = w² = z_k², then w = ±z_k
  If u = w² = (conj z_k)², then w = ±(conj z_k)

FINAL RESULT:
  The roots of Q_off_K(w) over ℂ are:
    R = {±z_1, ±(conj z_1), ±z_2, ±(conj z_2), ..., ±z_K, ±(conj z_K)}
  
  Cardinality: |R| = 4K (with multiplicities, since each z_k is distinct and nonreal)

EIGENVALUE INTERPRETATION:
  The characteristic polynomial of the matrix pencil satisfies:
    det(w·1 - G_off0⁻¹·K_off0) = Q_off_K(w)
  
  Therefore, the EIGENVALUES of G_off0⁻¹·K_off0 are exactly these ±z_k, ±conj(z_k).

KEY INSIGHT FOR PROVING THE CHARACTERISTIC POLYNOMIAL IDENTITY:
  To prove det(w·1 - G_off0⁻¹·K_off0) = Q_off_K(w), it suffices to show:
  1. Both sides are polynomials of degree 4K
  2. Both have the same roots (with multiplicities)
  3. Both have the same leading coefficient
  
  The root structure ±z_k, ±conj(z_k) provides the PRIMARY INVARIANT that
  connects the two sides. The proof would leverage this structure via a
  detailed analysis of how the interpolation basis (defining G_off0 and K_off0)
  encodes this eigenvalue structure implicitly.
-/

#check Q_off_K

