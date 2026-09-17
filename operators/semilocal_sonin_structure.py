#!/usr/bin/env python3
"""§cy371: Semilocal Sonin projection structure — does C-C's archimedean bound survive with primes?

RESEARCH QUESTION: Connes §1228 gives the SEMILOCAL trace formula
  -∑_{v∈S} W_v(f) = log(TW)f(1) + Tr(ϑ(f)(1−P_T^S−P̂_W^S))
where S is a finite set of places (arch + finitely many primes).

The archimedean case (cy370, formalized) shows:
  W_∞ ≥ Tr(ϑ(g) 𝔖 ϑ(g)^H) ≥ 0   via Sonin projection 𝔖 on narrow support.

QUESTION: Does a similar projection-lower-bound survive for SEMILOCAL W?
  That is: can we write
    ∑_{p∈S} W_p(f) = Tr(ϑ(f) (semilocal-projection) ϑ(f)^H) + boundary
  where (semilocal-projection) is manifestly nonnegative?

HONEST ANSWER (this cycle):
  The semilocal formula DOES have projection structure (1−P_T^S−P̂_W^S), BUT:
  • The projections P_T^S, P̂_W^S are defined on the SEMILOCAL ADELE CLASS SPACE
    Y_S = A_S / Γ_S (product of local fields), not on individual primes.
  • The positivity structure (projection = PSD) holds at EACH PLACE v ∈ S
    INDEPENDENTLY — it is NOT a single global projection, but a PRODUCT of local projections.
  • CENTRAL FACT (Connes §1236-1261): the module map |·|_S is MULTIPLICATIVE across
    primes: Mod_S(u) = ∏_{p∈S} |u_p| is a group homomorphism R^*_+ × K_S → R^*_+.
  • CONSEQUENCE: the projection (1−P_T^S−P̂_W^S) on Y_S FACTORS across places —
    each prime p contributes a LOCAL projection P_p^locally that is also an orthogonal
    projection on L²(Q_p/Z_p^×) (via the same interval projection trick on |·|_p).

PROBE: numerically verify that
  1. archimedean W_∞ ≥ Tr(ϑ 𝔖 ϑ^H) holds (known)
  2. for EACH PRIME p in S, the local W_p ALSO admits a Sonin-type lower bound
  3. the sum structure is PRESERVED: W_S = W_∞ + ∑_{p∈S} W_p still admits a global bound
  4. NO CANCELLATION happens between arch and primes that would break positivity.

EXPECTED FINDING:
  If each W_p is individually nonneg (via local Sonin on Q_p), then:
    ∑_{p∈S} W_p ≥ ∑_p Tr(local_{p}) ≥ 0   (sum of nonneg = nonneg)
  ⟹ semilocal Sonin bound EXTENDS cleanly (no new math needed; just bookkeeping).

WALL (if reached):
  If one prime p ∈ S has a NEGATIVE eigenmode (not caught by local Sonin), then
  the global bound BREAKS, and we need NEW structure (deficiency indices? indefinite Sonin?).
"""

import numpy as np
import sys
sys.path.insert(0, '/Users/scones/projects/riemann/operators')

# For now: document the STRUCTURE, not compute (data not readily accessible).
# The actual test requires:
#   (a) the local-field measure dμ_p on each Q_p (from ccm_discrepancy.py but extended to p-local)
#   (b) the Sonin projection on L²(local field) (requires Fourier on Q_p / Bruhat-Schwartz space)
#   (c) the scaling operator ϑ on each local field
# These are noncommutative-geometry objects; not a finite-matrix problem.

print("""
HONEST STATUS (Report-tier):

C-C §1228 FORMULA (semilocal trace formula) is PROVEN by Connes (Co-zeta 1999).
It takes the form:
  -∑_{v∈S} W_v(f) = log(TW)f(1) + Tr(ϑ(f)(1−P_T^S−P̂_W^S))

where P_T^S, P̂_W^S are projections on the semilocal adele class space Y_S.

KEY STRUCTURAL OBSERVATION (§1236-1261):
  • The semilocal space Y_S = A_S / Γ_S is a PRODUCT of local fields
  • The module Mod_S is MULTIPLICATIVE: |λ|_S = ∏_{p∈S} |λ_p|
  • Each LOCAL place p ∈ S has its own ORTHOGONAL PROJECTION P_p on L²(Q_p)
  • The GLOBAL projection (1−P_T^S−P̂_W^S) is the PRODUCT of local projections

CONSEQUENCE for POSITIVITY:
  Since each LOCAL projection P_p satisfies:
    - P_p is an orthogonal projection (self-adjoint idempotent) on L²(Q_p)
    - Therefore P_p is PSD
    - Conjugation by ϑ_p (local scaling) preserves PSD
    - Therefore Tr(ϑ_p P_p ϑ_p^*) ≥ 0 for EACH p

  And since Mod_S is MULTIPLICATIVE:
    ∑_{p∈S} W_p(f) can be written as a SUM of local contributions,
    each of which admits a Sonin-type lower bound

  The GLOBAL FORM ∑_p Tr(local_p) has NO CANCELLATION ⟹ remains PSD.

FORMAL CLAIM (to be verified):
  The SEMILOCAL W admits a manifestly-nonnegative lower bound of the SAME SHAPE
  as the archimedean case, by taking the SUM of local Sonin bounds.

  ∑_{p∈S} W_p(f) ≥ ∑_{p∈S} Tr(ϑ_p 𝔖_p ϑ_p^*) ≥ 0

  where 𝔖_p is the Sonin-like projection on L²(Q_p).

OPEN QUESTION (where the research happens):
  The C-C semilocal paper does NOT EXPLICITLY CONSTRUCT the local Sonin 𝔖_p
  on each Q_p. The ARCHIMEDEAN 𝔖 is well-understood (Burnol, Connes §1195).
  The PRIME SONIN spaces may have a DIFFERENT STRUCTURE due to the p-adic
  topology (e.g., the interval [-λ,λ] ⊂ ℝ has no direct p-adic analog).

  NEEDED: a LOCAL SONIN construction on p-adic fields, and a PROOF that it
  gives a nonnegative lower bound analogous to the archimedean case.

NEXT BRICK (if research confirms):
  Formalize the semilocal decomposition:
    theorem semilocal_weil_nonneg_of_local_sonin_bounds
      (S : Finset Prime) (f : test_function)
      (h_arch : arch_weil_nonneg_of_sonin_bound V_∞ S_∞ f)
      (h_primes : ∀ p ∈ S, local_weil_p_nonneg_of_sonin_bound V_p S_p f) :
      semilocal_weil_nonneg S f

  This would extend cy370 from arch-only to arch+finitely-many-primes.
  The CRUX is proving h_primes for each prime (requires p-adic analysis).
""")

print("\n=== ANTI-LOOP CHECK ===")
print("Probe status: REPORT-tier only (structural analysis, no new computation).")
print("No new computation = no new loop risk; we're reading Connes + documenting what's needed.")
print("\nTo push to HANDLE-tier (formalized brick):")
print("  1. Prove (or cite from literature) that Sonin space exists on Q_p")
print("  2. Establish local_weil_p ≥ Tr(ϑ_p 𝔖_p ϑ_p^*) for each p")
print("  3. Formalize the sum structure (no cancellation)")
print("  4. Bank as semilocal_sonin_bound.lean")
print("\nCRITICAL: if a prime p breaks the projection structure (e.g., deficiency indices ≠ 0),")
print("then the simple sum bound DOES NOT WORK, and we need INDEFINITE projection theory.")
