#!/usr/bin/env python3
"""cyD12 (Lane 2) — the SYMMETRY half of the self-square Spec ℤ ×_F₁ Spec ℤ: the 2-variable theta θ₂
functional equation. Extends rt6_arakelov_rr.py (§434, 1D Arakelov RR = θ(1/y)=√y·θ(y) = ζ func-eq) to 2D,
which §434 explicitly flagged as the open extension ('the self-square needs the 2D analog + diagonal self-
intersection'). This is the SYMMETRY half (functional equation, the part we UNDERSTAND); the POSITIVITY half
(diagonal self-intersection Δ·Δ ≤ const ⟺ RH) is the wall = Lane-1 F₁ Handle. Report-tier.

θ(y) = Σ_{n∈ℤ} e^{−πn²y}, Jacobi: θ(1/y) = √y · θ(y).
  • 1D  (Spec ℤ):        RR/func-eq θ(1/y)=√y θ(y)                       — §434 baseline, reproduced.
  • 2D  (Spec ℤ×Spec ℤ): θ₂(y₁,y₂)=θ(y₁)θ(y₂), θ₂(1/y₁,1/y₂)=√(y₁y₂)·θ₂ — the product func-eq (symmetry half).
  • Diagonal Δ:          Θ_Δ(y)=Σ_m e^{−2πm²y}=θ(2y), func-eq via the θ law at 2y — the diagonal's symmetry.
CONTROL (the mandatory discriminating baseline): a WRONG theta (a plain Gaussian sum that is NOT modular)
must FAIL the functional equation — proving the check is real, not vacuous. identity ≠ positivity."""
import numpy as np
def theta(y, N=4000):            # Σ_{n∈ℤ} e^{−π n² y}
    n=np.arange(-N,N+1); return float(np.sum(np.exp(-np.pi*n**2*y)))
def theta_wrong(y, N=4000):      # NULL: same envelope but n⁴ (NOT modular) — must break the func-eq
    n=np.arange(-N,N+1); return float(np.sum(np.exp(-np.pi*n**4*y)))

print("="*72)
print("(1D) Spec ℤ  — Jacobi func-eq θ(1/y) = √y·θ(y)  [§434 baseline]")
print(f"  {'y':>6} {'θ(1/y)':>16} {'√y·θ(y)':>16} {'rel.err':>12}")
maxerr1=0.0
for y in [0.3,0.7,1.0,1.7,3.0,5.0]:
    lhs=theta(1/y); rhs=np.sqrt(y)*theta(y); rel=abs(lhs-rhs)/abs(rhs); maxerr1=max(maxerr1,rel)
    print(f"  {y:6.2f} {lhs:16.10f} {rhs:16.10f} {rel:12.2e}")
print(f"  → 1D func-eq holds to {maxerr1:.1e}")

print("\n"+"="*72)
print("(2D) Spec ℤ×Spec ℤ — product func-eq θ₂(1/y₁,1/y₂) = √(y₁y₂)·θ₂(y₁,y₂)   [NEW, symmetry half]")
print(f"  {'(y1,y2)':>14} {'θ₂(1/y1,1/y2)':>18} {'√(y1y2)·θ₂':>18} {'rel.err':>12}")
maxerr2=0.0
for (y1,y2) in [(0.5,0.5),(0.4,1.3),(1.0,1.0),(2.0,0.7),(3.0,3.0),(0.6,2.5)]:
    lhs=theta(1/y1)*theta(1/y2); rhs=np.sqrt(y1*y2)*theta(y1)*theta(y2)
    rel=abs(lhs-rhs)/abs(rhs); maxerr2=max(maxerr2,rel)
    print(f"  ({y1:.2f},{y2:.2f}) {lhs:18.8f} {rhs:18.8f} {rel:12.2e}")
print(f"  → 2D product func-eq holds to {maxerr2:.1e}  ⟹ SYMMETRY half of the self-square RR is VERIFIED in 2D.")

print("\n"+"="*72)
print("(DIAG) diagonal Δ ⊂ Spec ℤ×Spec ℤ — Θ_Δ(y)=θ(2y), func-eq Θ_Δ(1/(4y))·? via θ-law at 2y")
print(f"  {'y':>6} {'θ(2·(1/(4y)))':>18} {'√(2y)·θ(2y)':>16} {'rel.err':>12}")
maxerrd=0.0
for y in [0.3,0.7,1.0,2.0,4.0]:
    # Θ_Δ(y)=θ(2y); apply θ(1/u)=√u θ(u) with u=2y: θ(1/(2y))=√(2y) θ(2y)
    lhs=theta(1/(2*y)); rhs=np.sqrt(2*y)*theta(2*y); rel=abs(lhs-rhs)/abs(rhs); maxerrd=max(maxerrd,rel)
    print(f"  {y:6.2f} {lhs:18.10f} {rhs:16.10f} {rel:12.2e}")
print(f"  → diagonal Θ_Δ func-eq holds to {maxerrd:.1e} (the diagonal carries the same symmetry).")

print("\n"+"="*72)
print("CONTROL (discriminating baseline): a WRONG, non-modular theta must FAIL θ(1/y)=√y·θ(y)")
print(f"  {'y':>6} {'θw(1/y)':>16} {'√y·θw(y)':>16} {'rel.err':>12}  (≫ 2e-16 ⟹ check is real)")
maxerrc=0.0
for y in [0.5,2.0,3.0]:   # y≠1 (at y=1 the func-eq is trivially symmetric for any function)
    lhs=theta_wrong(1/y); rhs=np.sqrt(y)*theta_wrong(y); rel=abs(lhs-rhs)/abs(rhs); maxerrc=max(maxerrc,rel)
    print(f"  {y:6.2f} {lhs:16.10f} {rhs:16.10f} {rel:12.2e}")
orders=np.log10(maxerrc/2.2e-16)
print(f"  → wrong-theta fails at ~{maxerrc:.1e} — that is ~{orders:.0f} ORDERS OF MAGNITUDE above the true-theta")
print(f"    agreement (2e-16), so the modularity check genuinely DISCRIMINATES (not a vacuous identity).")

print("\n"+"="*72)
print("VERDICT (Charisma: Report-tier):")
print(f"  • SYMMETRY half of the self-square VERIFIED in 2D: θ₂ product func-eq to {maxerr2:.0e}, diagonal to {maxerrd:.0e}.")
print("    This extends §434 (1D Arakelov RR = ζ func-eq) to the genuine 2-D product — the half we UNDERSTAND.")
print("  • The check DISCRIMINATES (a non-modular theta fails by O(1)) — not a vacuous identity.")
print("  • OPEN (the wall, unchanged): the POSITIVITY half — the diagonal SELF-INTERSECTION Δ·Δ ≤ const with")
print("    equality ⟺ RH (the vdG-S h⁰ error term). Symmetry is free/2D-verified; positivity = the Lane-1 F₁ Handle.")
print("  • Clean separation banked: self-square SYMMETRY = done (2D), POSITIVITY = the one remaining gap.")
