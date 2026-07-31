#!/usr/bin/env python3
"""ROSETTA — instrument the PROOF of RH-for-curves (the answer key), not its downstream shadow.
Weil's mechanism:  Hodge index on the surface C x C  ==>  Castelnuovo-Severi:  for a correspondence D
of bidegree (a,b),  D^2 <= 2ab  (the primitive part is intersection-NEGATIVE-definite: the (1,rho-1)
signature).  Applied to Gamma_{F^r} (graph of r-th Frobenius power) and the diagonal Delta:
    N_r := Gamma_{F^r} . Delta = #C(F_{q^r}) = q^r + 1 - sum_i alpha_i^r
    Hodge index / Castelnuovo-Severi  ==>  |sum_i alpha_i^r| = |q^r+1 - N_r| <= 2g * q^{r/2}   (ALL r).
That bound, holding for ALL r, FORCES |alpha_i| = sqrt(q) (any |alpha|>sqrt q makes sum grow like
q^{r(1/2+eps)}, breaking the bound for large r). THAT is the geometric positivity proving RH-for-curves.

We (1) verify the bound on real curves, (2) demonstrate the FORCING by a counterfactual (push one alpha
off the circle -> the Hodge bound breaks), (3) map the analogy to the number field and locate the GAP.
"""
import numpy as np

def alphas_elliptic(a_p, p):
    d = a_p*a_p - 4*p
    al = complex(a_p/2, np.sqrt(-d)/2) if d < 0 else complex((a_p+np.sqrt(d))/2, 0)
    return np.array([al, np.conj(al) if d < 0 else p/al]), 1   # (roots, genus)

def alphas_from_angles(thetas, q):
    """genus g = len(thetas): alpha = sqrt(q) e^{i theta}, conj pairs (RH-for-curves)."""
    a = [np.sqrt(q)*np.exp(1j*t) for t in thetas] + [np.sqrt(q)*np.exp(-1j*t) for t in thetas]
    return np.array(a), len(thetas)

def hodge_bound_check(alphas, q, g, rmax=12):
    print(f"    r :  |sum alpha^r|      2g q^(r/2)     ratio   (ratio<1 = Hodge bound holds)")
    ok = True
    for r in range(1, rmax+1):
        S = abs(np.sum(alphas**r)); B = 2*g*q**(r/2)
        ok = ok and S <= B*(1+1e-9)
        print(f"    {r:>2}: {S:>14.3f}  {B:>14.3f}  {S/B:>8.3f}")
    return ok

if __name__ == "__main__":
    print("="*72)
    print("(1) VERIFY the Hodge-index / Castelnuovo-Severi bound on real curves (must hold):")
    print("="*72)
    al, g = alphas_elliptic(3, 11)
    print(f"\n  E/F11 (a_p=3, g=1): alphas={np.round(al,3)}, |alpha|={abs(al[0]):.4f}=sqrt(11)={np.sqrt(11):.4f}")
    print(f"  Hodge bound |sum alpha^r| <= 2g q^(r/2):")
    print(f"  bound holds: {hodge_bound_check(al, 11, g, 8)}")

    print("\n" + "="*72)
    print("(2) THE FORCING — counterfactual: push one alpha OFF the circle to |alpha|=sqrt(q)*1.15.")
    print("    If RH failed, the Hodge bound WOULD BREAK. This is why geometry => RH-for-curves.")
    print("="*72)
    q = 11; g = 1
    th = np.angle(al[0])
    al_bad = np.array([1.15*np.sqrt(q)*np.exp(1j*th), (1/1.15)*np.sqrt(q)*np.exp(-1j*th)])  # keep product=q
    print(f"\n  perturbed alphas |.|={abs(al_bad[0]):.3f},{abs(al_bad[1]):.3f} (product still q={abs(np.prod(al_bad)):.1f}):")
    print(f"  Hodge bound holds: {hodge_bound_check(al_bad, q, g, 12)}")
    print("  => an off-circle zero VIOLATES the geometric bound at large r. Geometry forbids it. THAT is the proof.")

    print("\n" + "="*72)
    print("(3) THE NUMBER-FIELD ANALOG and the GAP")
    print("="*72)
    print("""
  Function field                          Number field (Spec Z)
  --------------------------------------  ------------------------------------------
  curve C / F_q                           Spec Z  (the 'curve')
  surface C x C                           Spec Z x_{F_1} Spec Z   <-- DOES NOT EXIST (no F_1 geometry)
  intersection form on NS(C x C)          ??? (no surface => no intersection form)
  Hodge index thm: signature (1, rho-1)   ??? (the positivity has NO geometric source)
  => Castelnuovo-Severi |N_r-(q^r+1)|      => the analog is the WEIL EXPLICIT FORMULA bound, which is
     <= 2g q^{r/2}  (a THEOREM)               EXACTLY our a=1 razor -- and it is NOT a theorem (= RH).
  Frobenius eigenvalues alpha_i on circle Riemann zeros gamma on the line

  THE GAP, precisely: in the function field the positivity that forces the zeros onto the circle is
  the INTERSECTION FORM of a REAL algebraic surface (C x C), and its (1, rho-1) signature is the Hodge
  index THEOREM. In the number field there is no surface Spec Z x_{F_1} Spec Z, so the identical
  positivity (Weil's) has no geometric source -- it is the a=1 razor we mapped five ways, unproven
  precisely because the geometry that would supply it does not (yet) exist.
  => The Connes-Consani F_1 program / Deninger's dynamical cohomology / the arithmetic-surface Hodge
     index are ALL attempts to BUILD that missing surface. That is where a real proof must come from,
     and it is exactly the object the session's operator/positivity lanes kept bottoming out against.
""")
