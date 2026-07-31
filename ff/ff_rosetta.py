#!/usr/bin/env python3
"""FUNCTION-FIELD ROSETTA (Lane 2 of the superprism) — instrument the SOLVED case of RH.
For a curve C/F_q, the zeta function is Z(T)=P(T)/((1-T)(1-qT)), P(T)=prod(1-alpha_i T), deg P = 2g,
and RH (a THEOREM, Weil 1948) says |alpha_i| = sqrt(q) exactly — all "zeros" on the circle |T|=q^{-1/2},
i.e. the analog of Re(s)=1/2 under T=q^{-s}.

We build the EXACT function-field analogs of our number-field instruments and read what a WORKING proof
looks like through them:
  - Weil explicit-formula positivity (finite sum, no truncation): does an a=1 razor edge appear, or is there
    genuine positive margin (because here RH is TRUE with the Hodge-index mechanism behind it)?
  - The prime side = sum over closed points (Frobenius), archimedean side replaced by the (1-T)(1-qT) poles.

Start: the simplest nontrivial curves — elliptic curves E/F_p (g=1, P(T)=1-a_p T+p T^2, |alpha|=sqrt p iff
|a_p|<=2 sqrt p = Hasse, ALWAYS true). Then higher genus. Everything is a FINITE polynomial computation.

Key object: N_r = #C(F_{q^r}) = q^r + 1 - sum alpha_i^r  (the "prime counts" / Frobenius power sums).
Explicit formula test function f on the alpha_i (points on |T|=q^{-1/2} under RH): the Weil pairing.
"""
import numpy as np

def ecount_primes(a_p, p, rmax=40):
    """N_r = p^r + 1 - (alpha^r + alphabar^r), alpha,alphabar roots of T^2 - a_p T + p (recip of P)."""
    # Frobenius eigenvalues: roots of x^2 - a_p x + p = 0 (these are the alpha_i with |alpha|=sqrt p)
    disc = a_p*a_p - 4*p
    if disc < 0:
        alpha = complex(a_p/2, np.sqrt(-disc)/2)
    else:
        alpha = complex((a_p+np.sqrt(disc))/2, 0)
    ab = p/alpha if alpha != 0 else 0
    tr = np.array([alpha**r + ab**r for r in range(1, rmax+1)])
    return np.real(tr)   # power sums s_r = sum alpha_i^r (real; conj pair)

def weil_form_ff(alphas, q, testcoef):
    """Function-field Weil explicit-formula quadratic form on a test 'function' given by Fourier coefs
    c_r (r=1..R) on the circle. The zero side is sum over alpha_i of |sum_r c_r alpha_i^{-r}|^2-type.
    Under RH |alpha_i|=sqrt q, alpha_i^{-r} = q^{-r/2} e^{-i r theta_i}. Zero side = sum_i |F(theta_i)|^2 >=0
    AUTOMATICALLY (it's a sum of squares) — that's the point: in the solved case positivity is MANIFEST.
    We instrument the analog of OUR number-field Weil form and check the margin structure."""
    R = len(testcoef)
    # angles of the alpha_i on the circle |T|=sqrt q  (Frobenius angles)
    thetas = np.angle(alphas)
    F = np.array([sum(testcoef[r]*np.exp(1j*(r+1)*th) for r in range(R)) for th in thetas])
    return np.sum(np.abs(F)**2)

if __name__ == "__main__":
    print("FUNCTION-FIELD ROSETTA — elliptic curves E/F_p (g=1), the SOLVED RH case.")
    print("Frobenius eigenvalues alpha, alphabar with |alpha|=sqrt(p) (Hasse/RH-for-curves, a THEOREM).\n")
    for (a_p, p, name) in [(0,5,"y^2=x^3+x /F5"),(-1,7,"y^2=x^3+1 /F7"),(3,11,"E/F11"),(-2,13,"E/F13")]:
        disc = a_p*a_p - 4*p
        alpha = complex(a_p/2, np.sqrt(-disc)/2) if disc < 0 else complex((a_p+np.sqrt(disc))/2,0)
        theta = np.angle(alpha)
        s = ecount_primes(a_p, p, 8)
        print(f"{name}: a_p={a_p}, |alpha|=sqrt({p})={np.sqrt(p):.4f}, |alpha|_actual={abs(alpha):.4f}, "
              f"Frob angle={theta:.4f}")
        print(f"   #E(F_p^r) - (p^r+1) = -s_r for r=1..8: {np.round(-s,2)}")
        print(f"   RH check: |alpha|=sqrt(p)? {abs(abs(alpha)-np.sqrt(p))<1e-9}")
    print("\nNEXT: build the exact Weil-form analog + look for the a=1 razor (or genuine margin) in the solved case.")
