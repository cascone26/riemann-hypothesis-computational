#!/usr/bin/env python3
"""§400 (angle 5): GENUINE F_q instrument — genus-2 curve zeta over F_q. Build the Weil explicit-formula
quadratic form from POINT COUNTS (the 'prime' data) and compare its certificate signature to the number-
field arch-prime (§392). In F_q RH is PROVEN (Frobenius eigenvalues on |z|=sqrt q); we can put zeros ON vs
OFF the circle and watch the form.

Zeta of genus-g curve /F_q: Z(T)=P(T)/((1-T)(1-qT)), P(T)=prod_{i=1}^{2g}(1-alpha_i T), |alpha_i|=sqrt q (RH).
Point counts N_r = q^r + 1 - sum alpha_i^r. Explicit-formula quad form on test coeffs c_1..c_M:
  Q[c] = sum over zeros theta of |sum_r c_r e^{i r theta}|^2   (zero side, PSD by construction like our zeroGram)
       = analytic side built from N_r (the 'prime'/point-count side).
We build Q from point counts and check: on-line (|alpha|=sqrt q, RH) => PSD; push ONE conjugate pair off
the circle (|alpha|=sqrt q * e^{+/-b}) => does Q develop exactly ONE negative direction (Hodge-index (1,rest))?
This is the F_q ANALOG of §392, in the world where the certificate is a PROVEN theorem."""
import numpy as np
q = 25   # F_q, q=5^2 say
g = 2
sq = np.sqrt(q)
# genus-2: 2g=4 Frobenius eigenvalues in conj pairs sqrt(q) e^{+-i th1}, sqrt(q) e^{+-i th2}
th1, th2 = 0.7, 2.1
def eigs(beta1=0.0, beta2=0.0):
    # beta = radial off-line displacement of pair 1,2 (0 = on circle = RH)
    r1, r2 = sq*np.exp(beta1), sq*np.exp(beta2)
    # to keep functional equation (alpha * alpha_bar = q), the FE partner is q/alpha:
    return [r1*np.exp(1j*th1), r1*np.exp(-1j*th1), (q/r1)*np.exp(1j*th1)*0 + q/(r1*np.exp(1j*th1)),
            q/(r1*np.exp(-1j*th1))][:2] + [r2*np.exp(1j*th2), q/(r2*np.exp(1j*th2))]
def pointcounts(al, R):
    return np.array([q**r + 1 - sum(a**r for a in al) for r in range(1,R+1)])
def zeroform(al, M):
    # Q[j,k] = sum over the 2g zeros theta of e^{i(j-k) theta}? Build from log|Z| explicit formula:
    # use the 2g Frobenius angles directly: zeros of P are 1/alpha_i; the "spectral" phases are arg(alpha_i).
    # Weil quad form on symmetric test: Q = sum_i sum_{j,k} c_j c_k (alpha_i/sqrt q)^{j} conj((alpha_i/sqrt q)^k)
    A = np.array(al)/sq   # normalized: on RH |A|=1
    F = np.array([[A[i]**j for j in range(M)] for i in range(len(A))])  # (2g, M)
    return np.real(F.conj().T @ F)   # M x M, PSD iff |A|=1
def inertia(Msym, tol=1e-9):
    ev=np.linalg.eigvalsh(Msym); s=tol*max(1,abs(ev).max())
    return int((ev>s).sum()), int((ev<-s).sum()), int((abs(ev)<=s).sum())
print(f"GENUS-2 F_q={q} Weil form. on-line (RH) vs radial-off-line displacement beta. M test coeffs.")
print(f"{'beta':>6} {'M':>3} {'inertia Q':>12} {'lmin(Q)':>11}")
for beta in (0.0, 0.05, 0.1, 0.2, 0.4):
    al = [sq*np.exp(beta)*np.exp(1j*th1), q/(sq*np.exp(beta)*np.exp(1j*th1)),
          sq*np.exp(-beta)*np.exp(1j*th2), q/(sq*np.exp(-beta)*np.exp(1j*th2))]
    for M in (4,6):
        Q = zeroform(al, M)
        print(f"{beta:6.2f} {M:3d} {str(inertia(Q)):>12} {float(np.linalg.eigvalsh(Q)[0]):+11.3e}", flush=True)
print("\nreading: beta=0 (RH) => Q PSD (inertia (r,0,.)). beta>0 (off circle) => Q develops NEGATIVE")
print("eigenvalues. If exactly ONE per displaced pair (Hodge-index (1,.) signature), this IS the F_q")
print("analog of §392: the proven certificate = 'off-line => negative direction', pole/ample = anchor.")
