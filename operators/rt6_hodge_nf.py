#!/usr/bin/env python3
"""§405 (handoff i): number-field Hodge-index hunt. In F_q the certificate lives on NS(CxC) with the
intersection form signature (1,rho-1); the ample class (fibers/diagonal) anchors it. On the number field
we have B_full = pole + arch - prime (PSD = RH). QUESTION: does B_full carry an intersection-form-like
INDEFINITE companion whose signature is the real certificate? Build the ROSATI-analog on the number-field
'eigenvalues' = the zeros themselves: A[r,s] = sum_gamma (rho_gamma)^r * (1-rho_gamma-bar)^s style using the
FE partner rho -> 1-rho. On the critical line rho=1/2+i*gamma, FE partner 1-rho = 1/2-i*gamma = conj(rho),
so A is a Hermitian Gram (PSD). Test the number-field Rosati form directly: A[r,s]=sum_g rho^r (1-rho)^s with
rho=1/2+i*gamma (finitely many zeros), on-line vs a hypothetical off-line gamma -> does it go indefinite like
F_q? THIS is the direct number-field analog of the proven §401 F_q object."""
import numpy as np
import l2_edge_sweep as L
Z = L.ZEROS  # gamma values
def rosati_nf(betas, M):
    # rho_k = 1/2 + beta_k + i*gamma_k ; FE partner = 1 - rho_k = 1/2 - beta_k - i*gamma_k
    rho = 0.5 + betas + 1j*Z
    fe  = 1.0 - rho
    A = np.zeros((M,M), dtype=complex)
    for r in range(M):
        for s in range(M):
            A[r,s] = np.sum(rho**r * fe**s)
    return np.real((A + A.conj().T)/2)
def inertia(Ms, tol=1e-9):
    ev=np.linalg.eigvalsh(Ms); t=tol*max(1,abs(ev).max())
    return int((ev>t).sum()), int((ev<-t).sum()), int((abs(ev)<=t).sum())
print("NUMBER-FIELD Rosati/Weil form A[r,s]=sum_g rho^r (1-rho)^s, rho=1/2+beta+i*gamma. RH: beta=0.")
print("(direct analog of proven F_q §401: on-line PSD, off-line indefinite?)")
print(f"{'beta':>7} {'M':>3} {'inertia A':>12} {'lmin':>11}")
for beta in (0.0, 0.02, 0.05, 0.1):
    b = np.zeros(len(Z)); b[0] = beta; b[1] = -beta   # displace 2 zeros off-line (keep FE closure via 1-rho)
    for M in (4, 6, 8):
        A = rosati_nf(b, M)
        print(f"{beta:7.3f} {M:3d} {str(inertia(A)):>12} {float(np.linalg.eigvalsh(A)[0]):+11.3e}", flush=True)
print("\nreading: if beta=0 PSD and beta>0 develops NEGATIVE eigenvalues (like F_q §401), the number-field")
print("Rosati form IS the certificate and off-line zeros break its PSD-ness DIRECTLY (Weil positivity in")
print("zero coordinates, not prime coordinates). This is well-posed (unlike §378: we move zeros in the ZERO")
print("form, not the prime form). The gap to a PROOF = why must this specific form be PSD (=the Hodge index).")
