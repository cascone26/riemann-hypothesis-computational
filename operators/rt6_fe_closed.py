#!/usr/bin/env python3
"""§446 (CORRECT FE-closed Rosati — decisive). §444/§445 failed because the q-factor is diagonal conjugation
(no inertia change) and my eigenvalue set wasn't FE-closed. §401 detected because its set was CLOSED under the FE
involution. Build it right for ζ: eigenvalues e^{ρ}, FE partner e^{1-ρ}. Rosati bilinear form
B(u,v)=Σ_pairs [ (Σ_r u_r e^{ρr})(Σ_s v_s e^{(1-ρ)s}) + FE-swap ], as a matrix A[r,s]=Σ_ρ [e^{ρr}e^{(1-ρ)s}+e^{(1-ρ)r}e^{ρs}].
On RH (1-ρ=ρ̄) this is a real Gram (PSD); off-line test whether it goes indefinite. This is the honest test of
whether the FE-closed number-field Rosati detects off-line."""
import numpy as np
G=np.load('zeros150.npy')
def rosati_fe(beta,M,K=40):
    g=G[:K]; rho=0.5+beta+1j*g; fe=1-rho     # ρ and FE partner 1-ρ
    A=np.zeros((M,M),dtype=complex)
    r=np.arange(M)
    for k in range(K):
        er=np.exp(rho[k]*r); ef=np.exp(fe[k]*r)
        A+=np.outer(er,ef)+np.outer(ef,er)
    return np.real((A+A.conj().T)/2)
def inertia(A,tol=1e-7):
    e=np.linalg.eigvalsh(A); m=max(1,abs(e).max())
    return int((e>tol*m).sum()),int((abs(e)<=tol*m).sum()),int((e<-tol*m).sum())
print("§446 FE-closed Rosati (eigenvalues e^ρ, partner e^{1-ρ}). detect off-line? K=40 zeros.")
print(f"{'β':>7} {'M':>3} {'inertia(+,0,−)':>15} {'λmin':>11}")
for beta in [0.0,0.05,0.1,0.2,0.3]:
    for M in (6,10):
        A=rosati_fe(beta,M); print(f"{beta:7.3f} {M:3d} {str(inertia(A)):>15} {np.linalg.eigvalsh(A)[0]:+11.3e}",flush=True)
print("""
reading: β=0 PSD (Gram, since 1-ρ=ρ̄); β>0 developing NEGATIVE eigenvalues => the FE-closed number-field Rosati
DOES detect off-line (the correct mechanism, honest fix of §444/445). Still PSD for all β => even the FE-closed
form degenerates for ℚ (deeper obstruction) — then the archimedean gamma-factor is genuinely needed in the form.
""")
