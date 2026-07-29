#!/usr/bin/env python3
"""§448 VERIFY: is §446/§447's off-line detection REAL or a conditioning artifact? §446's form = §444's up to
diag(e^{(r+s)/2}) conjugation (inertia-invariant), but that factor blows to e^12 at M=24 → ill-conditioned. Build
the WELL-CONDITIONED FE-closed Rosati (factor out e^{(r+s)/2}, i.e. use unit-modulus α=e^{iγ}, partner e^{-iγ},
with the β as radial e^{β}): the honest inertia. Also report cond number. If it STAYS PSD for all β → §446/447
were numerical artifacts and the number-field Rosati (this form) does NOT detect off-line (= §444 was right)."""
import numpy as np
G=np.load('zeros150.npy')
def rosati_cond(beta,M,K):   # §446 as-built (ill-conditioned, e^{ρr})
    g=G[:K]; rho=0.5+beta+1j*g; fe=1-rho; r=np.arange(M); A=np.zeros((M,M),dtype=complex)
    for k in range(K):
        er=np.exp(rho[k]*r); ef=np.exp(fe[k]*r); A+=np.outer(er,ef)+np.outer(ef,er)
    return np.real((A+A.conj().T)/2)
def rosati_wellcond(beta,M,K):  # factor out e^{(r+s)/2}: use λ=e^{β+iγ}, partner e^{-β-iγ} (unit-ish, bounded)
    g=G[:K]; r=np.arange(M); A=np.zeros((M,M),dtype=complex)
    for k in range(K):
        lam=np.exp(beta+1j*g[k]); par=np.exp(-beta-1j*g[k])   # |lam|=e^β, bounded ~1
        er=lam**r; ef=par**r; A+=np.outer(er,ef)+np.outer(ef,er)
    return np.real((A+A.conj().T)/2)
def report(f,beta,M,K):
    A=f(beta,M,K); e=np.linalg.eigvalsh(A); m=max(1,abs(e).max())
    nneg=int((e<-1e-9*m).sum()); cond=abs(e).max()/max(abs(e).min(),1e-300)
    return nneg,e[0],cond
print("§448 well-conditioned FE-closed Rosati vs the ill-conditioned §446 form. K=40.")
print(f"{'β':>6} {'M':>3} | {'ill: #neg':>9} {'cond':>9} | {'well: #neg':>10} {'well λmin':>11} {'well cond':>10}")
for beta in [0.0,0.1,0.2,0.3,0.5]:
    for M in (10,20):
        ni,ei,ci=report(rosati_cond,beta,M,40)
        nw,ew,cw=report(rosati_wellcond,beta,M,40)
        print(f"{beta:6.2f} {M:3d} | {ni:9d} {ci:9.1e} | {nw:10d} {ew:+11.3e} {cw:10.1e}",flush=True)
print("""
reading: if 'well: #neg'=0 for ALL β while 'ill: #neg'>0 at large β with huge cond => §446/447 detection was a
CONDITIONING ARTIFACT (the e^{(r+s)/2} blowup), and the honest FE-closed number-field Rosati does NOT detect
off-line (stays PSD, = §444). Then the pure-zero Rosati genuinely CANNOT see off-line — the archimedean Γ-factor
is ESSENTIAL (not just a sharpener). If 'well: #neg'>0 too => detection is real.
""")
