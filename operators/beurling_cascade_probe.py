#!/usr/bin/env python3
"""beurling_cascade_probe.py -- §253: does the near-null TELESCOPING CASCADE (§176) halt for Beurling primes?

§176/§176b (real primes, dps=150): the Weil-form near-null ψ has arch expectation +E cancelled by the prime sum
Σh(n) to ~59 orders, and it TELESCOPES: residual after prime-power n ≈ |h(next)| (each prime cancels the previous
residual). §L2 prediction: a Beurling comb (no additive/Poisson coherence) should HALT — residual plateaus, no clean
telescope. This float probe reads the FIRST FEW cascade steps (arch→n2→n3→n4→n5, real primes reach log10 −8.6 there,
float-visible) for REAL vs BEURLING combs. If Beurling halts early (residual doesn't track |h(next)|), the cascade
IS the integer-fingerprint = the concrete 'what a proof must consume'.

Reuses the connes_decisive tau structure (arch digamma + W02 pole - prime q-kernel) but with a CONFIGURABLE comb.
Float caveat: telescope only visible to ~1e-13; that covers arch..n5. Deep tail needs dps=150 (cvs_gate, Olympus).
"""
import numpy as np, mpmath as mp
from scipy.special import digamma as _digamma
mp.mp.dps = 30
LAMBDA_SQ = 14
L = float(2*mp.log(mp.sqrt(LAMBDA_SQ)))
TWO_PI = 2.0*np.pi; LOGPI = float(np.log(np.pi))

def omega_vec(r): return np.real(_digamma(0.25+0.5j*r))-LOGPI
def e_hat_vec(r,n):
    s=TWO_PI*n/L; d=s-r; out=np.empty_like(d,dtype=complex); sm=np.abs(d)<1e-12
    out[sm]=np.sqrt(L)+0j; dd=d[~sm]
    out[~sm]=(1/np.sqrt(L))*np.exp(0.5j*dd*L)*2*np.sin(0.5*dd*L)/dd; return out
def q_kernel(n,m,y):
    if n==m: return 2*(1-abs(y)/L)*np.cos(TWO_PI*n*y/L)
    return (np.sin(TWO_PI*m*y/L)-np.sin(TWO_PI*n*y/L))/(np.pi*(n-m))
def w02_mat(idx):
    D=len(idx); W=np.zeros((D,D)); sh=np.sinh(L/4)**2
    for ii,n in enumerate(idx):
        for jj,m in enumerate(idx):
            num=L**2-16*np.pi**2*m*n
            den=(L**2+16*np.pi**2*m**2)*(L**2+16*np.pi**2*n**2)
            W[ii,jj]=32*L*sh*num/den
    return W
def arch_pole(idx,R=60.0,G=24000):
    grid=np.linspace(-R,R,G); h=grid[1]-grid[0]; w=np.full(G,h); w[0]=w[-1]=0.5*h
    Om=omega_vec(grid); D=len(idx); E=np.empty((D,G),dtype=complex)
    for k,n in enumerate(idx): E[k]=e_hat_vec(grid,n)
    A=E*(w*Om); arch=(A@E.conj().T).real/(2*np.pi)
    return arch+w02_mat(idx)
def Vn_mat(idx, weight, pos):
    """single-'prime' matrix: weight * q_kernel at log-position 'pos' (pos=log q for real)."""
    D=len(idx); P=np.zeros((D,D))
    for ii,n in enumerate(idx):
        for jj in range(ii,D):
            m=idx[jj]; P[ii,jj]=P[jj,ii]=weight*q_kernel(n,m,pos)
    return P

def cascade(comb, idx, R=60.0, G=24000):
    """comb = list of (label, weight, logpos). Build tau=archpole - sum Vn, get near-null psi, telescope."""
    AP=arch_pole(idx,R,G)
    Vs=[(lab, Vn_mat(idx,w,p)) for (lab,w,p) in comb]
    tau=AP - sum(V for _,V in Vs); tau=(tau+tau.T)/2
    ev,U=np.linalg.eigh(tau); psi=U[:,0]
    E_arch=float(psi@AP@psi)
    hs=[(lab, float(psi@V@psi)) for lab,V in Vs]
    # telescope residual R_m = E_arch - sum_{<=m} h  (h defined with +sign since tau=AP - sumV, <psi|V|psi>)
    print(f"   arch+pole expectation = {E_arch:+.6e}   lambda_min={ev[0]:.3e}")
    Rm=E_arch
    print(f"   {'step':>8} {'h(step)':>13} {'R_after':>13} {'|R_after|/|h_next|':>18}")
    for i,(lab,h) in enumerate(hs):
        Rm-=h
        nxt=abs(hs[i+1][1]) if i+1<len(hs) else float('nan')
        ratio=abs(Rm)/nxt if nxt==nxt and nxt>0 else float('nan')
        print(f"   {lab:>8} {h:+13.5e} {Rm:+13.5e} {ratio:18.3f}")
    return E_arch, hs

def vm_real(kmax):
    combs=[]
    for k in range(2,kmax+1):
        n,p,f,ok=k,0,2,True
        while f*f<=n:
            if n%f==0:
                if p==0:p=f
                elif p!=f:ok=False;break
                while n%f==0:n//=f
            f+=1
        if n>1:
            if p==0:p=n
            elif p!=n:ok=False
        if ok and p>0: combs.append((f"n{k}", float(mp.log(p))/np.sqrt(k), np.log(k)))
    return combs

if __name__=="__main__":
    N=24; idx=list(range(-N,N+1))
    Lmax=int(np.floor(L))  # positions must have log q < L
    print("== REAL primes (telescopes to ~1e-59 at dps150; float sees arch..n5) ==")
    real=[c for c in vm_real(Lmax) if c[2]<L]
    cascade(real, idx)
    # Beurling B1: {2} only -> zeta_B has zeros OFF the 1/2 line (on Re=0). weight log2/sqrt(2^k) at k*log2
    print("\n== BEURLING B1: powers of 2 only (RH-false: zeros off critical line) ==")
    b1=[(f"2^{k}", float(mp.log(2))/np.sqrt(2**k), k*np.log(2)) for k in range(1,int(L/np.log(2))+1)]
    cascade(b1, idx)
    # Beurling B2: shift every real prime position by +0.15 (breaks PNT/additive coherence), same weights
    print("\n== BEURLING B2: real weights at SHIFTED positions (log q + 0.15) -- breaks additive structure ==")
    b2=[(lab,w,p+0.15) for (lab,w,p) in real if p+0.15<L]
    cascade(b2, idx)
    print("\nREAD: if REAL telescopes (ratio ~1.0 each step) but BEURLING ratios scatter/blow up => cascade halts")
    print("for Beurling = the integer-fingerprint. If both telescope => cascade is NOT the discriminator (float).")
