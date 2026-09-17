#!/usr/bin/env python3
"""§961 WHY κ≤1? — the COHERENCE/ALIGNMENT mechanism (connects §944 Sidon to the one negative square).

Target (§960): prove κ(arch−prime) ≤ 1. The pole gives κ≥1; RH is the upper bound κ≤1. MECHANISM hypothesis:
the prime form's "excess over arch" D = prime − arch has exactly ONE positive eigenvalue (=κ) because the
individual prime terms' contributions ALIGN coherently into a single direction (the arithmetic/Sidon coherence
of {log n}, §944). Perturbing the primes breaks the alignment → the excess spreads into several directions → κ>1.
TEST: for REAL vs perturbed primes, (a) #pos(D)=κ; (b) the ALIGNMENT of the per-prime excess directions —
build the 'excess operator' E_n = c_n C_n − (arch share), find how concentrated D's positive mass is
(participation ratio / effective rank of the positive part). High coherence (eff-rank→1) for real primes,
spreading for perturbed, would be the structural reason κ=1 — a candidate path to PROVING κ≤1 from {log n}'s
arithmetic. HANDLE (own code, honest — report the numbers, don't force a story)."""
import numpy as np
from scipy.special import digamma

def vonmangoldt(N):
    lam = np.zeros(N+1); isc = np.zeros(N+1, dtype=bool)
    for p in range(2, N+1):
        if not isc[p]:
            isc[2*p::p] = True; lp = np.log(p); pk = p
            while pk <= N: lam[pk] = lp; pk *= p
    return lam
NMAX=8000; LAM=vonmangoldt(NMAX); NS=np.arange(2,NMAX+1); LOGN=np.log(NS); LAMv=LAM[2:]

def build(a,b,W,N,freqs,weights,dr=0.004):
    RG=np.arange(-(b+120),(b+120),dr)
    OM=np.real(digamma(0.25+0.5j*RG))-np.log(np.pi)
    CEN=np.linspace(a,b,N)
    def fej(x):
        y=W*x;o=np.ones_like(y);nz=np.abs(y)>1e-12;o[nz]=(np.sin(y[nz])/y[nz])**2;return o
    P=np.array([fej(RG-c)+fej(RG+c) for c in CEN]);G=(P@P.T)*dr
    arch=(P*OM[None,:])@P.T*dr/(2*np.pi)
    # per-frequency prime pieces C_n (each the compression of cos(r*freq)); prime=Σ w_n C_n
    pieces=[]
    for f,w in zip(freqs,weights):
        cu=np.cos(f*RG); Cn=np.zeros((N,N))
        for j in range(N):
            Cn[j,:]=2.0*((P[j][None,:]*P)@cu)*dr/(2*np.pi)
        pieces.append((0.5*(Cn+Cn.T), w))
    prime=sum(w*Cn for Cn,w in pieces)
    return arch, prime, pieces, G

def whiten(G):
    gw,gv=np.linalg.eigh(G);k=gw>gw.max()*1e-12;return gv[:,k]/np.sqrt(gw[k])

def analyze(a,b,W,N,freqs,weights,label):
    arch,prime,pieces,G=build(a,b,W,N,freqs,weights)
    B=whiten(G)
    D = B.T@(prime-arch)@B                       # excess of prime over arch; #pos(D)=κ(arch-prime)
    ev,V=np.linalg.eigh(0.5*(D+D.T))
    s=max(1,abs(ev).max()); kap=int((ev>1e-9*s).sum())
    # coherence: for the TOP positive eigenvector u of D, per-prime contribution c_n (u^T C_n u)
    u=V[:,-1]
    contribs=np.array([w*(u@ (B.T@Cn@B) @u) for Cn,w in pieces])
    total=contribs.sum()
    # participation ratio of the (positive) contributions: PR=(Σx)²/Σx² ; low PR=few dominate, high=coherent spread
    pos=contribs[contribs>0]
    PR = (pos.sum()**2/ (pos**2).sum()) if pos.size else 0
    # coherence ratio: actual total vs incoherent (sqrt-sum) expectation
    coh = total/ (np.sqrt((contribs**2).sum())+1e-30)
    print(f"{label:>22} | κ={kap} | top-eig={ev[-1]:+.3e} | #primes={len(freqs)} | PR(+contribs)={PR:5.1f} | coherence(tot/rms)={coh:+6.2f}")
    return kap

rng=np.random.default_rng(961)
a,b,fourW=14,44,5.0; W=fourW/4; N=round(0.85*(b-a)*2*W/np.pi)
mask=LOGN<fourW; ln=LOGN[mask]; lw=LAMv[mask]/np.sqrt(NS[mask]); m=len(ln)
print("="*104)
print(f"§961 coherence mechanism for κ. window[{a},{b}] 4W={fourW} N={N}, {m} prime-power terms.")
print("="*104)
analyze(a,b,W,N,ln,lw,"REAL primes")
for d in [1e-3,1e-2,3e-2]:
    analyze(a,b,W,N,ln*(1+d*rng.uniform(-1,1,m)),lw,f"perturbed δ={d}")
analyze(a,b,W,N,rng.uniform(ln.min(),ln.max(),m),lw,"scrambled freqs")
print("="*104)
print("READ: if REAL primes give κ=1 with HIGH coherence (tot/rms≈1, contributions constructively aligned into")
print("ONE direction) while perturbed give κ>1 with LOW coherence (contributions spread/cancel), then κ=1 is the")
print("arithmetic COHERENCE of {log n} (§944 Sidon) forcing the prime-excess into a single direction — a candidate")
print("structural route to PROVING κ≤1. Honest: report whatever the coherence numbers actually show.")
print("="*104)
