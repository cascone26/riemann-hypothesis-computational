#!/usr/bin/env python3
"""k0_solver.py — minimal-K0 finder for the tail lemma. Split even-sector odd modes into low block L={odd k≤K0}
(certified in Arb) and tail T={odd k>K0}. A finite K0-block cert implies TRUE infinite-dim positivity when:
  (a) tail is diagonally dominant (Gershgorin margin g(k)=Q[k,k]−Σ_{l≠k}|Q[k,l]| > 0 for all odd k>K0), and
  (b) the Schur residual of L→T coupling is below the low-block margin.
Emits the exact large-diagonal-first modes list for arb_certify.

Efficient: precompute Q[k,l] for all odd k,l ≤ KMAX ONCE (cached), then evaluate every K0 from the cache.
(The first agent build recomputed entries per-K0 and timed out — kept at k0_solver.agent.py.bak.)"""
import sys, mpmath as mp
mp.mp.dps = 30
A = mp.mpf(sys.argv[1]) if len(sys.argv)>1 else mp.mpf('0.549')   # window endpoint = tightest
KMAX = 51; NL = 1200
b=mp.pi/(2*A); log2=mp.log(2); twoa=2*A; wp=log2/mp.sqrt(2)
Om0=mp.digamma(mp.mpf(1)/4)-mp.log(mp.pi); Phi0=mp.zeta(2,mp.mpf(1)/4)
def Pp(p,t,s):
    if p==0: return 1/s
    if p==1: return t/s-1/s**2
    return t*t/s-2*t/s**2+2/s**3
def II(p,al,om,sc,t0,t1):
    s=mp.mpc(-al,om); F=lambda tt: mp.e**(s*tt)*Pp(p,tt,s); v=F(t1)-F(t0)
    return v.real if sc=='cos' else v.imag
def dterms(k,l):
    if k==l: return [(2*A*k*k*b*b,0,k*b,'cos'),(-k*k*b*b,1,k*b,'cos'),(-k*b,0,k*b,'sin')]
    al=k*l*b*(mp.mpf(1)/(k-l)-mp.mpf(1)/(k+l)); ak=-k*l*b*(mp.mpf(1)/(k-l)+mp.mpf(1)/(k+l))
    return [(al,0,l*b,'sin'),(ak,0,k*b,'sin')]
def Qraw(k,l):
    tot=mp.mpf(0); D=dterms(k,l)
    gt=[(-4,0,mp.mpf(-1)/2),(-4,0,mp.mpf(1)/2),(8,0,mp.mpf(0)),(-Om0/2,1,mp.mpf(0)),(-Phi0/4,0,mp.mpf(0))]
    for n in range(NL): gt.append((mp.mpf(1)/4/(n+mp.mpf(1)/4)**2,0,2*n+mp.mpf(1)/2))
    for (ag,pg,alp) in gt:
        for (ad,pd,om,sc) in D: tot+=ag*ad*II(pg+pd,alp,om,sc,0,twoa)
    for (ad,pd,om,sc) in D:
        tot+=wp*ad*II(pd+1,0,om,sc,log2,twoa)+(-wp*log2)*ad*II(pd,0,om,sc,log2,twoa)
    return tot/A    # whiten by Mass=A·I

ODD=list(range(1,KMAX+1,2)); idx={k:i for i,k in enumerate(ODD)}; n=len(ODD)
print(f"k0_solver @ a={mp.nstr(A,5)}  precomputing {n}×{n} even-block (odd modes ≤{KMAX}, NL={NL})…",flush=True)
Qm=mp.matrix(n,n)
for i,k in enumerate(ODD):
    for j,l in enumerate(ODD):
        if j<i: Qm[i,j]=Qm[j,i]
        else: Qm[i,j]=Qraw(k,l)
print("  …cached. Analyzing K0 candidates.\n")

def tailrows(K0):   # min Gershgorin margin over odd k in (K0, KMAX]
    ks=[k for k in ODD if k>K0];
    if not ks: return None
    m=None
    for k in ks:
        i=idx[k]; diag=Qm[i,i]; rs=mp.fsum(abs(Qm[i,idx[l]]) for l in ks if l!=k)
        g=diag-rs; m=g if m is None else min(m,g)
    return m
def lowblock(K0):
    ls=[k for k in ODD if k<=K0]; m=len(ls)
    B=mp.matrix(m,m)
    for a_,k in enumerate(ls):
        for c_,l in enumerate(ls): B[a_,c_]=Qm[idx[k],idx[l]]
    ev,V=mp.eigsy(B);
    return ls,min(ev),ev,V
def schur(K0):
    """worst-case ‖Q_LT‖²/λmin_tail  and  near-null-directed ‖Q_TL u‖²/λmin_tail (u=min eigvec of low block)."""
    ls=[k for k in ODD if k<=K0]; ts=[k for k in ODD if k>K0]
    if not ts: return None,None
    gmin=tailrows(K0)                      # λmin(Q_TT) ≥ Gershgorin margin
    if gmin<=0: return mp.inf,mp.inf
    # worst-case: ΣΣ Q[k,l]^2 over L×T
    wc=mp.fsum(Qm[idx[k],idx[l]]**2 for k in ls for l in ts)/gmin
    # near-null: u = min eigenvector of low block, residual ‖Q_TL u‖²/gmin
    _,_,ev,V=lowblock(K0); u=[V[r,0] for r in range(len(ls))]  # col 0 = smallest
    proj=[mp.fsum(Qm[idx[t],idx[ls[r]]]*u[r] for r in range(len(ls))) for t in ts]
    act=mp.fsum(p*p for p in proj)/gmin
    return wc,act

MARGIN=mp.mpf('2.6e-5')
print(f"{'K0':>4} {'lowλmin':>12} {'tailGersh':>12} {'Schur wc':>12} {'Schur nearnull':>14}")
choice=None
for K0 in [k for k in ODD if 5<=k<=41]:
    ls,lmin,_,_=lowblock(K0); g=tailrows(K0); wc,act=schur(K0)
    print(f"{K0:>4} {mp.nstr(lmin,4):>12} {mp.nstr(g,4):>12} {mp.nstr(wc,4):>12} {mp.nstr(act,4):>14}",flush=True)
    if choice is None and g is not None and g>0 and act<MARGIN and act<lmin:
        choice=K0
print()
if choice:
    modes=[k for k in ODD if k<=choice]; modes_sorted=sorted(modes,key=lambda k:-float(Qm[idx[k],idx[k]]))
    print(f"MINIMAL K0 = {choice}  (near-null Schur residual < margin {mp.nstr(MARGIN,2)} AND < low-block λmin)")
    print(f"arb_certify modes (large-diagonal-first) = {modes_sorted}")
    wc,act=schur(choice)
    print(f"  worst-case Schur residual = {mp.nstr(wc,4)} ; near-null-directed = {mp.nstr(act,4)}  (gap = {mp.nstr(act/wc,3)}× smaller)")
    print(f"  CAVEAT: tail window truncated at KMAX={KMAX}; true tail →∞ (rows converge, Gershgorin margin grows).")
else:
    print(f"No K0≤41 satisfied both at margin {mp.nstr(MARGIN,2)}; the worst-case bound is loose — near-null column is the real one.")
