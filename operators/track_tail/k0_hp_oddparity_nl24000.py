#!/usr/bin/env python3
"""k0_hp_oddparity_nl24000.py — ODD-PARITY-SECTOR twin of rigor/k0_hp_nl24000.py.
Modes = EVEN sine indices {2,4,...,KMAX} (LEDGER §192: ODD-v = even-k sines). Same kernel, same
column definitions as the even table (Gershgorin tail margin as delta_T proxy; Schur wc =
|Q_LT|_F^2/delta; near-null = |Q_TL u|^2/delta, u = min eigvec of low block), float64 (adequate:
odd margins ~1e-5, three orders above the even sector's 1e-8; NL=24000 resolves 1e-7 per §227).
Adds: headroom column (lowlmin/nearnull, per Isa), direct tail-window lambda_min alongside the
Gershgorin proxy, and lambda_min of the K0=24 certified block (per Rie).
Written for the piece-3 repair (see ODD_SECTOR_GAP.md); deposits results/Qodd52_a0549_nl24000_float.npy"""
import math, cmath, numpy as np
from scipy.special import digamma, zeta as hurwitz
A = 0.549; b = math.pi/(2*A); log2 = math.log(2); twoa = 2*A
wp = log2/math.sqrt(2); Om0 = digamma(0.25) - math.log(math.pi); Phi0 = hurwitz(2, 0.25)
KMAX = 52; NL = 24000
def Pp(p,t,s):
    if p==0: return 1/s
    if p==1: return t/s-1/s**2
    return t*t/s-2*t/s**2+2/s**3
def II(p,al,om,sc,t0,t1):
    s = complex(-al,om); F = lambda tt: cmath.exp(s*tt)*Pp(p,tt,s)
    v = F(t1)-F(t0)
    return v.real if sc=='cos' else v.imag
gt = [(-4,0,-0.5),(-4,0,0.5),(8,0,0.0),(-Om0/2,1,0.0),(-Phi0/4,0,0.0)]
for n in range(NL): gt.append((0.25/(n+0.25)**2,0,2*n+0.5))
def dterms(k,l):
    if k==l: return [(2*A*k*k*b*b,0,k*b,'cos'),(-k*k*b*b,1,k*b,'cos'),(-k*b,0,k*b,'sin')]
    al = k*l*b*(1.0/(k-l)-1.0/(k+l)); ak = -k*l*b*(1.0/(k-l)+1.0/(k+l))
    return [(al,0,l*b,'sin'),(ak,0,k*b,'sin')]
def Qraw(k,l):
    tot = 0.0
    D = dterms(k,l)
    for (ag,pg,alp) in gt:
        for (ad,pd,om,sc) in D: tot += ag*ad*II(pg+pd,alp,om,sc,0,twoa)
    for (ad,pd,om,sc) in D:
        tot += wp*ad*II(pd+1,0,om,sc,log2,twoa)+(-wp*log2)*ad*II(pd,0,om,sc,log2,twoa)
    return tot/A
MODES = list(range(2,KMAX+1,2)); idx = {k:i for i,k in enumerate(MODES)}; n = len(MODES)
print(f"k0_hp ODD-PARITY twin @ a={A}  precomputing {n}x{n} odd-sector block (EVEN modes <= {KMAX}, NL={NL})...", flush=True)
Q = np.zeros((n,n))
for i,k in enumerate(MODES):
    for j,l in enumerate(MODES):
        if j<i: Q[i,j]=Q[j,i]
        else: Q[i,j]=Qraw(k,l)
np.save('/Users/scones/projects/riemann/results/Qodd52_a0549_nl24000_float.npy', Q)
print("  ...cached + saved. Analyzing K0 candidates.\n", flush=True)
def tailrows(K0):
    ks=[k for k in MODES if k>K0]
    if not ks: return None
    return min(Q[idx[k],idx[k]] - sum(abs(Q[idx[k],idx[l]]) for l in ks if l!=k) for k in ks)
def lowblock(K0):
    ls=[k for k in MODES if k<=K0]
    B = Q[np.ix_([idx[k] for k in ls],[idx[k] for k in ls])]
    ev, V = np.linalg.eigh(B)
    return ls, ev[0], ev, V
def schur(K0):
    ls=[k for k in MODES if k<=K0]; ts=[k for k in MODES if k>K0]
    if not ts: return None,None,None
    g = tailrows(K0)
    tw = Q[np.ix_([idx[t] for t in ts],[idx[t] for t in ts])]
    lam_tail = np.linalg.eigh(tw)[0][0]
    if g<=0: return float('inf'),float('inf'),lam_tail
    LT = Q[np.ix_([idx[t] for t in ts],[idx[k] for k in ls])]
    wc = float(np.sum(LT**2))/g
    _,_,ev,V = lowblock(K0); u = V[:,0]
    act = float(np.sum((LT@u)**2))/g
    return wc, act, lam_tail
print(f"{'K0':>4} {'lowlmin':>12} {'tailGersh':>12} {'tail-lmin':>12} {'Schur wc':>12} {'Schur nearnull':>14} {'headroom':>10}")
for K0 in [k for k in MODES if 4<=k<=40]:
    ls,lmin,_,_ = lowblock(K0); g = tailrows(K0); wc,act,lt = schur(K0)
    hr = lmin/act if (act and act>0 and act!=float('inf')) else float('nan')
    print(f"{K0:>4} {lmin:>12.4g} {g:>12.4g} {lt:>12.4g} {wc:>12.4g} {act:>14.4g} {hr:>10.3g}", flush=True)
ls, lmin24, ev24, _ = lowblock(24)
print(f"\nRIE ASK 1: lambda_min(certified K0=24 low block, modes {ls}) = {lmin24:.6g}")
print(f"  full spectrum: {[f'{e:.4g}' for e in ev24]}")
wc24, act24, lt24 = schur(24)
g24 = tailrows(24)
print(f"RIE ASK 2: odd Schur residual at K0=24: worst-case {wc24:.4g}, near-null-directed {act24:.4g}")
print(f"RIE ASK 3: delta_T proxies at K0=24: tail Gershgorin {g24:.4g}, direct tail-window lambda_min {lt24:.4g}")
print(f"ISA HEADROOM at K0=24: lowlmin/nearnull = {lmin24/act24:.4g}")
print(f"CAVEAT: tail window truncated at KMAX={KMAX}; infinite tail composes via the mirrored 4-prime architecture (d0 law verified both parities).")
