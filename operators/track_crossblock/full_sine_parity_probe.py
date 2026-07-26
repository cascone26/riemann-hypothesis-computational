#!/usr/bin/env python3
"""§209 angle-3 probe: build the p=2-only Weil form on the FULL sine basis (k=1..K, BOTH parities)
at a=0.549 (=½log3), using the exact screw entry engine copied verbatim from
operators/dbn_dynamics/arb_certify.py (only NLERCH lowered for a float structural probe; prec kept high).

Questions:
 1) Does the form parity-decouple? Measure the odd<->even CROSS-block norm vs the diagonal blocks.
    (§186 claims parity-decoupled. If cross-block ≈ 0, combining Groskin's periodic/even content with
     screw's odd-sine content is CLEAN — angle-3 in LEDGER §209.)
 2) Is the FULL K×K sine-basis form PD (positivity on the whole sine subspace, not just odd modes)?
 3) Report eigenvalues of full / odd-block / even-block.

NOTE: float/interval-midpoint probe (REPORT tier). The rigorous cert is arb_certify.py at NLERCH=8000.
"""
import sys
from flint import arb, acb, ctx
ctx.prec = 400

A = arb('0.549'); B = arb.pi()/(2*A); LOG2 = arb(2).log(); TWOA = 2*A
OM0 = (arb(1)/4).digamma() - arb.pi().log()
PHI0 = acb(2).zeta(acb(arb(1)/4)).real
WP = arb(2).log()/arb(2).sqrt()
NLERCH = 3000  # lowered from 8000 (float structural probe, not a rigorous cert)

def Pp(p, t, s):
    if p == 0: return 1/s
    if p == 1: return t/s - 1/s**2
    return t*t/s - 2*t/s**2 + 2/s**3
def II(p, alpha, omega, sc, t0, t1):
    s = acb(-alpha, omega)
    def F(tt):
        t = acb(tt)
        return (s*t).exp()*Pp(p, t, s)
    v = F(t1) - F(t0)
    return v.real if sc == 'cos' else v.imag
def dterms(k, l):
    kb = arb(k)*B; lb = arb(l)*B
    if k == l:
        return [(2*A*arb(k)**2*B**2, 0, kb, 'cos'), (-arb(k)**2*B**2, 1, kb, 'cos'), (-arb(k)*B, 0, kb, 'sin')]
    al = arb(k)*arb(l)*B*(arb(1)/(k-l)-arb(1)/(k+l)); ak = -arb(k)*arb(l)*B*(arb(1)/(k-l)+arb(1)/(k+l))
    return [(al, 0, arb(l)*B, 'sin'), (ak, 0, arb(k)*B, 'sin')]
def Dmax(k, l):
    tot = arb(0)
    for (ad, pd, om, sc) in dterms(k, l): tot += abs(ad)*(TWOA if pd == 1 else arb(1))
    return tot
def Q(k, l):
    tot = arb(0); D = dterms(k, l)
    gt = [(arb(-4), 0, arb(-1)/2), (arb(-4), 0, arb(1)/2), (arb(8), 0, arb(0)), (-OM0/2, 1, arb(0)), (-PHI0/4, 0, arb(0))]
    for n in range(NLERCH): gt.append((arb(1)/4/(arb(n)+arb(1)/4)**2, 0, 2*arb(n)+arb(1)/2))
    for (ag, pg, alp) in gt:
        for (ad, pd, om, sc) in D: tot += ag*ad*II(pg+pd, alp, om, sc, arb(0), TWOA)
    N = NLERCH; dm = Dmax(k, l)
    tailbound = dm/8 * (arb(1)/(2*(arb(N)-1)**2))
    tot += tailbound*arb(0, 1)
    for (ad, pd, om, sc) in D:
        tot += WP*ad*II(pd+1, arb(0), om, sc, LOG2, TWOA) + (-WP*LOG2)*ad*II(pd, arb(0), om, sc, LOG2, TWOA)
    return tot

import numpy as np
K = 8
modes = list(range(1, K+1))   # 1..8, BOTH parities
n = len(modes)
print(f"Building full sine-basis Weil form, modes={modes}, a=0.549 (p=2 only), NLERCH={NLERCH}, prec={ctx.prec}")
Mf = np.zeros((n, n))
for i, k in enumerate(modes):
    for j, l in enumerate(modes):
        if j < i:
            Mf[i, j] = Mf[j, i]; continue
        v = Q(k, l)
        Mf[i, j] = float(v.mid() if hasattr(v, 'mid') else v)
        Mf[j, i] = Mf[i, j]
    print(f"  row {k} done", flush=True)

# parity split: sine mode k is ODD-parity index if k odd, EVEN if k even
odd_idx = [i for i, k in enumerate(modes) if k % 2 == 1]
even_idx = [i for i, k in enumerate(modes) if k % 2 == 0]
Odd = Mf[np.ix_(odd_idx, odd_idx)]
Even = Mf[np.ix_(even_idx, even_idx)]
Cross = Mf[np.ix_(odd_idx, even_idx)]

def blknorm(X): return float(np.linalg.norm(X))
print("\n=== BLOCK STRUCTURE ===")
print(f" odd-block  ({len(odd_idx)}x{len(odd_idx)}) Frobenius norm = {blknorm(Odd):.6e}")
print(f" even-block ({len(even_idx)}x{len(even_idx)}) Frobenius norm = {blknorm(Even):.6e}")
print(f" CROSS-block(odd,even) Frobenius norm = {blknorm(Cross):.6e}")
print(f" cross / sqrt(odd*even) diag-scale ratio = {blknorm(Cross)/np.sqrt(blknorm(Odd)*blknorm(Even)):.6e}")
print(f" max |cross entry| = {np.max(np.abs(Cross)):.6e}   max|diag-block entry| = {max(np.max(np.abs(Odd)),np.max(np.abs(Even))):.6e}")

print("\n=== EIGENVALUES (float midpoint; REPORT tier) ===")
print(" FULL   :", np.sort(np.linalg.eigvalsh(Mf)))
print(" ODD    :", np.sort(np.linalg.eigvalsh(Odd)))
print(" EVEN   :", np.sort(np.linalg.eigvalsh(Even)))
print(f"\n FULL min-eig = {np.min(np.linalg.eigvalsh(Mf)):.6e}   PD? {np.all(np.linalg.eigvalsh(Mf)>0)}")
print(f" ODD  min-eig = {np.min(np.linalg.eigvalsh(Odd)):.6e}   PD? {np.all(np.linalg.eigvalsh(Odd)>0)}")
print(f" EVEN min-eig = {np.min(np.linalg.eigvalsh(Even)):.6e}   PD? {np.all(np.linalg.eigvalsh(Even)>0)}")
