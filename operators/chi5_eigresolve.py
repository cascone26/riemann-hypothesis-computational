#!/usr/bin/env python3
"""cy348 — RESOLVE cy346-347: compute the ACTUAL eigenvalues of the χ_5 even-char Weil matrix (not the LDLᵀ
inertia count). If the 3 'negatives' are ~1e-9 (resolution-floor near-nulls) => true n_neg=0 (GRH), the
even-char discriminator just needs higher prec. If O(0.01+) => a genuine even-char n_neg (real finding/bug)."""
import sys, os, numpy as np
sys.path.insert(0, os.path.join('..','external','groskin_2607.02828','code'))
from flint import arb, ctx
from arb_ldlt_certify import prime_powers_up_to
from dh_certified_tau import build_tau_general

def kj_inline(nu, cond, c, prec):
    ctx.prec = prec; L = arb(c).log(); nu = arb(nu); s = arb(0); k = 0; thr = arb(2)**(-(prec-8))
    while True:
        t = (-2*(k+nu)*L).exp()/(k+nu); s += t
        if t < thr and k > 2: break
        k += 1
    return -nu.digamma() + (arb.pi()/cond).log() - s

def chi5(n):
    r = n % 5
    return 0 if r == 0 else (1 if r in (1,4) else -1)

def chi5_prime_side(c, prec):
    ctx.prec = prec; w = []; pos = []
    for (q, p) in prime_powers_up_to(c):
        cq = chi5(q)
        if cq == 0: continue
        w.append(arb(cq)*arb(p).log()*(arb(q)**arb('-0.5'))); pos.append(arb(q).log())
    return w, pos

prec, N, c = 2000, 64, 25
ctx.prec = prec
kj = kj_inline('0.25', 5, c, prec)
w, pos = chi5_prime_side(c, prec)
A, DIM = build_tau_general(c, N, prec, nu=0.25, include_pole=False, kj_const=float(kj.mid()),
                           prime_weights=w, prime_positions=pos)
# extract arb_mat -> float symmetric matrix
M = np.zeros((DIM, DIM))
for i in range(DIM):
    for j in range(DIM):
        M[i, j] = float(A[i, j].mid())
M = 0.5*(M + M.T)
ev = np.linalg.eigvalsh(M)
neg = ev[ev < 0]
print(f"χ_5 even, c={c}, N={N}, DIM={DIM}: {len(neg)} negative eigenvalues (float, from Arb mids at prec{prec})")
print(f"  smallest 8 eigenvalues: {np.array2string(ev[:8], formatter={'float':lambda x:f'{x:+.3e}'})}")
print(f"  most-negative = {ev[0]:+.3e}   |most-neg| = {abs(ev[0]):.3e}")
if len(neg):
    print(f"  all negatives: {np.array2string(neg, formatter={'float':lambda x:f'{x:+.2e}'})}")
print(f"\nVERDICT: |most-neg| = {abs(ev[0]):.2e} => "
      + ("RESOLUTION-FLOOR near-null artifact (true n_neg=0, GRH); even-char discriminator needs higher prec"
         if abs(ev[0]) < 1e-6 else
         "GENUINE O(1) negative — a real even-char n_neg (not a floor artifact); needs deeper look"))
