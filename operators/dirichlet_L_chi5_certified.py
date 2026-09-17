#!/usr/bin/env python3
"""
cycle 346 — extend the AF-realization discriminator to a NEW primitive Dirichlet L: L(s, chi_5), the EVEN
quadratic character mod 5 (Legendre (n/5), conductor 5, nu=1/4). Different parity (EVEN) and conductor than
the certified chi_4 case (odd, cond 4). GRH holds ⟹ PREDICT n_neg=0 (multiplicative Euler product forces PSD).
A 0 broadens the Handle-tier validation of "multiplicativity ⟹ n_neg=0 ⟹ clean AF 2/3" to an even/cond-5 L.
"""
import sys, os
sys.path.insert(0, os.path.join('..','external','groskin_2607.02828','code'))
from flint import arb, ctx
from arb_ldlt_certify import prime_powers_up_to
from dh_exact import kj_exact
from dh_certified_tau import build_tau_general, build_arb_tau, margin

def chi5(n):
    r = n % 5
    if r == 0: return 0
    return 1 if r in (1,4) else -1     # Legendre (n/5): QR={1,4} -> +1, {2,3} -> -1 ; EVEN (chi(-1)=chi(4)=+1)

def Lchi5_prime_side(c, prec):
    ctx.prec = prec
    pdata = prime_powers_up_to(c)
    w = []; pos = []
    for (q, p) in pdata:
        cq = chi5(q)                    # completely multiplicative: chi(p^k)=chi(p)^k
        if cq == 0: continue            # p=5 drops
        w.append(arb(cq)*arb(p).log()*(arb(q)**arb('-0.5')))
        pos.append(arb(q).log())
    return w, pos

prec, N = 800, 64
print(f"CERTIFIED n_neg for L(s,chi_5) (EVEN quadratic, cond 5, nu=1/4, Euler product, GRH holds) vs zeta, N={N}:")
print(f"{'c':>4} | {'zeta n_neg':>10} | {'L(chi5) n_neg':>13} | {'(k+J) used':>12}")
for c in [25, 50, 100]:
    ctx.prec = prec
    nz, _, _ = margin(*build_arb_tau(c, N, prec))
    kj = float(kj_exact('0.25', 5, c, prec).mid())    # nu=1/4 (even chi), conductor 5
    w, pos = Lchi5_prime_side(c, prec)
    A, DIM = build_tau_general(c, N, prec, nu=0.25, include_pole=False, kj_const=kj,
                               prime_weights=w, prime_positions=pos)
    nn, pd, undet = margin(A, DIM)
    print(f"{c:>4} | {nz:>10} | {nn:>13} | {kj:>12.6f}", flush=True)
print("PREDICT n_neg=0 (GRH). If 0 across c ⟹ multiplicativity forces PSD for an EVEN cond-5 primitive L too,")
print("broadening the AF-realization validation (multiplicativity ⟹ n_neg=0 ⟹ clean 2/3) beyond chi_4.")
