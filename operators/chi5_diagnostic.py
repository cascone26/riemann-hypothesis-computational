#!/usr/bin/env python3
"""cy347 — diagnose chi_5 n_neg>0: artifact (tiny negatives, even-char near-nullity) or bug/real (O(1))?
Replicates kj_exact inline (avoids dh_exact's unguarded import side-effect). Prints margin magnitude."""
import sys, os
sys.path.insert(0, os.path.join('..','external','groskin_2607.02828','code'))
from flint import arb, ctx
from arb_ldlt_certify import prime_powers_up_to
from dh_certified_tau import build_tau_general, build_arb_tau, margin   # __main__-guarded, no side effect

def kj_inline(nu, cond, c, prec):
    ctx.prec = prec; L = arb(c).log(); nu = arb(nu)
    s = arb(0); k = 0; thr = arb(2)**(-(prec-8))
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

prec, N = 800, 64
print("chi_5 (even, cond5, nu=1/4) n_neg + MARGIN MAGNITUDE (artifact if |mlow|~1e-8, bug/real if O(1)):")
print(f"{'c':>4} | {'zeta n_neg':>10} | {'chi5 n_neg':>10} | {'chi5 mlow (smallest pivot lb)':>30}")
for c in [25, 50]:
    ctx.prec = prec
    nz, _, _ = margin(*build_arb_tau(c, N, prec))
    kj = kj_inline('0.25', 5, c, prec)
    w, pos = chi5_prime_side(c, prec)
    A, DIM = build_tau_general(c, N, prec, nu=0.25, include_pole=False, kj_const=float(kj.mid()),
                               prime_weights=w, prime_positions=pos)
    nn, ok, mlow = margin(A, DIM)
    print(f"{c:>4} | {nz:>10} | {nn:>10} | {float(mlow.mid() if hasattr(mlow,'mid') else mlow):>30.3e}", flush=True)
print("\nverdict: |mlow|~1e-8 => even-char finite-truncation near-null artifact (discriminator needs higher N/prec")
print("for even chars, NOT a GRH violation). |mlow|~O(0.01+) => genuine n_neg (bug in my even-char setup or real).")
