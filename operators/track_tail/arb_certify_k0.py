#!/usr/bin/env python3
"""arb_certify_k0.py — generalized §199.3c Arb certification for K0×K0 tail-lemma block.
Extends operators/dbn_dynamics/arb_certify.py (which is hardcoded 4×4, K0 unspecified).
Certifies the low block {odd modes 1,3,...,K0} rigorously PD in Arb intervals (LDL^T pivots).

STATUS/PROVENANCE (corrected 2026-07-24, roundtable session — supersedes the two stale claims below):
  - STALE CLAIM 1 (removed): "K0=7 ... proves even-sector infinite-dim Weil positivity via
    tail_delta_k0_lemma.md δ(K0)=0.4552". REFUTED at §226: at K0=7 the near-null Schur residual
    (1.783e-7) EXCEEDS the low-block margin (1.652e-7) — composition FAILS at K0=7. True minimal
    K0 is 9 (§227); current banking target is K0=23 (Schur headroom ~9,500x, tail-margin floor
    ~1.22 — see results/k0_hp_nl24000_out.txt).
  - STALE CLAIM 2 (reworded): "K0=13 fails at NL=8000" was true but PRECISION-FLOOR-SPECIFIC
    (mode-1 pivot radius ±1.19e-7 exceeds the ~7e-8 margin at NL=8000). At NL=50000 modes {1..17}
    certified at radius ±2.7e-9 (§227, arb_certify_multiprime). Do NOT read that note as a
    refutation of the block's PD.
  - tail_delta_k0_lemma.md (v1) is DEPRECATED -> see operators/track_tail/tail_lemma_v2_nearband_identity.md.
    This tool certifies the LOW BLOCK ONLY. Infinite-dim composition additionally requires the v2
    tail lemma (obligations 1-4, OPEN) + an Arb-interval Schur composition. The success banner
    below states exactly this and no more.

Usage:
  python3 arb_certify_k0.py [--k0 23] [--nl 50000] [--a 0.549]
"""
import sys, argparse
from flint import arb, acb, ctx
ctx.prec = 1400

parser = argparse.ArgumentParser()
parser.add_argument('--k0', type=int, default=7, help='K0 (largest mode in low block; default=7 is minimal operational)')
parser.add_argument('--nl', type=int, default=8000, help='NLERCH (Lerch series terms)')
parser.add_argument('--a', type=str, default='0.549', help='window parameter a')
args = parser.parse_args()

K0 = args.k0
NLERCH = args.nl
A = arb(args.a)
B = arb.pi() / (2*A)
LOG2 = arb(2).log()
TWOA = 2*A

OM0 = (arb(1)/4).digamma() - arb.pi().log()
PHI0 = acb(2).zeta(acb(arb(1)/4)).real
WP = arb(2).log() / arb(2).sqrt()

def Pp(p, t, s):
    if p == 0: return 1/s
    if p == 1: return t/s - 1/s**2
    return t*t/s - 2*t/s**2 + 2/s**3

def II(p, alpha, omega, sc, t0, t1):
    """∫_{t0}^{t1} t^p e^{-alpha t}{cos/sin}(omega t) dt  (rigorous acb)."""
    s = acb(-alpha, omega)
    def F(tt):
        t = acb(tt)
        return (s*t).exp() * Pp(p, t, s)
    v = F(t1) - F(t0)
    return v.real if sc == 'cos' else v.imag

def dterms(k, l):
    """D_kl closed-form terms (as amplitude, power-of-t, frequency, cos/sin)."""
    kb, lb = arb(k)*B, arb(l)*B
    if k == l:
        return [(2*A*arb(k)**2*B**2, 0, kb, 'cos'),
                (-arb(k)**2*B**2, 1, kb, 'cos'),
                (-arb(k)*B, 0, kb, 'sin')]
    al = arb(k)*arb(l)*B * (arb(1)/(k-l) - arb(1)/(k+l))
    ak = -arb(k)*arb(l)*B * (arb(1)/(k-l) + arb(1)/(k+l))
    return [(al, 0, arb(l)*B, 'sin'), (ak, 0, arb(k)*B, 'sin')]

def Dmax(k, l):
    """Bound max|D_kl(t)| on [0,2a] for rigorous Lerch tail."""
    tot = arb(0)
    for (ad, pd, om, sc) in dterms(k, l):
        tot += abs(ad) * (TWOA if pd == 1 else arb(1))
    return tot

def Q(k, l):
    """Compute Q[k,l] rigorously in arb intervals with Lerch tail bound."""
    tot = arb(0)
    D = dterms(k, l)

    # g smooth part: boundary + arch_core + (-Phi0/4) + Lerch(NLERCH terms)
    gt = [(arb(-4), 0, arb(-1)/2),
          (arb(-4), 0, arb(1)/2),
          (arb(8), 0, arb(0)),
          (-OM0/2, 1, arb(0)),
          (-PHI0/4, 0, arb(0))]
    for n in range(NLERCH):
        gt.append((arb(1)/4/(arb(n)+arb(1)/4)**2, 0, 2*arb(n)+arb(1)/2))

    # Integrate g·D_kl over [0,2a]
    for (ag, pg, alp) in gt:
        for (ad, pd, om, sc) in D:
            tot += ag * ad * II(pg + pd, alp, om, sc, arb(0), TWOA)

    # Lerch rigorous tail (n >= NLERCH): |term| <= (1/4)(n+1/4)^-2 * Dmax/(2n+1/2); Σ_{n>=N} n^-3 <= 1/(2(N-1)^2)
    N = NLERCH
    dm = Dmax(k, l)
    tailbound = dm/8 * (arb(1)/(2*(arb(N)-1)**2))
    tot += tailbound * arb(0, 1)  # arb(0,1) = [-1,1] encodes ±tailbound

    # Prime part on [log2, 2a]
    for (ad, pd, om, sc) in D:
        tot += WP*ad*II(pd+1, arb(0), om, sc, LOG2, TWOA) + (-WP*LOG2)*ad*II(pd, arb(0), om, sc, LOG2, TWOA)

    return tot

# Build K0 block with odd modes 1,3,5,…,K0
modes_all = [k for k in range(1, K0+1, 2)]  # odd only
# Sort by diagonal (large-first) to avoid LDL^T interval amplification
diag_vals = [(Q(k, k), k) for k in modes_all]
modes = [k for _, k in sorted(diag_vals, key=lambda x: -float(x[0].mid()))]

print(f"§202-tail RIGOROUS Arb certification, K0={K0} (even sector, odd modes), NLERCH={NLERCH}")
print(f"  a={A}, prec={ctx.prec}, modes (large-diagonal-first): {modes}")
print()

n = len(modes)
M = [[Q(k, l) for l in modes] for k in modes]
idx = {k: i for i, k in enumerate(modes)}

# LDL^T decomposition with rigorous interval pivots
d = [None] * n
L = [[arb(0)] * n for _ in range(n)]
ok = True
piv = []

for i in range(n):
    s = M[i][i]
    for kk in range(i):
        s = s - L[i][kk]**2 * d[kk]
    d[i] = s
    cert = (s > 0)
    piv.append((modes[i], s.str(6), cert))
    if not cert:
        ok = False
    for j in range(i+1, n):
        t = M[j][i]
        for kk in range(i):
            t = t - L[j][kk] * L[i][kk] * d[kk]
        L[j][i] = t / d[i]

print(f"LDL^T pivots (large-diagonal-first order):")
for k, sv, c in piv:
    print(f"  mode {k:2d}: {sv}  rigorously>0? {c}")

print()
if ok:
    print(f"*** CERTIFIED PD (all {n} pivots lower-bound >0) ***")
    print(f"*** K0={K0} LOW BLOCK rigorously PD (Arb, NLERCH={NLERCH}). ***")
    print(f"*** NOT a positivity proof by itself: infinite-dim composition still requires ***")
    print(f"*** the v2 tail lemma (near-band identity + row-sum bound, OPEN) + interval Schur. ***")
    sys.exit(0)
else:
    print("NOT CERTIFIED — some pivot(s) lower-bound ≤ 0")
    sys.exit(1)
