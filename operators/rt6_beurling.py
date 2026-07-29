#!/usr/bin/env python3
"""§391 (angle 9): is collapse-at-exact-RvM GENERIC among pair-matched Lambda-systems or prime-special?
Surrogates: (a) IAAFT on the dLambda field (preserves marginal+spectrum=pair stats); (b) Poisson-random
primes (1-point only). For each surrogate build B=pole+arch-P_surr, find collapse fill n* (W/ovl<0.02,
§373 detector) on [14,44] w=2. Compare to true n*=8 and RvM_eff=8.75."""
import numpy as np
import l2_edge_sweep as L
rng = np.random.default_rng(42)
w = 2.0
NM = 20000
lam = np.zeros(NM + 1); comp = np.zeros(NM + 1, dtype=bool)
for p in range(2, NM + 1):
    if not comp[p]:
        comp[2*p::p] = True; lp = np.log(p); pk = p
        while pk <= NM: lam[pk] = lp; pk *= p
L.NS = np.arange(2, NM + 1); L.LOGN = np.log(L.NS)
TRUE_LAMW = lam[2:] / np.sqrt(L.NS)

def collapse_fill(lamw):
    L.LAMW = lamw
    for n in range(4, 30):
        cs = np.linspace(14, 44, n)
        f = L.build_forms(cs, w)
        mo = float(np.linalg.eigvalsh(f['overlap'])[0]); mf = float(np.linalg.eigvalsh(f['full'])[0])
        if mo > 1e-300 and mf/mo < 0.02:
            return n
    return None

def iaaft(x, it=60):
    xs = np.sort(x); amp = np.abs(np.fft.rfft(x)); y = rng.permutation(x)
    for _ in range(it):
        Y = np.fft.rfft(y); y = np.fft.irfft(amp*np.exp(1j*np.angle(Y)), n=len(x))
        y = xs[np.argsort(np.argsort(y))]
    return y

print("TRUE Lambda:", collapse_fill(TRUE_LAMW), " (RvM_eff=8.75)")
dl = lam[2:] - 1.0
fills_i = []
for k in range(8):
    surr = iaaft(dl)
    lw = np.maximum(surr + 1.0, 0.0) / np.sqrt(L.NS)
    fills_i.append(collapse_fill(lw))
print("IAAFT pair-matched surrogates n*:", fills_i)
fills_p = []
for k in range(8):
    fake = np.zeros(NM + 1)
    mask = rng.random(NM + 1) < 1/np.log(np.maximum(np.arange(NM+1), 3))
    fake[mask] = np.log(np.maximum(np.arange(NM+1)[mask], 2))
    fills_p.append(collapse_fill(fake[2:]/np.sqrt(L.NS)))
print("Poisson 1-point surrogates n*:", fills_p)
print("verdict: surrogate n* == 8 mostly => collapse fill is PAIR-GENERIC (upgrade path);")
print("n* shifted/absent => prime-special (deep-arithmetic corner).")
