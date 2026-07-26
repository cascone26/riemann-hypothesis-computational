#!/usr/bin/env python3
"""§342a: W-SWEEP of the adaptive chain — first face of the uniformity wall.
For w in {1.5, 2, 2.5, 3, 4}: does reach-to-edge persist, and how do the chain ingredients scale?
Expected structure: prime coupling e^{-(log n)^2 w^2/4} strengthens as w shrinks => core grows.
Uniformity question in measurable form: core size + tail norm + reach as functions of w.
NMAX extended for small w."""
import numpy as np
import l2_edge_sweep as L
from band_sweep_cert import gap_bound

# extend prime cutoff for small w
NMAX2 = 3000
lam = np.zeros(NMAX2 + 1)
comp = np.zeros(NMAX2 + 1, dtype=bool)
for p in range(2, NMAX2 + 1):
    if not comp[p]:
        comp[2 * p::p] = True
        lp = np.log(p)
        pk = p
        while pk <= NMAX2:
            lam[pk] = lp
            pk *= p
L.NS = np.arange(2, NMAX2 + 1)
L.LOGN = np.log(L.NS)
L.LAMW = lam[2:] / np.sqrt(L.NS)

def prime_terms_w(cs, w, thresh=1e-15):
    out = {}
    for i, n in enumerate(L.NS):
        if L.LAMW[i] == 0.0:
            continue
        u = L.LOGN[i]
        scale = L.LAMW[i] * np.exp(-u * u * w * w / 4)
        if scale < thresh:
            continue
        M = np.zeros((len(cs), len(cs)))
        for j in range(len(cs)):
            for k in range(j, len(cs)):
                v = 2.0 * L.LAMW[i] * L.g_pair(cs[j], cs[k], w, u)
                M[j, k] = M[k, j] = v
        out[int(n)] = M
    return out

print("W-SWEEP: adaptive chain (gap k=3, growing core) vs w.  band [14,44].")
print(f"{'w':>4} {'s':>5} {'n':>3} {'#pp':>4} {'true':>10} {'chain':>10} {'R':>3} {'core_max_n':>10}")
for w in (1.5, 2.0, 2.5, 3.0, 4.0):
    # rank/Nyquist edge estimate: zeros in [14-2w, 44+2w]
    lo, hi = 14 - 2 * w, 44 + 2 * w
    nz = int(np.sum((L.ZEROS >= lo) & (L.ZEROS <= hi)))
    reach = None
    lastR = None
    for s in np.arange(6.0, 1.19, -0.1):
        cs = L.probes(L.BAND, s)
        if len(cs) < 4:
            continue
        f = L.build_forms(cs, w)
        true_m = float(np.linalg.eigvalsh(f['full'])[0])
        if true_m < 1e-13:
            print(f"{w:4.1f}  rank limit at s={s:.1f}; last certified s={reach}, R={lastR[0] if lastR else '-'}, "
                  f"core max n={lastR[2] if lastR else '-'}; Nyquist edge ~{30/nz:.2f}, Nz={nz}")
            break
        terms = prime_terms_w(cs, w)
        order = sorted(terms, key=lambda n: -np.linalg.norm(terms[n], 2))
        A = f['archpole']; P = f['prime']
        best = None
        for R in range(1, min(16, len(order) + 1)):
            Pres = sum(terms[n] for n in order[:R])
            b = gap_bound(A - Pres, P - Pres, k=3)
            if best is None or b > best[1]:
                best = (R, b, max(order[:R]))
            if b > 0:
                break
        if best[1] > 0:
            reach = s
            lastR = best
        else:
            # first failure: report and stop descending
            print(f"{w:4.1f} {s:5.1f} {len(cs):3d} {len(terms):4d} {true_m:+10.3e} {best[1]:+10.3e} "
                  f"{best[0]:3d}  FAIL (reach was {reach}, R={lastR[0] if lastR else '-'}"
                  f", core max n={lastR[2] if lastR else '-'}; Nyquist edge ~{30/nz:.2f}, Nz={nz})")
            break
    else:
        print(f"{w:4.1f}  ---  reached rank limit; last certified s={reach}, R={lastR[0]}, "
              f"core max n={lastR[2]}; edge ~{30/nz:.2f}")
        continue
