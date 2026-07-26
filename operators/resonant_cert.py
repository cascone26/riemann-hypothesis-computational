#!/usr/bin/env python3
"""§338 experiment: RESONANT-EXACT certificate chains.

Split P = P_res + P_off:  P_res = the Bragg-resonant prime-power term(s) (treated EXACTLY inside
the known part M = A - P_res, still zero-free/explicit),  P_off = all other prime powers (bounded).

Chains (rigorous lower bounds on lambda_min(A - P), all ingredients zero-free):
  c3norm : lambda_min(A - P_res) - lambda_max(P_off)        [resonance exact, remainder by norm]
  c3gap  : 2-block gap bound of M = A - P_res vs Q = P_off  [resonance exact + directional]
Baselines from §337: gersh (5.7), norm (5.5), directional k=3 (4.9).
Selection of resonant terms: top-R prime powers n by ||P_n||_2 (adaptive; check vs Bragg law 2pi k/log n).
"""
import numpy as np
import l2_edge_sweep as L

w = 3.0

def prime_terms(cs):
    """dict n -> P_n matrix, for prime powers n<=NMAX."""
    out = {}
    for i, n in enumerate(L.NS):
        if L.LAMW[i] == 0.0:
            continue
        u = L.LOGN[i]
        M = np.zeros((len(cs), len(cs)))
        for j in range(len(cs)):
            for k in range(j, len(cs)):
                v = 2.0 * L.LAMW[i] * L.g_pair(cs[j], cs[k], w, u)
                M[j, k] = M[k, j] = v
        if np.max(np.abs(M)) > 1e-14:
            out[int(n)] = M
    return out

def gap_bound(M, Q, k=3):
    """rigorous: lambda_min(M - Q) >= min eig [[lmin(L_b - Q_bb), -||Q_off||], [-||Q_off||, l_{k+1}(M) - lmax(Q_rest)]]"""
    n = M.shape[0]
    k = min(k, n - 1)
    lM, UM = np.linalg.eigh(M)
    Ub, Ur = UM[:, :k], UM[:, k:]
    top = float(np.linalg.eigvalsh(np.diag(lM[:k]) - Ub.T @ Q @ Ub)[0])
    q = float(np.linalg.norm(Ub.T @ Q @ Ur, 2))
    pr = float(np.linalg.eigvalsh(Ur.T @ Q @ Ur)[-1])
    M2 = np.array([[top, -q], [-q, lM[k] - pr]])
    return float(np.linalg.eigvalsh(M2)[0])

print("RESONANT-EXACT CHAINS (band [14,44], w=3). R = # resonant terms treated exactly.")
print(f"{'s':>5} {'n':>3} {'true':>10} {'n*(R=1)':>8} {'lmax(P)':>8} {'lmax(Poff1)':>11} "
      f"{'c3norm R1':>10} {'c3norm R2':>10} {'c3gap R1':>10} {'c3gap R2':>10} {'c3gap R3':>10}")
rows = []
for s in np.arange(6.5, 3.39, -0.1):
    cs = L.probes(L.BAND, s)
    f = L.build_forms(cs, w)
    A = f['archpole']; P = f['prime']
    true_m = float(np.linalg.eigvalsh(f['full'])[0])
    terms = prime_terms(cs)
    order = sorted(terms, key=lambda n: -np.linalg.norm(terms[n], 2))
    row = {'s': s, 'true': true_m, 'nstar': order[0],
           'lmaxP': float(np.linalg.eigvalsh(P)[-1])}
    for R in (1, 2, 3):
        res = order[:R]
        Pres = sum(terms[n] for n in res)
        Poff = P - Pres
        M = A - Pres
        row[f'c3norm{R}'] = float(np.linalg.eigvalsh(M)[0] - np.linalg.eigvalsh(Poff)[-1])
        row[f'c3gap{R}'] = gap_bound(M, Poff, k=3)
        if R == 1:
            row['lmaxPoff1'] = float(np.linalg.eigvalsh(Poff)[-1])
    rows.append(row)
    print(f"{s:5.1f} {len(cs):3d} {true_m:+10.3e} {row['nstar']:8d} {row['lmaxP']:8.4f} "
          f"{row['lmaxPoff1']:11.4f} {row['c3norm1']:+10.3e} {row['c3norm2']:+10.3e} "
          f"{row['c3gap1']:+10.3e} {row['c3gap2']:+10.3e} {row['c3gap3']:+10.3e}", flush=True)

print("\nreach (deepest s certified positive):")
for key, label in [('c3norm1', 'norm, R=1 exact'), ('c3norm2', 'norm, R=2 exact'),
                   ('c3gap1', 'gap k=3 + R=1 exact'), ('c3gap2', 'gap k=3 + R=2 exact'),
                   ('c3gap3', 'gap k=3 + R=3 exact')]:
    ok = [r['s'] for r in rows if r[key] > 0]
    print(f"  {label:22s}: s = {min(ok) if ok else '---':}")
print("\n(§337 baselines: gersh 5.7 | norm 5.5 | gap k=3 4.9)")
print("Bragg check — n* selected vs 2pi k / log n resonance prediction:")
for r in rows[::4]:
    s = r['s']; n = r['nstar']
    pred = 2 * np.pi / np.log(n)
    print(f"  s={s:.1f}: n*={n}  (2pi/log n* = {pred:.2f}; harmonics {pred/2:.2f}, {pred*2:.2f})")
