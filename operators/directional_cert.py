#!/usr/bin/env python3
"""§337 experiment: DENSITY + DIRECTIONAL certificate chains for the Weil form.

Exact split (zero-free on both sides): B_full = A - P,  A = arch + pole,  P = prime.
Fluctuation identity: B_zero - A = -P exactly (explicit formula), so "directional control of the
fluctuation on the minimizer" = "directional control of the prime form".

Chains (all ingredients computable WITHOUT zeros; each is a rigorous lower bound on lambda_min(B_full)):
  chain0: Gershgorin on B_full (the §285-286 crude cert; saturates ~2 widths)
  chain1: lambda_min(A) - lambda_max(P)                      [operator-norm / §277-shape]
  chain2: 2-block gap bound: decompose on A's bottom eigenvector u1:
          lambda_min(A-P) >= min eig [[l1 - p11, -q], [-q, l2 - pr]]
          l1,l2 = two lowest eigs of A; p11 = <u1,P u1>; q = ||P u1 - p11 u1||;
          pr = lambda_max(P restricted to u1-perp).            [directional refinement]
  chain2k: same with a k-dim bottom block of A (k=2,3): min eig of [[L_b - P_bb, -Q],[-Q, l_{k+1} - pr]]
          done rigorously as: lambda_min of [[lambda_min(L_b - P_bb), -||P_off||_2],[-||P_off||_2, l_{k+1}-pr]].

Cancellation measurement: per-n equidistribution gain on the A-minimizer:
  gain(n) = |sum_jk v_j v_k e^{-D^2/4w^2} cos(u_n m_jk)| / sum_jk |v_j v_k| e^{-D^2/4w^2}
(the factor a future a-priori bound must supply; triangle bound = gain 1).
"""
import numpy as np
import l2_edge_sweep as L

w = 3.0

def chains(s):
    cs = L.probes(L.BAND, s)
    f = L.build_forms(cs, w)
    A = f['archpole']; P = f['prime']; B = f['full']
    n = len(cs)
    true_m = float(np.linalg.eigvalsh(B)[0])
    g0 = L.diag_margin(B)
    lA, UA = np.linalg.eigh(A)
    lP = np.linalg.eigvalsh(P)
    c1 = lA[0] - lP[-1]
    out = {'s': s, 'n': n, 'true': true_m, 'g0': g0, 'c1': c1,
           'lA1': lA[0], 'lPmax': lP[-1]}
    for k in (1, 2, 3):
        if k >= n: break
        Ub = UA[:, :k]; Ur = UA[:, k:]
        Pbb = Ub.T @ P @ Ub
        Poff = Ub.T @ P @ Ur
        Prr = Ur.T @ P @ Ur
        Lb = np.diag(lA[:k])
        top = float(np.linalg.eigvalsh(Lb - Pbb)[0])
        q = float(np.linalg.norm(Poff, 2))
        pr = float(np.linalg.eigvalsh(Prr)[-1])
        M2 = np.array([[top, -q], [-q, lA[k] - pr]])
        out[f'c2k{k}'] = float(np.linalg.eigvalsh(M2)[0])
    # cancellation gains on A-minimizer
    v = UA[:, 0]
    gains = {}
    for nn in (2, 3, 4, 5, 7, 8, 9):
        u = np.log(nn)
        m = (cs[None, :] + cs[:, None]) / 2.0
        d = cs[None, :] - cs[:, None]
        ker = np.exp(-d**2 / (4 * w**2))
        num = abs(float(np.sum(np.outer(v, v) * ker * np.cos(u * m))))
        den = float(np.sum(np.abs(np.outer(v, v)) * ker))
        gains[nn] = num / den
    out['gains'] = gains
    out['p11'] = float(v @ P @ v)
    return out

print(f"DIRECTIONAL CERTIFICATE CHAINS  (band [14,44], w=3; all chain ingredients zero-free)")
print(f"{'s':>5} {'n':>3} {'true':>10} {'gersh':>8} {'c1:normA-P':>10} {'c2 k=1':>9} {'c2 k=2':>9} {'c2 k=3':>9} {'p11':>8}")
rows = []
for s in np.arange(6.5, 3.39, -0.1):
    o = chains(s)
    rows.append(o)
    print(f"{o['s']:5.1f} {o['n']:3d} {o['true']:+10.3e} {o['g0']:+8.3f} {o['c1']:+10.3e} "
          f"{o.get('c2k1', float('nan')):+9.3e} {o.get('c2k2', float('nan')):+9.3e} "
          f"{o.get('c2k3', float('nan')):+9.3e} {o['p11']:+8.4f}", flush=True)

print("\nreach of each chain (deepest s with positive certified margin):")
for key, label in [('g0', 'chain0 Gershgorin'), ('c1', 'chain1 norm'),
                   ('c2k1', 'chain2 k=1'), ('c2k2', 'chain2 k=2'), ('c2k3', 'chain2 k=3')]:
    ok = [o['s'] for o in rows if o.get(key, -1) is not None and o.get(key, -1) > 0]
    print(f"  {label:18s}: reaches s = {min(ok) if ok else '---'}")

print("\nper-n cancellation gain on the A-minimizer (1.0 = no cancellation = triangle bound):")
print(f"{'s':>5} " + " ".join(f"n={nn:<4d}" for nn in (2, 3, 4, 5, 7, 8, 9)))
for o in rows[::3]:
    print(f"{o['s']:5.1f} " + " ".join(f"{o['gains'][nn]:6.3f}" for nn in (2, 3, 4, 5, 7, 8, 9)))
