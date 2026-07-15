#!/usr/bin/env python3
"""
P_offdiag_resonance.py — Which primes drive the largest off-diagonal P coupling
between A+B eigenmodes?

Method:
  1. Build A, B, P at N=64 (and N=96 for comparison).
  2. Diagonalize A+B: get V (columns = eigenvectors), eigenvalues a_j.
  3. Compute P_tilde = V^T P V  (P in the A+B eigenbasis).
  4. Find pairs (j,k), j<k with largest |P_tilde[j,k]|.
  5. For each such pair, decompose P_tilde[j,k] prime-by-prime:
       P_tilde[j,k] = Σ_p  w_p * (V_j^T · P_p · V_k)
     where P_p is the contribution of prime power p to P.
  6. Show: diagonal gap a_j - P_tilde[j,j] vs off-diagonal |P_tilde[j,k]|.
     The off-diagonal mixing "costs" gap/(off-diag)^2 in second-order perturbation theory.

Key question: is there a particular prime that creates an anomalously large resonance
between specific A+B modes?
"""

import json, os
import numpy as np
from scipy.special import digamma as _digamma
import mpmath as mp

mp.mp.dps = 30
TWO_PI = 2.0 * np.pi
LOGPI  = float(np.log(np.pi))
LAMBDA_SQ = 14
L  = float(mp.log(LAMBDA_SQ))
sh = float(np.sinh(L/4)**2)

OUTFILE = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                       "../results/P_offdiag_resonance.json")


def vm_table(kmax):
    tab = []
    for k in range(2, kmax + 1):
        n, p, f, ok = k, 0, 2, True
        while f * f <= n:
            if n % f == 0:
                if p == 0: p = f
                elif p != f: ok = False; break
                while n % f == 0: n //= f
            f += 1
        if n > 1:
            if p == 0: p = n
            elif p != n: ok = False
        if ok and p > 0:
            tab.append((k, float(mp.log(p))))
    return tab


def build_all(N):
    VM  = vm_table(LAMBDA_SQ)
    idx = list(range(-N, N+1))
    D   = len(idx)
    R   = TWO_PI * N / L + 80.0
    G   = 60000
    grid = np.linspace(-R, R, G)
    h    = grid[1] - grid[0]
    w    = np.full(G, h); w[0] = w[-1] = 0.5*h
    Om   = np.real(_digamma(0.25 + 0.5j * grid)) - LOGPI

    E = np.empty((D, G), dtype=complex)
    for ki, n in enumerate(idx):
        s  = TWO_PI * n / L
        d  = s - grid
        sm = np.abs(d) < 1e-12
        E[ki, sm]  = np.sqrt(L)
        dd         = d[~sm]
        E[ki, ~sm] = (1/np.sqrt(L)) * np.exp(0.5j*dd*L) * 2*np.sin(0.5*dd*L) / dd

    arch = ((E * (w * Om)) @ E.conj().T).real / (2*np.pi)
    ns   = np.array(idx, dtype=float)
    B    = (32*L*sh * (L**2 - 16*np.pi**2*np.outer(ns,ns)) /
            np.outer(L**2+16*np.pi**2*ns**2, L**2+16*np.pi**2*ns**2))

    log_k = np.array([np.log(k) for k, _ in VM])
    wts   = np.array([lk * k**(-0.5) for k, lk in VM])

    # Build P_per_prime: one matrix per prime power
    P_per_pp = {}   # key = (k, log_p) -> (D,D) matrix
    P = np.zeros((D, D))
    for pp_idx, (k, lk) in enumerate(VM):
        Ppp = np.zeros((D, D))
        for ii, n in enumerate(idx):
            for jj in range(ii, D):
                m = idx[jj]
                q = np.where(n == m,
                    2*(1-np.abs(log_k[pp_idx])/L)*np.cos(TWO_PI*n*log_k[pp_idx]/L),
                    (np.sin(TWO_PI*m*log_k[pp_idx]/L)-np.sin(TWO_PI*n*log_k[pp_idx]/L))/(np.pi*(n-m)+1e-300))
                val = float(wts[pp_idx] * q)
                Ppp[ii,jj] = Ppp[jj,ii] = val
        P_per_pp[(k, round(lk,6))] = Ppp
        P += Ppp

    tau = arch + B - P
    return (tau+tau.T)/2, arch, B, P, P_per_pp, idx, VM, wts


def main():
    for N in [64, 96]:
        print(f"\n{'='*70}")
        print(f"N = {N}")
        tau, arch, B, P, P_per_pp, idx, VM, wts = build_all(N)

        # Diagonalize A+B
        AB = arch + B
        ab_ev, V = np.linalg.eigh((AB + AB.T)/2)

        # P in A+B eigenbasis
        Pt = V.T @ P @ V   # P_tilde

        D = len(idx)
        # Diagonal of tau in A+B basis (a_j - P_tilde[j,j])
        diag_gap = ab_ev - np.diag(Pt)

        # Off-diagonal elements: find top 30
        offdiag = []
        for j in range(D):
            for k in range(j+1, D):
                offdiag.append((abs(Pt[j,k]), j, k, Pt[j,k]))
        offdiag.sort(reverse=True)
        top = offdiag[:30]

        print(f"\nTop 20 off-diagonal |P_tilde[j,k]|:")
        print(f"  {'rank':>4}  {'j':>4}  {'k':>4}  {'a_j':>9}  {'a_k':>9}  "
              f"{'P~[j,j]':>9}  {'P~[k,k]':>9}  {'gap_j':>9}  {'gap_k':>9}  "
              f"{'|P~[j,k]|':>10}  {'2nd-ord-cost':>13}")
        results = []
        for rank, (val, j, k, Ptjk) in enumerate(top[:20]):
            aj = ab_ev[j]; ak = ab_ev[k]
            gj = diag_gap[j]; gk = diag_gap[k]
            # 2nd-order perturbation cost: |P~[j,k]|^2 / |a_j - a_k|
            denom = abs(aj - ak)
            cost = val**2 / denom if denom > 1e-10 else float('inf')

            # Decompose by prime power
            pp_contribs = {}
            for (kp, lk), Ppp in P_per_pp.items():
                Pppt_jk = float(V[:,j] @ Ppp @ V[:,k])
                pp_contribs[(kp, lk)] = Pppt_jk
            # top prime contributors to this off-diagonal element
            top_pp = sorted(pp_contribs.items(), key=lambda x: abs(x[1]), reverse=True)[:5]

            print(f"  {rank+1:>4}  {j:>4}  {k:>4}  {aj:>9.4f}  {ak:>9.4f}  "
                  f"{Pt[j,j]:>9.4f}  {Pt[k,k]:>9.4f}  {gj:>9.4f}  {gk:>9.4f}  "
                  f"{val:>10.5f}  {cost:>13.5f}")
            pp_str = ", ".join(f"k={kp}(p={round(np.exp(lk))},c={c:.4f})"
                                for (kp,lk),c in top_pp[:3])
            print(f"       top primes: {pp_str}")
            print(f"       j-modes: {idx[j]}  k-modes: {idx[k]}")

            results.append({
                'rank': rank+1, 'j': j, 'k': k,
                'a_j': float(aj), 'a_k': float(ak),
                'diag_gap_j': float(gj), 'diag_gap_k': float(gk),
                'Pt_jk': float(Ptjk), 'abs_Pt_jk': float(val),
                'second_order_cost': float(cost),
                'top_prime_contribs': [
                    {'k': kp, 'log_p': lk, 'contribution': c}
                    for (kp,lk),c in top_pp
                ],
                'mode_j': idx[j], 'mode_k': idx[k]
            })

        # Summary: min eigenvalue of tau vs predicted 2nd-order correction
        ev_tau = np.linalg.eigvalsh(tau)
        print(f"\n  min eig(tau) = {ev_tau[0]:.4e}")
        print(f"  min diag gap = {diag_gap.min():.4f} (at j={diag_gap.argmin()})")
        print(f"  Total 2nd-order correction = {sum(r['second_order_cost'] for r in results):.4f}")

        # Are near-null A+B eigenvectors the ones with smallest diagonal gaps?
        # Show the 10 smallest diagonal gaps and their corresponding a_j
        gap_order = np.argsort(diag_gap)
        print(f"\n  10 smallest diagonal gaps (A+B - P_diag):")
        print(f"    {'j':>4}  {'a_j':>9}  {'P~[j,j]':>9}  {'diag_gap':>10}  {'mode':>5}")
        for j in gap_order[:10]:
            print(f"    {j:>4}  {ab_ev[j]:>9.4f}  {Pt[j,j]:>9.4f}  {diag_gap[j]:>10.5f}  {idx[j]:>5}")

        with open(OUTFILE.replace('.json', f'_N{N}.json'), 'w') as f:
            json.dump({'N': N, 'top_offdiag': results,
                       'min_eig_tau': float(ev_tau[0]),
                       'min_diag_gap': float(diag_gap.min())}, f, indent=2)

    print(f"\nDone.")


if __name__ == "__main__":
    main()
