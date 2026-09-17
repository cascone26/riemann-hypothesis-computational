#!/usr/bin/env python3
"""Honest analysis: Can Alpöge-Furman's rank-trace technique transfer to κ(M_σ) ≤ 1?

KEY FINDING from first attempt:
- Naive rank-trace bounds give κ ≤ 2, but actual κ = 11-15 at 4W=5,6
- This means either: (1) the technique needs deeper/tighter application, or (2) it doesn't
  transfer to our specific problem without additional structure

ANALYSIS:
The A-F technique has a specific PURPOSE: count simple/multiple/off-line zeros via matrix
decomposition. The MECHANISM: they decompose the Weil matrix G̃ into blocks by zero type,
where each zero's contribution has known inertia (simple-on-line→PSD, off-line→sig-(1,1)).

To transfer to our κ problem, we would need to:
(A) Decompose M_σ similarly by zero type
(B) Use Sylvester's law to bound inertia of each piece
(C) Sum to bound total κ

But step (A) requires KNOWING zero types (simple? on-line? off-line?), which is what
A-F compute FROM the matrix, using our very hard problem.

HONEST VERDICT:
The A-F rank-trace technique is designed to EXTRACT zero information from the matrix,
not to bound κ(M_σ) without that information. A direct transfer would be circular.

HOWEVER: There are two potential non-circular routes:
(1) Use A-F's own OUTPUTS (≥67% simple, ≥83% distinct) as inputs to bound off-line pairs,
    then use Groskin's result (κ = #off-line) to bound κ. But this gives κ ≤ ~0.33·N,
    not κ ≤ 1.
(2) Look for a STRUCTURAL reason why κ ≤ 1 that doesn't depend on knowing zeros
    individually. This is the Krein-Langer/Pontryagin route already identified.

Let me verify the κ=1 claim more carefully: maybe the issue is that 4W=5,6 are NOT the
bottom edge where κ=1 lives.
"""
import numpy as np
from scipy.special import digamma
from scipy.linalg import eigh
import matplotlib.pyplot as plt

def vonmangoldt(N):
    lam = np.zeros(N+1)
    isc = np.zeros(N+1, dtype=bool)
    for p in range(2, N+1):
        if not isc[p]:
            isc[2*p::p] = True
            lp = np.log(p)
            pk = p
            while pk <= N:
                lam[pk] = lp
                pk *= p
    return lam

def compute_kappa_vs_height(fourW=5.0, cen_lo_list=None, cen_hi_list=None,
                            frac=0.8, dr=0.004):
    """
    Compute κ(M_σ) for different window HEIGHTS to find where κ=1 occurs.

    The RH_CORE_LOCALIZATION.md claim is that κ=1 at the BOTTOM EDGE of the
    spectrum. We need to find which (cen_lo, cen_hi) range that is.
    """
    if cen_lo_list is None:
        cen_lo_list = [10, 12, 14, 16, 20, 30, 40, 50, 60]
    if cen_hi_list is None:
        cen_hi_list = [l + 30 for l in cen_lo_list]

    NMAX = 8000
    LAM = vonmangoldt(NMAX)
    NS = np.arange(2, NMAX+1)
    LOGN = np.log(NS)
    LAMv = LAM[2:]

    results = []

    for cen_lo, cen_hi in zip(cen_lo_list, cen_hi_list):
        W = fourW / 4
        L = cen_hi - cen_lo
        N = max(4, round(frac * L * 2*W / np.pi))

        RG = np.arange(-(cen_hi+120), (cen_hi+120), dr)
        OM = np.real(digamma((0.5+0.5j*RG)/2)) + np.log(1 / np.pi)
        CEN = np.linspace(cen_lo, cen_hi, N)

        def fej(x):
            y = W*x
            o = np.ones_like(y)
            nz = np.abs(y) > 1e-12
            o[nz] = (np.sin(y[nz]) / y[nz])**2
            return o

        P = np.array([fej(RG - c) + fej(RG + c) for c in CEN])
        G = (P @ P.T) * dr

        arch = (P * OM[None, :]) @ P.T * dr / (2*np.pi)
        mask = LOGN < fourW
        ln = LOGN[mask]
        chiv = np.ones(np.sum(mask))
        wts = LAMv[mask] * chiv / np.sqrt(NS[mask])

        cosu = np.cos(np.outer(ln, RG))
        prime = np.zeros((N, N))
        for j in range(N):
            gp = ((P[j][None, :] * P) @ cosu.T) * dr / (2*np.pi)
            prime[j, :] = 2.0 * (gp * wts[None, :]).sum(axis=1)

        M = arch - 0.5 * (prime + prime.T)

        gw, gv = eigh(G)
        k = gw > gw.max() * 1e-12
        B = gv[:, k] / np.sqrt(gw[k])

        ev = np.linalg.eigvalsh(B.T @ M @ B)
        kappa = int((ev < -1e-9*max(1, abs(ev).max())).sum())

        # Also compute the full form Q = arch + pole - prime (pole is rank ≤ 2)
        # For now, ignore pole and look at M_sigma alone
        Q_compact = B.T @ (arch - 0.5*(prime+prime.T)) @ B
        kappa_Q = int((np.linalg.eigvalsh(Q_compact) < -1e-9).sum())

        margin = min(ev)  # Smallest eigenvalue (most negative)
        n_pos_Q = np.sum(np.linalg.eigvalsh(Q_compact) > 1e-10)
        n_neg_Q = int((np.linalg.eigvalsh(Q_compact) < -1e-10).sum())

        height_label = f"T ∈ [{cen_lo}, {cen_hi}]"
        results.append({
            "cen_lo": cen_lo,
            "cen_hi": cen_hi,
            "height_label": height_label,
            "kappa": kappa,
            "kappa_Q": kappa_Q,
            "margin": margin,
            "N": N,
            "n_pos": n_pos_Q,
            "n_neg": n_neg_Q,
            "eigvals_bottom_5": sorted(np.linalg.eigvalsh(Q_compact))[:5],
        })

    return results

def report_kappa_landscape(fourW=5.0):
    """Survey κ(M_σ) across different heights to find κ=1 region."""
    print("\n" + "="*80)
    print(f"LANDSCAPE SURVEY: κ(M_σ) vs. HEIGHT at 4W={fourW}")
    print("="*80)
    print(f"{'Height':>20} {'κ':>4} {'Margin':>12} {'#pos':>4} {'#neg':>4} {'Bottom 3 eigvals':>50}")
    print("-"*80)

    results = compute_kappa_vs_height(fourW=fourW)

    for r in results:
        eigvals_str = ", ".join([f"{e:.3e}" for e in r["eigvals_bottom_5"][:3]])
        print(f"{r['height_label']:>20} {r['kappa']:>4} {r['margin']:>12.6e} "
              f"{r['n_pos']:>4} {r['n_neg']:>4}   {eigvals_str:>50}")

    return results

if __name__ == "__main__":
    print("\n" + "="*80)
    print("DIAGNOSTIC: WHERE IS κ = 1?")
    print("="*80)

    print("\nChecking 4W=5.0 (large bandwidth — likely interior, not bottom edge):")
    results_5 = report_kappa_landscape(fourW=5.0)

    print("\n\nChecking 4W=3.0 (smaller bandwidth — might be closer to bottom edge):")
    results_3 = report_kappa_landscape(fourW=3.0)

    print("\n\nChecking 4W=2.0 (even smaller — very bottom edge):")
    results_2 = report_kappa_landscape(fourW=2.0)

    print("\n" + "="*80)
    print("SUMMARY: Locating κ = 1")
    print("="*80)
    for fourW_val, results in [(5.0, results_5), (3.0, results_3), (2.0, results_2)]:
        min_kappa = min([r['kappa'] for r in results])
        kappa_1_cases = [r for r in results if r['kappa'] == 1]
        if kappa_1_cases:
            print(f"  4W={fourW_val}: κ=1 found at {[c['height_label'] for c in kappa_1_cases]}")
        else:
            print(f"  4W={fourW_val}: κ=1 NOT found; minimum κ={min_kappa}")
