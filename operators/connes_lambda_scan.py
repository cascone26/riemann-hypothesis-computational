"""
Connes Lambda Scan: How does ε_N depend on λ²?

With λ²=14 (only primes 2,3,5,7,11,13), ε_N plateaus at ~0.32.
Question: does larger λ² (more primes, different L) drive ε_N lower?

Also tests: diagonal-only vs full off-diagonal archimedean,
to identify where the systematic error comes from.

λ² values: 14, 30, 50, 100, 200
Fixed N = 80 (fast, captures the pattern)
"""

import mpmath
import numpy as np
import json
import os
import time

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

DATA_DIR    = os.path.join(os.path.dirname(os.path.dirname(__file__)), "data")
RESULTS_DIR = os.path.join(os.path.dirname(os.path.dirname(__file__)), "results")
PLOTS_DIR   = os.path.join(os.path.dirname(os.path.dirname(__file__)), "plots")

mpmath.mp.dps = 60
N_TRUNC = 80


def primes_up_to(n):
    sieve = list(range(n+1))
    for i in range(2, int(n**0.5)+1):
        if sieve[i]:
            for j in range(i*i, n+1, i):
                sieve[j] = 0
    return [x for x in sieve[2:] if x]


def von_mangoldt(k, primes):
    if k <= 1: return 0.0
    for p in primes:
        if p > k: break
        pk = p
        while pk <= k:
            if pk == k: return float(mpmath.log(p))
            pk *= p
    return 0.0


def build_matrix(lambda_sq, N, mode='full'):
    """
    mode: 'full'     = diagonal + off-diagonal archimedean
          'diagonal' = archimedean diagonal only (off-diagonal = 0)
          'no_arch'  = boundary - prime only (no archimedean)
    """
    L      = float(2 * mpmath.log(mpmath.sqrt(lambda_sq)))
    primes = primes_up_to(lambda_sq)
    dim    = 2*N + 1
    pi     = float(mpmath.pi)
    tau    = np.zeros((dim, dim))

    k_max  = int(float(mpmath.exp(mpmath.mpf(str(L))))) + 1

    for i in range(dim):
        n = i - N
        for j in range(i, dim):
            m = j - N
            val = 0.0

            # Archimedean term
            if mode in ('full', 'diagonal'):
                if n == m:
                    if n == 0:
                        arch = float(-mpmath.euler - mpmath.log(2))
                    else:
                        s    = pi * n / L
                        arch = float(mpmath.re(mpmath.digamma(mpmath.mpf('0.25') + 1j*s))) + float(mpmath.log(2))
                    val += arch
                elif mode == 'full' and n != m:
                    k_nm  = n - m
                    s_avg = pi * (abs(n) + abs(m)) / (2*L)
                    if s_avg < 0.01: s_avg = 0.01
                    psi_im = float(mpmath.im(mpmath.digamma(mpmath.mpf('0.25') + 1j*s_avg)))
                    val += ((-1)**abs(k_nm)) * psi_im / (pi * abs(k_nm))

            # Boundary term
            L_mp = mpmath.mpf(str(L))
            num  = 32*L_mp * mpmath.sinh(L_mp/4)**2 * (L_mp**2 - 16*pi**2*m*n)
            den  = (L_mp**2 + 16*pi**2*m**2) * (L_mp**2 + 16*pi**2*n**2)
            if abs(float(den)) > 1e-50:
                val += float(num / den)

            # Prime term (subtract)
            for k in range(2, k_max + 1):
                lam_k = von_mangoldt(k, primes)
                if lam_k == 0: continue
                y      = float(mpmath.log(k))
                weight = lam_k * k**(-0.5)
                if n != m:
                    q = (np.sin(2*pi*m*y/L) - np.sin(2*pi*n*y/L)) / (pi*(n-m))
                else:
                    q = 2*(1 - abs(y)/L) * np.cos(2*pi*n*y/L)
                val -= weight * q

            tau[i, j] = val
            tau[j, i] = val

    evals = np.linalg.eigvalsh(tau)
    return float(evals[0]), evals[:5].tolist(), L, len(primes)


def main():
    print("CONNES LAMBDA SCAN: ε_N vs λ²")
    print("=" * 60)
    print(f"Fixed N={N_TRUNC}")
    print()

    lambda_sqs = [14, 30, 50, 100, 200]
    results    = []

    print(f"{'λ²':>5}  {'L':>6}  {'#primes':>8}  {'ε_N (full)':>12}  {'ε_N (diag)':>12}  {'time':>6}")
    print("-" * 60)

    for lsq in lambda_sqs:
        t0 = time.time()
        eps_full, evals5_full, L, np_count = build_matrix(lsq, N_TRUNC, mode='full')
        t1 = time.time()
        eps_diag, evals5_diag, _, _        = build_matrix(lsq, N_TRUNC, mode='diagonal')
        t2 = time.time()

        print(f"{lsq:>5}  {L:>6.3f}  {np_count:>8d}  "
              f"{eps_full:>+12.6f}  {eps_diag:>+12.6f}  {t2-t0:>5.1f}s")

        results.append({
            "lambda_sq": lsq, "L": L, "n_primes": np_count,
            "epsilon_full": eps_full, "epsilon_diagonal": eps_diag,
            "top5_evals_full": evals5_full,
        })

    print()
    print("Analysis:")
    print("  'full' uses diagonal + off-diagonal archimedean approximation")
    print("  'diag' uses diagonal archimedean only (off-diagonal = 0)")
    print()

    # Find which λ² gives minimum ε_N
    best_full = min(results, key=lambda r: r["epsilon_full"])
    best_diag = min(results, key=lambda r: r["epsilon_diagonal"])
    print(f"Best ε_N (full):     λ²={best_full['lambda_sq']}, ε_N={best_full['epsilon_full']:.6f}")
    print(f"Best ε_N (diagonal): λ²={best_diag['lambda_sq']}, ε_N={best_diag['epsilon_diagonal']:.6f}")

    # ── Plot ────────────────────────────────────────────────────
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    lsqs      = [r["lambda_sq"]      for r in results]
    eps_full  = [r["epsilon_full"]   for r in results]
    eps_diag  = [r["epsilon_diagonal"] for r in results]
    Ls        = [r["L"]              for r in results]
    n_primes  = [r["n_primes"]       for r in results]

    ax = axes[0]
    ax.plot(lsqs, eps_full, 'b^-', markersize=8, linewidth=2, label='Full (diag + off-diag)')
    ax.plot(lsqs, eps_diag, 'rs-', markersize=8, linewidth=2, label='Diagonal only')
    ax.axhline(0, color='k', linestyle='--', linewidth=1, label='ε = 0 (RH target)')
    ax.set_xlabel('λ²')
    ax.set_ylabel('ε_N (minimum eigenvalue)')
    ax.set_title(f'ε_N vs λ² (N={N_TRUNC}): Does more primes help?')
    ax.legend()
    ax.grid(True, alpha=0.3)
    for i, (lsq, ef, ed) in enumerate(zip(lsqs, eps_full, eps_diag)):
        ax.annotate(f'{ef:.3f}', (lsq, ef), textcoords='offset points', xytext=(4,4), fontsize=7, color='blue')
        ax.annotate(f'{ed:.3f}', (lsq, ed), textcoords='offset points', xytext=(4,-10), fontsize=7, color='red')

    ax = axes[1]
    ax.plot(n_primes, eps_full, 'b^-', markersize=8, linewidth=2, label='Full')
    ax.plot(n_primes, eps_diag, 'rs-', markersize=8, linewidth=2, label='Diagonal only')
    ax.axhline(0, color='k', linestyle='--', linewidth=1)
    ax.set_xlabel('Number of primes included')
    ax.set_ylabel('ε_N')
    ax.set_title('ε_N vs #primes: Does more primes improve?')
    ax.legend()
    ax.grid(True, alpha=0.3)

    plt.suptitle(f'Connes Lambda Scan (N={N_TRUNC}) — effect of λ² on ε_N', fontsize=11)
    plt.tight_layout()
    out = os.path.join(PLOTS_DIR, "connes_lambda_scan.png")
    plt.savefig(out, dpi=150)
    plt.close()
    print(f"\nPlot saved: {out}")

    with open(os.path.join(RESULTS_DIR, "connes_lambda_scan.json"), "w") as f:
        json.dump({"N": N_TRUNC, "results": results}, f, indent=2)
    print("Results saved.")


if __name__ == "__main__":
    main()
