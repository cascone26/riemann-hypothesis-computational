"""
Connes ε_N Convergence Test

The key claim of Connes' spectral approach:
  ε_N → 0 as N → ∞   ⟺   RH is true

where ε_N = minimum eigenvalue of the Weil quadratic form τ(V_n, V_m).

We test this numerically: run N = 20, 40, 60, 80, 100, 120, 150, 200
and plot ε_N vs N. If ε_N is monotone decreasing toward 0, that's
numerical evidence consistent with RH.

Previous data points:
  N=30:  ε_N = 1.144  (v1, diagonal archimedean only)
  N=120: ε_N = 0.326  (v2, off-diagonal archimedean included)
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
LAMBDA_SQ = 14


def primes_up_to(n):
    sieve = list(range(n+1))
    for i in range(2, int(n**0.5)+1):
        if sieve[i]:
            for j in range(i*i, n+1, i):
                sieve[j] = 0
    return [x for x in sieve[2:] if x]


def von_mangoldt(k, primes):
    if k <= 1:
        return 0.0
    for p in primes:
        if p > k:
            break
        pk = p
        while pk <= k:
            if pk == k:
                return float(mpmath.log(p))
            pk *= p
    return 0.0


def W_archimedean(n, m, L):
    pi = float(mpmath.pi)
    if n == m:
        if n == 0:
            return float(-mpmath.euler - mpmath.log(2))
        s = pi * n / L
        val = float(mpmath.re(mpmath.digamma(mpmath.mpf('0.25') + 1j*s)))
        return val + float(mpmath.log(2))
    else:
        k = n - m
        s_avg = pi * (abs(n) + abs(m)) / (2*L)
        if s_avg < 0.01:
            s_avg = 0.01
        psi_im = float(mpmath.im(mpmath.digamma(mpmath.mpf('0.25') + 1j*s_avg)))
        return ((-1)**abs(k)) * psi_im / (pi * abs(k))


def W_boundary(n, m, L):
    pi = float(mpmath.pi)
    L_mp = mpmath.mpf(str(L))
    num = 32*L_mp * mpmath.sinh(L_mp/4)**2 * (L_mp**2 - 16*pi**2*m*n)
    den = (L_mp**2 + 16*pi**2*m**2) * (L_mp**2 + 16*pi**2*n**2)
    if abs(float(den)) < 1e-50:
        return 0.0
    return float(num / den)


def W_prime(n, m, L, primes):
    pi = float(mpmath.pi)
    k_max = int(float(mpmath.exp(mpmath.mpf(str(L))))) + 1
    total = 0.0
    for k in range(2, k_max + 1):
        lam_k = von_mangoldt(k, primes)
        if lam_k == 0:
            continue
        y = float(mpmath.log(k))
        weight = lam_k * k**(-0.5)
        if n != m:
            q = (np.sin(2*pi*m*y/L) - np.sin(2*pi*n*y/L)) / (pi * (n - m))
        else:
            q = 2*(1 - abs(y)/L) * np.cos(2*pi*n*y/L)
        total += weight * q
    return total


def compute_epsilon(N, verbose=True):
    """Build Weil matrix of size (2N+1)×(2N+1) and return minimum eigenvalue."""
    L      = float(2 * mpmath.log(mpmath.sqrt(LAMBDA_SQ)))
    primes = primes_up_to(LAMBDA_SQ)
    dim    = 2*N + 1
    tau    = np.zeros((dim, dim))
    t0     = time.time()

    for i in range(dim):
        n = i - N
        for j in range(i, dim):
            m = j - N
            val = (W_archimedean(n, m, L)
                   + W_boundary(n, m, L)
                   - W_prime(n, m, L, primes))
            tau[i, j] = val
            tau[j, i] = val

    build_time = time.time() - t0
    evals = np.linalg.eigvalsh(tau)
    eps   = float(evals[0])

    if verbose:
        print(f"  N={N:4d}: dim={dim:4d}, ε_N={eps:+.6f}  [{build_time:.1f}s]")

    return eps, build_time


def main():
    print("CONNES ε_N CONVERGENCE TEST")
    print("=" * 55)
    print(f"λ²={LAMBDA_SQ}, testing N = 20 to 200")
    print("Theory: ε_N → 0 as N → ∞ iff RH")
    print()

    # Previous results for reference
    known = {30: 1.144, 120: 0.326}

    N_values = [20, 40, 60, 80, 100, 120, 150, 180]
    results  = []

    for N in N_values:
        eps, t = compute_epsilon(N, verbose=True)
        results.append({"N": N, "dim": 2*N+1, "epsilon_N": eps, "build_time": t})

    # ── Analysis ────────────────────────────────────────────────
    print()
    print("=" * 55)
    print("CONVERGENCE ANALYSIS")
    print("=" * 55)
    Ns   = [r["N"]         for r in results]
    epss = [r["epsilon_N"] for r in results]

    print(f"{'N':>6}  {'ε_N':>12}  {'Δε':>10}")
    prev = None
    for r in results:
        delta = f"{r['epsilon_N']-prev:+.6f}" if prev is not None else "     ---"
        print(f"{r['N']:>6}  {r['epsilon_N']:>+12.6f}  {delta:>10}")
        prev = r["epsilon_N"]

    print()
    if epss[-1] < epss[0]:
        trend = "DECREASING"
    elif epss[-1] > epss[0]:
        trend = "INCREASING"
    else:
        trend = "FLAT"

    # Fit power law: ε_N ~ A · N^β
    log_Ns   = np.log(Ns)
    log_epss = np.log(np.abs(epss))
    coeffs   = np.polyfit(log_Ns, log_epss, 1)
    beta     = coeffs[0]
    A        = np.exp(coeffs[1])
    print(f"Power-law fit: ε_N ≈ {A:.4f} · N^{beta:.3f}")
    print(f"Trend: {trend}")

    # Extrapolate
    for N_ext in [500, 1000, 10000]:
        eps_ext = A * N_ext**beta
        print(f"  ε_N extrapolated at N={N_ext}: {eps_ext:.4f}")

    # ── Plot ────────────────────────────────────────────────────
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    ax = axes[0]
    ax.plot(Ns, epss, 'b^-', markersize=8, linewidth=2, label='Computed ε_N')
    # Add known reference points
    for nk, ek in known.items():
        if nk not in Ns:
            ax.scatter([nk], [ek], c='gray', s=80, marker='o', zorder=5,
                       label=f'Reference (N={nk})')
    ax.axhline(0, color='red', linestyle='--', linewidth=1.5, label='ε_N = 0 (RH target)')
    ax.set_xlabel('N (truncation size)')
    ax.set_ylabel('ε_N (minimum eigenvalue)')
    ax.set_title(f'Connes ε_N Convergence (λ²={LAMBDA_SQ})')
    ax.legend()
    ax.grid(True, alpha=0.3)

    ax = axes[1]
    ax.loglog(Ns, np.abs(epss), 'b^-', markersize=8, linewidth=2, label='|ε_N|')
    # Power-law fit line
    N_fit = np.array([Ns[0], Ns[-1]*2])
    ax.loglog(N_fit, A * N_fit**beta, 'r--', linewidth=1.5,
              label=f'Fit: {A:.3f}·N^{beta:.3f}')
    ax.axhline(0.01, color='green', linestyle=':', linewidth=1, label='0.01 threshold')
    ax.set_xlabel('N (log scale)')
    ax.set_ylabel('|ε_N| (log scale)')
    ax.set_title('Log-Log: Power Law Convergence?')
    ax.legend()
    ax.grid(True, alpha=0.3)

    plt.suptitle(
        f'Connes ε_N → 0? λ²={LAMBDA_SQ} | '
        f'Trend: {trend} | Power-law: N^{beta:.2f}',
        fontsize=11
    )
    plt.tight_layout()
    out = os.path.join(PLOTS_DIR, "connes_epsilon_convergence.png")
    plt.savefig(out, dpi=150)
    plt.close()
    print(f"\nPlot saved: {out}")

    with open(os.path.join(RESULTS_DIR, "connes_epsilon_convergence.json"), "w") as f:
        json.dump({"lambda_sq": LAMBDA_SQ, "results": results,
                   "trend": trend, "power_law_beta": float(beta),
                   "power_law_A": float(A)}, f, indent=2)
    print("Results saved.")


if __name__ == "__main__":
    main()
