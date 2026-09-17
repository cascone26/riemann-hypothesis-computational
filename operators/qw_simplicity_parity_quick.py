#!/usr/bin/env python3
"""
qw_simplicity_parity_quick.py  --  Fast version of the Connes §1174 probe.

Same as qw_simplicity_parity_probe.py but with reduced parameters for quick feedback:
  - Fewer λ values: [4, 10, 20]
  - Smaller N range: [10, 15, 20, 25]
  - Faster turnaround
"""

import sys, json, time
import mpmath as mp
import numpy as np

mp.mp.dps = 40

LAMBDA_SQS = [4, 10, 20]      # Smaller subset of λ² values
N_RANGE = [10, 15, 20, 25]    # Smaller N range

# ---------------------------------------------------------------------------
# von Mangoldt data
# ---------------------------------------------------------------------------
def vonmangoldt_table(kmax):
    tab = []
    for k in range(2, int(kmax) + 1):
        n = k
        p = 0
        f = 2
        ok = True
        while f * f <= n:
            if n % f == 0:
                if p == 0:
                    p = f
                elif p != f:
                    ok = False
                    break
                while n % f == 0:
                    n //= f
            f += 1
        if n > 1:
            if p == 0:
                p = n
            elif p != n:
                ok = False
        if ok and p > 0:
            tab.append((k, mp.log(p)))
    return tab

# ---------------------------------------------------------------------------
# Archimedean density
# ---------------------------------------------------------------------------
def rho(y):
    return mp.e ** (y / 2) / (mp.e ** y - mp.e ** (-y))

def rho_reg(y):
    if y == 0:
        return mp.mpf(1) / 4
    return rho(y) - 1 / (2 * y)

# ---------------------------------------------------------------------------
# Build tau
# ---------------------------------------------------------------------------
def build_tau_fast(lambda_sq, N):
    L = mp.log(lambda_sq)
    TWO_PI = 2 * mp.pi
    VM = vonmangoldt_table(lambda_sq)

    def _w(k):
        return TWO_PI * k / L

    def Ssin(k):
        if k == 0:
            return mp.mpf(0)
        w = _w(k)
        f = lambda y: mp.sin(w * y) * rho(y)
        return mp.quad(f, [0, L / 2, L])

    def Scos_diag(k):
        w = _w(k)
        g = lambda y: (1 - y / L) * mp.cos(w * y)
        reg_part = mp.quad(lambda y: g(y) * rho_reg(y), [0, L / 2, L])
        head_smooth2 = mp.quad(lambda y: (g(y) - 1) / (2 * y) if y != 0 else
                               ( -1/L ), [0, L/2, L])
        C_R = (mp.log(4 * mp.pi) + mp.euler) / 2
        return reg_part + head_smooth2 + C_R

    def Psin(k):
        s = mp.mpf(0)
        w = _w(k)
        for kk, Lam in VM:
            y = mp.log(kk)
            s += Lam * kk ** (mp.mpf(-1) / 2) * mp.sin(w * y)
        return s

    def Pcos_diag(k):
        s = mp.mpf(0)
        w = _w(k)
        for kk, Lam in VM:
            y = mp.log(kk)
            s += Lam * kk ** (mp.mpf(-1) / 2) * 2 * (1 - y / L) * mp.cos(w * y)
        return s

    def W02(n, m):
        num = (L ** 2 - 16 * mp.pi ** 2 * m * n)
        den = (L ** 2 + 16 * mp.pi ** 2 * m ** 2) * (L ** 2 + 16 * mp.pi ** 2 * n ** 2)
        return 32 * L * mp.sinh(L / 4) ** 2 * num / den

    idx = list(range(-N, N + 1))
    D = 2 * N + 1

    b = {}
    a = {}
    for k in range(0, 2 * N + 1):
        bk = -(Ssin(k) - Psin(k)) / mp.pi
        b[k] = bk
        b[-k] = -bk
    for k in range(0, N + 1):
        ak = (Scos_diag(k) - Pcos_diag(k)) + W02(k, k)
        a[k] = ak
        a[-k] = ak

    tau = mp.zeros(D, D)
    for ii, n in enumerate(idx):
        for jj, m in enumerate(idx):
            if n == m:
                tau[ii, jj] = a[n]
            else:
                tau[ii, jj] = (b[n] - b[m]) / (n - m) + W02(n, m)

    for ii in range(D):
        for jj in range(ii + 1, D):
            v = (tau[ii, jj] + tau[jj, ii]) / 2
            tau[ii, jj] = v
            tau[jj, ii] = v

    return tau, idx

# ---------------------------------------------------------------------------
# Analyze
# ---------------------------------------------------------------------------
def analyze_smallest_eigenvalue(lambda_sq, N):
    t0 = time.time()
    tau, idx = build_tau_fast(lambda_sq, N)
    D = 2 * N + 1

    A = np.array([[float(tau[i, j]) for j in range(D)] for i in range(D)],
                 dtype=float)
    A = (A + A.T) / 2

    evals, evecs = np.linalg.eigh(A)

    ε1 = float(evals[0])
    ε2 = float(evals[1]) if len(evals) > 1 else ε1 + 1.0
    gap = ε2 - ε1

    v = evecs[:, 0]
    even = np.array([(v[i] + v[D - 1 - i]) / 2 for i in range(D)])
    parity_loss = np.linalg.norm(v - even) / np.linalg.norm(v) if np.linalg.norm(v) > 0 else 1.0

    t_elapsed = time.time() - t0

    return {
        "lambda_sq": lambda_sq,
        "N": N,
        "D": D,
        "epsilon_1": ε1,
        "epsilon_2": ε2,
        "gap": gap,
        "gap_relative": gap / abs(ε1) if ε1 != 0 else 1.0,
        "parity_loss_norm": parity_loss,
        "time_seconds": t_elapsed
    }

# ---------------------------------------------------------------------------
def run_quick_sweep():
    results = {
        "metadata": {
            "title": "Connes 2026 §1174: Simplicity & parity (QUICK)",
            "target": "Is ε(λ) simple with even eigenvector?",
        },
        "runs": []
    }

    print("\n" + "="*80)
    print("QW_λ SIMPLICITY & PARITY (QUICK PROBE)")
    print("="*80)
    print(f"{'λ²':>6} {'N':>3} {'D':>4} {'ε₁':>14} {'ε₂':>14} {'gap':>14} {'gap/|ε₁|':>10} {'parity_loss':>12} {'time':>8}")
    print("-"*80)

    for lambda_sq in LAMBDA_SQS:
        for N in N_RANGE:
            try:
                result = analyze_smallest_eigenvalue(lambda_sq, N)
                results["runs"].append(result)

                lsq = f"{lambda_sq}"
                n_str = f"{N}"
                d_str = f"{result['D']}"
                eps1_str = f"{result['epsilon_1']:.6e}"
                eps2_str = f"{result['epsilon_2']:.6e}"
                gap_str = f"{result['gap']:.6e}"
                gap_rel = f"{result['gap_relative']:.4e}"
                parity = f"{result['parity_loss_norm']:.4e}"
                time_str = f"{result['time_seconds']:.1f}s"

                print(f"{lsq:>6} {n_str:>3} {d_str:>4} {eps1_str:>14} {eps2_str:>14} "
                      f"{gap_str:>14} {gap_rel:>10} {parity:>12} {time_str:>8}")
            except Exception as e:
                print(f"{lambda_sq:>6} {N:>3}  ERROR: {str(e)[:60]}")

    print("\n" + "="*80)
    print("REPORT-TIER SUMMARY:")
    print("-"*80)

    gaps = [r["gap"] for r in results["runs"]]
    parities = [r["parity_loss_norm"] for r in results["runs"]]

    print(f"Simplicity (gap ε₂ - ε₁):")
    print(f"  Min: {min(gaps):.6e},  Max: {max(gaps):.6e},  Mean: {np.mean(gaps):.6e}")
    print(f"  All > 0? {all(g > 0 for g in gaps)} ✓" if all(g > 0 for g in gaps) else f"  Some ≤ 0 ✗")

    print(f"\nParity (eigenvector loss after even projection):")
    print(f"  Min: {min(parities):.6e},  Max: {max(parities):.6e},  Mean: {np.mean(parities):.6e}")

    print("\n" + "="*80 + "\n")

    return results

# ---------------------------------------------------------------------------
if __name__ == "__main__":
    results = run_quick_sweep()
    op = "/Users/scones/projects/riemann/results/qw_simplicity_quick.json"
    with open(op, "w") as fh:
        json.dump(results, fh, indent=2)
    print(f"Results saved: {op}\n")
