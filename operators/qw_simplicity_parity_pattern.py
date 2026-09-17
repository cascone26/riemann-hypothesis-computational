#!/usr/bin/env python3
"""
qw_simplicity_parity_pattern.py  --  Map the pattern of eigenvector parity.

FINDING: The smallest eigenvalue changes parity (even ↔ odd) as N varies.
This is a crucial observation for Connes §1174.

Goal: Map which (λ, N) pairs have EVEN vs ODD smallest eigenvector.
"""

import sys, json
import mpmath as mp
import numpy as np

mp.mp.dps = 40

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

def rho(y):
    return mp.e ** (y / 2) / (mp.e ** y - mp.e ** (-y))

def rho_reg(y):
    if y == 0:
        return mp.mpf(1) / 4
    return rho(y) - 1 / (2 * y)

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

def check_parity(lambda_sq, N):
    """
    Returns:
      - parity: "EVEN", "ODD", or "MIXED"
      - even_norm: ||even part||/||v||
      - odd_norm: ||odd part||/||v||
    """
    tau, idx = build_tau_fast(lambda_sq, N)
    D = 2 * N + 1

    A = np.array([[float(tau[i, j]) for j in range(D)] for i in range(D)],
                 dtype=float)
    A = (A + A.T) / 2

    evals, evecs = np.linalg.eigh(A)
    v = evecs[:, 0]

    even = np.array([(v[i] + v[D - 1 - i]) / 2 for i in range(D)])
    odd = np.array([(v[i] - v[D - 1 - i]) / 2 for i in range(D)])

    norm_v = np.linalg.norm(v)
    norm_even = np.linalg.norm(even) / norm_v if norm_v > 0 else 0
    norm_odd = np.linalg.norm(odd) / norm_v if norm_v > 0 else 0

    threshold = 0.1
    if norm_even > 1 - threshold and norm_odd < threshold:
        parity = "EVEN"
    elif norm_odd > 1 - threshold and norm_even < threshold:
        parity = "ODD"
    else:
        parity = "MIXED"

    return parity, norm_even, norm_odd

def run_pattern_scan():
    print("\n" + "="*80)
    print("PARITY PATTERN SCAN: Which (λ, N) have EVEN vs ODD smallest eigenvector?")
    print("="*80)

    lambda_sqs = [4, 6, 8, 10, 14, 20]
    n_range = list(range(8, 36, 2))  # N = 8, 10, 12, ..., 34

    results = {}

    for lambda_sq in lambda_sqs:
        print(f"\nλ² = {lambda_sq}:")
        print(f"  N:     {' '.join(f'{n:3d}' for n in n_range)}")
        print(f"  Par:   ", end="")

        results[lambda_sq] = {}
        parity_seq = []

        for N in n_range:
            try:
                parity, even_norm, odd_norm = check_parity(lambda_sq, N)
                results[lambda_sq][N] = {
                    "parity": parity,
                    "even_norm": even_norm,
                    "odd_norm": odd_norm
                }
                # Symbol: E = EVEN, O = ODD, M = MIXED
                symbol = "E" if parity == "EVEN" else ("O" if parity == "ODD" else "M")
                parity_seq.append(symbol)
                print(f"{symbol:3s} ", end="")
            except Exception as e:
                print(f"!   ", end="")
                parity_seq.append("!")

        print()

        # Analyze pattern
        pattern_str = "".join(parity_seq)
        print(f"  Pattern: {pattern_str}")

        # Count alternations
        alternations = sum(1 for i in range(len(pattern_str)-1) if pattern_str[i] != pattern_str[i+1])
        print(f"  Alternations: {alternations} (strict alternation would be {len(pattern_str)-1})")

    print("\n" + "="*80)
    print("REPORT-TIER INTERPRETATION:")
    print("-"*80)
    print("• EVEN eigenvector (E): smallest eigenvalue is symmetric v[i] = v[mirror]")
    print("• ODD eigenvector (O): smallest eigenvalue is antisymmetric v[i] = -v[mirror]")
    print("• MIXED: neither EVEN nor ODD dominates")
    print()
    print("KEY FINDING:")
    print("• If the smallest eigenvalue alternates between EVEN and ODD with N,")
    print("  then it is NOT unique (simple) in the usual sense.")
    print("• The smallest eigenvalue eigenspace may be 2-dimensional (degenerate):")
    print("  one even eigenvector and one odd eigenvector, both with the same ε.")
    print()
    print("This would REFUTE the hypothesis that ε(λ) is SIMPLE with unique parity.")
    print("="*80 + "\n")

    return results

if __name__ == "__main__":
    results = run_pattern_scan()

    # Save detailed results
    op = "/Users/scones/projects/riemann/results/qw_parity_pattern.json"
    # Convert results to serializable format
    serializable = {}
    for lsq, data in results.items():
        serializable[str(lsq)] = {str(n): v for n, v in data.items()}

    with open(op, "w") as fh:
        json.dump(serializable, fh, indent=2)
    print(f"Detailed results saved: {op}\n")
