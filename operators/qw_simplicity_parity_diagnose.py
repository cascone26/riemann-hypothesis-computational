#!/usr/bin/env python3
"""
qw_simplicity_parity_diagnose.py  --  Diagnostic version for problem cases.

Investigate why parity_loss = 1.0 for some (λ, N) pairs.
Check eigenvector structure and even/odd decomposition.
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

def diagnose(lambda_sq, N):
    """Detailed diagnostic of parity for one (λ, N) pair."""
    print(f"\n{'='*70}")
    print(f"DIAGNOSIS: λ² = {lambda_sq}, N = {N}")
    print(f"{'='*70}")

    tau, idx = build_tau_fast(lambda_sq, N)
    D = 2 * N + 1

    A = np.array([[float(tau[i, j]) for j in range(D)] for i in range(D)],
                 dtype=float)
    A = (A + A.T) / 2

    evals, evecs = np.linalg.eigh(A)

    ε1 = float(evals[0])
    ε2 = float(evals[1]) if len(evals) > 1 else ε1 + 1.0

    v = evecs[:, 0]
    print(f"Smallest eigenvalue ε₁ = {ε1:.10e}")
    print(f"Second eigenvalue ε₂ = {ε2:.10e}")
    print(f"Gap = {ε2 - ε1:.10e}")

    # Analyze parity structure
    print(f"\nEigenvector structure (first 5 entries):")
    for i in range(min(5, D)):
        mirror_i = D - 1 - i
        val_i = v[i]
        val_mirror = v[mirror_i]
        diff = abs(val_i - val_mirror)
        ratio = diff / (abs(val_i) + 1e-30)
        print(f"  v[{i:2d}] = {val_i:+.8e},  v[{mirror_i:2d}] = {val_mirror:+.8e},  diff = {diff:.2e},  ratio = {ratio:.2e}")

    # Check even/odd decomposition
    even = np.array([(v[i] + v[D - 1 - i]) / 2 for i in range(D)])
    odd = np.array([(v[i] - v[D - 1 - i]) / 2 for i in range(D)])

    norm_v = np.linalg.norm(v)
    norm_even = np.linalg.norm(even)
    norm_odd = np.linalg.norm(odd)

    print(f"\nParity decomposition:")
    print(f"  ||v||          = {norm_v:.10e}")
    print(f"  ||even part||  = {norm_even:.10e}")
    print(f"  ||odd part||   = {norm_odd:.10e}")
    print(f"  even/||v||     = {norm_even/norm_v:.10e}")
    print(f"  odd/||v||      = {norm_odd/norm_v:.10e}")

    # Check reconstruction
    reconstructed = even + odd
    recon_error = np.linalg.norm(v - reconstructed) / norm_v
    print(f"  reconstruction error (even+odd vs v): {recon_error:.2e}")

    # Loss metric check
    loss_metric = np.linalg.norm(v - even) / norm_v
    print(f"\nParity loss metric ||v - even||/||v|| = {loss_metric:.10e}")

    # Direct evenness check
    is_even_max_diff = max([abs(v[i] - v[D - 1 - i]) for i in range(D)])
    print(f"Max |v[i] - v[mirror]| = {is_even_max_diff:.10e}")

    # If parity loss is 1.0, diagnose why
    if loss_metric > 0.99:
        print(f"\n⚠️  PARITY LOSS NEAR 1.0 — DIAGNOSTIC:")
        if norm_v > 0:
            print(f"   ||even|| vs ||odd||: {norm_even:.2e} vs {norm_odd:.2e}")
            if norm_odd > 0.5 * norm_even:
                print(f"   → Eigenvector has significant ODD component!")
                print(f"   → This suggests the eigenvector is NOT even.")
            elif norm_even < 1e-10:
                print(f"   → Even part is nearly zero!")
                print(f"   → Normalization issue?")
        else:
            print(f"   → ||v|| = 0; normalization failed?")

    return {
        "lambda_sq": lambda_sq,
        "N": N,
        "epsilon_1": ε1,
        "gap": ε2 - ε1,
        "parity_loss": loss_metric,
        "norm_even_ratio": norm_even / norm_v,
        "norm_odd_ratio": norm_odd / norm_v,
        "is_even_max_diff": is_even_max_diff,
    }

if __name__ == "__main__":
    # Diagnose the problem cases: parity = 1.0
    problem_cases = [
        (4, 15),   # parity = 1.0
        (10, 25),  # parity = 1.0
    ]

    # Also diagnose some good cases for comparison
    good_cases = [
        (4, 10),   # parity = 2.2e-15
        (10, 20),  # parity = 3.6e-14
    ]

    print("\n" + "="*70)
    print("DIAGNOSING PROBLEM CASES (parity_loss ≈ 1.0)")
    print("="*70)
    for lsq, n in problem_cases:
        result = diagnose(lsq, n)

    print("\n" + "="*70)
    print("GOOD CASES FOR COMPARISON (parity_loss ≈ 0)")
    print("="*70)
    for lsq, n in good_cases:
        result = diagnose(lsq, n)

    print("\n" + "="*70 + "\n")
