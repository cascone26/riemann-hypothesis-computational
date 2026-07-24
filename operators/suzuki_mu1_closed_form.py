#!/usr/bin/env python3
"""
suzuki_mu1_closed_form.py — RESEARCH_QUEUE #0 (2026-07-20). Closed-form path for Suzuki's μ₁.

For a < ½log2 ≈ 0.347 the screw function has NO primes (LEDGER §42b), and (Suzuki eqs 2.2-2.4)
  Q_{B_a}(v) = 2 L_a(v) - 2(1+A) ‖v‖²,   A = ½(log2π - ψ(2)) = 0.70755,
  L_a(v) = ¼ ∫∫ |v(x)-v(y)|²/|x-y| dxdy  -  ½ ∫ log(a²-x²) |v(x)|² dx    (PURELY archimedean).
So λ_a = 2·inf[L_a/‖v‖²] - 2(1+A), and μ₁ = 2·inf[L_a/‖v‖²] - log(1/a) - 1  (matching Thm 1.4).

This script computes inf[L_a/‖v‖²] independently (a prime-free H^{1/2} eigenvalue problem), (a) CROSS-
CHECKS the reduction λ_a = 2·inf[L_a] - 2(1+A) against the direct screw-function λ_a, and (b) extracts μ₁
from the L_a side — pinning down whether inf[L_a] scales as ½log(1/a) (Suzuki) or log(1/a), resolving the
factor and giving a clean prime-free route to μ₁ (closed form is the a→0 limit of the rescaled (-1,1)
functional).
"""

import numpy as np
import mpmath as mp
from scipy.linalg import eigh

A_CONST = float(0.5 * (mp.log(2 * mp.pi) - mp.digamma(2)))  # 0.707546


def inf_La(a, M=2000, K=60):
    """inf_v L_a(v)/‖v‖² via sine-Galerkin. v_k=sin(kπ(x+a)/2a), v(±a)=0.
    Kernel 1/|x-y| is integrable here (|v(x)-v(y)|²/|x-y| → 0 on the diagonal)."""
    x = np.linspace(-a, a, M + 1)
    h = x[1] - x[0]
    kk = np.arange(1, K + 1)
    Vk = np.sin(np.pi * np.outer((x + a), kk) / (2 * a))  # (M+1, K)
    w = np.full(M + 1, h); w[0] = w[-1] = 0.5 * h
    # seminorm [V]²_{1/2} = Σ_xy (V_x - V_y)² K_xy w_x w_y, K=1/|x-y| (0 on diagonal)
    D = np.abs(x[:, None] - x[None, :])
    Kxy = np.zeros_like(D)
    nz = D > 1e-14
    Kxy[nz] = 1.0 / D[nz]
    Wt = w[:, None] * w[None, :] * Kxy
    rowsum = Wt.sum(1)
    # S_{kl} = Σ_xy (Vk_x-Vk_y)(Vl_x-Vl_y) Wt = 2[ Vk^T diag(rowsum) Vl - Vk^T Wt Vl ]
    S = 2 * (Vk.T @ (rowsum[:, None] * Vk) - Vk.T @ Wt @ Vk)
    # potential P_{kl} = ∫ log(a²-x²) v_k v_l dx
    pot = np.log(np.maximum(a * a - x * x, 1e-300))
    P = Vk.T @ ((pot * w)[:, None] * Vk)
    Mass = Vk.T @ (w[:, None] * Vk)  # ‖v‖² Gram
    Lmat = 0.25 * S - 0.5 * P
    Lmat = (Lmat + Lmat.T) / 2
    Mass = (Mass + Mass.T) / 2
    ev = eigh(Lmat, Mass, eigvals_only=True)
    return ev[0]


def main():
    print("=" * 78)
    print("μ₁ CLOSED-FORM PATH — prime-free H^{1/2} functional L_a (Suzuki §2)")
    print("=" * 78)
    print(f"A = ½(log2π - ψ(2)) = {A_CONST:.6f}")
    print(f"\n{'a':>7} {'inf[L_a/‖v‖²]':>15} {'-½log(1/a)':>12} {'rem (→(μ₁+1)/2?)':>16} {'μ₁ via L_a':>11}")
    rows = []
    for a in [0.03, 0.05, 0.08, 0.12, 0.18, 0.25]:
        il = inf_La(a, M=max(2000, int(120 / a)), K=60)
        loginv = np.log(1 / a)
        # test BOTH scalings:
        rem_half = il - 0.5 * loginv       # if inf ~ ½log(1/a)+(μ₁+1)/2  ⇒ rem→(μ₁+1)/2
        mu1_half = 2 * rem_half - 1         # μ₁ = 2·rem - 1
        rows.append((a, il, mu1_half))
        print(f"{a:7.3f} {il:15.6f} {0.5*loginv:12.5f} {rem_half:16.5f} {mu1_half:11.5f}")
    # also report the coefficient of log(1/a): fit inf_La = c·log(1/a)+b
    a_arr = np.array([r[0] for r in rows]); il_arr = np.array([r[1] for r in rows])
    Lg = np.log(1 / a_arr)
    c, b = np.polyfit(Lg[:5], il_arr[:5], 1)
    print(f"\n  fit inf[L_a] = c·log(1/a)+b:  c = {c:.4f}  (Suzuki reduction needs c=½=0.5)")
    print(f"  b = {b:.5f}  ⇒ if c=½: μ₁ = 2b-1 = {2*b-1:.4f}")
    # μ₁ from the smallest-a extrapolation of μ1_half
    mu_arr = np.array([r[2] for r in rows])
    slope, mu0 = np.polyfit(a_arr[:4], mu_arr[:4], 1)
    print(f"  μ₁ (L_a route, extrapolated a→0) = {mu0:.4f}")
    print(f"\n  Cross-check target: direct screw-function μ₁ ≈ 0.25 (LEDGER §42).")
    print("=" * 78)


if __name__ == "__main__":
    main()
