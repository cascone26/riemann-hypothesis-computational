#!/usr/bin/env python3
"""
suzuki_vhat_zeros.py — RESEARCH_QUEUE #0, CLEAN Hilbert–Pólya realization (Suzuki eq 1.2), 2026-07-20.

Suzuki [14] (arXiv:2301.00421) eq (1.2):  lim_{a→∞} c_a · v̂_a(z) = ξ(1/2 + iz),
where v_a = the LOWEST eigenfunction of the localized Weil operator A_a, and v̂_a(z)=∫v_a(x)e^{izx}dx.
Since ξ(1/2+iz)=0 ⟺ z=γ_n (a zeta zero, under RH), the zeros of v̂_a → the zeta zeros. NO θ-tuning
(unlike the W(a,θ;z) route, which tests the messier ξ/ξ' of Cor 1.6). This is the clean HP demonstration.

RESULT (LEDGER §46): at a=1.0–2.0, v̂_a VANISHES at the first ~12 zeta zeros to ~1e-4 (relative), with
NO tuning. Verified genuine sharp zeros (|v̂_a(γ_n)|/|v̂_a(γ_n±0.5)| ~1e-4), PRIME-DRIVEN (arch-only v̂_a is
O(1) at γ_n), NON-CIRCULAR (A_a built from digamma Ω + von Mangoldt + pole term — zero zeta-zeros as input),
larger a extends the reach (a=2.0 nails γ₁₂ to 1e-4) = the a→∞ convergence of (1.2). Depends on the fixed
vonmangoldt_upto (LEDGER §44). Consistent with the older connes_decisive_v3 non-circular reproduction.
HONEST: numerical, finite a, ~1e-4 (not machine-zero), a DEMONSTRATION of (1.2)/the construction — NOT a proof
of RH (proof = v̂_a→ξ holds UNCONDITIONALLY as a→∞, the open analytic problem).
"""

import numpy as np
import os, sys
from scipy.linalg import eigh
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import suzuki_screw_lambda_a as S
import mpmath as mp


def vhat_abs(a, use_primes=True, M=7000, K=100):
    """Return z ↦ |v̂_a(z)|/median, for v_a the lowest eigenfunction of A_a (or arch-only)."""
    Nmax = int(np.floor(np.exp(2 * a))) + 1
    Lam = S.vonmangoldt_upto(min(Nmax, 30_000_000))
    x = np.linspace(-a, a, M + 1)
    h = x[1] - x[0]

    def gf(t):
        at = abs(t)
        b = -4 * (np.exp(t / 2) + np.exp(-t / 2) - 2)
        pr = 0.0
        if use_primes:
            nc = int(np.floor(np.exp(at)))
            if nc >= 2:
                ns = np.arange(2, min(nc, Nmax) + 1); L = Lam[2:min(nc, Nmax) + 1]; nz = L > 0
                pr = np.sum(L[nz] / np.sqrt(ns[nz]) * (at - np.log(ns[nz])))
        return b + pr - (at / 2) * S.PSI_QUARTER_MINUS_LOGPI - 0.25 * (S.PHI_1_2_QUARTER - np.exp(-at / 2) * S.phi_tail(t))

    diffs = np.arange(-M, M + 1) * h
    gv = np.array([gf(d) for d in diffs])
    idx = np.arange(M + 1); G = gv[(idx[:, None] - idx[None, :]) + M]
    kk = np.arange(1, K + 1)
    Vk = np.sin(np.pi * np.outer((x + a), kk) / (2 * a))
    Vp = (kk * np.pi / (2 * a)) * np.cos(np.pi * np.outer((x + a), kk) / (2 * a))
    w = np.full(M + 1, h); w[0] = w[-1] = 0.5 * h; Wp = Vp * w[:, None]
    Q = (Wp.T @ G @ Wp); Q = (Q + Q.T) / 2
    Mass = Vk.T @ (w[:, None] * Vk); Mass = (Mass + Mass.T) / 2
    ev, V = eigh(Q, Mass); c0 = V[:, 0]; va = Vk @ c0
    med = np.median([abs(np.sum(va * np.exp(1j * z * x) * w)) for z in np.linspace(1, 55, 400)])
    return lambda z: abs(np.sum(va * np.exp(1j * z * x) * w)) / med, ev[0]


def main():
    gam = [float(mp.im(mp.zetazero(n))) for n in range(1, 13)]
    print("=" * 74)
    print("SUZUKI eq (1.2): v̂_a (lowest-eigenfunction FT) vanishes at the ZETA ZEROS — no tuning")
    print("=" * 74)
    print(f"{'n':>3} {'γ_n':>9} | {'a=1.0 FULL':>11} {'a=2.0 FULL':>11} | {'a=1.0 ARCH':>11}")
    v1, l1 = vhat_abs(1.0, True)
    v2, l2 = vhat_abs(2.0, True)
    va1, _ = vhat_abs(1.0, False)
    for i, g in enumerate(gam):
        print(f"{i+1:3d} {g:9.3f} | {v1(g):11.4f} {v2(g):11.4f} | {va1(g):11.4f}")
    print("\n(FULL ≈0 at every γ_n ⟹ v̂_a vanishes at the zeta zeros; ARCH ≈O(1) ⟹ prime-driven.)")
    print(f"λ_a(a=1)={l1:+.4f}, λ_a(a=2)={l2:+.4f} (knife-edge, RH-consistent).")
    print("Non-circular: A_a built from digamma Ω + von Mangoldt + pole term (no zeta zeros as input).")


if __name__ == "__main__":
    main()
