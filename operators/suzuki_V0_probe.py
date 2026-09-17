#!/usr/bin/env python3
"""Frontier probe (2026-08-26): does Suzuki's V(0)=L²(0,∞)∩K·L²(0,∞) ≠ 0? (RH is FALSE if V(0)=0, Rmk 5.2).

NON-CIRCULAR route: the whole de Branges construction (E, Θ, K, V(0)) flows from the SCREW FUNCTION
g=gξ = A+B-P (archimedean digamma Ω + pole + von Mangoldt primes) — computable from the explicit formula
WITHOUT the zeta zeros. V(0)≠0 ⟺ the ground state v_a of the localized Weil/screw operator A_a converges to a
NONZERO L²(0,∞) function as a→∞ (the a→∞ limit = the K-eigenfunction). So we track v_a (built exactly as in
suzuki_vhat_zeros, verified: v̂_a→ξ zeros) across growing a and test:
  (1) does the NORMALIZED ground state v_a(x) on (0,a) CONVERGE in shape as a grows? (overlap on common support)
  (2) is its mass CONCENTRATED (genuine L²(0,∞), decaying) vs SPREADING (would fail L²)?
  (3) does the eigenvalue λ_a (knife-edge) stay ≥0?
CONVERGE + CONCENTRATED ⟹ V(0)≠0 numerically, RH-consistent, and characterizes the K-eigenfunction shape.
"""
import numpy as np, os, sys
from scipy.linalg import eigh
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import suzuki_screw_lambda_a as S


def ground_state(a, M=6000, K=120):
    """Ground eigenfunction v_a(x) of the localized Weil/screw operator on [-a,a], and λ_a. (as suzuki_vhat_zeros)"""
    Nmax = int(np.floor(np.exp(2 * a))) + 1
    Lam = S.vonmangoldt_upto(min(Nmax, 30_000_000))
    x = np.linspace(-a, a, M + 1); h = x[1] - x[0]

    def gf(t):
        at = abs(t)
        b = -4 * (np.exp(t / 2) + np.exp(-t / 2) - 2)
        pr = 0.0
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
    # sign-fix and L2-normalize on [-a,a]
    va = va * np.sign(va[np.argmax(np.abs(va))])
    va = va / np.sqrt(np.sum(va**2 * w))
    return x, va, ev[0]


def main():
    print("=" * 78)
    print("FRONTIER PROBE: V(0) = L²(0,∞) ∩ K·L²(0,∞) ≠ 0 ?  (RH FALSE if V(0)=0, Suzuki Rmk 5.2)")
    print("Non-circular: ground state v_a of the screw operator; V(0)≠0 ⟺ v_a → nonzero L²(0,∞) as a→∞")
    print("=" * 78)
    aa = [1.0, 1.5, 2.0, 2.5, 3.0]
    data = {}
    for a in aa:
        x, va, lam = ground_state(a)
        # restrict to (0, a): the L²(0,∞) part
        pos = x > 0
        xp = x[pos]; vp = va[pos]
        w = (x[1] - x[0])
        norm_pos = np.sqrt(np.sum(vp**2) * w)           # L²(0,a) mass of the ground state
        # concentration: fraction of positive-side mass in inner [0, a/2] vs outer [a/2, a]
        inner = xp <= a / 2
        frac_inner = np.sum(vp[inner]**2) / max(np.sum(vp**2), 1e-30)
        data[a] = (xp, vp, norm_pos, frac_inner, lam)
        print(f"  a={a:.1f}: λ_a={lam:+.4f}  ‖v_a‖_{{L²(0,a)}}={norm_pos:.4f}  "
              f"mass-frac in [0,a/2]={frac_inner:.3f}  (concentrated if →1, spreading if →0.5)")
    # convergence: overlap of normalized v_a on common support (0, min a = 1.0)
    print("\n  Shape convergence (|⟨v_a, v_a'⟩| on common support [0, 1.0], normalized):")
    ref_a = aa[-1]
    xr, vr, *_ = data[ref_a]
    for a in aa[:-1]:
        xp, vp, *_ = data[a]
        # interpolate both onto common grid on [0, 1.0]
        cg = np.linspace(0.02, 0.98, 300)
        fa = np.interp(cg, xp, vp); fr = np.interp(cg, xr, vr)
        ov = abs(np.dot(fa, fr)) / (np.linalg.norm(fa) * np.linalg.norm(fr) + 1e-30)
        print(f"    a={a:.1f} vs a={ref_a:.1f}: overlap on [0,1] = {ov:.4f}")
    print("\n  VERDICT: convergent shape (overlap→1) + concentrated mass (frac→1) ⟹ v_a → nonzero L²(0,∞) limit")
    print("           ⟹ V(0)≠0 numerically ⟹ RH-consistent, and this limit IS the K-eigenfunction (Rmk 5.2).")


if __name__ == "__main__":
    main()
