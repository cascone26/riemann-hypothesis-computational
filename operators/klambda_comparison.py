#!/usr/bin/env python3
"""
CYCLE 127 — Direct comparison of ξ_λ (Weil-CS ground eigenvector) to k_λ (Connes' PSWF
educated guess). This tests Connes' GAP (2): "ξ_λ ≈ k_λ" — the main remaining obstacle.

From CCM [2511.22755] eq (line 1232/1266):
  h(u) = (π/2) u² (2πu² - 3) e^{-πu²}               [Hermite-like Gaussian]
  k_λ(u) = u^{1/2} Σ_{n≤λ/u} h(nu),  u ∈ [λ^{-1}, λ]  [Connes' cE map applied to h_λ≈h]

In log-coordinate x = log u, a = log λ:
  k_λ(e^x) = e^{x/2} Σ_{n=1}^{floor(e^{a-x})} h(n·e^x)

Compare to ξ_λ = θ(x) (ground eigenvector of Weil-CS form on [-a,a] in log-coords).

Metric: after finding the best scalar match c* = argmin_c ‖θ - c·k_λ‖,
  residual = ‖θ - c*·k_λ‖ / ‖k_λ‖

If residual → 0 as a→∞, the approximation ξ_λ ≈ k_λ is confirmed.

Note: k_λ defined this way uses h (not h_λ), the λ→∞ limit. The O(λ^{-1/2-α}) error from
using h vs h_λ is Connes' PROVEN rate (CCM Lemma hermfact1), separate from the ξ_λ≈k_λ step.
"""

import numpy as np
import sys, os
sys.path.insert(0, os.path.dirname(__file__))
from suzuki_screw_lambda_a import screw_g, vonmangoldt_upto

def h_hermite(u):
    """h(u) = (π/2) u² (2πu² - 3) e^{-πu²}  [CCM line 1232]"""
    return (np.pi / 2) * u**2 * (2 * np.pi * u**2 - 3) * np.exp(-np.pi * u**2)

def k_lambda_logcoord(a, xgrid):
    """
    Compute k_λ in log-coordinate: k_λ(e^x) = e^{x/2} Σ_{n=1}^{floor(e^{a-x})} h(n·e^x)
    xgrid: array of x values in [-a, a]
    Uses h (large-λ limit); the approximation error h_λ→h is O(λ^{-1/2-α}) (proven).
    """
    lam = np.exp(a)
    k = np.zeros(len(xgrid))
    for i, x in enumerate(xgrid):
        u = np.exp(x)
        nmax = int(np.floor(lam / u))  # n·u ≤ λ ⟺ n ≤ λ/u = e^{a-x}
        if nmax < 1:
            k[i] = 0.0
            continue
        ns = np.arange(1, nmax + 1)
        k[i] = u**0.5 * np.sum(h_hermite(ns * u))
    return k

def weil_min_eigvec(a, M):
    """Ground eigenvector θ(x) of the plain Weil-CS form on [-a,a], L²-normalized."""
    Nmax = min(int(np.floor(np.exp(2*a)))+1, 30_000_000)
    Lam = vonmangoldt_upto(Nmax)
    x = np.linspace(-a, a, M+1); h = x[1]-x[0]
    diffs = np.arange(-M, M+1)*h
    gv = np.array([screw_g(d, Lam, Nmax) for d in diffs])
    idx = np.arange(M+1); G = gv[(idx[:,None]-idx[None,:])+M]
    w = np.sqrt(np.full(M+1, h)); w[0]=w[-1]=np.sqrt(0.5*h)
    Gs = w[:,None]*G*w[None,:]; Gs=(Gs+Gs.T)/2
    ev, V = np.linalg.eigh(Gs)
    theta = V[:,0]/w
    theta /= np.sqrt(np.sum(theta*theta)*h)
    return x, h, theta, ev[0]

def best_scalar_match(theta, k, h):
    """c* = argmin_c ‖θ - c·k‖² = <θ,k>/<k,k>."""
    num = np.sum(theta * k) * h
    den = np.sum(k * k) * h
    if abs(den) < 1e-20:
        return 0.0
    return float(num / den)

print("="*70)
print("ξ_λ vs k_λ COMPARISON — testing Connes' gap (2): ξ_λ ≈ k_λ ?")
print("k_λ(e^x) = e^{x/2} Σ_{n=1}^{floor(e^{a-x})} h(n·e^x), h = Hermite/Gaussian")
print("="*70)
print(f"\n{'a':>4} {'M':>5} {'λ=e^a':>7} | {'‖k_λ‖':>7} {'c*':>8} {'‖θ-c*k‖':>9} {'residual':>9} {'cos(θ,k)':>9}")
print("-"*75)

results = []
for a in [2.0, 3.0, 4.0, 5.0, 6.0]:
    M = int(400 * a)
    xg, h_step, theta, ev0 = weil_min_eigvec(a, M)

    # k_λ on same grid
    k = k_lambda_logcoord(a, xg)

    norm_k = np.sqrt(np.sum(k * k) * h_step)
    if norm_k < 1e-15:
        print(f"{a:>4.1f} {M:>5d} {np.exp(a):>7.1f} | k_λ is ZERO on this grid")
        continue

    # Best scalar match
    c_star = best_scalar_match(theta, k, h_step)
    residual_vec = theta - c_star * k
    res_norm = np.sqrt(np.sum(residual_vec**2) * h_step)
    relative_res = res_norm / norm_k

    # Cosine similarity (angle between theta and k)
    cos_sim = np.sum(theta * k) * h_step / norm_k  # norm_theta = 1

    print(f"{a:>4.1f} {M:>5d} {np.exp(a):>7.1f} | {norm_k:>7.4f} {c_star:>8.4f} {res_norm:>9.4f} {relative_res:>9.4f} {cos_sim:>9.4f}")
    results.append((a, np.exp(a), relative_res, cos_sim))

print()
if len(results) >= 3:
    avals = np.array([r[0] for r in results])
    resids = np.array([r[2] for r in results])
    cosims = np.array([r[3] for r in results])

    print("Monotone check on relative residual ‖θ - c*k‖/‖k‖:")
    for i, r in enumerate(results):
        tag = "✓ " if i == 0 or resids[i] < resids[i-1] else "↑ "
        print(f"  {tag}a={r[0]:.1f}: residual={r[2]:.4f}  cos(θ,k)={r[3]:.4f}")

    print()
    if np.all(resids[1:] < resids[:-1]):
        print("  RESIDUAL MONOTONE DECREASING → consistent with ξ_λ→k_λ")
    else:
        print("  NON-MONOTONE — convergence to k_λ not confirmed by this diagnostic")

    print()
    print("Rate fit: residual vs a")
    log_r = np.log(resids + 1e-15)
    A = np.column_stack([np.ones(len(avals)), avals])
    c, _, _, _ = np.linalg.lstsq(A, log_r, rcond=None)
    print(f"  residual ~ exp({c[0]:.3f}) * e^({c[1]:.3f}*a)  → rate p={-c[1]:.3f}")
    if -c[1] > 0.1:
        print(f"  Exponential decay confirmed. Compare Connes' rate cλ^{{-1/2-α}} = ce^{{-(1/2+α)a}}: p={-c[1]:.3f} vs 0.5+α")
    else:
        print(f"  Rate p={-c[1]:.3f}: {'decay' if -c[1]>0 else 'growth — not converging with this metric'}")

print()
print("NOTE: If residual→0, Connes' gap (2) is confirmed numerically. The k_λ used here is the")
print("large-λ limit (h, not h_λ); the O(λ^{-1/2-α}) error from using h vs h_λ is already PROVEN.")
