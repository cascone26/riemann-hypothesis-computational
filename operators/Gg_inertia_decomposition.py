#!/usr/bin/env python3
"""
CYCLE 127 — G_g inertia decomposition.

Tests whether G_g(t,u) = g(t-u) - g(t) - g(-u) + g(0) is nonneg-def (PSD) for:
  (A) full g (arch + boundary + prime) — should be PSD iff RH holds
  (B) arch+boundary only (no primes) — UNCONDITIONAL test
  (C) prime only — the RH-conditional piece

If (B) is PSD unconditionally, then #neg≤1 for the plain form is ONLY blocked by
the prime term (and that block is exactly what RH resolves). If (B) is NOT PSD,
then even the arch+boundary kernel is not in G_∞, and the screw membership is a
global property of the whole assembly.

This directly answers NEXT LEAD (A) from cycle 126.
"""

import numpy as np
import sys, os
sys.path.insert(0, os.path.dirname(__file__))
from suzuki_screw_lambda_a import screw_g, screw_g_no_primes, vonmangoldt_upto

def screw_g_prime_only(t, Lam, Nmax):
    """Only the prime sum: g_prime(t) = Σ Λ(n)/√n (|t|-log n) for n≤exp|t|."""
    at = abs(t)
    ncut = int(np.floor(np.exp(at)))
    if ncut < 2:
        return 0.0
    ns = np.arange(2, min(ncut, Nmax) + 1)
    L = Lam[2:min(ncut, Nmax) + 1]
    nz = L > 0
    if not np.any(nz):
        return 0.0
    return float(np.sum(L[nz] / np.sqrt(ns[nz]) * (at - np.log(ns[nz]))))

def build_Gg_matrix(a, N, g_func, Lam=None, Nmax=None):
    """Discretize G_g(t_i, t_j) = g(t_i-t_j) - g(t_i) - g(-t_j) + g(0) on N points in [-a,a]."""
    tg = np.linspace(-a, a, N)
    g0 = g_func(0.0, Lam, Nmax) if Lam is not None else g_func(0.0)

    def geval(t):
        return g_func(t, Lam, Nmax) if Lam is not None else g_func(t)

    # g values at grid points
    gvals = np.array([geval(ti) for ti in tg])
    gvals_neg = np.array([geval(-ti) for ti in tg])

    # G_g matrix: G[i,j] = g(ti-tj) - g(ti) - g(-tj) + g(0)
    G = np.zeros((N, N))
    for i in range(N):
        for j in range(N):
            G[i, j] = geval(tg[i] - tg[j]) - gvals[i] - gvals_neg[j] + g0
    G = 0.5 * (G + G.T)  # symmetrize
    return G

def count_neg(ev, tol=1e-9):
    scale = max(1.0, np.abs(ev).max())
    return int(np.sum(ev < -tol * scale))

def inertia_report(ev, tol=1e-9):
    scale = max(1.0, np.abs(ev).max())
    nneg = int(np.sum(ev < -tol * scale))
    npos = int(np.sum(ev > tol * scale))
    nzero = len(ev) - nneg - npos
    min_ev = ev.min()
    max_neg = ev[ev < 0].max() if nneg > 0 else float('nan')
    return nneg, npos, nzero, min_ev, max_neg

print("="*80)
print("G_g INERTIA DECOMPOSITION: full / arch+bdry-only / prime-only")
print("G_g(t,u) = g(t-u) - g(t) - g(-u) + g(0); PSD ⟺ g ∈ G_∞ (screw function class)")
print("Full g ∈ G_∞ ⟺ RH (Suzuki Thm 1.2 Su22).")
print("="*80)

# Load Mangoldt table for primes up to exp(10) ≈ 22026
Nmax_table = 30000
Lam = vonmangoldt_upto(Nmax_table)

N = 40  # grid points in [-a,a]; keep small for speed

for a in [2.0, 3.0, 4.0, 5.0, 6.0]:
    print(f"\na = {a:.1f}  (N={N} grid points on [-a,a])")

    # --- FULL g ---
    G_full = build_Gg_matrix(a, N, screw_g, Lam, Nmax_table)
    ev_full = np.linalg.eigvalsh(G_full)
    nn_f, np_f, nz_f, minev_f, maxneg_f = inertia_report(ev_full)
    print(f"  FULL g:        #neg={nn_f}  #zero={nz_f}  #pos={np_f}  "
          f"min_ev={minev_f:+.3e}  {'PSD ✓' if nn_f==0 else f'FAILS: largest-neg={maxneg_f:+.3e}'}")

    # --- ARCH+BOUNDARY only (no primes) ---
    def g_arch_bdry(t, *_):
        return screw_g_no_primes(t, None, None)
    G_arch = build_Gg_matrix(a, N, lambda t, L, M: screw_g_no_primes(t, L, M), Lam, Nmax_table)
    ev_arch = np.linalg.eigvalsh(G_arch)
    nn_a, np_a, nz_a, minev_a, maxneg_a = inertia_report(ev_arch)
    print(f"  ARCH+BDRY g:   #neg={nn_a}  #zero={nz_a}  #pos={np_a}  "
          f"min_ev={minev_a:+.3e}  {'PSD ✓' if nn_a==0 else f'FAILS: largest-neg={maxneg_a:+.3e}'}")

    # --- PRIME only ---
    def g_prime(t, L, M):
        return screw_g_prime_only(t, L, M)
    G_prime = build_Gg_matrix(a, N, g_prime, Lam, Nmax_table)
    ev_prime = np.linalg.eigvalsh(G_prime)
    nn_p, np_p, nz_p, minev_p, maxneg_p = inertia_report(ev_prime)
    print(f"  PRIME-ONLY g:  #neg={nn_p}  #zero={nz_p}  #pos={np_p}  "
          f"min_ev={minev_p:+.3e}  {'PSD ✓' if nn_p==0 else f'FAILS: largest-neg={maxneg_p:+.3e}'}")

print()
print("="*80)
print("INTERPRETATION:")
print("  If ARCH+BDRY is PSD (G_∞) unconditionally → the RH-content is entirely in PRIME term.")
print("  If ARCH+BDRY FAILS PSD → even the arch assembly is not a screw fn; G_∞ is global property.")
print("  If PRIME-ONLY FAILS PSD → that's the RH-conditional piece, consistent with Suzuki Thm 1.2.")
print("  FULL g PSD everywhere → consistent with RH for primes up to exp(2a).")
print("="*80)
