#!/usr/bin/env python3
"""
weil_cancellation_analysis.py — Decompose tau_N into A+B-P for the minimum
eigenvector and study HOW the near-cancellation A(g)+B(g)≈P(g) arises.

For the minimum eigenvector g_min of tau_N (the direction where tau is
tightest), we have:
    A(g_min) + B(g_min) - P(g_min) = lambda_min ≈ 0^+

So A+B and P nearly cancel for this specific function.  Studying the ratio
A/(A+B), B/(A+B), P/(A+B) vs N reveals the structure of the cancellation.

If the cancellation has an analytic explanation tied to the digamma symbol Omega,
it might suggest the proof path for uniform positivity.

Also probes:
- Omega decomposition: Omega(r) term-by-term vs r
- The r-distribution of the integrand |g_hat(1/2+ir)|^2 * Omega(r)
  weighted by the minimum eigenvector (where is most of A(g_min) coming from?)
"""

import json, os, sys
import numpy as np
from scipy.special import digamma as _digamma
import mpmath as mp

mp.mp.dps = 30

LAMBDA_SQ = 14
L = float(2 * mp.log(mp.sqrt(LAMBDA_SQ)))
TWO_PI = 2.0 * np.pi
LOGPI = float(np.log(np.pi))
OUTDIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "../results")
OUTFILE = os.path.join(OUTDIR, "weil_cancellation_analysis.json")


def vm_table(kmax):
    tab = []
    for k in range(2, kmax + 1):
        n, p, f, ok = k, 0, 2, True
        while f * f <= n:
            if n % f == 0:
                if p == 0: p = f
                elif p != f: ok = False; break
                while n % f == 0: n //= f
            f += 1
        if n > 1:
            if p == 0: p = n
            elif p != n: ok = False
        if ok and p > 0:
            tab.append((k, float(mp.log(p))))
    return tab

VM = vm_table(LAMBDA_SQ)


def omega_vec(r):
    return np.real(_digamma(0.25 + 0.5j * r)) - LOGPI


def ehat_vec(r, n):
    s = TWO_PI * n / L
    d = s - r
    out = np.empty_like(d, dtype=complex)
    sm = np.abs(d) < 1e-12
    out[sm] = np.sqrt(L) + 0j
    dd = d[~sm]
    out[~sm] = (1.0/np.sqrt(L)) * np.exp(0.5j*dd*L) * 2.0*np.sin(0.5*dd*L) / dd
    return out


def build_tau_and_components(N):
    idx = list(range(-N, N+1))
    D = len(idx)
    R = TWO_PI * N / L + 80.0
    G = max(60_000, 3000 * N)
    grid = np.linspace(-R, R, G)
    h = grid[1] - grid[0]
    w = np.full(G, h); w[0] = w[-1] = 0.5*h
    Om = omega_vec(grid)

    E = np.empty((D, G), dtype=complex)
    for k, n in enumerate(idx):
        E[k] = ehat_vec(grid, n)

    # Archimedean term
    A_mat = (E * (w * Om) @ E.conj().T).real / (2*np.pi)

    # Boundary/pole term
    sh = np.sinh(L/4)**2
    ns = np.array(idx, dtype=float)
    num = L**2 - 16*np.pi**2 * np.outer(ns, ns)
    den = (L**2 + 16*np.pi**2 * ns**2)
    B_mat = 32*L*sh * num / np.outer(den, den)

    # Prime term — y argument to q is log(k), weight is Lambda(k)*k^{-1/2}
    log_k = np.array([np.log(k) for k, _ in VM])
    weights = np.array([lk * k**(-0.5) for k, lk in VM])
    P_mat = np.zeros((D, D))
    for ii, n in enumerate(idx):
        for jj in range(ii, D):
            m = idx[jj]
            q = np.where(n == m,
                         2*(1-np.abs(log_k)/L)*np.cos(TWO_PI*n*log_k/L),
                         (np.sin(TWO_PI*m*log_k/L)-np.sin(TWO_PI*n*log_k/L))/(np.pi*(n-m)+1e-300))
            P_mat[ii,jj] = P_mat[jj,ii] = float(weights @ q)

    tau = (A_mat + B_mat - P_mat)
    tau = (tau + tau.T) / 2
    return tau, A_mat, B_mat, P_mat, idx, grid, w, Om, E


def ghat_on_grid(v, idx, grid):
    E = np.empty((len(idx), len(grid)), dtype=complex)
    for k, n in enumerate(idx):
        E[k] = ehat_vec(grid, n)
    return v @ E   # (G,) complex


def analyze_N(N):
    print(f"\n=== N={N} ===", flush=True)
    tau, A_mat, B_mat, P_mat, idx, grid, w, Om, E = build_tau_and_components(N)

    ev, evec = np.linalg.eigh(tau)
    lam_min = float(ev[0])
    v_min = evec[:, 0]  # minimum eigenvector

    # Decompose quadratic forms for the minimum eigenvector
    A_val = float(v_min @ A_mat @ v_min)
    B_val = float(v_min @ B_mat @ v_min)
    P_val = float(v_min @ P_mat @ v_min)
    total = A_val + B_val - P_val

    print(f"  lambda_min = {lam_min:.6e}", flush=True)
    print(f"  A(g_min)   = {A_val:.6e}", flush=True)
    print(f"  B(g_min)   = {B_val:.6e}", flush=True)
    print(f"  P(g_min)   = {P_val:.6e}", flush=True)
    print(f"  A+B-P      = {total:.6e}  (should match lambda_min)", flush=True)
    print(f"  A/(A+B)    = {A_val/(A_val+B_val+1e-300):.6f}", flush=True)
    print(f"  B/(A+B)    = {B_val/(A_val+B_val+1e-300):.6f}", flush=True)
    print(f"  P/(A+B)    = {P_val/(A_val+B_val+1e-300):.6f}", flush=True)

    # Spatial distribution of integrand: Omega(r) * |g_hat(1/2+ir)|^2
    g_hat = v_min @ E   # (G,) complex
    integrand = Om * np.abs(g_hat)**2
    density = w * integrand / (2*np.pi)   # contribution to A per grid point

    # Where does A(g_min) come from? Find the r-distribution
    # Split into positive (Omega>0) and negative (Omega<0) contributions
    pos_mask = Om > 0
    A_pos = float(density[pos_mask].sum())
    A_neg = float(density[~pos_mask].sum())
    r_Om_zero = float(np.abs(grid[np.argmin(np.abs(Om))]))
    print(f"  Omega zero-crossing at |r| ≈ {r_Om_zero:.3f}", flush=True)
    print(f"  A positive part: {A_pos:.6e}", flush=True)
    print(f"  A negative part: {A_neg:.6e}  (sum = {A_pos+A_neg:.6e})", flush=True)

    # Where is |g_hat|^2 concentrated?
    power = w * np.abs(g_hat)**2
    total_power = power.sum()
    peak_r = float(np.abs(grid[np.argmax(np.abs(g_hat)**2)]))
    power_in_neg = float(power[~pos_mask].sum()) / (total_power + 1e-300)
    print(f"  |g_hat|^2 peak at |r| ≈ {peak_r:.3f}", flush=True)
    print(f"  Fraction of |g_hat|^2 in Omega<0 region: {power_in_neg:.4f}", flush=True)

    # Omega decomposition: term-by-term Poisson components at r=0,5,10,50
    probe_r = np.array([0.0, 5.0, 10.0, 20.0, 50.0, 100.0])
    omega_probe = omega_vec(probe_r)

    # Sensitivity to a (perturbation): compute d(lambda_min)/d(a) numerically at a=0.25
    da = 0.001
    def omega_a(r, a):
        return np.real(_digamma(a + 0.5j*r)) - LOGPI
    Om_plus  = omega_a(grid, 0.25 + da)
    Om_minus = omega_a(grid, 0.25 - da)
    dA_da_vec = (Om_plus - Om_minus) / (2*da) * np.abs(g_hat)**2
    dA_da = float((w * dA_da_vec / (2*np.pi)).sum())
    print(f"  dA(g_min)/da|_{{a=1/4}} ≈ {dA_da:.6e}", flush=True)

    return {
        "N": N, "lambda_min": lam_min,
        "A_val": A_val, "B_val": B_val, "P_val": P_val, "total": total,
        "A_pos": A_pos, "A_neg": A_neg,
        "omega_zero_crossing_r": r_Om_zero,
        "g_hat_peak_r": peak_r,
        "power_fraction_in_omega_neg": power_in_neg,
        "omega_at_probe_r": dict(zip(probe_r.tolist(), omega_probe.tolist())),
        "dA_da_at_a025": dA_da,
        "5_smallest_ev": ev[:5].tolist(),
    }


def main():
    print("WEIL CANCELLATION ANALYSIS", flush=True)
    print(f"lambda^2={LAMBDA_SQ}, L={L:.4f}", flush=True)

    N_values = [16, 24, 32, 48, 64]
    results = []

    for N in N_values:
        r = analyze_N(N)
        results.append(r)
        with open(OUTFILE, "w") as f:
            json.dump({"lambda_sq": LAMBDA_SQ, "L": L, "results": results}, f, indent=2)
        print(f"  Saved.", flush=True)

    # Cross-N analysis
    print("\n=== CROSS-N PATTERNS ===", flush=True)
    print(f"{'N':>4}  {'lam_min':>12}  {'A':>12}  {'B':>12}  {'P':>12}  "
          f"{'A/(A+B)':>8}  {'dA/da':>12}", flush=True)
    for r in results:
        print(f"  {r['N']:>4}  {r['lambda_min']:>12.4e}  {r['A_val']:>12.4e}  "
              f"{r['B_val']:>12.4e}  {r['P_val']:>12.4e}  "
              f"{r['A_val']/(r['A_val']+r['B_val']+1e-300):>8.5f}  "
              f"{r['dA_da_at_a025']:>12.4e}", flush=True)

    print(f"\nResults: {OUTFILE}", flush=True)


if __name__ == "__main__":
    main()
