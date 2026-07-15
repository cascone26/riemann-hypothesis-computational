#!/usr/bin/env python3
"""
fm_cauchy.py — Compare F_even (Groskin's Cauchy formula) vs our sinc F
applied to the same float64 eigenvector.

F_even(τ) = Re[e^{iτL/2} * Σ_k c_k * (e^{-iτL}-1)/(i*(2πk/L - τ))] / √L
           = -2√L * sin(τL/2) * Σ_k c_k / (2πk - τL)   [for real τ]

Our F(τ) = √L * Σ_n c_n * e^{iα_n/2} * sinc(α_n/2), α_n = 2πn - τL

Both have zeros at Riemann zeros. F_even is sharper (rational, not sinc).
"""
import sys, json, os, time
import numpy as np
from scipy.special import digamma as _digamma
import scipy.linalg
import mpmath as mp
from scipy.signal import argrelmin

mp.mp.dps = 50
TWO_PI = 2.0 * np.pi
LOGPI  = float(np.log(np.pi))
G      = 60000

KNOWN_ZEROS = [
    14.134725141734693, 21.022039638771555, 25.010857580145688,
    30.424876125859513, 32.935061587739189, 37.586178158825671,
    40.918719012147495, 43.327073280914999, 48.005150881167159,
    49.773832477672302,
]

OUTFILE = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                       "../results/fm_cauchy.json")


def vm_table(kmax):
    tab = []
    for k in range(2, int(kmax) + 1):
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


def build_tau(lambda_sq, N):
    L  = float(mp.log(lambda_sq))
    sh = float(np.sinh(L / 4) ** 2)
    VM = vm_table(lambda_sq)
    idx = np.arange(-N, N + 1, dtype=float)
    ns  = idx
    R   = TWO_PI * N / L + 80.0
    grid = np.linspace(-R, R, G)
    h = grid[1] - grid[0]; w = np.full(G, h); w[0] = w[-1] = 0.5 * h
    Om = np.real(_digamma(0.25 + 0.5j * grid)) - LOGPI
    d = (TWO_PI * ns / L)[:, None] - grid[None, :]
    tiny = np.abs(d) < 1e-12
    E = np.where(tiny, np.sqrt(L),
                 (1 / np.sqrt(L)) * np.exp(0.5j * d * L)
                 * 2 * np.sin(0.5 * d * L) / (d + 1e-300 * tiny))
    arch = ((E * (w * Om)) @ E.conj().T).real / (2.0 * np.pi)
    B = (32.0 * L * sh * (L ** 2 - 16.0 * np.pi ** 2 * np.outer(ns, ns)) /
         np.outer(L ** 2 + 16.0 * np.pi ** 2 * ns ** 2,
                  L ** 2 + 16.0 * np.pi ** 2 * ns ** 2))
    if VM:
        lkpp = np.array([np.log(k) for k, _ in VM])
        wts  = np.array([lk * k ** (-0.5) for k, lk in VM])
        ph   = TWO_PI * np.outer(ns, lkpp) / L
        S    = np.sin(ph) @ wts
        dn   = ns[:, None] - ns[None, :]
        with np.errstate(divide='ignore', invalid='ignore'):
            P = (S[None, :] - S[:, None]) / (np.pi * dn)
        P[dn == 0] = 0.0
        np.fill_diagonal(P, np.cos(ph) @ (2 * wts * (1 - lkpp / L)))
    else:
        P = np.zeros((len(idx), len(idx)))
    tau = arch + B - P
    return (tau + tau.T) / 2.0, idx.astype(int)


def min_kinetic_eigvec(ev, evec, idx, thresh=1e-2):
    mask = np.abs(ev) < thresh
    nn_evec = evec[:, mask]
    K = int(mask.sum())
    if K == 0:
        return evec[:, 0], 0
    n2 = (idx ** 2).astype(float)
    K_proj = (nn_evec * n2[:, None]).T @ nn_evec
    K_proj = (K_proj + K_proj.T) / 2
    kev, kevec = np.linalg.eigh(K_proj)
    c = nn_evec @ kevec[:, 0]
    return c / np.linalg.norm(c), K


def f_sinc(c, idx, L, z_grid):
    """Our Fourier-Mellin (sinc formula). Complex-valued on real z."""
    alpha = TWO_PI * idx[:, None] - L * z_grid[None, :]
    half  = 0.5 * alpha
    with np.errstate(divide='ignore', invalid='ignore'):
        kern = np.where(np.abs(half) < 1e-10, 1.0, np.sin(half) / half)
    return np.sqrt(L) * np.sum(c[:, None] * np.exp(1j * half) * kern, axis=0)


def f_cauchy(c, idx, L, z_grid):
    """Groskin's F_even (Cauchy formula). Real-valued on real z."""
    alpha = TWO_PI * idx[:, None] - L * z_grid[None, :]
    # F_even = -2√L * sin(τL/2) * Σ c_k / α_k
    # where α_k = 2πk - τL, τ = z
    sin_term = np.sin(0.5 * L * z_grid)  # sin(zL/2) for each z
    with np.errstate(divide='ignore', invalid='ignore'):
        cauchy_sum = np.where(np.abs(alpha) < 1e-10,
                              0.0,  # skip if singular (resonant)
                              c[:, None] / alpha)
    cauchy_total = np.sum(cauchy_sum, axis=0)
    return -2.0 * np.sqrt(L) * sin_term * cauchy_total


def analyze(c, idx, L, z_grid, label):
    F_sinc = np.abs(f_sinc(c, idx, L, z_grid))
    F_cau  = np.abs(f_cauchy(c, idx, L, z_grid))
    Fmax_s = F_sinc.max()
    Fmax_c = F_cau.max()

    print(f"\n  [{label}]")
    print(f"  {'k':>3}  {'gamma_k':>10}  {'|F_sinc|/max':>14}  {'|F_cauchy|/max':>16}  {'min_s':>6}  {'min_c':>6}")
    checks = []
    for ki, gk in enumerate(KNOWN_ZEROS[:8]):
        iz = np.argmin(np.abs(z_grid - gk))
        fs = float(F_sinc[iz]) / Fmax_s
        fc = float(F_cau[iz]) / Fmax_c
        lo, hi = max(0, iz-200), min(len(F_sinc)-1, iz+200)
        ms = '*' if F_sinc[iz] == F_sinc[lo:hi+1].min() else ' '
        mc = '*' if F_cau[iz] == F_cau[lo:hi+1].min() else ' '
        print(f"  {ki+1:>3}  {gk:>10.5f}  {fs:>14.4e}  {fc:>16.4e}  {ms:>6}  {mc:>6}")
        checks.append({'gamma': gk, 'F_sinc_rel': fs, 'F_cauchy_rel': fc, 'sinc_is_min': ms=='*', 'cauchy_is_min': mc=='*'})
    return checks


def refine_cauchy(c, idx, L, gk, dps=40):
    """Find zero of F_cauchy near gk using mpmath findroot."""
    mp.mp.dps = dps + 5
    c_mp = [mp.mpf(float(x)) for x in c]
    idx_mp = [mp.mpf(int(n)) for n in idx]
    L_mp = mp.mpf(L)
    sqrt2L = mp.sqrt(L_mp)
    pi2 = 2 * mp.pi

    def F_cau(z):
        z = mp.mpf(z)
        alpha = [pi2 * n - L_mp * z for n in idx_mp]
        total = sum(cn / a if abs(a) > 1e-100 else 0 for cn, a in zip(c_mp, alpha))
        return -2 * sqrt2L * mp.sin(L_mp * z / 2) * total

    try:
        mp.mp.dps = dps + 5
        root = mp.findroot(F_cau, mp.mpf(gk), tol=mp.mpf(10)**(-(dps-5)))
        err = float(root) - gk
        return float(root), err
    except Exception as e:
        return float(gk), float('nan')


def main():
    lambda_sq = 50.0
    N = 192

    L = float(np.log(lambda_sq))
    print(f"Building tau at lambda_sq={lambda_sq}, N={N}...", flush=True)
    t0 = time.time()
    tau, idx = build_tau(lambda_sq, N)
    ev, evec = scipy.linalg.eigh(tau)
    t1 = time.time()
    print(f"  build+eigh: {t1-t0:.1f}s  lam_min={ev[0]:.4e}", flush=True)

    c, K = min_kinetic_eigvec(ev, evec, idx, thresh=1e-2)
    print(f"  K(1e-2)={K}", flush=True)

    z_grid = np.linspace(5.0, 55.0, 100000)
    checks = analyze(c, idx, L, z_grid, f"lambda_sq={lambda_sq}")

    # Refine F_cauchy zeros with mpmath
    print(f"\n  Newton refinement of F_cauchy zeros (dps=30):")
    print(f"  {'k':>3}  {'gamma_true':>20}  {'z_cauchy':>20}  {'error':>14}  {'digits':>8}")
    refined = []
    for ki, gk in enumerate(KNOWN_ZEROS[:8]):
        z_ref, err = refine_cauchy(c, idx, L, gk, dps=30)
        digits = -np.log10(abs(err/gk)) if abs(err) > 0 else 30
        print(f"  {ki+1:>3}  {gk:>20.15f}  {z_ref:>20.15f}  {err:>+14.8e}  {digits:>8.3f}", flush=True)
        refined.append({'k': ki+1, 'gamma_true': gk, 'z_cauchy': z_ref, 'error': err, 'digits': digits})

    out = {'lambda_sq': lambda_sq, 'N': N, 'K': K, 'zero_checks': checks, 'refined': refined}
    with open(OUTFILE, 'w') as f:
        json.dump(out, f, indent=2)
    print(f"\nSaved: {OUTFILE}")


if __name__ == "__main__":
    main()
