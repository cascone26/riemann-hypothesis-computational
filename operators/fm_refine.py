#!/usr/bin/env python3
"""
fm_refine.py — Newton refinement of Fourier-Mellin zero candidates.

Takes the min_kinetic eigenvector at lambda_sq=50, N=192 and uses mpmath
findroot to refine zero locations beyond grid resolution.

Also sweeps lambda_sq=50,100,200 to see how precision scales.

Shows how many significant digits we recover at each lambda_sq.
"""
import sys, json, os, time
import numpy as np
from scipy.special import digamma as _digamma
import scipy.linalg
import mpmath as mp

TWO_PI = 2.0 * np.pi
LOGPI  = float(np.log(np.pi))
G      = 60000

KNOWN_ZEROS = [
    14.134725141734693790, 21.022039638771554993, 25.010857580145688764,
    30.424876125859513210, 32.935061587739189690, 37.586178158825671257,
    40.918719012147495187, 43.327073280914999519,
]

OUTFILE = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                       "../results/fm_refine.json")


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


def make_fm_func(c_np, idx_np, L, dps=50):
    """Return mpmath-precision Fourier-Mellin function F(z) for complex z."""
    mp.mp.dps = dps + 10
    c_mp  = [mp.mpf(float(x)) for x in c_np]
    idx_mp = [mp.mpf(int(n)) for n in idx_np]
    sqL   = mp.sqrt(L)
    L_mp  = mp.mpf(L)
    pi2   = 2 * mp.pi

    def F(z):
        z = mp.mpc(z)
        s = mp.mpc(0)
        for cn, n in zip(c_mp, idx_mp):
            alpha = pi2 * n - L_mp * z
            half  = alpha / 2
            if mp.fabs(half) < mp.mpf('1e-20'):
                kern = mp.mpc(1)
            else:
                kern = mp.sin(half) / half
            s += cn * mp.exp(mp.mpc(0, 1) * half) * kern
        return sqL * s

    return F


def coarse_scan(c, idx, L, Nz=300000):
    """Fast float64 scan for local minima."""
    alpha = TWO_PI * idx[:, None] - L * np.linspace(5.0, 55.0, Nz)[None, :]
    half  = 0.5 * alpha
    with np.errstate(divide='ignore', invalid='ignore'):
        kern = np.where(np.abs(half) < 1e-10, 1.0, np.sin(half) / half)
    F_abs = np.abs(np.sqrt(L) * np.sum(c[:, None] * np.exp(1j * half) * kern, axis=0))
    z_grid = np.linspace(5.0, 55.0, Nz)
    return z_grid, F_abs


def run_case(lambda_sq, N, dps_refine=25):
    L = float(mp.log(lambda_sq))
    print(f"\n{'='*65}")
    print(f"lambda_sq={lambda_sq}, N={N}, L={L:.4f}", flush=True)

    t0 = time.time()
    tau, idx = build_tau(lambda_sq, N)
    ev, evec = scipy.linalg.eigh(tau)
    t1 = time.time()
    print(f"  build+eigh: {t1-t0:.1f}s  lam_min={ev[0]:.4e}", flush=True)

    c, K = min_kinetic_eigvec(ev, evec, idx, thresh=1e-2)
    print(f"  K(1e-2)={K}", flush=True)

    z_grid, F_abs = coarse_scan(c, idx, L)
    F = make_fm_func(c, idx, L, dps=dps_refine)

    Fmax = F_abs.max()

    zero_results = []
    print(f"  {'k':>3}  {'gamma_true':>20}  {'z_refined':>22}  {'error':>14}  {'digits':>8}  {'|F0|/max':>12}")
    for ki, gk in enumerate(KNOWN_ZEROS):
        if gk > lambda_sq:
            continue

        # Find grid minimum near gk
        iz = np.argmin(np.abs(z_grid - gk))
        wlo = max(0, iz - 1000)
        whi = min(len(F_abs) - 1, iz + 1000)
        local_min_iz = wlo + np.argmin(F_abs[wlo:whi])
        z0 = float(z_grid[local_min_iz])
        f0_rel = float(F_abs[local_min_iz]) / float(Fmax)

        # Newton refinement with mpmath
        try:
            mp.mp.dps = dps_refine + 10
            z_ref = mp.findroot(F, mp.mpf(z0), solver='muller', tol=mp.mpf(10)**(-dps_refine))
            err = float(z_ref - gk)
            # Count matching significant digits
            if abs(err) > 0:
                digits = -float(mp.log10(abs(err / gk)))
            else:
                digits = dps_refine
            z_str = mp.nstr(z_ref, dps_refine)
        except Exception as e:
            z_ref = z0
            err = z0 - gk
            digits = 0 if abs(err) > 1e-10 else 6
            z_str = str(z0)

        print(f"  {ki+1:>3}  {gk:>20.15f}  {z_str:>22}  {err:>+14.8f}  {digits:>8.2f}  {f0_rel:>12.4e}", flush=True)

        zero_results.append({
            'k': ki + 1,
            'gamma_true': gk,
            'z_refined': str(z_ref),
            'error': err,
            'digits': digits,
            'F_abs_rel': f0_rel,
        })

    return {
        'lambda_sq': lambda_sq, 'N': N, 'K': K,
        'lam_min': float(ev[0]),
        'zeros': zero_results,
    }


def main():
    cases = [
        (50, 192, 15),
        (100, 256, 15),
        (50, 300, 15),
    ]

    all_results = []
    for lsq, N, dps in cases:
        try:
            r = run_case(lsq, N, dps_refine=dps)
            all_results.append(r)
        except Exception as e:
            print(f"ERROR lsq={lsq} N={N}: {e}")
            import traceback; traceback.print_exc()

    with open(OUTFILE, 'w') as f:
        json.dump(all_results, f, indent=2)
    print(f"\nSaved: {OUTFILE}")


if __name__ == "__main__":
    import scipy
    main()
