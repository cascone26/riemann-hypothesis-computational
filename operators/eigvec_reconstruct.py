#!/usr/bin/env python3
"""
eigvec_reconstruct.py — Reconstruct near-null eigenvectors as actual functions.

The eigenvectors c = [c_n : n=-N..N] represent:
  f(t) = (1/sqrt(L)) * Σ_n c_n * exp(2πint/L)    t ∈ [0, L]
where t = log(x) for x ∈ [1, lambda^2], L = log(lambda^2).

Prime powers k ≤ lambda^2 live at t = log(k) — special points.

Questions:
  1. What do these functions look like at the prime-power positions?
  2. Symmetric vs antisymmetric eigenvectors — different shapes?
  3. Recognizable special functions?
"""
import json, os
import numpy as np
import scipy.linalg
from scipy.special import digamma as _digamma
import mpmath as mp

mp.mp.dps = 30
TWO_PI = 2.0 * np.pi
LOGPI  = float(np.log(np.pi))
LAMBDA_SQ = 14
THRESH = 1e-3
OUTFILE = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                       "../results/eigvec_reconstruct.json")
N = 96   # use N=96: 24 near-null eigenvectors confirmed


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


def build_tau(N, lambda_sq=14):
    L  = float(mp.log(lambda_sq))
    sh = float(np.sinh(L/4)**2)
    VM = vm_table(lambda_sq)
    idx = np.arange(-N, N+1, dtype=float)
    D   = len(idx)
    R   = TWO_PI * N / L + 80.0
    G   = max(4000, int(np.ceil(4*R)))
    grid = np.linspace(-R, R, G)
    h = grid[1]-grid[0]; w = np.full(G,h); w[0]=w[-1]=0.5*h
    Om = np.real(_digamma(0.25+0.5j*grid)) - LOGPI
    svals = TWO_PI*idx/L
    d = svals[:,None]-grid[None,:]; tiny=(np.abs(d)<1e-12)
    E = np.where(tiny, np.sqrt(L),
                 (1/np.sqrt(L))*np.exp(0.5j*d*L)*2*np.sin(0.5*d*L)/(d+1e-300*tiny))
    arch = ((E*(w*Om))@E.conj().T).real/(2*np.pi)
    ns = idx
    B = (32*L*sh*(L**2-16*np.pi**2*np.outer(ns,ns))/
         np.outer(L**2+16*np.pi**2*ns**2, L**2+16*np.pi**2*ns**2))
    log_k_pp = np.array([np.log(k) for k,_ in VM])   # log(prime_power) for phases
    wts      = np.array([lk*k**(-0.5) for k,lk in VM])
    phases = TWO_PI*np.outer(ns,log_k_pp)/L
    S = np.sin(phases)@wts
    dn = ns[:,None]-ns[None,:]    # dn[i,j] = ns[i]-ns[j]
    with np.errstate(divide='ignore', invalid='ignore'):
        P = (S[None,:]-S[:,None])/(np.pi*dn)
    P[dn==0]=0.0
    np.fill_diagonal(P, np.cos(phases)@(2*wts*(1-log_k_pp/L)))
    tau = arch+B-P; return (tau+tau.T)/2, idx.astype(int), L


def reconstruct(c, idx, L, t_pts):
    """f(t) = (1/sqrt(L)) * Σ_n c_n * exp(2πint/L)."""
    n_arr = idx.astype(float)
    f = np.array([np.sum(c * np.exp(2j*np.pi*n_arr*t/L)) / np.sqrt(L)
                  for t in t_pts])
    return np.real(f)


def sym_frac(c, idx):
    D = len(c); vs = np.zeros(D); va = np.zeros(D)
    for i in range(D):
        j = D-1-i
        vs[i] = (c[i]+c[j])/2; va[i] = (c[i]-c[j])/2
    return float(np.dot(vs,vs))


def main():
    tau, idx, L = build_tau(N, LAMBDA_SQ)
    ev, evec = scipy.linalg.eigh(tau)
    mask = np.abs(ev) < THRESH
    nn_evs  = ev[mask]
    nn_evec = evec[:, mask]
    dim = int(mask.sum())
    print(f"N={N}, lambda_sq={LAMBDA_SQ}, near-null dim={dim}, L={L:.4f}")
    if dim == 0:
        print(f"lam_min={ev.min():.6e}, no near-null cluster — check tau build")
        return
    print(f"eigenvalues range: [{nn_evs.min():.3e}, {nn_evs.max():.3e}]")

    VM = vm_table(LAMBDA_SQ)
    t_primes = np.array([np.log(k) for k, _ in VM])
    pp_labels = [f"log({k})" for k, _ in VM]
    t_eval = np.concatenate([[0.0], t_primes, [L]])
    eval_labels = ['t=0'] + pp_labels + ['t=L']
    t_uniform = np.linspace(0, L, 200)

    results = []
    sym_fracs = [sym_frac(nn_evec[:,k], idx) for k in range(dim)]
    sym_labels = ['SYM' if s > 0.5 else 'ANTI' for s in sym_fracs]

    print(f"\n{'k':>4}  {'eig':>12}  {'type':>5}  {'sym_f':>6}  "
          + '  '.join(f"{lbl:>8}" for lbl in eval_labels[:6]))
    for k in range(dim):
        c  = nn_evec[:, k]
        sf = sym_fracs[k]
        fv = reconstruct(c, idx, L, t_eval)
        vals = '  '.join(f"{v:>8.4f}" for v in fv[:6])
        print(f"{k:>4}  {nn_evs[k]:>12.5e}  {sym_labels[k]:>5}  {sf:>6.3f}  {vals}")

        f_profile = reconstruct(c, idx, L, t_uniform)
        # amplitude at each prime power (normalized by max)
        amp = np.max(np.abs(f_profile)) + 1e-12
        pp_vals = reconstruct(c, idx, L, t_primes)
        results.append({
            'k': k, 'eigenvalue': float(nn_evs[k]),
            'sym_frac': sf, 'kind': sym_labels[k],
            'f_at_eval': fv.tolist(), 'eval_labels': eval_labels,
            'f_at_primes': pp_vals.tolist(), 'pp_labels': pp_labels,
            'f_profile': f_profile.tolist(), 't_uniform': t_uniform.tolist(),
            'f_norm_max': float(amp),
            'dominant_modes': [
                {'n': int(idx[i]), 'coeff': float(c[i]), 'weight': float(c[i]**2)}
                for i in np.argsort(c**2)[-6:][::-1]
            ]
        })

    # Pattern analysis
    print(f"\n=== Function values at prime powers (rows=eigvec, cols=primes) ===")
    print(f"{'k':>3} {'type':>5}  " + '  '.join(f"{lbl:>9}" for lbl in pp_labels))
    for k in range(dim):
        pp_vals = results[k]['f_at_primes']
        vals = '  '.join(f"{v:>9.4f}" for v in pp_vals)
        print(f"{k:>3} {sym_labels[k]:>5}  {vals}")

    # Legendre polynomial comparison
    print(f"\n=== Legendre correlation (first 8 symmetric eigvecs) ===")
    from numpy.polynomial import legendre as Pleg
    t_norm = t_uniform / L  # [0,1]
    sym_ks = [k for k in range(dim) if sym_fracs[k] > 0.5]
    for k in sym_ks[:8]:
        profile = np.array(results[k]['f_profile'])
        profile = profile / (np.max(np.abs(profile))+1e-10)
        corrs = []
        for deg in range(8):
            leg = Pleg.legval(2*t_norm-1, [0]*deg+[1])
            leg = leg / (np.sqrt(np.mean(leg**2))+1e-10)
            profile_n = profile / (np.sqrt(np.mean(profile**2))+1e-10)
            corrs.append(float(abs(np.mean(profile_n*leg))))
        best_deg = int(np.argmax(corrs))
        print(f"  eigvec {k} (eig={nn_evs[k]:.3e}): best Legendre P_{best_deg}, corr={max(corrs):.4f}")

    # Indicator function check: is f(t) ≈ constant on [0, L]?
    print(f"\n=== Smoothness: std(f) / mean(|f|) for each eigvec ===")
    for k in range(min(8, dim)):
        profile = np.array(results[k]['f_profile'])
        ratio = np.std(profile) / (np.mean(np.abs(profile))+1e-10)
        print(f"  eigvec {k} ({sym_labels[k]}): std/mean = {ratio:.3f}  "
              f"[1=rough, <0.5=smooth, 0=constant]")

    with open(OUTFILE, 'w') as f:
        json.dump(results, f, indent=2)
    print(f"\nSaved: {OUTFILE}")


if __name__ == "__main__":
    main()
