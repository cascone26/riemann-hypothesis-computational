#!/usr/bin/env python3
"""
lambda_sq_scan.py — Scan many lambda^2 values, track near-null dim.

KEY QUESTION: does near-null dim = 2 * (# prime powers <= lambda^2)?
Or some other function of lambda^2?

Tests: lambda^2 in {4,5,6,7,8,9,10,11,12,13,14,16,18,20,25,36,49,64,100}
At each: N=48, count near-null dim and number of prime powers.
"""
import sys, json, os, time
import numpy as np
from scipy.special import digamma as _digamma
import scipy.linalg
import mpmath as mp

mp.mp.dps = 30
TWO_PI = 2.0 * np.pi
LOGPI  = float(np.log(np.pi))
THRESH = 1e-3
OUTFILE = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                       "../results/lambda_sq_scan.json")


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


def build_and_scan(lambda_sq, N=48):
    L  = float(mp.log(lambda_sq))
    sh = float(np.sinh(L/4)**2)
    VM = vm_table(lambda_sq)
    n_pp = len(VM)

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

    if len(VM) == 0:
        P = np.zeros((D,D))
    else:
        log_k_pp = np.array([np.log(k) for k,_ in VM])  # log(prime_power) for phases
        wts      = np.array([lk*k**(-0.5) for k,lk in VM])
        phases = TWO_PI*np.outer(ns,log_k_pp)/L
        S = np.sin(phases)@wts
        dn = ns[:,None]-ns[None,:]   # dn[i,j] = ns[i]-ns[j]
        with np.errstate(divide='ignore', invalid='ignore'):
            P = (S[None,:]-S[:,None])/(np.pi*dn)
        P[dn==0]=0.0
        np.fill_diagonal(P, np.cos(phases)@(2*wts*(1-log_k_pp/L)))

    tau = arch+B-P; tau=(tau+tau.T)/2
    ev = scipy.linalg.eigvalsh(tau)
    mask = np.abs(ev) < THRESH
    dim  = int(mask.sum())
    lam_min = float(ev[0])
    above = ev[~mask & (ev > 0)]
    gap   = float(above.min() - ev[mask].max()) if (dim > 0 and len(above) > 0) else float('nan')
    return {
        'lambda_sq': lambda_sq, 'L': round(L,6), 'N': N,
        'n_prime_powers': n_pp,
        'near_null_dim': dim,
        'dim_per_pp': round(dim/n_pp, 4) if n_pp > 0 else None,
        'dim_minus_2pp': dim - 2*n_pp,
        'lam_min': lam_min, 'spectral_gap': gap,
        'prime_powers': [k for k,_ in VM]
    }


def main():
    # Only scan integer lambda_sq values (non-integer would give non-integer kmax)
    lambdas = [4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 16, 18, 20, 25, 36, 49, 64, 100]
    N = 48
    results = []
    print(f"lambda^2 scan | N={N} | threshold={THRESH}")
    print(f"{'lam_sq':>7}  {'n_pp':>4}  {'dim':>4}  {'dim/n_pp':>9}  {'dim-2*pp':>9}  "
          f"{'lam_min':>12}  {'gap':>10}")
    for lsq in lambdas:
        t0 = time.time()
        try:
            r = build_and_scan(lsq, N)
            print(f"{lsq:>7}  {r['n_prime_powers']:>4}  {r['near_null_dim']:>4}  "
                  f"{r['dim_per_pp']:>9}  {r['dim_minus_2pp']:>9}  "
                  f"{r['lam_min']:>12.4e}  {r['spectral_gap']:>10.4e}  ({time.time()-t0:.1f}s)")
            sys.stdout.flush()
            results.append(r)
        except Exception as e:
            print(f"lam_sq={lsq} ERROR: {e}")

    # Summary: is dim always = 2*n_pp + constant?
    print("\n=== Summary: near_null_dim vs 2*n_prime_powers ===")
    for r in results:
        print(f"  lam_sq={r['lambda_sq']:>4}: dim={r['near_null_dim']}, 2*pp={2*r['n_prime_powers']}, "
              f"diff={r['dim_minus_2pp']}, dim/pp={r['dim_per_pp']}")

    with open(OUTFILE, 'w') as f:
        json.dump(results, f, indent=2)
    print(f"\nSaved: {OUTFILE}")


if __name__ == "__main__":
    main()
