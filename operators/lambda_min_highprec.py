#!/usr/bin/env python3
"""
lambda_min_highprec.py — Is lambda_min genuinely 0 or just near-zero?

Two approaches:
  1. Full mpmath tau at 100 dps for small N (N=8,12,16,20,24) — exact answer.
  2. Rayleigh quotient: float64 min eigenvector + mpmath B (exact closed form) + float64 arch/P.

Key: at 100 dps, 0.0 vs 1e-10 are distinguishable.
"""
import json, os, time
import numpy as np
import mpmath as mp
from scipy.special import digamma as _digamma

DPS = 100
LAMBDA_SQ = 14
OUTFILE = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                       "../results/lambda_min_highprec.json")


def vm_table_int(kmax):
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
            tab.append((k, p))   # (prime power, prime base)
    return tab


def build_tau_loop_f64(N, lambda_sq=14):
    """Correct loop-based float64 builder (matches original verified code)."""
    TWO_PI = 2.0 * np.pi
    L  = float(mp.log(lambda_sq))
    sh = float(np.sinh(L/4)**2)
    VM = [(k, float(np.log(p))) for k, p in vm_table_int(lambda_sq)]
    idx = list(range(-N, N+1))
    D   = 2*N+1
    ns  = np.array(idx, dtype=float)
    G   = 60000   # match original scripts — coarser G aliases arch at high N
    R   = TWO_PI * N / L + 80.0
    grid = np.linspace(-R, R, G)
    h = grid[1]-grid[0]; w = np.full(G,h); w[0]=w[-1]=0.5*h
    Om = np.real(_digamma(0.25+0.5j*grid)) - float(np.log(np.pi))
    d = (TWO_PI*ns/L)[:,None]-grid[None,:]
    tiny=(np.abs(d)<1e-12)
    E = np.where(tiny, np.sqrt(L),
                 (1/np.sqrt(L))*np.exp(0.5j*d*L)*2*np.sin(0.5*d*L)/(d+1e-300*tiny))
    arch = ((E*(w*Om))@E.conj().T).real/(TWO_PI)
    B = (32*L*sh*(L**2-16*np.pi**2*np.outer(ns,ns))/
         np.outer(L**2+16*np.pi**2*ns**2, L**2+16*np.pi**2*ns**2))
    log_k_pp = np.array([np.log(k) for k,_ in VM])  # log(prime_power) for phases
    wts      = np.array([lk*k**(-0.5) for k,lk in VM])
    P = np.zeros((D,D))
    for ii, n in enumerate(idx):
        for jj in range(ii, D):
            m = idx[jj]
            q = np.where(n==m,
                2*(1-log_k_pp/L)*np.cos(TWO_PI*n*log_k_pp/L),
                (np.sin(TWO_PI*m*log_k_pp/L)-np.sin(TWO_PI*n*log_k_pp/L))/(np.pi*(n-m)+1e-300))
            P[ii,jj]=P[jj,ii]=float(wts@q)
    tau=arch+B-P; return (tau+tau.T)/2, arch, B, P, idx


def build_tau_mpmath(N, lambda_sq, dps=100):
    """Full mpmath tau at dps precision. Slow, use only for small N."""
    mp.mp.dps = dps
    L   = mp.log(lambda_sq)
    sh  = mp.sinh(L/4)**2
    pi  = mp.pi
    # (prime_power, prime_base, log_prime_power, log_prime_base)
    VM  = [(k, p, mp.log(k), mp.log(p)) for k, p in vm_table_int(lambda_sq)]
    idx = list(range(-N, N+1))
    D   = 2*N+1
    ns  = [mp.mpf(n) for n in idx]

    # arch: computed in float64 then converted (accurate to 1e-14, much faster)
    tau_f, arch_f, _, _, _ = build_tau_loop_f64(N, lambda_sq)
    arch_mp = mp.matrix([[mp.mpf(float(arch_f[i,j])) for j in range(D)] for i in range(D)])

    # B: exact in mpmath
    B_mp = mp.matrix(D, D)
    for i in range(D):
        di = L**2 + 16*pi**2*ns[i]**2
        for j in range(D):
            dj = L**2 + 16*pi**2*ns[j]**2
            B_mp[i,j] = 32*L*sh*(L**2 - 16*pi**2*ns[i]*ns[j]) / (di*dj)

    # P: exact in mpmath (CORRECT SIGN: denominator = n_row - n_col)
    # VM = (k, p, lk_pp=log(prime_power), lk=log(prime_base))
    P_mp = mp.matrix(D, D)
    for i in range(D):
        ni = ns[i]
        for j in range(D):
            nj = ns[j]
            if i == j:
                P_mp[i,j] = sum(2*lk*mp.mpf(k)**mp.mpf(-0.5)
                                *(1-lk_pp/L)*mp.cos(2*pi*ni*lk_pp/L)
                                for k,p,lk_pp,lk in VM)
            else:
                P_mp[i,j] = sum(lk*mp.mpf(k)**mp.mpf(-0.5)
                                *(mp.sin(2*pi*nj*lk_pp/L)-mp.sin(2*pi*ni*lk_pp/L))
                                /(pi*(ni-nj))
                                for k,p,lk_pp,lk in VM)

    tau_mp = arch_mp + B_mp - P_mp
    return tau_mp


def rayleigh_mp(N, lambda_sq, dps=100):
    """
    Rayleigh quotient at dps precision.
    Get min eigenvector from float64, compute v^T tau_mp v in mpmath.
    B is closed-form in mpmath; arch+P from float64.
    """
    mp.mp.dps = dps
    tau_f, arch_f, B_f, P_f, idx = build_tau_loop_f64(N, lambda_sq)
    ev, evec = np.linalg.eigh(tau_f)
    v = evec[:, 0]    # float64 min eigenvector
    D = len(v)

    L   = mp.log(lambda_sq)
    sh  = mp.sinh(L/4)**2
    pi  = mp.pi
    vm  = mp.matrix([mp.mpf(float(x)) for x in v])
    ns  = [mp.mpf(n) for n in idx]

    # B Rayleigh quotient — exact: v^T B v = 32Lsh*(L²*s0²-16π²*s1²)
    denom = [L**2 + 16*pi**2*ns[i]**2 for i in range(D)]
    s0 = sum(vm[i]/denom[i] for i in range(D))
    s1 = sum(ns[i]*vm[i]/denom[i] for i in range(D))
    Bvv_mp = 32*L*sh*(L**2*s0*s0 - 16*pi**2*s1*s1)

    # arch and P from float64 (accurate to ~1e-14)
    Avv = float(v @ arch_f @ v)
    Pvv = float(v @ P_f @ v)

    # Total Rayleigh quotient
    rq = mp.mpf(Avv) + Bvv_mp - mp.mpf(Pvv)
    return {
        'N': N, 'lambda_sq': lambda_sq,
        'lam_min_f64': float(ev[0]),
        'Avv_f64': Avv, 'Bvv_mp': float(Bvv_mp), 'Pvv_f64': Pvv,
        'rq_hybrid': float(rq),
        'Bvv_mp_str': mp.nstr(Bvv_mp, 40),
        'rq_str': mp.nstr(rq, 40),
    }


def main():
    mp.mp.dps = DPS
    results = []

    # Full mpmath at small N
    print(f"=== Full {DPS}-dps tau eigenvalue (small N) ===")
    for N in [8, 12, 16, 20, 24]:
        t0 = time.time()
        try:
            tau_mp = build_tau_mpmath(N, LAMBDA_SQ, dps=DPS)
            ev_mp = mp.eigh(tau_mp, eigvals_only=True)
            lam_min = min(ev_mp)
            dt = time.time()-t0
            s = mp.nstr(lam_min, 30)
            print(f"  N={N:>3}: lam_min = {s}  ({dt:.1f}s)")
            results.append({'N': N, 'method': f'full_mpmath_{DPS}dps',
                            'lam_min': s, 'time': round(dt,2)})
        except Exception as e:
            print(f"  N={N} ERROR: {e}")
            import traceback; traceback.print_exc()

    # Rayleigh quotient at larger N
    print(f"\n=== Hybrid Rayleigh quotient (B exact at {DPS}dps, arch+P float64) ===")
    print(f"{'N':>4}  {'lam_min_f64':>14}  {'rq_hybrid':>14}  {'Avv':>10}  {'Bvv_mp':>14}  {'Pvv':>10}")
    for N in [16, 24, 32, 48, 64, 96, 128]:
        t0 = time.time()
        try:
            r = rayleigh_mp(N, LAMBDA_SQ, dps=DPS)
            dt = time.time()-t0
            print(f"{N:>4}  {r['lam_min_f64']:>14.6e}  {r['rq_hybrid']:>14.6e}  "
                  f"{r['Avv_f64']:>10.6f}  {r['Bvv_mp']:>14.6e}  {r['Pvv_f64']:>10.6f}  ({dt:.1f}s)")
            r['time'] = round(dt,2); results.append(r)
        except Exception as e:
            print(f"N={N} ERROR: {e}")
            import traceback; traceback.print_exc()

    with open(OUTFILE, 'w') as f:
        json.dump(results, f, indent=2, default=str)
    print(f"\nSaved: {OUTFILE}")


if __name__ == "__main__":
    main()
