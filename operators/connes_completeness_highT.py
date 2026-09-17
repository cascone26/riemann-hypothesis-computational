#!/usr/bin/env python3
"""
connes_completeness_highT.py — §1209 follow-up: does the Connes/Weil spectral construction stay
SPECTRALLY COMPLETE at high T? Fill-in ratio = (#true zeros matched < tol) / (#true zeros in window),
swept over increasing prime reach LAMBDA_SQ, dimension N, and window height T.

Same non-circular construction as connes_decisive_v3 (arch = digamma Omega, prime = von Mangoldt,
boundary = W_{0,2}; NO zeta zeros as input). If fill-in -> 1 as (LAMBDA_SQ,N) grow to resolve height T,
the operator captures ALL zeros in-window (completeness). If it plateaus < 1, incompleteness signal.
"""
import sys, time
import numpy as np
import mpmath as mp
from scipy.special import digamma as _digamma

mp.mp.dps = 30
TWO_PI = 2.0 * np.pi
LOGPI = float(np.log(np.pi))

# true zeros (independent reference) from the 10k dataset
ZTRUE = np.load("operators/zeros10k.npy")

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

def omega_vec(r):
    return np.real(_digamma(0.25 + 0.5j * r)) - LOGPI

def e_hat_vec(r, n, L):
    s = TWO_PI * n / L
    d = s - r
    out = np.empty_like(d, dtype=complex)
    small = np.abs(d) < 1e-12
    out[small] = np.sqrt(L) + 0j
    dd = d[~small]
    out[~small] = (1.0/np.sqrt(L)) * np.exp(0.5j*dd*L) * 2.0*np.sin(0.5*dd*L)/dd
    return out

def q_kernel(n, m, y, L):
    if n == m:
        return 2.0*(1.0 - abs(y)/L)*np.cos(TWO_PI*n*y/L)
    return (np.sin(TWO_PI*m*y/L) - np.sin(TWO_PI*n*y/L))/(np.pi*(n-m))

def prime_mat(idx, VM, L):
    D = len(idx); P = np.zeros((D, D))
    for ii, n in enumerate(idx):
        for jj in range(ii, D):
            m = idx[jj]; s = 0.0
            for k, Lam in VM:
                s += Lam * k**(-0.5) * q_kernel(n, m, np.log(k), L)
            P[ii, jj] = P[jj, ii] = s
    return P

def w02_mat(idx, L):
    D = len(idx); W = np.zeros((D, D)); sh = np.sinh(L/4)**2
    for ii, n in enumerate(idx):
        for jj, m in enumerate(idx):
            num = (L**2 - 16*np.pi**2*m*n)
            den = (L**2 + 16*np.pi**2*m**2)*(L**2 + 16*np.pi**2*n**2)
            W[ii, jj] = 32*L*sh*num/den
    return W

def build_tau(N, L, VM):
    idx = list(range(-N, N+1)); D = len(idx)
    R = TWO_PI*N/L + 80.0
    G = max(60000, 3000*N)
    grid = np.linspace(-R, R, G); h = grid[1]-grid[0]
    w = np.full(G, h); w[0]=w[-1]=0.5*h
    Om = omega_vec(grid)
    E = np.empty((D, G), dtype=complex)
    for k, n in enumerate(idx): E[k] = e_hat_vec(grid, n, L)
    A = E*(w*Om)
    arch = (A @ E.conj().T).real/(2.0*np.pi)
    tau = arch + w02_mat(idx, L) - prime_mat(idx, VM, L)
    return (tau+tau.T)/2, idx

def xi_hat_factory(xi, idx, L):
    poles = np.array([TWO_PI*j/L for j in idx]); xi = np.array(xi, float)
    def xh(z):
        d = z - poles; d[np.abs(d)<1e-12] = 1e-12
        return (2.0/np.sqrt(L))*np.sin(z*L/2.0)*np.sum(xi/d)
    return xh

def find_zeros(xh, lo, hi, ng):
    xs = np.linspace(lo, hi, ng); vs = np.array([xh(x) for x in xs]); out = []
    for i in range(1, len(xs)):
        if vs[i-1] == 0 or (vs[i-1] < 0) != (vs[i] < 0):
            a, b = xs[i-1], xs[i]
            for _ in range(80):
                mid = 0.5*(a+b)
                if (xh(a)<0) != (xh(mid)<0): b = mid
                else: a = mid
            r = 0.5*(a+b)
            if not out or abs(r-out[-1]) > 1e-4: out.append(r)
    return out

def run(LAMBDA_SQ, N, lo, hi, tol=0.01):
    L = float(np.log(LAMBDA_SQ))
    VM = vm_table(LAMBDA_SQ)
    t0 = time.time()
    tau, idx = build_tau(N, L, VM)
    ev, evec = np.linalg.eigh(tau)
    k0 = int(np.argmin(np.abs(ev)))
    D = 2*N+1; v = evec[:, k0]
    even = np.array([(v[i]+v[D-1-i])/2 for i in range(D)])
    s = even.sum()/np.sqrt(L)
    if abs(s) > 1e-30: even = even/s
    xh = xi_hat_factory(even, idx, L)
    got = find_zeros(xh, lo, hi, max(4000, int((hi-lo)*80)))
    true_in = ZTRUE[(ZTRUE >= lo) & (ZTRUE <= hi)]
    matched = 0; errs = []
    for tz in true_in:
        if got:
            g = min(got, key=lambda z: abs(z-tz)); e = abs(g-tz)
            if e < tol: matched += 1; errs.append(e)
    ratio = matched/len(true_in) if len(true_in) else float('nan')
    med = float(np.median(errs)) if errs else float('nan')
    return dict(LAMBDA_SQ=LAMBDA_SQ, N=N, L=L, nprimes=len(VM), window=(lo,hi),
                expected=len(true_in), matched=matched, ratio=ratio,
                med_err=med, eps_min=float(ev[0]), secs=time.time()-t0)

if __name__ == "__main__":
    # sweep: climb the window height T, grow (LAMBDA_SQ, N) to try to resolve it
    configs = [
        # (LAMBDA_SQ, N, lo, hi)
        (14,  32,  10,  60),
        (50,  48,  10,  90),
        (100, 64,  10, 120),
        (200, 96,  10, 160),
        (400, 128, 10, 200),
    ]
    print(f"{'LSQ':>5} {'N':>4} {'nP':>5} {'window':>12} {'exp':>4} {'match':>5} {'ratio':>6} {'medErr':>9} {'eps_min':>11} {'s':>6}")
    for (lsq, N, lo, hi) in configs:
        r = run(lsq, N, lo, hi)
        print(f"{r['LAMBDA_SQ']:>5} {r['N']:>4} {r['nprimes']:>5} "
              f"[{lo:>4},{hi:>4}] {r['expected']:>4} {r['matched']:>5} "
              f"{r['ratio']:>6.3f} {r['med_err']:>9.2e} {r['eps_min']:>11.3e} {r['secs']:>5.1f}",
              flush=True)
