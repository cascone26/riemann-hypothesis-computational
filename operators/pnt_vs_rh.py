#!/usr/bin/env python3
"""
pnt_vs_rh.py — Does sc → 0 require RH or just PNT?

CRITICAL THEORETICAL QUESTION:
  sc(σ) = ||V_null^T τ V_null||_F / ||V_null^T (A+B) V_null||_F

We observe sc(0.5) → 0 as L^{-1.9}. Is this:
  (a) Unconditional consequence of PNT (prime distribution ~ N/logN)?
  (b) A consequence of RH (zeros on critical line)?
  (c) Something about the operator structure?

TEST STRATEGY:
Construct a FAKE operator τ_fake where:
1. We replace the prime sum P with a different function that follows PNT density
   but does NOT encode actual prime locations.
2. Specifically: P_random uses primes shifted by a random offset (still PNT-distributed).

If sc(τ_fake) → 0 at the same rate → sc convergence is PNT-only, NOT RH-sensitive.
If sc(τ_fake) → constant → sc convergence requires actual prime positions (which
   encode RH via their zeros).

This would definitively answer: is our sub-cancel result meaningful for RH?

IMPLEMENTATION:
- Build τ_true: actual prime sum
- Build τ_shifted: primes × e^{i×random_offset} (randomized positions)
- Compare sc for both

If τ_shifted gives same sc → PNT only.
If τ_shifted gives different sc → actual zero locations matter.
"""
import time, json, random
import numpy as np
from scipy.special import digamma as _digamma
import scipy.linalg
import mpmath as mp

mp.mp.dps = 30
TWO_PI = 2.0 * np.pi
LOGPI  = float(np.log(np.pi))
G      = 50000
THRESH = 0.02
OUTFILE = '/Users/scones/projects/riemann/results/pnt_vs_rh.out'

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

def build_tau_variants(lsq, N, shift=0.0):
    """
    Build τ with prime sum where phases are shifted by 'shift'.
    shift=0.0: true τ
    shift=φ: each log(k) replaced by log(k) + φ/N
    shift='random': each prime contribution gets independent random phase
    """
    L  = float(mp.log(lsq))
    sh = float(np.sinh(L / 4) ** 2)
    VM = vm_table(lsq)
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

    def make_P(phase_func):
        if VM:
            lkpp = np.array([np.log(k) for k, _ in VM])
            wts  = np.array([lk * k ** (-0.5) for k, lk in VM])
            # phase_func(lkpp, index) → modified phases
            lkpp_mod = phase_func(lkpp)
            ph   = TWO_PI * np.outer(ns, lkpp_mod) / L
            S    = np.sin(ph) @ wts
            dn   = ns[:, None] - ns[None, :]
            with np.errstate(divide='ignore', invalid='ignore'):
                P = (S[None, :] - S[:, None]) / (np.pi * dn)
            P[dn == 0] = 0.0
            np.fill_diagonal(P, np.cos(ph) @ (2 * wts * (1 - lkpp_mod / L)))
            return P
        return np.zeros((len(idx), len(idx)))

    # True P: no modification
    P_true = make_P(lambda lk: lk)

    # Phase-shifted P: shift all log(k) by a constant offset (preserves PNT structure)
    phase_offset = shift * L / (2 * np.pi)  # shift by 'shift' in log space
    P_shift = make_P(lambda lk: lk + phase_offset)

    # Random-phase P: each prime gets a random log shift (destroys correlations)
    rng = np.random.RandomState(42)
    random_offsets = rng.uniform(-np.pi, np.pi, size=len(VM)) * L / (2 * np.pi)
    P_rand = make_P(lambda lk: lk + random_offsets[:len(lk)])

    AB = arch + B

    results = {}
    for name, P in [('true', P_true), ('shifted', P_shift), ('random', P_rand)]:
        tau = (AB - P); tau = (tau + tau.T) / 2.0
        ev, V = scipy.linalg.eigh(tau)
        lam0 = float(ev[0])
        null_mask = np.abs(ev) < THRESH
        dim = int(null_mask.sum())
        if dim > 0:
            Vn = V[:, null_mask]
            ABp = Vn.T @ AB @ Vn
            Pp  = Vn.T @ P  @ Vn
            sc = float(np.linalg.norm(ABp - Pp, 'fro') / (np.linalg.norm(ABp, 'fro') + 1e-100))
        else:
            sc = float('nan')
        results[name] = {'sc': sc, 'lam0': lam0, 'dim': dim}

    return results, L

def main():
    configs = [(50, 250), (100, 500)]

    print(f"PNT vs RH test: does sc → 0 require actual prime positions?")
    print(f"shift = π (half-period phase offset)")
    print()
    with open(OUTFILE, 'w') as f:
        f.write(f"pnt_vs_rh {time.strftime('%Y-%m-%d %H:%M:%S')}\n\n")

    for lsq, N in configs:
        t0 = time.time()
        print(f"--- lsq={lsq}, N={N} ---")
        results, L = build_tau_variants(lsq, N, shift=np.pi)
        sc_pred = 3.073e-3 / L ** 1.9151

        print(f"  Predicted sc: {sc_pred:.4e}")
        for name, r in results.items():
            ratio = r['sc'] / sc_pred if sc_pred > 0 and not np.isnan(r['sc']) else float('nan')
            print(f"  {name:8s}: sc={r['sc']:.4e}  ratio={ratio:.3f}  lam0={r['lam0']:.4e}  dim={r['dim']}")

        print(f"  [{time.time()-t0:.0f}s]")
        print()
        with open(OUTFILE, 'a') as f:
            f.write(f"lsq={lsq} N={N} L={L:.4f}\n")
            f.write(f"  pred_sc={sc_pred:.4e}\n")
            for name, r in results.items():
                ratio = r['sc'] / sc_pred if sc_pred > 0 and not np.isnan(r['sc']) else float('nan')
                f.write(f"  {name:8s}: sc={r['sc']:.4e}  ratio={ratio:.3f}  lam0={r['lam0']:.4e}\n")
            f.write(f"  [{time.time()-t0:.0f}s]\n\n")

    print(f"INTERPRETATION:")
    print(f"  'true' ratio ~ 1.0: expected (baseline)")
    print(f"  'shifted' ratio ~ 1.0: PNT structure preserved, sc not prime-position-sensitive")
    print(f"  'shifted' ratio >> 1.0: prime phases matter → sc IS sensitive to actual positions")
    print(f"  'random' ratio >> 1.0: completely non-prime sc much larger → encoding is real")
    print(f"  'random' ratio ~ 1.0: sc is not sensitive to any structure beyond PNT density")

if __name__ == "__main__":
    main()
