#!/usr/bin/env python3
"""
eigvec_modes.py — What Fourier modes dominate each near-null eigenvector?

The 24 near-null eigenvectors c^(k) = [c_n : n=-N..N] at lambda^2=14, N=96.
Each eigenvector represents f(t) = (1/sqrt(L)) * sum_n c_n * exp(2*pi*i*n*t/L).

Key questions:
  1. Are dominant modes correlated with prime powers (n ~ log(p)/L * N)?
  2. Do symmetric modes concentrate at different n than antisymmetric?
  3. Is there a pattern in |c_n|^2 (power spectrum)?
  4. Do modes cluster at specific "prime-resonant" frequencies?
  5. What is |c_n|^2 summed over the whole 24-dim subspace?
"""
import json, os
import numpy as np
import scipy.linalg
from scipy.special import digamma as _digamma
import mpmath as mp

mp.mp.dps = 50
TWO_PI = 2.0 * np.pi
LOGPI  = float(np.log(np.pi))
G      = 60000
N      = 96
LAMBDA_SQ = 14
THRESH = 1e-3
OUTFILE = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                       "../results/eigvec_modes.json")


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
    D, ns = len(idx), idx
    R = TWO_PI * N / L + 80.0
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
        P = np.zeros((D, D))
    tau = arch + B - P
    return (tau + tau.T) / 2.0, idx.astype(int)


def sym_class(c, idx):
    D = len(c); vs = np.zeros(D); va = np.zeros(D)
    for i in range(D):
        j = D - 1 - i
        vs[i] = (c[i] + c[j]) / 2
        va[i] = (c[i] - c[j]) / 2
    sf = float(np.dot(vs, vs))
    return sf, 'SYM' if sf > 0.5 else 'ANTI'


def main():
    L = float(mp.log(LAMBDA_SQ))
    VM = vm_table(LAMBDA_SQ)
    prime_powers = [k for k, _ in VM]
    # "prime-resonant" Fourier indices: n such that 2*pi*n/L * log(k)/L is near integer
    # Actually: prime power k sits at log(k), which is at Fourier index 2*pi*n/L ~ log(k)/L * 2*pi*N
    # Better: n_k = round(N * log(k) / L) would be the "resonant" mode for prime k
    resonant_n = {k: round(N * float(mp.log(k)) / L) for k, _ in VM}
    print(f"Prime resonant Fourier indices (N={N}, L={L:.4f}):")
    for k, _ in VM:
        print(f"  k={k:>3}: log(k)/L={float(mp.log(k))/L:.4f}, resonant n={resonant_n[k]}")

    tau, idx = build_tau(LAMBDA_SQ, N)
    ev, evec = scipy.linalg.eigh(tau)
    mask = np.abs(ev) < THRESH
    nn_evs  = ev[mask]
    nn_evec = evec[:, mask]
    dim = int(mask.sum())
    print(f"\nN={N}, lambda_sq={LAMBDA_SQ}: dim={dim}, thresh={THRESH}")

    # Power spectrum of each eigenvector
    nn_power = nn_evec ** 2   # (D, 24) — power at each mode n for each eigvec

    # Subspace total power at each Fourier index
    subspace_power = np.sum(nn_power, axis=1)  # (D,) — sum over all 24 modes

    # For each near-null eigvec: dominant modes and prime-resonance
    records = []
    print(f"\n{'k':>3} {'eig':>12} {'type':>5}  {'top5 modes (n, power)':}")
    for k in range(dim):
        c  = nn_evec[:, k]
        sf, kind = sym_class(c, idx)
        power = c ** 2
        top5_idx = np.argsort(power)[-5:][::-1]
        top5 = [(int(idx[i]), float(power[i])) for i in top5_idx]

        # Check if dominant mode matches a prime resonance
        dominant_n = idx[np.argmax(power)]
        closest_prime = min(prime_powers, key=lambda pp: abs(resonant_n[pp] - dominant_n))
        closest_dist  = abs(resonant_n[closest_prime] - dominant_n)

        print(f"{k:>3} {nn_evs[k]:>12.4e} {kind:>5}  "
              + ' '.join(f"n={n}:{p:.3f}" for n, p in top5[:3])
              + f"  [nearest prime: k={closest_prime}, dist={closest_dist}]")
        records.append({
            'k': k, 'eigenvalue': float(nn_evs[k]),
            'sym_frac': float(sf), 'kind': kind,
            'top10_modes': [(int(idx[i]), float(power[i]))
                            for i in np.argsort(power)[-10:][::-1]],
            'power_spectrum': power.tolist(),
            'dominant_n': int(dominant_n),
            'closest_prime_k': int(closest_prime),
            'closest_prime_dist': int(closest_dist),
        })

    # Subspace power spectrum (sum over all 24 eigvecs)
    print(f"\n=== Subspace power spectrum (sum |c_n|^2 over all {dim} near-null modes) ===")
    print("Top 20 Fourier indices by subspace power:")
    top20_idx = np.argsort(subspace_power)[-20:][::-1]
    for i in top20_idx:
        n_val = idx[i]
        sp    = subspace_power[i]
        # Is this n close to a prime resonance?
        closest_prime = min(prime_powers, key=lambda pp: abs(resonant_n[pp] - n_val))
        dist = abs(resonant_n[closest_prime] - n_val)
        print(f"  n={n_val:>5}: power={sp:.4f}  "
              f"[nearest prime resonance: k={closest_prime}, dist={dist}]")

    # Prime resonance check: power at exactly resonant indices
    print(f"\n=== Subspace power at prime-resonant Fourier indices ===")
    for k, lk in VM:
        nr = resonant_n[k]
        idx_pos = np.where(idx == nr)[0]
        idx_neg = np.where(idx == -nr)[0]
        p_pos = float(subspace_power[idx_pos[0]]) if len(idx_pos) > 0 else 0
        p_neg = float(subspace_power[idx_neg[0]]) if len(idx_neg) > 0 else 0
        print(f"  k={k:>3} (log(k)/L={float(mp.log(k))/L:.4f}): "
              f"n=+{nr}: power={p_pos:.4f}  n=-{nr}: power={p_neg:.4f}")

    with open(OUTFILE, 'w') as f:
        json.dump({
            'N': N, 'lambda_sq': LAMBDA_SQ, 'dim': dim, 'L': L,
            'prime_powers': prime_powers,
            'resonant_n': {str(k): v for k, v in resonant_n.items()},
            'eigvecs': records,
            'subspace_power': subspace_power.tolist(),
            'idx': idx.tolist(),
        }, f, indent=2)
    print(f"\nSaved: {OUTFILE}")


if __name__ == "__main__":
    main()
