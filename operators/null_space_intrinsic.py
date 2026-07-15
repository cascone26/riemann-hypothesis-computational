#!/usr/bin/env python3
"""
null_space_intrinsic.py — Test if the near-null subspace of τ is intrinsic.

HYPOTHESIS: The near-null eigenvectors of τ(σ=0.5) are approximately
the same as those of τ(σ=0.5+ε), even though the eigenvalues change.

If true: the null space is "pinned" by the functional equation symmetry,
independent of the precise zero locations. This would be a new structural
result, not circular.

TEST: Compute the principal angles between V_null(σ=0.5) and V_null(σ=0.5+ε)
for ε = 0.01, 0.05, 0.1, 0.2. If principal angles are small (cosines ≈ 1),
the subspaces are nearly identical → null space is intrinsic.

lsq=50, N=250 (fast reference)
lsq=100, N=500 (medium)
"""
import time, json
import numpy as np
from scipy.special import digamma as _digamma
import scipy.linalg
import mpmath as mp

mp.mp.dps = 30
TWO_PI = 2.0 * np.pi
LOGPI  = float(np.log(np.pi))
G      = 40000
THRESH = 0.02
OUTFILE = '/Users/scones/projects/riemann/results/null_space_intrinsic.out'

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

def build_tau(lsq, N, sigma=0.5):
    L  = float(mp.log(lsq))
    sh = float(np.sinh(L / 4) ** 2)
    VM = vm_table(lsq)
    idx = np.arange(-N, N + 1, dtype=float)
    ns  = idx

    R   = TWO_PI * N / L + 80.0
    grid = np.linspace(-R, R, G)
    h = grid[1] - grid[0]; w = np.full(G, h); w[0] = w[-1] = 0.5 * h

    # A matrix: use sigma-shifted digamma kernel
    # At sigma=0.5: standard Omega(xi) = Re psi(1/4 + i*xi/2) - log pi
    # At sigma=0.5+eps: shift the digamma argument
    sigma_arg = sigma  # 0.5 corresponds to 1/4 + i*xi/2
    # Standard: A_{mn} = (1/2pi) int Omega(xi) E_m*(xi) E_n(xi) dxi
    # where Omega = Re psi(1/4 + i*xi/2) - log pi
    # For sigma deviation: the digamma shifts the spectrum but A still uses 1/4
    # We implement the sensitivity test: perturb each zero's Lorentzian contribution
    # For simplicity: use the full digamma (sigma deviation shifts eigenvalues only slightly)
    # The EXACT way: modify the prime sum weight from 1/k^(1/2) to 1/k^sigma
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
        # KEY: use sigma in the prime weight: Λ(k) / k^sigma instead of k^0.5
        wts  = np.array([lk * k ** (-sigma) for k, lk in VM])
        ph   = TWO_PI * np.outer(ns, lkpp) / L
        S    = np.sin(ph) @ wts
        dn   = ns[:, None] - ns[None, :]
        with np.errstate(divide='ignore', invalid='ignore'):
            P = (S[None, :] - S[:, None]) / (np.pi * dn)
        P[dn == 0] = 0.0
        np.fill_diagonal(P, np.cos(ph) @ (2 * wts * (1 - lkpp / L)))
    else:
        P = np.zeros((len(idx), len(idx)))

    AB = arch + B
    tau = (AB - P); tau = (tau + tau.T) / 2.0
    return tau, AB

def null_subspace(tau, thresh=THRESH):
    ev, V = scipy.linalg.eigh(tau)
    null_mask = np.abs(ev) < thresh
    dim = int(null_mask.sum())
    lam_min = float(ev[0])
    return V[:, null_mask], ev[null_mask], dim, lam_min

def principal_angles(A, B):
    """Compute cosines of principal angles between column spaces of A and B."""
    Qa, _ = np.linalg.qr(A)
    Qb, _ = np.linalg.qr(B)
    M = Qa.T @ Qb
    sv = np.linalg.svd(M, compute_uv=False)
    return np.clip(sv, 0, 1)

def main():
    test_cases = [
        (50, 250),
        (100, 500),
    ]
    sigmas = [0.5, 0.505, 0.51, 0.52, 0.55, 0.6]

    print(f"Null space intrinsicness test")
    print(f"Testing: do eigenvectors of V_null change when σ moves from 0.5?")
    print()

    with open(OUTFILE, 'w') as f:
        f.write(f"null_space_intrinsic {time.strftime('%Y-%m-%d %H:%M:%S')}\n\n")

    for lsq, N in test_cases:
        print(f"lsq={lsq}, N={N}")
        t0 = time.time()

        # Reference: sigma=0.5
        tau0, AB0 = build_tau(lsq, N, sigma=0.5)
        Vn0, ev0, dim0, lam0 = null_subspace(tau0)
        print(f"  σ=0.500: dim={dim0} lam_min={lam0:.4e}")

        with open(OUTFILE, 'a') as f:
            f.write(f"lsq={lsq} N={N} dim(σ=0.5)={dim0} lam_min(0.5)={lam0:.4e}\n")
            f.write(f"{'sigma':>8} {'dim':>5} {'lam_min':>12} {'mean_cos':>10} {'min_cos':>10} {'overlap':>10}\n")
            f.write("-" * 65 + "\n")

        for sigma in sigmas:
            t1 = time.time()
            tau_s, AB_s = build_tau(lsq, N, sigma=sigma)
            Vn_s, ev_s, dim_s, lam_s = null_subspace(tau_s)
            dt = time.time() - t1

            # Principal angles between V_null(0.5) and V_null(sigma)
            # Use the smaller subspace to avoid dimension mismatch
            d_min = min(dim0, dim_s)
            if d_min > 0:
                cosines = principal_angles(Vn0[:, :d_min], Vn_s[:, :d_min])
                mean_cos = float(np.mean(cosines))
                min_cos = float(np.min(cosines))
                # Subspace overlap = normalized sum of squared cosines
                overlap = float(np.sum(cosines**2)) / d_min
            else:
                mean_cos = min_cos = overlap = float('nan')

            status = "SAME" if min_cos > 0.99 else ("CLOSE" if min_cos > 0.95 else "DIFFERENT")
            print(f"  σ={sigma:.3f}: dim={dim_s:3d} lam_min={lam_s:.4e}  mean_cos={mean_cos:.4f} min_cos={min_cos:.4f}  → {status} ({dt:.0f}s)")

            with open(OUTFILE, 'a') as f:
                f.write(f"  {sigma:>8.4f} {dim_s:>5} {lam_s:>12.4e} {mean_cos:>10.4f} {min_cos:>10.4f} {overlap:>10.4f}  {status}\n")

        print(f"  [lsq={lsq} done in {time.time()-t0:.0f}s]")
        print()
        with open(OUTFILE, 'a') as f:
            f.write(f"  [done in {time.time()-t0:.0f}s]\n\n")

    print(f"INTERPRETATION:")
    print(f"  mean_cos ≈ 1.0: null space is intrinsic (good for proof)")
    print(f"  mean_cos << 1.0: null space changes with σ (circular dependency confirmed)")
    print(f"  See {OUTFILE}")

if __name__ == "__main__":
    main()
