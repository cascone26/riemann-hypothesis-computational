#!/usr/bin/env python3
"""
null_space_coverage.py — Test the null-space covering hypothesis.

HYPOTHESIS: The deepest resonances occur when the number of aligned zeros
(Δ_k < threshold) matches the dimension D of the null space of A+B.

At a resonant lsq:
- A+B has a null space of dimension D (eigenvalues < THRESH)
- Each resonating zero γ_k (mode n_k near integer) "covers" one null direction
- When D zeros are aligned, the null space is fully covered → perfect cancellation
- When < D zeros are aligned, partial coverage → less deep resonance

Test by:
1. Computing dim(null space) at each deep resonance lsq
2. Counting aligned zeros (Δ < 0.08)
3. Testing if deep resonance depth correlates with coverage ratio

Also measure: WHICH null-space modes are "covered" (which modes n_k are near integer).
"""
import numpy as np
from scipy.special import digamma as _digamma
import scipy.linalg
import mpmath as mp

mp.mp.dps = 30
TWO_PI = 2.0 * np.pi
LOGPI  = float(np.log(np.pi))
G      = 50000
THRESH = 0.02   # null-space eigenvalue threshold
OUTFILE = '/Users/scones/projects/riemann/results/null_space_coverage.out'

GAMMAS_20 = [
    14.134725, 21.022040, 25.010858, 30.424876, 32.935062,
    37.586178, 40.918720, 43.327073, 48.005151, 49.773832,
    52.970321, 56.446247, 59.347044, 60.831779, 65.112544,
    67.079810, 69.546402, 72.067157, 75.704691, 77.144840,
]

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

def build_tau(lsq, N):
    """Build full τ matrix and return (τ, A+B, P, eigenvalues, eigenvectors)."""
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
    AB = arch + B

    if VM:
        lkpp = np.array([float(mp.log(k)) for k, _ in VM])
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

    tau = (AB - P); tau = (tau + tau.T) / 2.0
    ev, V = scipy.linalg.eigh(tau)
    return ev, V, AB, P, L, ns

def aligned_zeros(lsq, delta_thresh=0.08):
    L = float(mp.log(lsq))
    aligned = []
    for k, g in enumerate(GAMMAS_20):
        nk = g * L / TWO_PI
        f = nk % 1.0
        d = min(f, 1.0 - f)
        if d < delta_thresh:
            aligned.append((k+1, int(round(nk)), d, g))
    return aligned

def analyze(lsq, N, label=""):
    """Full null-space coverage analysis at given lsq, N."""
    ev, V, AB, P, L, ns = build_tau(lsq, N)
    lam_min = float(ev[0])

    # Null space of A+B
    ev_AB, V_AB = scipy.linalg.eigh((AB + AB.T)/2)
    null_AB_mask = np.abs(ev_AB) < THRESH
    D_null = int(null_AB_mask.sum())
    null_modes = ns[null_AB_mask.nonzero()[0]] if D_null > 0 else []

    # Null space of tau
    null_tau_mask = np.abs(ev) < THRESH
    D_null_tau = int(null_tau_mask.sum())

    # Which modes are resonant (near integer γ alignment)
    az = aligned_zeros(lsq, 0.08)

    # Mode coverage: for each aligned zero, find which mode n_k it aligns with
    # and check if that mode is in the null space of A+B
    coverage = []
    for k_idx, m, delta, gamma in az:
        # Is mode m in the null space of A+B?
        in_null = m in null_modes if len(null_modes) > 0 else False
        coverage.append((k_idx, m, delta, gamma, in_null))

    print(f"\n{'='*60}")
    print(f"lsq={lsq}, N={N}: {label}")
    print(f"{'='*60}")
    print(f"  lam_min(τ) = {lam_min:.4e}")
    print(f"  D_null(A+B) = {D_null} modes in null space (|ev|<{THRESH})")
    if len(null_modes) > 0:
        print(f"  Null modes n: {sorted(null_modes.tolist())}")
    print(f"  D_null(τ) = {D_null_tau} (near-zero τ eigenvalues)")
    print(f"  Aligned zeros (Δ<0.08, 20 zeros): {len(az)}")
    for k_idx, m, delta, gamma, in_null in coverage:
        flag = "IN NULL" if in_null else "not in null"
        print(f"    γ{k_idx}={gamma:.3f}: mode n={m}, Δ={delta:.4f} [{flag}]")

    coverage_score = sum(1 for _, _, _, _, in_null in coverage if in_null)
    print(f"  Coverage: {coverage_score}/{D_null} null directions covered by aligned zeros")
    if D_null > 0:
        print(f"  Coverage ratio: {coverage_score/D_null:.3f}")

    return {
        'lsq': lsq, 'N': N, 'lam_min': lam_min,
        'D_null_AB': D_null, 'D_null_tau': D_null_tau,
        'n_aligned': len(az), 'coverage': coverage_score,
        'ratio': coverage_score/D_null if D_null > 0 else 0
    }

def main():
    N = 150  # Smaller N for speed (full eigh needed)

    print("Null-space coverage analysis")
    print(f"N={N}, THRESH={THRESH}")
    print("HYPOTHESIS: deep resonance ⟺ (# aligned zeros) ≥ dim(null space of A+B)")

    # Test at key resonance points and contrast cases
    test_cases = [
        # (lsq, label, is_deep)
        (625, "γ₁₀(n=51) deep", True),
        (605, "γ₁₁(n=54,Δ=0.001) single NOT deep", False),  # High Δ single
        (715, "γ₂(n=22,Δ=0.011) deep", True),
        (710, "near γ₂ mode but NOT deep", False),
        (780, "γ₁+γ₁₅ dual deep", True),
        (800, "high score but NOT deep", False),
    ]

    results = []
    with open(OUTFILE, 'w') as f:
        f.write("Null-space coverage analysis\n\n")

    for lsq, label, is_deep in test_cases:
        r = analyze(lsq, N, label)
        r['is_deep'] = is_deep
        results.append(r)

        with open(OUTFILE, 'a') as f:
            f.write(f"lsq={lsq} ({'DEEP' if is_deep else 'not deep'}): "
                    f"D_null={r['D_null_AB']}, aligned={r['n_aligned']}, "
                    f"coverage={r['coverage']}/{r['D_null_AB']} "
                    f"(ratio={r['ratio']:.3f}), lam_min={r['lam_min']:.4e}\n")

    print(f"\n\nSUMMARY:")
    print(f"{'lsq':>5} {'deep?':>6} {'D_null':>7} {'aligned':>8} {'covered':>8} {'ratio':>7}")
    print("-" * 50)
    for r in results:
        print(f"{r['lsq']:5d} {'YES' if r['is_deep'] else 'NO':>6} {r['D_null_AB']:>7} "
              f"{r['n_aligned']:>8} {r['coverage']:>8} {r['ratio']:>7.3f}")

    with open(OUTFILE, 'a') as f:
        f.write("\nSUMMARY\n")
        f.write(f"{'lsq':>5} {'deep':>5} {'D_null':>7} {'aligned':>8} {'covered':>8} {'ratio':>7}\n")
        for r in results:
            f.write(f"{r['lsq']:5d} {'Y' if r['is_deep'] else 'N':>5} {r['D_null_AB']:>7} "
                    f"{r['n_aligned']:>8} {r['coverage']:>8} {r['ratio']:>7.3f}\n")

    print(f"\nResults written to {OUTFILE}")

if __name__ == "__main__":
    main()
