#!/usr/bin/env python3
"""ANGLE A: Investigate NEAR-NULL CLUSTER arithmetic identity. The near-null positive eigenvalues
(~[1e-8, 1e-4]) of whitened M_σ should be band-limited functions that VANISH at the window's ζ-zeros.
Test: (1) build M_σ, whiten, get eigendecomposition; (2) identify near-null cluster; (3) for each
near-null eigenvector, reconstruct the band-limited function φ(r)=Σ coef·P_center(r); (4) evaluate
|φ(γ)|² at the window's ζ-zeros vs total mass; (5) report counts, vanishing ratios, #near-null ≈ #zeros?"""
import numpy as np
from scipy.special import digamma
from scipy.linalg import eigh
import mpmath as mp

def vonmangoldt(N):
    lam = np.zeros(N+1); isc = np.zeros(N+1, dtype=bool)
    for p in range(2, N+1):
        if not isc[p]:
            isc[2*p::p] = True; lp = np.log(p); pk = p
            while pk <= N: lam[pk] = lp; pk *= p
    return lam

NMAX = 4000; LAM = vonmangoldt(NMAX); NS = np.arange(2, NMAX+1)
LOGN = np.log(NS); LAMW = LAM[2:]/np.sqrt(NS)

# Cache ζ-zeros via mpmath
mp.mp.dps = 20
ZT = []
k = 1
while True:
    g = float(mp.im(mp.zetazero(k))); ZT.append(g)
    if g > 360: break
    k += 1
ZT = np.array(ZT)

def zeros_in(a, b):
    return ZT[(ZT >= a) & (ZT <= b)]

def analyze_corner(a, b, fourW, frac=0.65, dr=0.004):
    """Build M_σ, whiten, eigendecompose; identify and analyze near-null cluster."""
    W = fourW/4
    L = b - a
    Ncap = L * 2 * W / np.pi
    N = max(4, int(round(frac * Ncap)))

    RGRID = np.arange(-(b+120), (b+120), dr)
    OMEGA = np.real(digamma(0.25 + 0.5j*RGRID)) - np.log(np.pi)
    CEN = np.linspace(a, b, N)

    def fej(x):
        y = W*x
        out = np.ones_like(y)
        nz = np.abs(y) > 1e-12
        out[nz] = (np.sin(y[nz])/y[nz])**2
        return out

    # Basis functions P_c(r) = fej(r-c) + fej(r+c)
    P = np.array([fej(RGRID-c) + fej(RGRID+c) for c in CEN])

    # Gram matrix
    G = (P @ P.T) * dr

    # Architecture term
    arch = (P * OMEGA[None, :]) @ P.T * dr / (2*np.pi)

    # Prime term
    mask = LOGN < fourW
    ln = LOGN[mask]
    lw = LAMW[mask]
    cosu = np.cos(np.outer(ln, RGRID))
    prime = np.zeros((N, N))
    for j in range(N):
        gp = ((P[j][None, :] * P) @ cosu.T) * dr / (2*np.pi)
        prime[j, :] = 2.0 * (gp * lw[None, :]).sum(axis=1)
    prime = 0.5 * (prime + prime.T)

    # Pole term
    def fc(c):
        y = W*(0.5j - c)
        y2 = W*(0.5j + c)
        return (np.sin(y)/y)**2 + (np.sin(y2)/y2)**2

    fp = np.array([fc(c) for c in CEN])
    pole = 2.0 * np.real(np.outer(fp, fp))

    # Whiten by G^{-1/2}
    gw, gv = eigh(G)
    B = gv @ np.diag(1/np.sqrt(gw)) @ gv.T

    # M_σ = arch - prime (whitened)
    Msig = arch - prime
    Msig_w = B @ Msig @ B

    # Eigendecompose whitened M_σ
    ev, V = eigh(Msig_w)

    # Identify near-null cluster (positive, small eigenvalues)
    # Threshold: eigenvalues in [1e-8, 1e-4] (after whitening)
    nn_mask = (ev > 1e-8) & (ev < 1e-4)
    nn_indices = np.where(nn_mask)[0]
    nn_eigenvalues = ev[nn_indices]
    nn_eigenvectors = V[:, nn_indices]  # columns of V (in whitened coordinates)

    # Also get the negative eigenvalue for context
    neg_mask = ev < 0
    neg_eigenvalues = ev[neg_mask]

    # Get window's ζ-zeros
    zeta_zeros = zeros_in(a, b)
    num_zeros = len(zeta_zeros)

    # For each near-null eigenvector (in whitened coords), reconstruct the band-limited function
    # in the ORIGINAL basis and evaluate at ζ-zeros
    # Unwhiten: eigenvector in original coords = G^{1/2} @ V_w[:, i]
    Gsqrt = gv @ np.diag(np.sqrt(gw)) @ gv.T

    vanishing_data = []
    for i, (eig_val, v_w) in enumerate(zip(nn_eigenvalues, nn_eigenvectors.T)):
        # Unwhiten to original basis
        v_orig = Gsqrt @ v_w

        # Reconstruct the band-limited function φ(r) = Σ v_orig[j] * P_j(r)
        # Evaluate at ζ-zeros and elsewhere
        phi_at_zeros = np.array([np.sum(v_orig * P[:, np.argmin(np.abs(RGRID - z))]) for z in zeta_zeros])

        # More careful: interpolate via P dot product
        # Actually, we have P as shape (N, len(RGRID)), so P[j, k] = P_j(RGRID[k])
        # φ(r) = Σ v_orig[j] P[j, idx(r)]
        # Better: evaluate φ at the grid points closest to the zeros
        phi_vals_full = (v_orig[:, None] * P).sum(axis=0)  # shape (len(RGRID),)

        # Find indices closest to zeros
        zero_indices = [np.argmin(np.abs(RGRID - z)) for z in zeta_zeros]
        phi_at_zeros = phi_vals_full[zero_indices]

        # Total mass of the function
        total_mass = np.sum(phi_vals_full**2) * dr

        # Mass at zeros (sum of squares at zero locations)
        mass_at_zeros = np.sum(phi_at_zeros**2) * (dr if len(zero_indices) > 0 else 1)

        # Vanishing ratio: how much of the mass concentrates at the zeros?
        vanishing_ratio = mass_at_zeros / total_mass if total_mass > 0 else 0

        vanishing_data.append({
            'eig_val': eig_val,
            'total_mass': total_mass,
            'mass_at_zeros': mass_at_zeros,
            'vanishing_ratio': vanishing_ratio,
            'phi_at_zeros': phi_at_zeros
        })

    return {
        'corner': (a, b),
        'fourW': fourW,
        'N': N,
        'Ncap': Ncap,
        'num_zeros': num_zeros,
        'num_near_null': len(nn_eigenvalues),
        'near_null_eigenvalues': nn_eigenvalues,
        'near_null_cluster_size': len(nn_eigenvalues),
        'vanishing_data': vanishing_data,
        'neg_eigenvalues': neg_eigenvalues,
        'zeta_zeros': zeta_zeros,
        'CEN': CEN,
        'P': P,
        'RGRID': RGRID,
        'V': V,
        'ev': ev
    }

# Run on corner configs from briefing
corners = [
    (14, 44, 4.0),
    (14, 34, 4.0),
    (14, 44, 5.0)
]

print("="*140)
print("ANGLE A: NEAR-NULL CLUSTER ARITHMETIC IDENTITY")
print("Test: band-limited near-null eigenvectors vanish at window ζ-zeros")
print("="*140)

all_results = []
for a, b, fourW in corners:
    result = analyze_corner(a, b, fourW)
    all_results.append(result)

    print(f"\n{'='*140}")
    print(f"CORNER: [{a}, {b}], 4W={fourW}")
    print(f"{'='*140}")
    print(f"Capacity N_cap = {result['Ncap']:.1f};  Actual N = {result['N']}")
    print(f"Number of ζ-zeros in window: {result['num_zeros']}")
    print(f"Number of NEGATIVE eigenvalues in M_σ: {len(result['neg_eigenvalues'])}")
    print(f"  Negative eigenvalues: {np.array2string(result['neg_eigenvalues'], precision=4)}")
    print(f"\nNEAR-NULL CLUSTER (eigenvalues in [1e-8, 1e-4]):")
    print(f"  Count: {result['num_near_null']}")
    if result['num_near_null'] > 0:
        print(f"  Eigenvalues: {np.array2string(result['near_null_eigenvalues'], precision=4)}")
        print(f"\nVANISHING ANALYSIS (does mass concentrate at zeros?):")
        print(f"  {'idx':>3} {'λ':>12} {'total_mass':>14} {'mass@zeros':>14} {'ratio':>10}")
        for i, data in enumerate(result['vanishing_data']):
            print(f"  {i:>3} {data['eig_val']:>12.4e} {data['total_mass']:>14.4e} {data['mass_at_zeros']:>14.4e} {data['vanishing_ratio']:>10.4f}")

        # Summary stats
        vanishing_ratios = np.array([d['vanishing_ratio'] for d in result['vanishing_data']])
        print(f"\n  VANISHING RATIO SUMMARY:")
        print(f"    Min: {vanishing_ratios.min():.4f}")
        print(f"    Mean: {vanishing_ratios.mean():.4f}")
        print(f"    Max: {vanishing_ratios.max():.4f}")
        print(f"\n  Hypothesis Test: Near-null modes vanish strongly at zeros?")
        strong_vanish = (vanishing_ratios > 0.1).sum()
        print(f"    {strong_vanish} / {len(vanishing_ratios)} modes with vanishing_ratio > 0.1 → ", end="")
        if strong_vanish / len(vanishing_ratios) > 0.5:
            print("LIKELY (>50%)")
        else:
            print("UNLIKELY (<50%)")
    else:
        print("  [No near-null modes detected in this corner]")

print(f"\n{'='*140}")
print("SUMMARY TABLE: #near-null vs #zeros")
print(f"{'='*140}")
print(f"{'corner':>12} {'4W':>4} {'#zeros':>7} {'#neg':>4} {'#near-null':>11} {'ratio':>8}")
for result in all_results:
    a, b = result['corner']
    ratio_str = f"{result['num_near_null']}/{result['num_zeros']}" if result['num_zeros'] > 0 else "—"
    print(f"{f'[{a},{b}]':>12} {result['fourW']:>4.1f} {result['num_zeros']:>7} {len(result['neg_eigenvalues']):>4} {result['num_near_null']:>11} {ratio_str:>8}")

print(f"{'='*140}")
print("HYPOTHESIS: If #near-null ≈ #zeros AND vanishing_ratios >> 0, then near-null cluster = band-limited")
print("functions concentrated at the window's ζ-zeros, confirming §855 (arithmetic identity in zero spacing).")
print(f"{'='*140}")
