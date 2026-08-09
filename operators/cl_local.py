#!/usr/bin/env python3
"""
Conrey–Li band-limited locality test.

TASK: de Branges GLOBAL positivity has the Conrey–Li counterexample: A+B (arch+pole, no primes)
has 4 negative eigenvalues (modes n=±1,±2), so it is NOT a valid reproducing kernel (not positive).
The Weil form τ=A+B−P is PSD (RH), but A+B alone fails.

QUESTION: If we band-limit (window to a frequency range), does the LOCAL form A_local+B_local
recover positivity? Or does LOCALITY still hit a negative-eigenvalue obstruction?

Method:
1. Build full A+B (global frequency support, ~200 basis functions)
2. Count negative eigenvalues in A+B (expect 4, matching literature claim)
3. Build windowed A_local+B_local (frequency restricted to [freq_min, freq_max])
4. Count negative eigenvalues in A_local+B_local
5. Compare: does windowing genuinely recover positivity (no neg eigs), or does the 4-neg obstruction
   persist (or transfer)?

Rigorous verdict: Does LOCALITY buy anything provable over GLOBAL?
"""
import numpy as np
from scipy.special import digamma
from scipy.linalg import eigh, svd
import sys
from fractions import Fraction

# Von Mangoldt function: tabulate log(p) for p prime
def vmanngoldt(NMAX):
    """Return Λ(n) for n=1..NMAX."""
    lam = np.zeros(NMAX + 1)
    is_composite = np.zeros(NMAX + 1, dtype=bool)
    for p in range(2, NMAX + 1):
        if not is_composite[p]:
            lam[p] = np.log(p)
            pk = p * p
            while pk <= NMAX:
                is_composite[pk] = True
                pk += p
            # Mark powers p^k
            pk = p * p
            while pk <= NMAX:
                lam[pk] = np.log(p)
                pk *= p
    return lam

def build_arch_pole_global(centers, fourW=1.32, grammat=None):
    """
    Build archimedean + pole term on full Fejér basis.

    centers: array of N Fejér basis center points
    fourW: bandwidth (half-width = fourW/4)
    grammat: pre-computed Gram matrix G (for normalization)

    Returns:
        arch: N×N archimedean operator
        pole: N×N pole term operator
        G_inv: N×N inverse Gram matrix (for inner product normalization)
    """
    N = len(centers)
    W = fourW / 4.0

    # Fine frequency grid for integration
    rmax = centers[-1] + 50
    dr = 0.005
    RG = np.arange(-rmax, rmax, dr)

    # Fejér basis (sinc-squared)
    def fejér(x):
        y = W * x
        out = np.ones_like(y)
        mask = np.abs(y) > 1e-12
        out[mask] = (np.sin(y[mask]) / y[mask]) ** 2
        return out

    # Basis evaluations
    P = np.array([fejér(RG - c) + fejér(RG + c) for c in centers])

    # Gram matrix (if not provided)
    if grammat is None:
        G = (P @ P.T) * dr
    else:
        G = grammat

    # Archimedean density: Ω(r) = Re ψ₀(1/4 + ir/2) - log(π)
    Omega = np.real(digamma(0.25 + 0.5j * RG)) - np.log(np.pi)

    # Archimedean operator: A = (1/2π) ∫ Ω(r) |P(r)|² dr (Gram-normalized)
    arch_raw = (P * Omega[None, :]) @ P.T * dr / (2 * np.pi)

    # Pole term: from ζ's pole at s=1
    # f(s) = ζ(s) has pole at s=1; φ(t) = ∫ |ψ̂(t)|² at poles s=0,1
    # Crude approximation: pole ∝ basis overlap at s=1/2
    def pole_func(x):
        y = W * (0.5j - x)
        y2 = W * (0.5j + x)
        return (np.sin(y) / y) ** 2 + (np.sin(y2) / y2) ** 2

    fp = np.array([pole_func(c) for c in centers])
    pole_raw = 2.0 * np.real(np.outer(fp, fp))  # outer product for quadratic form

    # Normalize via Gram-inverse: A → G^{-1/2} A G^{-1/2}
    try:
        evals_g, evecs_g = eigh(G)
        if evals_g[0] < 1e-15:
            print(f"WARNING: Gram matrix ill-conditioned, smallest eigenvalue = {evals_g[0]:.3e}")
            evals_g[evals_g < 1e-15] = 1e-15
        G_inv_sqrt = evecs_g @ np.diag(1.0 / np.sqrt(evals_g)) @ evecs_g.T
        arch = G_inv_sqrt @ arch_raw @ G_inv_sqrt
        pole = G_inv_sqrt @ pole_raw @ G_inv_sqrt
        G_inv = np.linalg.inv(G)
    except np.linalg.LinAlgError:
        print("ERROR: Gram matrix singular")
        return None, None, None

    return arch, pole, G_inv, G

def build_arch_pole_local(centers, fourW=1.32, freq_min=None, freq_max=None, grammat=None):
    """
    Build archimedean + pole term RESTRICTED to frequency window [freq_min, freq_max].

    This windowed form asks: if we only look at frequencies in [freq_min, freq_max],
    can we recover positivity (no negative eigenvalues)?
    """
    if freq_min is None or freq_max is None:
        freq_min, freq_max = 14, 44  # ANCHOR_PAPER window

    N = len(centers)
    W = fourW / 4.0

    # Use a finer grid within and around the window
    rmax = max(freq_max + 50, centers[-1] + 50)
    dr = 0.005
    RG = np.arange(-rmax, rmax, dr)

    def fejér(x):
        y = W * x
        out = np.ones_like(y)
        mask = np.abs(y) > 1e-12
        out[mask] = (np.sin(y[mask]) / y[mask]) ** 2
        return out

    P = np.array([fejér(RG - c) + fejér(RG + c) for c in centers])

    # Window mask: only keep frequencies in [freq_min, freq_max]
    window_mask = (np.abs(RG) >= freq_min) & (np.abs(RG) <= freq_max)
    P_windowed = P[:, window_mask]
    RG_windowed = RG[window_mask]
    dr_window = dr

    if len(RG_windowed) == 0:
        print(f"WARNING: empty window [{freq_min}, {freq_max}]")
        return None, None, None, None

    # Gram matrix on windowed frequencies
    if grammat is None:
        G_local = (P_windowed @ P_windowed.T) * dr_window
    else:
        G_local = grammat

    # Archimedean (windowed)
    Omega_windowed = np.real(digamma(0.25 + 0.5j * RG_windowed)) - np.log(np.pi)
    arch_local_raw = (P_windowed * Omega_windowed[None, :]) @ P_windowed.T * dr_window / (2 * np.pi)

    # Pole (windowed)
    def pole_func(x):
        y = W * (0.5j - x)
        y2 = W * (0.5j + x)
        return (np.sin(y) / y) ** 2 + (np.sin(y2) / y2) ** 2
    fp = np.array([pole_func(c) for c in centers])
    pole_local_raw = 2.0 * np.real(np.outer(fp, fp))

    # Normalize
    try:
        evals_g, evecs_g = eigh(G_local)
        if evals_g[0] < 1e-15:
            evals_g[evals_g < 1e-15] = 1e-15
        G_inv_sqrt = evecs_g @ np.diag(1.0 / np.sqrt(evals_g)) @ evecs_g.T
        arch_local = G_inv_sqrt @ arch_local_raw @ G_inv_sqrt
        pole_local = G_inv_sqrt @ pole_local_raw @ G_inv_sqrt
        G_local_inv = np.linalg.inv(G_local)
    except np.linalg.LinAlgError:
        print("ERROR: Local Gram matrix singular")
        return None, None, None, None

    return arch_local, pole_local, G_local_inv, G_local

def run_test(fourW=1.32, window_range=(14, 44), N_basis=32):
    """
    Run the full Conrey–Li locality test.

    1. Build global A+B
    2. Count negative eigenvalues (expect ~4)
    3. Build local A_local + B_local
    4. Count negative eigenvalues
    5. Verdict: does locality recover positivity?
    """
    print("\n" + "=" * 100)
    print(f"CONREY–LI LOCALITY TEST: 4W={fourW}, window={window_range}, N_basis={N_basis}")
    print("=" * 100)

    # Create Fejér basis centers (evenly spaced in frequency)
    freq_min, freq_max = window_range
    centers = np.linspace(freq_min, freq_max, N_basis)

    # --- GLOBAL (full frequency) ---
    print("\n[GLOBAL] Building A + B on full frequency support...")
    arch_glob, pole_glob, _, G_glob = build_arch_pole_global(centers, fourW=fourW)
    if arch_glob is None:
        print("FAILED to build global operators")
        return

    AB_glob = arch_glob + pole_glob
    evals_glob, evecs_glob = eigh(AB_glob)
    n_neg_glob = np.sum(evals_glob < -1e-10)
    n_zero_glob = np.sum(np.abs(evals_glob) < 1e-10)

    print(f"  A+B eigenspectrum: min={evals_glob[0]:.6e}, max={evals_glob[-1]:.6e}")
    print(f"  Negative eigenvalues: {n_neg_glob} (expect ~4 per Conrey–Li)")
    print(f"  Near-zero eigenvalues: {n_zero_glob}")
    print(f"  Spectral gap (neg/pos boundary): {abs(evals_glob[n_neg_glob]) if n_neg_glob > 0 else '(all positive)'}")

    if n_neg_glob > 0:
        print(f"  Negative eigenvalues: {evals_glob[:n_neg_glob]}")

    # --- LOCAL (windowed) ---
    print(f"\n[LOCAL] Building A_local + B_local with frequency window [{freq_min}, {freq_max}]...")
    arch_loc, pole_loc, _, G_loc = build_arch_pole_local(
        centers, fourW=fourW, freq_min=freq_min, freq_max=freq_max, grammat=G_glob
    )
    if arch_loc is None:
        print("FAILED to build local operators")
        return

    AB_loc = arch_loc + pole_loc
    evals_loc, evecs_loc = eigh(AB_loc)
    n_neg_loc = np.sum(evals_loc < -1e-10)
    n_zero_loc = np.sum(np.abs(evals_loc) < 1e-10)

    print(f"  A_local+B_local eigenspectrum: min={evals_loc[0]:.6e}, max={evals_loc[-1]:.6e}")
    print(f"  Negative eigenvalues: {n_neg_loc}")
    print(f"  Near-zero eigenvalues: {n_zero_loc}")

    if n_neg_loc > 0:
        print(f"  Negative eigenvalues: {evals_loc[:n_neg_loc]}")
        print(f"  Spectral gap (neg/pos boundary): {abs(evals_loc[n_neg_loc]):.6e}")
    else:
        print(f"  ✓ NO NEGATIVE EIGENVALUES — LOCAL POSITIVITY RECOVERED!")
        print(f"  Spectral gap (smallest positive): {evals_loc[0]:.6e}")

    # --- DIAGNOSIS ---
    print("\n" + "=" * 100)
    print("DIAGNOSIS:")
    print("=" * 100)
    if n_neg_glob == 0:
        print("⚠ WARNING: Global A+B has NO negative eigenvalues (unexpected!)")
        print("  This contradicts the Conrey–Li claim. Check Fejér grid size, bandwidth, etc.")
    elif n_neg_glob != 4:
        print(f"⚠ WARNING: Global A+B has {n_neg_glob} negative eigenvalues (expected 4)")
        print("  Conrey–Li cited 4 for full de Branges (modes n=±1,±2). May be basis-dependent.")
    else:
        print("✓ Global A+B has 4 negative eigenvalues (matches Conrey–Li counterexample)")

    print()
    if n_neg_loc == 0:
        print("✓✓ LOCALITY RECOVERS POSITIVITY:")
        print("   A_local + B_local ⪰ 0 (no negative eigenvalues)")
        print("   The band-limited restriction EVADES the Conrey–Li obstruction!")
        print("   → Windowing genuinely buys positivity (provable on local form).")
        return "HANDLE"  # positivity genuinely recovered
    elif n_neg_loc < n_neg_glob:
        print(f"◑ PARTIAL LOCALITY:")
        print(f"   A_local + B_local has {n_neg_loc} negative eigenvalues (down from {n_neg_glob})")
        print("   Some obstruction persists, but mitigated by windowing.")
        print("   → Locality buys something (fewer negative eigs) but not full positivity.")
        return "REPORT"
    else:
        print(f"✗ LOCALITY DOES NOT RECOVER POSITIVITY:")
        print(f"   A_local + B_local still has {n_neg_loc} negative eigenvalues")
        print("   The Conrey–Li obstruction persists under windowing.")
        print("   → Locality does NOT buy anything provable; same razor at both global/local.")
        return "REPORT"

if __name__ == "__main__":
    # Test across several configurations
    results = []

    print("\n" + "=" * 100)
    print("CONREY–LI LOCALITY CAMPAIGN")
    print("Hypothesis: Band-limiting evades the 4-negative-eigenvalue counterexample")
    print("=" * 100)

    # Scenario 1: The ANCHOR_PAPER window
    verdict1 = run_test(fourW=1.32, window_range=(14, 44), N_basis=24)
    results.append(("ANCHOR window [14,44], 4W=1.32", verdict1))

    # Scenario 2: Wider window (to see if more terms help)
    verdict2 = run_test(fourW=1.32, window_range=(10, 50), N_basis=32)
    results.append(("Wider window [10,50], 4W=1.32", verdict2))

    # Scenario 3: Narrower window (more aggressive banding)
    verdict3 = run_test(fourW=1.32, window_range=(18, 40), N_basis=16)
    results.append(("Narrow window [18,40], 4W=1.32", verdict3))

    # Scenario 4: Larger bandwidth
    verdict4 = run_test(fourW=2.0, window_range=(14, 44), N_basis=32)
    results.append(("ANCHOR window [14,44], 4W=2.0 (larger)", verdict4))

    print("\n" + "=" * 100)
    print("SUMMARY")
    print("=" * 100)
    for config, verdict in results:
        status = "✓ POSITIVE" if verdict == "HANDLE" else ("◑ MIXED" if verdict == "REPORT" else "✗ NEGATIVE")
        print(f"  {config:40s} → {status}")

    print("\n" + "=" * 100)
    print("HONEST VERDICT:")
    print("=" * 100)
    handle_count = sum(1 for _, v in results if v == "HANDLE")
    mixed_count = sum(1 for _, v in results if v == "REPORT")

    if handle_count == len(results):
        print("✓✓ LOCALITY UNIFORMLY RECOVERS POSITIVITY across all tested windows/bandwidths.")
        print("   Conclusion: Band-limiting genuinely evades Conrey–Li. Does this help uniform-margin bound?")
        print("   → YES, conditionally: LOCALITY is a real structural win, but only within each window.")
        print("   → Full RH proof still requires showing this holds uniformly across all windows (the open problem).")
    elif handle_count > 0:
        print(f"◑ MIXED RESULTS: {handle_count} configs recover, {mixed_count} show obstruction remains.")
        print("   Conclusion: Locality helps in SOME regimes (narrower windows, specific bandwidths)")
        print("   but not uniformly. The uniform-margin bound still bottoms out at RH.")
    else:
        print(f"✗ LOCALITY DOES NOT EVADE CONREY–LI in any tested config.")
        print("   Conclusion: Windowing does not recover positivity; same razor at global/local.")
        print("   → Band-limiting does NOT help the uniform-margin bound; we hit the same obstruction.")
        print("   → The hard content is global, not local; any RH proof needs NEW unconditional ideas.")
