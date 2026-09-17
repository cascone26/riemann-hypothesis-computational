#!/usr/bin/env python3
"""variational_setup_v2.py -- Corrected decomposition of the Weil form.

CORRECTED LENS: The Weil form Q integrates with a MINUS sign (integration by parts):
  Q = ∫∫ (−g'') v_k v_l

So the arch pieces (Qcore, Qtail) from build_Q carry the minus sign. To extract the
ARCH and PRIME quadratic forms for the generalized eigenvalue problem, we need to
FLIP THE SIGNS:

  arch = −(Qcore + Qtail)  [positive-definite, per §965]
  prime = −Qprime          [flip the minus to get positive]

Then Q = Qpole − arch − prime = Qpole − (arch + prime), or equivalently:
  prime + arch = −Q + Qpole

The generalized eigenvalue λ of (prime, arch) tells us the ratio by which prime
exceeds arch. RH ⟺ κ = #{λ > 1} ≤ 1.
"""

import numpy as np
from scipy.linalg import eigh
import sys
sys.path.insert(0, '/Users/scones/projects/riemann/operators')
from suzuki_exact_Q import build_Q

def extremal_analysis_corrected(a, K=80, NG=320, verbose=True):
    """Variational analysis with CORRECTED arch/prime sign decomposition."""

    Q, pieces = build_Q(a, K, NG, return_pieces=True)

    # CORRECTED DECOMPOSITION (flip signs due to integration by parts)
    Qcore = pieces["core"]     # −2c·a·I  (negative diagonal)
    Qtail = pieces["tail"]     # positive from quadrature
    Qprime_signed = pieces["prime"]  # comes with minus sign in the definition
    Qpole = pieces["pole"]

    # Flip to get the positive-definite versions
    arch = -(Qcore + Qtail)     # This should be PSD
    prime = -Qprime_signed      # This should be PSD

    # Verify: Q should equal Qpole - arch - prime
    Q_check = Qpole - arch - prime
    Q_error = np.linalg.norm(Q - Q_check, 'fro') / (np.linalg.norm(Q, 'fro') + 1e-14)

    if verbose:
        print(f"\n{'='*80}")
        print(f"WINDOW [{-a:+.2f}, {a:+.2f}] K={K}, NG={NG}, dim={K}")
        print(f"{'='*80}")
        print(f"Decomposition check |Q - (Qpole - arch - prime)|_F / |Q|_F = {Q_error:.2e}")

    # Verify arch is PSD
    arch_eigs = np.linalg.eigvalsh(arch)
    arch_neg = np.sum(arch_eigs < -1e-13)
    arch_pos = np.sum(arch_eigs > 1e-13)

    if verbose:
        print(f"\nARCH inertia (should be PSD per §965):")
        print(f"  neg={arch_neg}, pos={arch_pos}")
        print(f"  min eig={arch_eigs[0]:.3e}, max eig={arch_eigs[-1]:.3e}")
        if arch_neg > 0:
            print(f"  ⚠ ARCH has {arch_neg} negative eigenvalues — structure incomplete?")

    # Verify prime is PSD
    prime_eigs = np.linalg.eigvalsh(prime)
    prime_neg = np.sum(prime_eigs < -1e-13)

    if verbose:
        print(f"\nPRIME inertia (should be PSD):")
        print(f"  neg={prime_neg}, pos={np.sum(prime_eigs > 1e-13)}")
        print(f"  min eig={prime_eigs[0]:.3e}, max eig={prime_eigs[-1]:.3e}")

    # Generalized spectrum: prime x = λ arch x
    # Use scipy's eigh to solve the generalized problem
    try:
        lam, X = eigh(prime, arch)
    except np.linalg.LinAlgError:
        if verbose:
            print("⚠ Generalized eigenvalue solve failed — arch may be singular")
        return None

    # Count λ > 1
    kappa = np.sum(lam > 1.0 + 1e-10)

    if verbose:
        print(f"\nGENERALIZED SPECTRUM (prime x = λ arch x):")
        print(f"  κ = #{'{λ > 1}'} = {kappa}")
        print(f"  λ_max = {lam[-1]:.6e}  (margin = {lam[-1]-1:.6e})")
        if len(lam) >= 2:
            print(f"  λ_2   = {lam[-2]:.6e}  (margin to 1 = {1-lam[-2]:.6e})")
        if len(lam) >= 3:
            print(f"  λ_3   = {lam[-3]:.6e}")
        print(f"  min λ = {lam[0]:.6e}")

    # Pole structure (should be rank-1)
    pole_eigs = np.linalg.eigvalsh(Qpole)
    pole_rank_1_eig = pole_eigs[-1]
    U_pole = np.linalg.eigh(Qpole)[1][:, -1]

    if verbose:
        print(f"\nPOLE structure:")
        print(f"  Top eigenvalue = {pole_rank_1_eig:.6e}")
        print(f"  Other top eigs = {pole_eigs[-5:]}")

    # Containment test (§976): arch−prime restricted to pole⊥ should be PSD
    M = arch - prime  # Since both arch and prime are now positive, this is the difference
    u_norm = U_pole / (np.linalg.norm(U_pole) + 1e-14)
    P_perp = np.eye(K) - np.outer(u_norm, u_norm)
    M_perp = P_perp @ M @ P_perp
    eigs_M_perp = np.linalg.eigvalsh(M_perp)

    if verbose:
        print(f"\nGEOMETRIC CONTAINMENT (§976): M|pole⊥ with M = arch − prime")
        print(f"  min eig(M|pole⊥) = {eigs_M_perp[0]:.3e}  (should be ≥0 for RH)")
        print(f"  κ on pole⊥ = {np.sum(eigs_M_perp < -1e-13)}")

    # Extremal analysis
    x_top = X[:, -1]
    x_2nd = X[:, -2] if len(lam) >= 2 else None

    align_top = abs(np.dot(U_pole, x_top))
    align_2nd = abs(np.dot(U_pole, x_2nd)) if x_2nd is not None else np.nan

    if verbose:
        print(f"\nEXTREMAL alignment:")
        print(f"  |⟨pole, x_top⟩| = {align_top:.6f}")
        if x_2nd is not None:
            print(f"  |⟨pole, x_2nd⟩| = {align_2nd:.6f}")

    return {
        'a': a, 'K': K, 'NG': NG,
        'lambda': lam, 'X': X, 'kappa': kappa,
        'arch': arch, 'prime': prime,
        'U_pole': U_pole,
        'x_top': x_top, 'x_2nd': x_2nd,
        'align_top': align_top, 'align_2nd': align_2nd,
        'M_perp_min': eigs_M_perp[0],
    }


def bandwidth_sweep_second_eigenvalue(a_list=[0.5, 1.0, 1.5, 2.0], K=80):
    """Track (1 - λ_2) across windows."""
    print("\n" + "="*90)
    print("BANDWIDTH SWEEP: λ_2 margin 1−λ_2 across window size a")
    print("="*90)
    print(f"{'a':>6} {'λ_max':>12} {'λ_max−1':>12} {'λ_2':>12} {'1−λ_2':>12} {'κ':>4}")
    print("-"*70)

    results = []
    for a in a_list:
        res = extremal_analysis_corrected(a, K=K, NG=320, verbose=False)
        if res is None:
            print(f"{a:6.2f} --- generalized eigenvalue solve failed")
            continue
        lam = res['lambda']
        lam_2 = lam[-2] if len(lam) >= 2 else np.nan
        print(f"{a:6.2f} {lam[-1]:12.6e} {lam[-1]-1:12.6e} {lam_2:12.6e} "
              f"{1-lam_2:12.6e} {res['kappa']:4d}")
        results.append(res)

    return results


def analyze_second_eigenvalue(a=1.5, K=80):
    """Detailed structure of λ_2."""
    res = extremal_analysis_corrected(a, K=K, NG=320, verbose=True)

    if res is None:
        print("Analysis failed.")
        return res

    print(f"\n{'='*80}")
    print(f"λ_2 DETAILED STRUCTURE")
    print(f"{'='*80}")

    lam = res['lambda']
    x_2nd = res['x_2nd']
    U_pole = res['U_pole']

    if x_2nd is None:
        print("Only one eigenvalue computed.")
        return res

    # Project x_2nd to pole⊥
    u_norm = U_pole / np.linalg.norm(U_pole)
    x_2nd_perp = x_2nd - np.dot(u_norm, x_2nd) * u_norm

    print(f"Alignment of x_2nd with pole direction: {abs(np.dot(U_pole, x_2nd)):.6f}")
    print(f"Norm of x_2nd projected to pole⊥: {np.linalg.norm(x_2nd_perp):.6e}")

    # Spectrum
    print(f"\nGeneralized spectrum tail:")
    for i in range(max(0, len(lam)-8), len(lam)):
        marker = " ← λ_top" if i == len(lam)-1 else (" ← λ_2nd" if i == len(lam)-2 else "")
        print(f"  λ[{i}] = {lam[i]:12.6e} {marker}")

    return res


if __name__ == "__main__":
    print("\n" + "="*90)
    print("VARIATIONAL ANALYSIS — Corrected Weil decomposition")
    print("="*90)

    # Detailed analysis
    res_detailed = analyze_second_eigenvalue(a=1.5, K=80)

    # Bandwidth sweep
    results = bandwidth_sweep_second_eigenvalue(a_list=[0.5, 0.8, 1.2, 1.5, 2.0], K=80)

    print("\n✓ Analysis complete.")
