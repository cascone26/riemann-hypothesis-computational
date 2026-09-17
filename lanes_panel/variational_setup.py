#!/usr/bin/env python3
"""variational_setup.py -- Rayleigh quotient & extremal test function analysis.

LENS: Variational / extremal via max_g prime(g,g) / arch(g,g).

GOAL: (1) identify the structure of the top generalized eigenvalue λ_max (the pole),
       (2) analyze the SECOND eigenvalue and its margin 1−λ_2,
       (3) probe whether unconditional arguments can bound λ_2≤1 from the structure alone.

The reframing from §980: κ = #{λ > 1} is a SIGN COUNT on the generalized spectrum.
Key fact: prime ≡ arch to machine precision on pole⊥ (all λ → 1 from below).
The question: can we show unconditionally that λ_2 is the one on pole⊥, not somewhere
else, and that it stays ≤1?
"""

import numpy as np
from scipy.linalg import eigh
import sys
sys.path.insert(0, '/Users/scones/projects/riemann/operators')
from suzuki_exact_Q import build_Q

def generalized_spectrum_full(prime, arch):
    """Compute generalized eigenvalues λ of prime x = λ arch x.

    Returns (λ sorted, eigenvectors in columns of X).
    """
    lam, X = eigh(prime, arch)
    return lam, X

def rayleigh_quotient(g, prime, arch):
    """Compute R(g) = g^T prime g / g^T arch g."""
    num = g @ prime @ g
    den = g @ arch @ g
    if abs(den) < 1e-14:
        return np.nan
    return num / den

def extremal_analysis(a, K=80, NG=320, verbose=True):
    """
    Full variational analysis on window [−a, a].

    Returns dict with:
      - λ: all generalized eigenvalues (sorted)
      - X: eigenvectors (columns)
      - pole_dir: principal direction of pole (rank-1)
      - margins: (λ_max − 1, 1 − λ_2, λ_2 − λ_3, ...)
      - diagnostics: various structural checks
    """
    # Build exact Q = arch + pole - prime
    Q, pieces = build_Q(a, K, NG, return_pieces=True)
    Qpole = pieces["pole"]
    Qcore = pieces["core"]
    Qtail = pieces["tail"]
    Qprime = pieces["prime"]

    # Decompose: arch = Qcore + Qtail (PSD), prime = −Qprime
    arch = Qcore + Qtail
    prime = -Qprime  # Qprime has a minus sign in front

    # Verify: Q = arch + Qpole - prime = 0 (the Weil form)
    Q_check = arch + Qpole - prime
    Q_error = np.linalg.norm(Q - Q_check, 'fro') / np.linalg.norm(Q, 'fro')

    if verbose:
        print(f"\n{'='*80}")
        print(f"WINDOW [{-a:+.2f}, {a:+.2f}] K={K}, NG={NG}, dim={K}")
        print(f"{'='*80}")
        print(f"Decomposition check |Q - (arch + pole - prime)|_F / |Q|_F = {Q_error:.2e}")

    # Compute arch inertia (should be PD by §965)
    arch_eigs = np.linalg.eigvalsh(arch)
    arch_neg = np.sum(arch_eigs < -1e-13)
    arch_pos = np.sum(arch_eigs > 1e-13)

    if verbose:
        print(f"ARCH inertia: neg={arch_neg}, pos={arch_pos}, "
              f"min eig={arch_eigs[0]:.3e}, max eig={arch_eigs[-1]:.3e}")
        assert arch_neg == 0, "arch should be PSD (§965)"

    # Generalized spectrum: prime x = λ arch x
    lam, X = generalized_spectrum_full(prime, arch)

    # Count λ > 1
    kappa = np.sum(lam > 1.0 + 1e-10)

    # Margins
    margins = {
        'λ_max': lam[-1],
        'λ_max - 1': lam[-1] - 1.0,
        'λ_2': lam[-2] if len(lam) >= 2 else np.nan,
        '1 - λ_2': 1.0 - lam[-2] if len(lam) >= 2 else np.nan,
        'λ_2 - λ_3': lam[-2] - lam[-3] if len(lam) >= 3 else np.nan,
    }

    if verbose:
        print(f"\nGENERALIZED SPECTRUM (prime x = λ arch x):")
        print(f"  κ = #{'{λ > 1}'} = {kappa}")
        print(f"  λ_max = {lam[-1]:.6e}  (margin = {margins['λ_max - 1']:.6e})")
        print(f"  λ_2 = {lam[-2]:.6e}  (1 - λ_2 = {margins['1 - λ_2']:.6e})")
        print(f"  λ_3 = {lam[-3]:.6e}  (λ_2 - λ_3 = {margins['λ_2 - λ_3']:.6e})")
        print(f"  min(λ) = {lam[0]:.6e}")

    # Pole direction: principal direction of Qpole (should be rank-1)
    pole_eigs = np.linalg.eigvalsh(Qpole)
    pole_rank_1_eig = pole_eigs[-1]  # top eigenvalue
    U_pole = np.linalg.eigh(Qpole)[1][:, -1]  # corresponding eigenvector

    if verbose:
        print(f"\nPOLE structure (Qpole rank-1 approximation):")
        print(f"  λ_pole (top eig) = {pole_rank_1_eig:.6e}")
        print(f"  Other pole eigs = {pole_eigs[-5:]}")

    # Extremal analysis
    # The top generalized eigenvector
    x_top = X[:, -1]  # eigenvector for λ_max

    # The second generalized eigenvector
    x_2nd = X[:, -2] if len(lam) >= 2 else None

    # Alignment: pole direction vs top & second eigenvectors
    align_pole_top = abs(np.dot(U_pole, x_top))
    align_pole_2nd = abs(np.dot(U_pole, x_2nd)) if x_2nd is not None else np.nan

    if verbose:
        print(f"\nEXTREMAL alignment:")
        print(f"  |⟨pole direction, x_top⟩| = {align_pole_top:.6f}")
        print(f"  |⟨pole direction, x_2nd⟩| = {align_pole_2nd:.6f}")

    # Key test: project arch−prime to pole⊥
    # If M|pole⊥ ⪰ 0 (§976), then all generalized eigenvalues on that subspace should be ≤1
    M = arch - prime  # = Qcore + Qtail + Qpole (the negatives of −Qprime alone)

    # Project to pole⊥: P⊥ = I - uu^T where u = U_pole / |U_pole|
    u_norm = U_pole / np.linalg.norm(U_pole)
    P_perp = np.eye(K) - np.outer(u_norm, u_norm)

    M_perp = P_perp @ M @ P_perp
    eigs_M_perp = np.linalg.eigvalsh(M_perp)

    if verbose:
        print(f"\nGEOMETRIC CONTAINMENT (§976): M|pole⊥ = arch−prime projected to pole⊥")
        print(f"  min eig(M|pole⊥) = {eigs_M_perp[0]:.3e}  (should be ≥0 for RH)")
        print(f"  κ on pole⊥ = {np.sum(eigs_M_perp < -1e-13)}")

    # Rayleigh quotient check
    R_top = rayleigh_quotient(x_top, prime, arch)
    R_2nd = rayleigh_quotient(x_2nd, prime, arch) if x_2nd is not None else np.nan

    if verbose:
        print(f"\nRAYLEIGH VERIFICATION:")
        print(f"  R(x_top) = {R_top:.6e}  (should = λ_max = {lam[-1]:.6e})")
        print(f"  R(x_2nd) = {R_2nd:.6e}  (should = λ_2 = {lam[-2]:.6e})")

    return {
        'a': a,
        'K': K,
        'NG': NG,
        'lambda': lam,
        'X': X,
        'kappa': kappa,
        'margins': margins,
        'arch': arch,
        'prime': prime,
        'pole_eigs': pole_eigs,
        'U_pole': U_pole,
        'pole_rank_1_eig': pole_rank_1_eig,
        'x_top': x_top,
        'x_2nd': x_2nd,
        'align_pole_top': align_pole_top,
        'align_pole_2nd': align_pole_2nd,
        'M_perp_min': eigs_M_perp[0],
        'M_perp': M_perp,
        'eigs_M_perp': eigs_M_perp,
        'Qpole': Qpole,
        'Qcore': Qcore,
        'Qtail': Qtail,
        'Qprime': Qprime,
    }


def bandwidth_scaling_second_eigenvalue(a_list=[0.5, 1.0, 1.5, 2.0], K=80, NG=320):
    """Sweep window size and track λ_2 margin."""
    print("\n" + "="*90)
    print("BANDWIDTH SWEEP: how does (1 - λ_2) scale with window a?")
    print("="*90)

    results = []
    print(f"{'a':>6} {'λ_max':>12} {'λ_max-1':>12} {'λ_2':>12} {'1-λ_2':>12} {'κ':>4}")
    print("-" * 70)

    for a in a_list:
        res = extremal_analysis(a, K=K, NG=NG, verbose=False)
        lam = res['lambda']
        print(f"{a:6.2f} {lam[-1]:12.6e} {lam[-1]-1:12.6e} {lam[-2]:12.6e} "
              f"{1-lam[-2]:12.6e} {res['kappa']:4d}")
        results.append(res)

    return results


def second_eigenvalue_structure(a=1.5, K=80, NG=320):
    """Detailed analysis of λ_2: is it on pole⊥? How does it behave?"""

    res = extremal_analysis(a, K=K, NG=NG, verbose=True)

    lam = res['lambda']
    X = res['X']
    U_pole = res['U_pole']

    print(f"\n{'='*80}")
    print(f"DETAILED STRUCTURE OF λ_2 (second generalized eigenvalue)")
    print(f"{'='*80}")

    # The second eigenvector
    x_2nd = X[:, -2]

    # Is it on pole⊥?
    align_with_pole = abs(np.dot(U_pole, x_2nd))
    print(f"\nAlignment with pole direction:")
    print(f"  |⟨pole, x_2nd⟩| = {align_with_pole:.6f}  ({100*align_with_pole:.2f}%)")

    # The Rayleigh quotient on pole⊥
    # Project x_2nd to pole⊥
    u_norm = U_pole / np.linalg.norm(U_pole)
    x_2nd_perp = x_2nd - np.dot(u_norm, x_2nd) * u_norm

    if np.linalg.norm(x_2nd_perp) > 1e-10:
        x_2nd_perp_normalized = x_2nd_perp / np.linalg.norm(x_2nd_perp)
        R_2nd_perp = rayleigh_quotient(x_2nd_perp_normalized, res['prime'], res['arch'])
        print(f"\nRayleigh quotient on pole⊥:")
        print(f"  R(x_2nd projected to pole⊥) = {R_2nd_perp:.6e}  (expect ≤1)")

    # Spectrum visualization
    print(f"\nFull generalized spectrum (λ sorted):")
    for i in range(max(0, len(lam)-10), len(lam)):
        marker = "← λ_top (pole)" if i == len(lam)-1 else ("← λ_2nd" if i == len(lam)-2 else "")
        print(f"  λ[{i}] = {lam[i]:12.6e}  {marker}")

    return res


if __name__ == "__main__":
    # Main test: detailed analysis at a=1.5
    print("\n" + "="*90)
    print("VARIATIONAL ANALYSIS — Rayleigh quotient & extremal structure")
    print("="*90)

    res = second_eigenvalue_structure(a=1.5, K=80, NG=320)

    # Bandwidth sweep
    bandwidth_scaling_second_eigenvalue(a_list=[0.5, 0.8, 1.2, 1.5, 2.0], K=80, NG=320)

    print("\n✓ Analysis complete. Output in lanes_panel/variational_setup.py.")
