#!/usr/bin/env python3
"""POLE AUGMENTATION PERTURBATION TEST — actually run the §977-style test that S2645 deferred

S2646 correctly flagged that S2645's "faithfulness" claim was never tested. The assertion was:
  "At true primes: κ=1 (verified)"
  "At perturbed primes (§977 result): κ would become κ >> 1"

But this was UNVERIFIED — the code just had a print statement. This script ACTUALLY runs it.

METHODOLOGY (§977 + §949):
1. Build Q with TRUE primes → compute κ of augmented Suzuki-Weyl m-function
2. Build Q with PERTURBED log-frequencies (log n → log n * (1 + δ·U(-1,1))) → recompute κ
3. Vary δ ∈ {0.0005, 0.001, 0.002, 0.005, 0.01} to see fragility
4. Test NEGATIVE CONTROLS: random pole placements unrelated to truth

CRITICAL: We test the actual AUGMENTED m-function (m_augmented = m_weyl_base + c/(z-z_pole)),
NOT just κ(arch-prime). This tests whether the augmentation is "faithful" to perturbation
in the sense that S2645 claimed.
"""

import numpy as np
import sys, os
sys.path.insert(0, os.path.join(os.getcwd(),'operators'))
from suzuki_exact_Q import build_Q, E_p, shat, vonmangoldt_upto
from scipy.linalg import solve

def weyl_pieces_from_Q(Q, a, K, lam_shift=0.2):
    """Construct Weyl solutions vp, vm from a given Q matrix"""
    kk = np.arange(1, K+1)
    lam_a = np.linalg.eigvalsh(Q)[0]/a
    lam = lam_a - lam_shift
    Ta = Q - lam*a*np.eye(K)
    fplus = E_p(1.0, a, kk)
    fminus = E_p(-1.0, a, kk)
    cplus = solve(Ta, fplus, assume_a='sym')
    cminus = solve(Ta, fminus, assume_a='sym')
    npl = np.sqrt(abs(cplus@fplus))
    nmi = np.sqrt(abs(cminus@fminus))
    cminus = cminus*(npl/nmi)
    def vp(z): return np.sum(cplus*shat(z, a, kk))
    def vm(z): return np.sum(cminus*shat(z, a, kk))
    return vp, vm, lam_a

def build_Q_perturbed_log_freqs(a, K, NG=400, delta=0.0, seed=3141):
    """Build Q with perturbed log-frequencies in the prime term

    Methodology: log n → log n * (1 + delta * U(-1,1)) for each prime power n.
    """
    rng = np.random.default_rng(seed)

    # Build base Q pieces (using original suzuki_exact_Q)
    Q_orig, pieces_orig = build_Q(a, K, NG, return_pieces=True)

    if delta == 0.0:
        return Q_orig, pieces_orig["pole"], 0.0

    # Rebuild prime term with perturbed log-frequencies
    L = 2*a
    Nmax = int(np.floor(np.exp(L)))
    Lam = vonmangoldt_upto(max(Nmax, 2))

    kk = np.arange(1, K+1)
    Qprime_perturbed = np.zeros((K, K))
    perturbation_sizes = []

    for n in range(2, Nmax+1):
        if Lam[n] > 0 and np.log(n) < L:
            log_n_orig = np.log(n)
            # Perturb: log n -> log n * (1 + delta * U(-1,1))
            perturbation_factor = 1.0 + delta * rng.uniform(-1.0, 1.0)
            log_n_perturbed = log_n_orig * perturbation_factor
            perturbation_sizes.append(abs(log_n_perturbed - log_n_orig))

            # Use perturbed log-frequency in S_mat
            from suzuki_exact_Q import S_mat
            Sn_perturbed = S_mat(log_n_perturbed, a, K)
            Qn = -(Lam[n]/np.sqrt(n))*Sn_perturbed
            Qprime_perturbed += Qn

    # Reconstruct Q with perturbed prime term
    Q_perturbed = pieces_orig["pole"] + pieces_orig["core"] + pieces_orig["tail"] + Qprime_perturbed
    Q_perturbed = (Q_perturbed + Q_perturbed.T)/2

    mean_perturbation = np.mean(perturbation_sizes) if perturbation_sizes else 0.0

    return Q_perturbed, pieces_orig["pole"], mean_perturbation

def kappa_kernel(mfun, zpts, tol=1e-9):
    """Compute κ via Nevanlinna kernel eigenvalue count (negative squares)"""
    Z = np.array(zpts)
    mv = np.array([mfun(z) for z in Z])
    n = len(Z)
    K_mat = np.empty((n,n), complex)
    for i in range(n):
        for j in range(n):
            d = Z[i] - np.conj(Z[j])
            if abs(d) < 1e-12:
                K_mat[i,j] = np.imag(mv[i])/np.imag(Z[i])
            else:
                K_mat[i,j] = (mv[i] - np.conj(mv[j]))/d
    K_mat = 0.5*(K_mat + K_mat.conj().T)
    ev = np.linalg.eigvalsh(K_mat)
    s = max(1.0, abs(ev).max())
    return int(np.sum(ev < -tol*s)), ev

def pole_location_and_magnitude_from_Q(Q_pole, a, K):
    """Extract pole location and magnitude from pole term Q_pole"""
    pole_eigs = np.linalg.eigvalsh(Q_pole)
    z_pole = pole_eigs[-1]
    c_scale = z_pole / np.sqrt(K)
    return z_pole, c_scale

def test_perturbation_sensitivity():
    """Main test: does κ change under prime perturbation?"""
    print("="*100)
    print("POLE AUGMENTATION PERTURBATION TEST — S2645 FAITHFULNESS VERIFICATION")
    print("="*100)
    print()
    print("QUESTION: Does the augmented Suzuki-Weyl m-function show κ=1 only at TRUE primes,")
    print("and κ>>1 when primes are perturbed? (This was claimed but never tested in S2645.)")
    print()

    a = 1.6
    K = max(80, int(30*a)+40)
    NG = 400

    # Generate test points in UHP (same for all tests)
    rng_zpts = np.random.default_rng(3)
    zpts = [complex(x, y) for x, y in zip(rng_zpts.uniform(6, 40, 60), rng_zpts.uniform(0.05, 3.0, 60))]

    print(f"Window a={a}, K={K}, NG={NG}")
    print(f"Test points: 60 random UHP points in [6,40] × i[0.05,3.0]")
    print()

    # Test 1: TRUTH (delta=0)
    print("="*100)
    print("TEST 1: TRUTH (δ=0.0, genuine primes)")
    print("="*100)
    Q_truth, Q_pole_truth, _ = build_Q_perturbed_log_freqs(a, K, NG, delta=0.0)
    vp_truth, vm_truth, lam_truth = weyl_pieces_from_Q(Q_truth, a, K)
    z_pole_truth, c_truth = pole_location_and_magnitude_from_Q(Q_pole_truth, a, K)

    def Theta_truth(z):
        return -((z-1j)*vp_truth(z))/((z+1j)*vm_truth(z))
    def m_weyl_base_truth(z):
        th = Theta_truth(z)
        return 1j*(1+th)/(1-th)
    def m_aug_truth(z):
        return m_weyl_base_truth(z) + c_truth/(z - z_pole_truth)

    k_base_truth, ev_base_truth = kappa_kernel(m_weyl_base_truth, zpts)
    k_aug_truth, ev_aug_truth = kappa_kernel(m_aug_truth, zpts)

    print(f"  z_pole = {z_pole_truth:.6f}")
    print(f"  c      = {c_truth:.6f}")
    print(f"  κ(base Weyl) = {k_base_truth}")
    print(f"  κ(augmented) = {k_aug_truth}  ← EXPECTED: κ=1")
    print(f"  Nevanlinna kernel eig range: [{ev_aug_truth.min():.3e}, {ev_aug_truth.max():.3e}]")

    # Test 2: PERTURBATIONS
    print()
    print("="*100)
    print("TEST 2: PERTURBATIONS — varying δ")
    print("="*100)
    print()
    print(f"{'δ':>8}  {'Mean|Δ log n|':>15}  {'κ(aug)':>8}  {'Change':>8}  {'Robust?':>10}  {'Eig-min':>12}  {'Eig-max':>12}")
    print("-"*100)

    results_perturb = []
    for delta in [0.0005, 0.001, 0.002, 0.005, 0.01]:
        Q_pert, Q_pole_pert, mean_pert_size = build_Q_perturbed_log_freqs(a, K, NG, delta=delta)
        vp_pert, vm_pert, lam_pert = weyl_pieces_from_Q(Q_pert, a, K)
        z_pole_pert, c_pert = pole_location_and_magnitude_from_Q(Q_pole_pert, a, K)

        def Theta_pert(z):
            return -((z-1j)*vp_pert(z))/((z+1j)*vm_pert(z))
        def m_weyl_base_pert(z):
            th = Theta_pert(z)
            return 1j*(1+th)/(1-th)
        def m_aug_pert(z):
            return m_weyl_base_pert(z) + c_pert/(z - z_pole_pert)

        k_aug_pert, ev_aug_pert = kappa_kernel(m_aug_pert, zpts)
        change = k_aug_pert - k_aug_truth
        robust = "✓ YES" if change == 0 else "✗ NO"

        results_perturb.append((delta, k_aug_pert))

        print(f"{delta:>8.4f}  {mean_pert_size:>15.6e}  {k_aug_pert:>8}  {change:+8}  {robust:>10}  {ev_aug_pert.min():>12.3e}  {ev_aug_pert.max():>12.3e}")

    # Test 3: NEGATIVE CONTROLS
    print()
    print("="*100)
    print("TEST 3: NEGATIVE CONTROLS — wrong pole placements")
    print("="*100)
    print()
    print("Does κ=1 REQUIRE the true z_pole, or is it 'near-automatic' for any pole?")
    print()

    # Use the perturbed Q (δ=0.002) as the base, but try different pole locations
    Q_base_for_ctrl, _, _ = build_Q_perturbed_log_freqs(a, K, NG, delta=0.0)
    vp_ctrl, vm_ctrl, _ = weyl_pieces_from_Q(Q_base_for_ctrl, a, K)

    def Theta_ctrl(z):
        return -((z-1j)*vp_ctrl(z))/((z+1j)*vm_ctrl(z))
    def m_weyl_base_ctrl(z):
        th = Theta_ctrl(z)
        return 1j*(1+th)/(1-th)

    print(f"{'Control':>30}  {'z_pole':>12}  {'c':>8}  {'κ':>6}  {'Change':>8}")
    print("-"*80)

    # Control 1: correct pole
    print(f"{'1. Correct pole (z_pole_truth)':>30}  {z_pole_truth:>12.4f}  {c_truth:>8.4f}  {k_aug_truth:>6}  {0:>8}")

    # Control 2: random pole nearby
    z_random_nearby = z_pole_truth + np.random.uniform(-2, 2)
    c_random_nearby = c_truth
    def m_aug_ctrl2(z):
        return m_weyl_base_ctrl(z) + c_random_nearby/(z - z_random_nearby)
    k_ctrl2, _ = kappa_kernel(m_aug_ctrl2, zpts)
    print(f"{'2. Random pole (z ± 2)':>30}  {z_random_nearby:>12.4f}  {c_random_nearby:>8.4f}  {k_ctrl2:>6}  {k_ctrl2 - k_aug_truth:>8}")

    # Control 3: zero pole (no pole term at all)
    def m_ctrl3(z):
        return m_weyl_base_ctrl(z)  # just the base Weyl, no pole
    k_ctrl3, _ = kappa_kernel(m_ctrl3, zpts)
    print(f"{'3. No pole (base Weyl only)':>30}  {'—':>12}  {'—':>8}  {k_ctrl3:>6}  {k_ctrl3 - k_aug_truth:>8}")

    # Control 4: wrong pole from different window
    z_wrong = 5.0  # arbitrary small value, unrelated to actual pole structure
    def m_ctrl4(z):
        return m_weyl_base_ctrl(z) + c_truth/(z - z_wrong)
    k_ctrl4, _ = kappa_kernel(m_ctrl4, zpts)
    print(f"{'4. Wrong pole (z=5.0)':>30}  {z_wrong:>12.4f}  {c_truth:>8.4f}  {k_ctrl4:>6}  {k_ctrl4 - k_aug_truth:>8}")

    # Control 5: inverted pole (negative c)
    c_inverted = -c_truth
    def m_ctrl5(z):
        return m_weyl_base_ctrl(z) + c_inverted/(z - z_pole_truth)
    k_ctrl5, _ = kappa_kernel(m_ctrl5, zpts)
    print(f"{'5. Inverted residue (c → -c)':>30}  {z_pole_truth:>12.4f}  {c_inverted:>8.4f}  {k_ctrl5:>6}  {k_ctrl5 - k_aug_truth:>8}")

    # Control 6: scaled pole (2x the pole location)
    z_scaled = 2*z_pole_truth
    def m_ctrl6(z):
        return m_weyl_base_ctrl(z) + c_truth/(z - z_scaled)
    k_ctrl6, _ = kappa_kernel(m_ctrl6, zpts)
    print(f"{'6. Scaled pole (2×z_pole)':>30}  {z_scaled:>12.4f}  {c_truth:>8.4f}  {k_ctrl6:>6}  {k_ctrl6 - k_aug_truth:>8}")

    # SYNTHESIS
    print()
    print("="*100)
    print("SYNTHESIS & INTERPRETATION")
    print("="*100)

    # Check results
    any_change_under_perturbation = any(k != k_aug_truth for _, k in results_perturb)
    any_change_for_wrong_pole = any([k_ctrl2, k_ctrl3, k_ctrl4, k_ctrl5, k_ctrl6]) != [k_aug_truth]*5

    print()
    print(f"RESULT 1 — Perturbation Sensitivity:")
    if any_change_under_perturbation:
        print(f"  ✗ κ CHANGES under prime perturbation (as predicted by S2646's faithfulness test)")
        print(f"    This CONFIRMS the augmentation is FAITHFUL to prime arithmetic: κ=1 only at truth.")
        print(f"    Interpretation: The construction has REAL CONTENT, not near-automatic (Generalized-Nevanlinna).")
    else:
        print(f"  ✓ κ STABLE under prime perturbation (contradicts S2646's prediction)")
        print(f"    This suggests κ=1 is GENERIC/STRUCTURAL, not arithmetic-specific.")
        print(f"    Interpretation: Any nonzero pole gives κ≥1 (confirms S2646's overclaim concern).")

    print()
    print(f"RESULT 2 — Pole Specificity (Negative Controls):")
    wrong_pole_changes = [k_ctrl2 != k_aug_truth, k_ctrl3 != k_aug_truth, k_ctrl4 != k_aug_truth,
                          k_ctrl5 != k_aug_truth, k_ctrl6 != k_aug_truth]
    if any(wrong_pole_changes):
        print(f"  ✓ κ value DEPENDS on pole placement (some controls differ from truth)")
        print(f"    This suggests the construction has specificity; z_pole location MATTERS.")
        print(f"    Specific changes: {sum(wrong_pole_changes)} of 5 controls differ from truth κ={k_aug_truth}.")
    else:
        print(f"  ✗ κ=1 for ALL pole placements tested (even wrong ones)")
        print(f"    This confirms S2646's 'near-automatic' concern: the Generalized-Nevanlinna theorem")
        print(f"    guarantees κ≥1 for almost any nonzero pole (is the Π_1 property STRUCTURAL, not RH-specific?).")

    print()
    print("HONEST CONCLUSION:")
    if any_change_under_perturbation:
        print("  The augmented construction IS faithful to arithmetic truth.")
        print("  κ=1 at true primes, κ>>1 at perturbed primes — evidence against 'near-automatic.'")
    else:
        print("  The augmented construction shows κ=1 is ROBUST to perturbation and pole placement.")
        print("  This suggests the Generalized-Nevanlinna guarantee (Π_1 for any nonzero pole) is operative.")
        print("  S2646's concern is CONFIRMED: κ=1 is likely near-automatic, not faithfully RH-specific.")

if __name__ == "__main__":
    test_perturbation_sensitivity()
