"""
Candidate (C): Canonical system with boundary conditions AT infinity.

Task: Build a de Branges / canonical-system Hamiltonian H(t) on [0,∞)
whose Weyl m-function directly encodes the explicit formula (arch+pole-prime)
rather than truncating at finite a and taking a→∞ limit.

The finite-a Suzuki construction (suzuki_weyl_kappa.py) has:
- H(t) supported on [0,a], truncated
- κ=0 at every finite a (by Thm 1.5, W-zeros are real)
- a→∞ limit is open (Cor 1.6, Suzuki's own papers state this "requires detailed
  arithmetic analysis which lies beyond their approach")

Candidate (C) approach:
(1) Understand H(t) for finite a (what is its explicit form?)
(2) Determine whether there exists a well-defined limit H_∞(t) on [0,∞)
(3) Classify H_∞ as limit-point or limit-circle at infinity
(4) If limit-point, compute/estimate its Weyl m-function pole structure

This is genuinely hard research. The field's own open problem per Suzuki's papers.
We push as far as honest technical progress allows.
"""

import numpy as np
import mpmath as mp
from scipy.linalg import solve, eigh
from scipy.interpolate import interp1d
import sys, os

sys.path.insert(0, os.path.join(os.getcwd(), 'operators'))

# ============================================================================
# PART 1: UNDERSTAND THE FINITE-a SUZUKI HAMILTONIAN
# ============================================================================

def build_suzuki_screw_Q(a, K, NG=400):
    """
    Build Suzuki's Q matrix (the Weil form truncated to window [a,b]).
    This is arch - prime in Fejér basis.
    
    From suzuki_exact_Q.py, this builds the exact Weil form.
    Now: what is the underlying Hamiltonian H(t)?
    
    In canonical-system theory, the Weil form structure comes from:
    H(t) = a 2x2 matrix (Hamiltonian), usually:
           H(t) = [[0, -ρ(t)], [σ(t), 0]]
    where ρ(t), σ(t) are the symbol/weight functions.
    
    For the explicit formula:
    - arch ↔ the archimedean density Ω(r)/2π
    - prime ↔ the von Mangoldt prime-power weights
    - pole ↔ the ζ(s) pole at s=1
    
    The "screw function" (Suzuki's W characteristic function) arises from
    solving the differential-equation system:
    
        y'' + λ(t)·y = 0   (Schrödinger-like on [0,a])
    
    with appropriate boundary conditions encoding the explicit formula.
    """
    
    # For now: import the existing build_Q and understand its structure
    from suzuki_exact_Q import build_Q
    Q = build_Q(a, K, NG)
    
    # Q is N×N symmetric; diagonal entries relate to arch, off-diagonal to prime.
    # The Hamiltonian H(t) in height t should encode this structure.
    
    return Q

def analyze_suzuki_frequency_structure(a, K, NG=400):
    """
    Analyze the explicit-formula structure in Suzuki's construction.
    
    Key insight: Fejér basis is {cos(2πn(t-a)/2L), sin(...)} for n=0,1,...,N-1
    on window [a, a+2L] of width 2L.
    
    Each basis function has frequency n/(2L) in height-coordinate t.
    The Hamiltonian H(t) should be built from compression of multiplication by:
    
        σ(t) = Ω(t)/2π - P(t)
    
    where:
    - Ω(t) = Re ψ₀(1/4 + it/2) - ln(π)  [archimedean]
    - P(t) = weighted von Mangoldt sum  [prime]
    
    At finite a: P(t) has compact support [say, 2...exp(a)], Ω(t) is smooth.
    """
    
    # Read the explicit formula components
    print("Analyzing Suzuki finite-a structure:")
    print(f"  a = {a}, K = {K} basis functions")
    print(f"  Window width ≈ 2L where 2L relates to bandwidth in screw formulation")
    print(f"  Highest frequency mode: n≈K, frequency ≈ K/(2L)")
    print()
    
    # The symbol σ(t) = Ω(t)/2π - P(t) is pointwise negative on ~30-40% of window
    # yet its compression to Fejér space has EXACTLY κ=1 negative eigenvalue.
    # This is a subtle feature of the Fourier compression.
    
    return {
        'description': 'Fejér bandlimited compression of explicit-formula symbol',
        'K': K,
        'a': a,
    }

# ============================================================================
# PART 2: CHARACTERIZE THE a→∞ LIMIT (Formal Analysis)
# ============================================================================

def formal_limiting_hamiltonian():
    """
    REPORT-tier formal analysis: what should H_∞(t) look like?
    
    As a→∞ in Suzuki's construction:
    1. The window [a, a+2L] slides to ∞
    2. The Fejér basis functions become increasingly localized in frequency
    3. The limiting object should be a canonical system on [0,∞) whose
       Weyl m-function encodes the FULL explicit formula (no truncation)
    
    Expected structure of H_∞(t):
    
    For a canonical system on [0,∞), the Hamiltonian H(t) is a 2×2 matrix:
    
        dX/dt = J H(t) X(t)
    
    where J = [[0, -1], [1, 0]] (symplectic structure) and H(t) ≥ 0 (positivity).
    
    For the explicit formula:
    
        H(t) = [[0, -ρ(t)], [σ(t), 0]]
    
    where σ(t), ρ(t) encode arch and prime weights over [0,∞).
    
    At finite a: σ(t) = 0 for t>a (truncated primes).
    At infinite a: σ(t) should be the FULL σ(t) = Ω(t)/2π - P(t) for ALL t ∈ [0,∞).
    
    The key question: does P(t) (the von Mangoldt sum) grow fast enough
    that the Weyl m-function develops a pole at s=1 (in height coordinates)?
    """
    
    print("="*80)
    print("FORMAL ANALYSIS: Limiting H_∞(t) structure")
    print("="*80)
    print()
    print("Canonical system differential equation:")
    print("  dX/dt = J H(t) X(t),  where J = [[0,-1],[1,0]]")
    print("  H(t) = [[0, -ρ(t)], [σ(t), 0]]  with H(t) ≥ 0 (PSD)")
    print()
    print("For the explicit formula:")
    print("  σ(t) = Ω(t)/2π - P(t)")
    print("  Ω(t) = Re ψ₀(1/4 + it/2) - ln(π)  [archimedean]")
    print("  P(t) = (1/π) Σ_{log n ≤ window} Λ(n)/√n · weight(log n)")
    print()
    print("Finite a (Suzuki): P(t) = 0 for t > a")
    print("Infinite a (Candidate C): P(t) = full sum for ALL t ∈ [0,∞)")
    print()
    print("Critical observation:")
    print("  At finite a: W's characteristic function has all real zeros (Thm 1.5)")
    print("  At a→∞: do zeros approach the critical line (RH ⟺ all real)?")
    print("           or do some leave the real line (violation)?")
    print()
    print("The Weyl m-function of H_∞ should:")
    print("  - Be a generalized-Nevanlinna function (real-symmetric, ℑ m(z)>0 in UHP)")
    print("  - Have pole-count κ ≤ 1 (Krein–Langer / Nevanlinna theorem)")
    print("  - κ = 1 iff RH (the pole at s=1 ↔ one negative square)")
    print()
    
    return {
        'form': 'canonical_system_2x2',
        'H_shape': '[[0, -rho(t)], [sigma(t), 0]]',
        'sigma_t': 'Omega(t)/2π - P(t)',
        'domain': '[0, ∞)',
        'truncation': 'none',
    }

# ============================================================================
# PART 3: LIMIT-POINT VS LIMIT-CIRCLE CLASSIFICATION
# ============================================================================

def weyl_classification_criterion():
    """
    For a canonical system on [0,∞), determine limit-point vs limit-circle.
    
    STANDARD CRITERION (de Branges / Kaltenbäck–Woracek):
    A canonical system with Hamiltonian H(t) on [0,∞) is:
    
    LIMIT-POINT at ∞ if:
        ∫₀^∞ tr H(t) dt < ∞  (trace-integrable)
    or more generally, if the Weyl m-function is UNIQUE (no free boundary data needed).
    
    LIMIT-CIRCLE at ∞ if:
        ∫₀^∞ tr H(t) dt = ∞  (non-integrable)
    and extra boundary conditions at ∞ are needed to pin down the m-function.
    
    For our H_∞(t):
        H(t) = [[0, -ρ(t)], [σ(t), 0]]
        tr H(t) = 0 + 0 = 0   (always!)
    
    So the trace criterion is inconclusive. We need a STRONGER test.
    
    STRONGER TEST (existence of square-integrable solutions):
    For z ∈ ℂ with ℑz ≠ 0, the system
        dX/dt = (z J H(t)) X
    is said to be in limit-point class if there is EXACTLY ONE solution X(t)
    with ∫₀^∞ |X(t)|² dt < ∞ (up to scalar multiple).
    
    If there are TWO linearly-independent L² solutions, it's limit-circle.
    
    PRACTICAL TEST: compute the Weyl m-function at z and check whether
    it's well-defined (unique) without boundary data at ∞.
    """
    
    print("="*80)
    print("WEYL CLASSIFICATION: Limit-Point vs Limit-Circle at ∞")
    print("="*80)
    print()
    print("Standard test (trace criterion):")
    print("  If ∫₀^∞ tr H(t) dt < ∞  ⟹  limit-point")
    print("  If ∫₀^∞ tr H(t) dt = ∞  ⟹  limit-circle (inconclusive)")
    print()
    print("For our H(t) = [[0, -ρ], [σ, 0]]:  tr H = 0 always.")
    print("  ⟹  trace test is INCONCLUSIVE.")
    print()
    print("Stronger test (square-integrable solutions):")
    print("  Count L² solutions to dX/dt = z·J·H(t)·X for z ∈ ℂ, ℑz>0")
    print("  If exactly one (up to scalar): limit-point")
    print("  If two or more: limit-circle")
    print()
    print("Numerical approach:")
    print("  Solve the screw system with increasing domain [0,a]")
    print("  Compute Weyl m-function m_a(z) for each a")
    print("  Check convergence as a→∞:")
    print("    - m_a(z) → m_∞(z) (unique limit): limit-point")
    print("    - m_a(z) oscillates or diverges: potential limit-circle")
    print()
    
    return {
        'criterion': 'square_integrable_solutions',
        'test': 'weyl_function_convergence_as_a_to_infinity',
    }

# ============================================================================
# PART 4: PRELIMINARY NUMERICAL INVESTIGATION
# ============================================================================

def suzuki_weyl_convergence_test(a_values=None, K_base=80):
    """
    NUMERICAL TEST (REPORT-tier):
    Track the Weyl m-function of Suzuki's construction as a grows.
    
    If m_a(z) converges to a well-defined m_∞(z), that's evidence
    of limit-point behavior (and supports the existence of H_∞).
    
    If m_a(z) diverges or oscillates wildly, that suggests limit-circle
    or non-standard behavior.
    """
    
    if a_values is None:
        a_values = np.array([1.5, 2.0, 2.5, 3.0, 3.5, 4.0])
    
    print("="*80)
    print("NUMERICAL TEST: Suzuki Weyl m-function convergence as a→∞")
    print("="*80)
    print()
    
    # Test point in UHP
    z_test = complex(15.0, 1.0)
    
    results = []
    
    for a in a_values:
        K = max(80, int(30*a) + 40)
        
        try:
            from suzuki_exact_Q import build_Q, E_p, shat
            from weyl_nevanlinna_kappa import weyl_pieces_correct
            
            # Try to get the Weyl pieces (solution couple)
            # This requires understanding the actual suzuki_weyl_kappa.py implementation
            
            Q = build_Q(a, K, 400)
            lam_a = np.linalg.eigvalsh(Q)[0] / a
            
            print(f"a = {a:4.1f}, K = {K:3d}: λ_min = {lam_a:.6e}", end="")
            print(f" | m_a(z) convergence tracking...")
            
            results.append({
                'a': a,
                'K': K,
                'lambda_min': lam_a,
            })
            
        except Exception as e:
            print(f"  [error: {str(e)[:50]}]")
    
    print()
    print("STATUS (REPORT): Convergence tracking initiated.")
    print("Full test requires running the actual Weyl construction at each a.")
    print()
    
    return results

# ============================================================================
# PART 5: THEORETICAL ROADMAP FOR CANDIDATE (C)
# ============================================================================

def roadmap_candidate_c():
    """
    Summary of Candidate (C) research plan and honest assessment.
    """
    
    print("\n" + "="*80)
    print("CANDIDATE (C) ROADMAP: Canonical System with Boundary Conditions at ∞")
    print("="*80)
    print()
    
    print("STEP 1: UNDERSTAND FINITE-a HAMILTONIAN [NEEDED]")
    print("  - Read Suzuki 1204.1827 / 2606.09096 to get explicit H_a(t)")
    print("  - Understand how primes are encoded in the Hamiltonian structure")
    print("  - Determine growth rates and limits as a→∞")
    print()
    
    print("STEP 2: CLASSIFY a→∞ LIMIT [OPEN RESEARCH]")
    print("  - Does H_a(t) converge to a well-defined H_∞(t) on [0,∞)?")
    print("  - Is the limit limit-point or limit-circle?")
    print("  - Test via Weyl m-function convergence: m_a(z) →? m_∞(z)")
    print()
    
    print("STEP 3: COMPUTE POLE-COUNT OF H_∞'S WEYL m-FUNCTION [IF LIMIT-POINT]")
    print("  - If limit-point: m_∞ is uniquely determined")
    print("  - Use Krein–Langer / Nevanlinna kernel to compute κ")
    print("  - κ should equal #poles of m_∞ in UHP")
    print("  - RH ⟺ κ = 1")
    print()
    
    print("STEP 4: FAITHFULNESS TEST [CRITICAL]")
    print("  - Does κ(H_∞) = 1 track RH under prime perturbation?")
    print("  - Perturb: P(t) → P(t) + δP(t) (arithmetic perturbation)")
    print("  - Check: does κ change?")
    print("  - If κ(H_∞) is ROBUST (κ stays 1 under perturbation),")
    print("    then H_∞ captures RH-specific arithmetic (non-circular).")
    print("  - If κ stays 1 even for wrong pole locations,")
    print("    then it's a generic Nevanlinna result (circular).")
    print()
    
    print("KNOWN BLOCKERS & HONEST ASSESSMENT:")
    print("  - Suzuki's own papers state a→∞ limit is open")
    print("  - Candidate (A) (inverse-spectral) hit a real category error")
    print("  - Candidate (B) (pole augmentation) is near-automatic by Nevanlinna theorem")
    print("  - Candidate (C) requires building H_∞ explicitly")
    print()
    
    print("ACHIEVABLE MILESTONES (realistic targets):")
    print("  ✓ REPORT: characterize H_a(t) explicitly from Suzuki papers")
    print("  ✓ REPORT: formulate H_∞(t) as a canonical system on [0,∞)")
    print("  ? REPORT→HANDLE: test limit-point via Weyl convergence")
    print("  ? HANDLE: if limit-point, compute κ(H_∞) via pole-count")
    print("  ? HANDLE: faithfulness test under perturbation")
    print()

# ============================================================================
# EXECUTION
# ============================================================================

if __name__ == "__main__":
    print("\n" + "="*80)
    print("CANDIDATE (C) RESEARCH INITIATIVE")
    print("Canonical System with Boundary Conditions at Infinity")
    print("="*80)
    print()
    
    # Part 1: Understand finite-a Suzuki structure
    suzuki_info = analyze_suzuki_frequency_structure(a=2.0, K=88)
    
    # Part 2: Formal analysis of limiting Hamiltonian
    h_inf_formal = formal_limiting_hamiltonian()
    
    # Part 3: Weyl classification theory
    classification = weyl_classification_criterion()
    
    # Part 4: Preliminary numerical test
    convergence_results = suzuki_weyl_convergence_test()
    
    # Part 5: Research roadmap
    roadmap_candidate_c()
    
    print("\n" + "="*80)
    print("NEXT STEP: Read Suzuki 1204.1827 & 2606.09096 to extract explicit H_a(t)")
    print("="*80)

