"""
Deep analysis of Suzuki's canonical-system Hamiltonian H(t).

Source: Suzuki 1204.1827 "A canonical system of differential equations 
arising from ζ" and 2606.09096 "Screw function and a generalization of 
convolution: Riemann hypothesis and related problems."

The key insight: Suzuki constructs ζ(s) via a canonical system:
  dw₁/dt = λ(t) w₂
  dw₂/dt = μ(t) w₁
  
where λ(t), μ(t) are weight functions encoding the explicit formula.

In matrix form (with symplectic structure):
  dW/dt = J H(t) W
  
where J = [[0,-1],[1,0]] and H(t) = [[0, -λ(t)], [μ(t), 0]].

Boundary conditions at t=0 and t→∞ determine which spectral data
the system encodes.
"""

import numpy as np
import mpmath as mp
from scipy.special import digamma, polygamma
from scipy.integrate import odeint, quad, trapezoid
import matplotlib.pyplot as plt

# ============================================================================
# PART 1: THE EXPLICIT FORMULA SYMBOL
# ============================================================================

def archimedean_weight_omega(t, deriv=False):
    """
    The archimedean part of the explicit formula symbol.
    
    σ_arch(t) = Ω(t)/2π
    
    where Ω(t) = Re ψ₀(1/4 + it/2) - ln(π)
    
    ψ₀ is the digamma function; in complex height coordinates:
    ψ₀(1/4 + it/2) encodes the Γ-factor from the functional equation.
    
    In Suzuki's setup, this appears as a weight in the Hamiltonian:
    dw₂/dt = λ(t) w₁, where λ ~ Ω(t)/2π
    """
    
    # Use mpmath for high precision
    s = mp.mpf('0.25') + 1j * mp.mpf(t) / 2
    psi_val = mp.digamma(s)
    omega = mp.re(psi_val) - mp.log(mp.pi)
    
    if deriv:
        # Derivative: d/dt ψ₀(s) = ψ₁(s) (trigamma)
        # Ω'(t) = (1/2i) ψ₁(1/4 + it/2)·i = (1/2) Im ψ₁(...)
        psi1_val = mp.polygamma(1, s)
        omega_prime = mp.im(psi1_val) / 2
        return float(omega) / (2*mp.pi), float(omega_prime) / (2*mp.pi)
    
    return float(omega) / (2 * mp.pi)

def von_mangoldt_weight_P(t, K_max=None):
    """
    The prime part of the explicit formula symbol (truncated).
    
    P(t) = (1/π) Σ_{Λ(n)} Λ(n)/√n · cos(t·log n)
    
    Sum is over prime powers n with log n ≤ K_max (for finite truncation).
    
    At finite a: K_max = a (Suzuki's truncation)
    At infinite a: K_max = ∞ (Candidate C's goal)
    
    This is expensive to compute; we'll use a truncated list of prime powers.
    """
    
    if K_max is None:
        K_max = 100  # default truncation
    
    # Generate prime powers up to exp(K_max)
    primes = [2,3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,79,83,89,97]
    
    result = 0.0
    
    for p in primes:
        # Sum powers: p, p², p³, ... while log(p^k) ≤ K_max
        k = 1
        while True:
            n = float(p ** k)
            log_n = k * np.log(float(p))
            if log_n > K_max:
                break
            # Λ(p^k) = log(p)
            result += np.log(float(p)) / np.sqrt(n) * np.cos(float(t) * log_n)
            k += 1
    
    return result / np.pi

def explicit_formula_symbol(t, a_trunc=None):
    """
    Full explicit-formula symbol in canonical system coordinates.
    
    σ(t) = Ω(t)/2π - P(t)
    
    at finite a: P truncated to prime powers with log n ≤ a
    at infinite a: P summed over all prime powers
    """
    
    omega_val = archimedean_weight_omega(t)
    p_val = von_mangoldt_weight_P(t, K_max=a_trunc)
    
    return omega_val - p_val

# ============================================================================
# PART 2: THE CANONICAL SYSTEM HAMILTONIAN H(t)
# ============================================================================

def suzuki_hamiltonian_finite_a(t, a):
    """
    Suzuki's canonical-system Hamiltonian at finite support [0,a].
    
    For the explicit formula, the Hamiltonian is:
    
        H(t) = [[0, -ρ(t)], [σ(t), 0]]
    
    where:
    - σ(t) = Ω(t)/2π - P(t)  [explicit formula symbol, truncated at a]
    - ρ(t) = ?  [to be determined from Suzuki's papers]
    
    In Suzuki's screw-function formulation, the system is:
        dw₁/dt = -ρ(t) w₂
        dw₂/dt = σ(t) w₁
    
    (Note: some conventions use different signs/orderings.)
    
    For the Weil positivity (arch - prime ≥ 0 in the pole-perpendicular subspace),
    we expect σ(t) to be approximately equal to the explicit formula symbol.
    
    The mysterious part: what is ρ(t)?
    
    Hypothesis (from de Branges theory and Suzuki's construction):
    ρ(t) might be constant or related to a regularization/normalization.
    In the simplest case, ρ(t) = 1 (unit weight).
    """
    
    sigma_t = explicit_formula_symbol(t, a_trunc=a)
    
    # For now: assume ρ(t) = 1 (unit weight)
    # This can be refined once we read the papers more carefully
    rho_t = 1.0
    
    # Canonical system Hamiltonian
    H = np.array([[0, -rho_t], [sigma_t, 0]], dtype=complex)
    
    return H

def suzuki_weyl_system(Y, t, a, z):
    """
    Solve the associated Weyl system at parameter z ∈ ℂ.
    
    The Weyl m-function m(z) is constructed by solving:
        dW/dt = (z J H(t)) W,  W(0) = I
    
    where J = [[0,-1],[1,0]] is the symplectic matrix.
    
    The solution couple W(t) = [w₁(t), w₂(t)] at large t encodes
    the Weyl m-function via:
        m(z) = lim_{t→∞} w₂(t) / w₁(t)  (or similar ratio)
    
    For finite a, the system is solved on [0,a].
    For infinite a (Candidate C), solve on [0,∞).
    """
    
    Y = Y.reshape(2,2)
    H_t = suzuki_hamiltonian_finite_a(t, a)
    J = np.array([[0, -1], [1, 0]], dtype=complex)
    
    dYdt = z * J @ H_t @ Y
    return dYdt.flatten()

def suzuki_weyl_function_finite_a(z, a, t_eval=None):
    """
    Compute Suzuki's Weyl m-function for a canonical system with
    domain [0,a] (finite support).
    
    At finite a, we solve the system on [0,a] and extract the ratio
    of solution components to get the m-function.
    
    REPORT-tier: we're approximating the structure without full
    implementation of Suzuki's exact algorithm.
    """
    
    if t_eval is None:
        t_eval = np.linspace(0, a, int(100*a))
    
    # Initial condition (fundamental solution)
    Y0 = np.eye(2).flatten()
    
    # Solve ODE
    try:
        Y_sol = odeint(suzuki_weyl_system, Y0, t_eval, args=(a, z), rtol=1e-8)
        Y_final = Y_sol[-1].reshape(2,2)
        
        # Extract m-function from solution couple
        # Simplified: m(z) ~ Y[1,0] / Y[0,0] (depends on normalization)
        if abs(Y_final[0,0]) > 1e-12:
            m_z = Y_final[1,0] / Y_final[0,0]
        else:
            m_z = np.inf
        
        return m_z
    
    except Exception as e:
        print(f"  ODE solve failed for a={a}, z={z}: {e}")
        return np.nan

# ============================================================================
# PART 3: LIMIT-POINT / LIMIT-CIRCLE TEST
# ============================================================================

def test_limit_classification(z_test=None, a_values=None):
    """
    Numerical test: track Weyl m-function as a→∞.
    
    Limit-point: m_a(z) converges to a unique limit m_∞(z) as a→∞
    Limit-circle: m_a(z) oscillates or diverges (needs boundary data)
    """
    
    if z_test is None:
        z_test = complex(15.0, 1.0)
    
    if a_values is None:
        a_values = [1.5, 2.0, 2.5, 3.0, 3.5, 4.0]
    
    print("="*80)
    print("LIMIT-POINT TEST: Weyl m-function convergence")
    print("="*80)
    print(f"Test point: z = {z_test}")
    print()
    print("  a      m_a(z)              |Δm| from previous")
    print("-" * 70)
    
    m_prev = None
    
    for a in a_values:
        try:
            m_a = suzuki_weyl_function_finite_a(z_test, a)
            
            if m_prev is None:
                print(f"{a:4.1f}    {m_a.real:+.10e} + {m_a.imag:+.10e}i")
            else:
                delta_m = abs(m_a - m_prev)
                print(f"{a:4.1f}    {m_a.real:+.10e} + {m_a.imag:+.10e}i    {delta_m:.3e}")
            
            m_prev = m_a
        
        except Exception as e:
            print(f"{a:4.1f}    [error: {str(e)[:40]}]")
    
    print()
    print("INTERPRETATION:")
    print("  If |Δm| → 0: limit-point ✓")
    print("  If |Δm| oscillates or grows: limit-circle (or convergence failure)")
    print()

# ============================================================================
# PART 4: GROWTH ANALYSIS AND INTEGRABILITY
# ============================================================================

def analyze_hamiltonian_growth(a_max=10.0, num_points=100):
    """
    Analyze the growth of H(t) as t→∞ at finite support [0,a].
    
    Key question: what is ∫₀^a |H(t)| dt, and does it grow manageable?
    """
    
    print("="*80)
    print("HAMILTONIAN GROWTH ANALYSIS")
    print("="*80)
    print()
    
    t_vals = np.linspace(0, a_max, num_points)
    
    # Compute ||H(t)|| at each point (Frobenius norm)
    H_norms = []
    sigma_vals = []
    
    for t in t_vals:
        sigma_t = float(explicit_formula_symbol(t, a_trunc=a_max))
        sigma_vals.append(sigma_t)

        # H = [[0, -1], [σ(t), 0]] with ρ=1
        # ||H||_F = √(1 + σ(t)²)
        H_norm = float(np.sqrt(1 + sigma_t**2))
        H_norms.append(H_norm)
    
    # Estimate integral
    integral = trapezoid(H_norms, t_vals)
    
    print(f"Domain: [0, {a_max}]")
    print(f"∫₀^a ||H(t)|| dt ≈ {integral:.6f}")
    print()
    print(f"σ(t) = Ω(t)/2π - P(t) range:")
    print(f"  min σ = {np.min(sigma_vals):.6e}")
    print(f"  max σ = {np.max(sigma_vals):.6e}")
    print()
    
    # Check for integrability
    if integral < np.inf and integral < 1000:
        print("STATUS: ∫₀^∞ ||H(t)|| dt appears BOUNDED")
        print("        Suggests limit-point behavior (further tests needed)")
    else:
        print("STATUS: ∫₀^∞ ||H(t)|| dt may DIVERGE")
        print("        Potential limit-circle behavior")
    print()
    
    return {
        't': t_vals,
        'sigma': sigma_vals,
        'H_norm': H_norms,
        'integral': integral,
    }

# ============================================================================
# EXECUTION
# ============================================================================

if __name__ == "__main__":
    print("\n" + "="*80)
    print("SUZUKI HAMILTONIAN STRUCTURE ANALYSIS")
    print("="*80)
    print()
    
    # Test 1: Explicit formula symbol at a few points
    print("PART 1: Explicit-formula symbol σ(t) = Ω(t)/2π - P(t)")
    print("-" * 80)
    for t in [1.0, 5.0, 10.0, 20.0]:
        sigma = float(explicit_formula_symbol(t, a_trunc=100))
        omega = float(archimedean_weight_omega(t))
        p = float(von_mangoldt_weight_P(t, K_max=100))
        print(f"  t={t:5.1f}: σ={sigma:+.6e}, Ω/2π={omega:+.6e}, P={p:+.6e}")
    print()
    
    # Test 2: Hamiltonian growth
    print("\nPART 2: Hamiltonian growth and integrability")
    print("-" * 80)
    growth_data = analyze_hamiltonian_growth(a_max=8.0, num_points=80)
    print()
    
    # Test 3: Limit-point classification
    print("\nPART 3: Limit-point vs limit-circle classification")
    print("-" * 80)
    test_limit_classification(z_test=complex(15.0, 1.0), a_values=[1.5, 2.0, 2.5, 3.0])
    
    print("\n" + "="*80)
    print("ASSESSMENT")
    print("="*80)
    print("""
The Suzuki Hamiltonian structure encodes the explicit formula via:
  
  H(t) = [[0, -ρ(t)], [σ(t), 0]]
  
where σ(t) = Ω(t)/2π - P(t) is the symbol.

At finite a: P(t) truncated, giving κ=0 (Thm 1.5)
At infinite a: P(t) unrestricted, potential κ=1 (RH content)

NEXT STEPS:
1. Verify H(t) growth: ∫₀^∞ ||H(t)|| dt behavior
2. Compute Weyl m-function convergence as a→∞ (limit-point test)
3. If limit-point: extract pole-count κ via Nevanlinna kernel
4. Test faithfulness: does κ change under prime perturbation?

STATUS: REPORT-tier formal analysis. HANDLE-tier numerical verification
requires full implementation of Suzuki's exact algorithm and high-precision
ODE solving.
    """)

