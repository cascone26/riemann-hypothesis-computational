#!/usr/bin/env python3
"""
Plan #1: Bost-Connes exactness check.

The Bost-Connes (1995) C*-dynamical system has partition function Z(β) = ζ(β).
Its Hamiltonian H_BC is defined by: the time evolution σ_t acts as
  σ_t(μ_n) = n^{it} μ_n    (n ∈ N)
so H_BC is the operator H_BC |n⟩ = log(n)|n⟩ in the natural basis {|n⟩}.

In the KMS condition, the Weil explicit formula APPEARS as the KMS 2-point function.
The arch+B term in τ_N is (provably) related to the KMS 2-point correlator of H_BC.

Question: Is τ_N the EXACT discretization of the B-C Hamiltonian,
or only an approximation?

This script:
A. Computes H_BC matrix in the same Fourier basis {e^{2πin/L}} used by τ_N.
   H_BC in Fourier basis: H_{BC,nm} = Σ_k log(k) ⟨n|k⟩⟨k|m⟩
   where |k⟩ are position eigenstates.
   But B-C lives on the completion of Q/Z — so we need to be careful about the basis.

B. Alternative: the B-C system's KMS state at β gives the 2-point function:
   ⟨μ_n* μ_m⟩_β = Z(β)^{-1} Tr(e^{-βH_BC} μ_n* μ_m) = n^{-β}/ζ(β) δ_{nm}
   This is diagonal! And the KMS condition:
   ⟨A σ_{iβ}(B)⟩_β = ⟨BA⟩_β
   is automatically satisfied when H_BC = log(n) on |n⟩.

C. The CvS Weil explicit formula: the arch+B-P = τ_N ⪰ 0 is the Weil positivity.
   Connes showed Weil positivity ↔ RH.
   The connection to B-C: the arch term IS the β=1/2 KMS 2-point function
   of a regularized version of H_BC.

D. Direct comparison: build the 'log(n) Hamiltonian' in truncated number-theory basis
   and compare to τ_N's arch+B component.
"""
import math, numpy as np, scipy.linalg, json
from scipy.special import digamma

TWO_PI = 2*math.pi
LOGPI  = math.log(math.pi)
GAMMA1 = 14.134725141734693

def build_tau_components(N, lsq, G=8000, K_BATCH=5000):
    L = 2*math.log(math.sqrt(lsq))
    ns = np.arange(-N, N+1)
    D = len(ns)
    R = TWO_PI*N/L + 60.0
    grid = np.linspace(-R, R, G)
    h = grid[1]-grid[0]; w=np.full(G,h); w[0]=w[-1]=0.5*h
    Om = np.real(digamma(0.25+0.5j*grid)) - LOGPI
    c = TWO_PI*ns[:,None]/L
    d = c - grid[None,:]
    sm = np.abs(d)<1e-12
    E = np.where(sm, np.sqrt(L)+0j,
                 (1/np.sqrt(L))*np.exp(0.5j*d*L)*2*np.sin(0.5*d*L)/d)
    arch = (E*(w*Om)@E.conj().T).real/TWO_PI
    sh = np.sinh(L/4)**2; ni=ns[:,None]; nj=ns[None,:]
    B = 32*L*sh*(L**2-16*math.pi**2*ni*nj)/((L**2+16*math.pi**2*ni**2)*(L**2+16*math.pi**2*nj**2))
    sieve=np.ones(lsq+1,dtype=bool); sieve[0]=sieve[1]=False
    for i in range(2,int(lsq**0.5)+1):
        if sieve[i]: sieve[i*i::i]=False
    lkpp,wts_p=[],[]
    for p in range(2,lsq+1):
        if not sieve[p]: continue
        pk=p
        while pk<=lsq:
            lkpp.append(math.log(pk)); wts_p.append(math.log(p)/math.sqrt(pk)); pk*=p
    lkpp=np.array(lkpp); wts_p=np.array(wts_p)
    K=len(lkpp); S=np.zeros(D); Sd=np.zeros(D)
    for k0 in range(0,K,K_BATCH):
        k1=min(k0+K_BATCH,K)
        ph=TWO_PI*np.outer(ns,lkpp[k0:k1])/L
        S+=np.sin(ph)@wts_p[k0:k1]
        Sd+=np.cos(ph)@(2*wts_p[k0:k1]*(1-lkpp[k0:k1]/L))
    Si=S[:,None]; Sj=S[None,:]
    dni=ns[:,None]-ns[None,:]
    with np.errstate(divide='ignore',invalid='ignore'):
        P=np.where(dni!=0,(Sj-Si)/(math.pi*dni),0.0)
    np.fill_diagonal(P,Sd)
    return arch, B, P, ns, L

def build_bc_hamiltonian_number_basis(N_max):
    """
    Build H_BC = log(n) diagonal in the natural number basis {|1⟩, |2⟩, ..., |N_max⟩}.
    This is the EXACT B-C Hamiltonian in its own basis.
    """
    ns = np.arange(1, N_max+1)
    H_BC = np.diag(np.log(ns.astype(float)))
    return H_BC, ns

def bc_kms_state(beta, N_max):
    """
    KMS state at inverse temperature β: ρ_β = e^{-βH_BC} / Z(β)
    In the natural basis: ρ_{β,nn} = n^{-β} / ζ(β)
    """
    ns = np.arange(1, N_max+1, dtype=float)
    diag = ns**(-beta)
    Z_beta = diag.sum()  # finite truncation of ζ(β)
    return diag / Z_beta, Z_beta

def bc_kms_twopoint(beta, n, m, N_max):
    """
    B-C KMS 2-point function: ⟨μ_n* μ_m⟩_β
    In the Abelian subalgebra: = n^{-β}/ζ(β) if n==m, else 0.
    """
    diag, Z = bc_kms_state(beta, N_max)
    if n == m and 1 <= n <= N_max:
        return diag[n-1]
    return 0.0

def fourier_to_number_change_of_basis(ns_fourier, L, N_max):
    """
    Change of basis: Fourier modes {e^{2πin/L}} → number basis {|k⟩}.
    ⟨k|n_fourier⟩ = e^{2πi*n_fourier*k/L} / sqrt(L)  (sampling at integer k)

    This is how the B-C natural basis relates to the CvS Fourier basis.
    The 'position' x here corresponds to log(k) in the B-C dictionary.
    """
    U = np.zeros((N_max, len(ns_fourier)), dtype=complex)
    for i, k in enumerate(range(1, N_max+1)):
        for j, n in enumerate(ns_fourier):
            U[i,j] = np.exp(2j*math.pi*n*k/L) / math.sqrt(L)
    return U

print("Bost-Connes Exactness Check")
print("="*60)
print()
print("Testing if τ_N components = B-C Hamiltonian in disguise")
print()

results = {}

for N, lsq in [(10, 50), (12, 50)]:
    label = f"N={N}, lsq={lsq}"
    print(f"\n── {label} ──")
    arch, B, P, ns, L = build_tau_components(N, lsq)
    D = 2*N+1
    tau = arch + B - P; tau = (tau+tau.T)/2

    AB = arch + B
    evs_tau = scipy.linalg.eigvalsh(tau)
    evs_AB  = scipy.linalg.eigvalsh(AB)
    evs_arch = scipy.linalg.eigvalsh(arch)

    print(f"  τ_N min_ev = {evs_tau[0]:.4e}")
    print(f"  arch+B min_ev = {evs_AB[0]:.4e}")
    print(f"  arch alone min_ev = {evs_arch[0]:.4e}")

    # B-C in number basis up to N_max ~ lsq
    N_max = min(lsq, 200)
    H_BC, ns_nat = build_bc_hamiltonian_number_basis(N_max)
    evs_BC = np.diag(H_BC)  # = log(1), log(2), ..., log(N_max)
    print(f"  H_BC eigenvalues (first 8): {evs_BC[:8].round(4)}")

    # KMS 2-point function at β=1/2 (critical strip center)
    diag_kms, Z_kms = bc_kms_state(0.5, N_max)
    print(f"  KMS(β=1/2) state (first 8): {diag_kms[:8].round(6)}")
    print(f"  Z(β=1/2) truncated = {Z_kms:.4f} (actual ζ(1/2) diverges — critical line)")

    # The arch kernel in the Weil basis:
    # arch_{nm} = (1/2π) ∫ Ω(r) ψ_n(r)* ψ_m(r) dr
    # where Ω(r) = Re(ψ(1/4+ir/2)) - log(π) is the Weil "distribution"
    # Compare diagonal of arch to B-C 2-point function
    arch_diag = np.diag(arch)
    print(f"  arch diagonal (first 5): {arch_diag[:5].round(6)}")
    print(f"  arch+B diagonal (first 5): {np.diag(AB)[:5].round(6)}")

    # KEY: is the arch diagonal = KMS 2-point function in some form?
    # In the B-C basis at β=1/2: the KMS state gives n^{-1/2}/ζ(1/2)
    # In our Fourier basis n=0,±1,...,±N: the arch_{nn} should be related to n^{-1/2}?
    # Test: is arch_{nn} ~ |n|^{-α} for large |n|?
    pos_ns = ns[ns > 0]
    pos_arch = arch_diag[ns > 0]
    if len(pos_ns) >= 3:
        log_n  = np.log(pos_ns.astype(float))
        log_arch = np.log(np.abs(pos_arch))
        coeffs = np.polyfit(log_n, log_arch, 1)
        print(f"  arch diagonal power law: arch_nn ~ n^{coeffs[0]:.4f} (B-C predicts n^(-β) at β=1/2: exponent=-0.5)")
        if abs(coeffs[0] + 0.5) < 0.15:
            print(f"  *** MATCH: arch diagonal ~ n^(-1/2) as predicted by B-C KMS at β=1/2! ***")

    # Does arch+B diagonal have B-C structure?
    ab_diag = np.diag(AB)
    pos_ab  = ab_diag[ns > 0]
    if len(pos_ns) >= 3:
        log_ab = np.log(np.abs(pos_ab))
        coeffs_ab = np.polyfit(log_n, log_ab, 1)
        print(f"  arch+B diagonal power law: (arch+B)_nn ~ n^{coeffs_ab[0]:.4f}")

    # Structural comparison: is P = H_BC projected onto our basis?
    # The P term comes from prime powers: P_{nm} ↔ Σ_p log(p)/sqrt(p) oscillatory
    # In B-C language, the primes are the atoms of the system
    # The "prime projection" should relate to the B-C Hecke operators T_p
    print(f"\n  P matrix: Frobenius norm = {np.linalg.norm(P,'fro'):.4f}")
    print(f"  arch+B matrix: Frobenius norm = {np.linalg.norm(AB,'fro'):.4f}")
    print(f"  Ratio ||P||/||arch+B|| = {np.linalg.norm(P,'fro')/np.linalg.norm(AB,'fro'):.4f}")

    # Is there a scaling law: as lsq → ∞, does τ_N approach the B-C Hamiltonian?
    # The B-C result: the KMS state at β>1 is unique and = ζ(β)^{-1} Σ n^{-β}
    # As the truncation N_max → ∞: the unique KMS state is well-defined
    # τ_N → 0 means the boundary of the "feasible region" for the KMS condition
    # τ_N ⪰ 0 means the KMS condition CAN be satisfied
    print()
    print(f"  B-C EXACTNESS VERDICT:")
    print(f"  - arch diagonal ~ n^{coeffs[0]:.3f} vs B-C prediction n^(-0.5)")
    discrepancy = abs(coeffs[0] + 0.5)
    if discrepancy < 0.1:
        print(f"  CLOSE MATCH (discrepancy {discrepancy:.4f}) — arch IS approximately the B-C KMS kernel")
    elif discrepancy < 0.25:
        print(f"  PARTIAL MATCH (discrepancy {discrepancy:.4f}) — arch is related to B-C but with corrections")
    else:
        print(f"  NO MATCH (discrepancy {discrepancy:.4f}) — arch is NOT simply the B-C KMS kernel")
    print(f"  NOTE: CvS uses arch in Weil basis (t=Im(ρ) real axis), not direct B-C number basis.")
    print(f"  A non-trivial change of basis could reconcile them.")

    results[label] = {
        "tau_min": float(evs_tau[0]),
        "arch_diag_power": float(coeffs[0]),
        "archB_diag_power": float(coeffs_ab[0]),
        "bc_kms_match": bool(discrepancy < 0.25),
    }

print()
print("=== SYNTHESIS: B-C Exactness ===")
print()
print("If arch_{nn} ~ n^{-1/2}: arch IS the B-C KMS 2-point function at β=1/2.")
print("This would mean: τ_N = (B-C KMS kernel) + boundary - prime correction.")
print("The KMS condition at β=1/2 is the CRITICAL KMS state (β=1 is phase transition,")
print("β<1 gives multiple states — β=1/2 sits at the CRITICAL POINT of this system).")
print()
print("If τ_N IS the B-C Hamiltonian evaluated at the critical KMS state β=1/2:")
print("  τ_N ⪰ 0 ↔ KMS condition holds at β=1/2 ↔ ζ has no zeros off line")
print("This is the EXACT circular dependency, but formulated in C*-algebraic terms.")
print("The B-C axioms don't prove it — they DEFINE it.")
print()
print("KEY OPEN QUESTION: Is there a C*-algebraic UNIQUENESS theorem that forces")
print("the critical KMS state to be positive? If yes, τ_N ⪰ 0 is not an assumption.")

with open("results/bost_connes_exact.json","w") as fh:
    json.dump(results, fh, indent=2)
print("\nSaved → results/bost_connes_exact.json")
