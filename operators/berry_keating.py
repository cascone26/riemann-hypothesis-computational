#!/usr/bin/env python3
"""
Plan #4: Berry-Keating operator H = xp + px in the Weil basis.

τ_N eigenvalues are Poisson (integrable). τ_N is NOT the operator whose
eigenvalues are the ζ zeros. The Berry-Keating conjecture says that operator
is H = xp + px on L²(R₊). Build the finite-dim truncation in the SAME
Fourier basis used by τ_N and:
A. Compute eigenvalues → compare to ζ zeros (GUE check)
B. Compute τ_N acting on Berry-Keating eigenvectors → is τ_N a function of H_BK?
C. Compare H_BK to arch+B+P term by term — is τ_N = |H_BK|² or similar?

If τ_N = f(H_BK) for a simple f, then:
  τ_N ⪰ 0 ↔ f(H_BK) ⪰ 0
and we need to show f(λ) ≥ 0 for all λ = B-K eigenvalues.
That could break the circle if f ≥ 0 is provable independently.
"""
import math, numpy as np, scipy.linalg, json

TWO_PI = 2*math.pi
GAMMA1 = 14.134725141734693

# Known low Riemann zeros (imaginary parts)
KNOWN_ZEROS = [14.1347, 21.0220, 25.0109, 30.4249, 32.9351,
               37.5862, 40.9187, 43.3271, 48.0052, 49.7738,
               52.9703, 56.4462, 59.3470, 60.8318, 65.1125,
               67.0798, 69.5465, 72.0672, 75.7047, 77.1448]

def build_tau_N(N, lsq, G=8000, G_BATCH=4000, K_BATCH=5000):
    """Build full τ_N = arch + B - P."""
    L = 2*math.log(math.sqrt(lsq))
    ns = np.arange(-N, N+1)
    D = len(ns)
    R = TWO_PI*N/L + 60.0
    grid = np.linspace(-R, R, G)
    h = grid[1]-grid[0]
    w = np.full(G, h); w[0]=w[-1]=0.5*h

    from scipy.special import digamma
    Om = np.real(digamma(0.25+0.5j*grid)) - math.log(math.pi)
    c = TWO_PI*ns[:,None]/L
    d = c - grid[None,:]
    sm = np.abs(d)<1e-12
    E = np.where(sm, np.sqrt(L)+0j,
                 (1/np.sqrt(L))*np.exp(0.5j*d*L)*2*np.sin(0.5*d*L)/d)
    arch = (E*(w*Om) @ E.conj().T).real / TWO_PI
    sh = np.sinh(L/4)**2
    ni = ns[:,None]; nj = ns[None,:]
    num = L**2 - 16*np.pi**2*ni*nj
    den = (L**2+16*np.pi**2*ni**2)*(L**2+16*np.pi**2*nj**2)
    B = 32*L*sh*num/den

    sieve = np.ones(lsq+1, dtype=bool); sieve[0]=sieve[1]=False
    for i in range(2, int(lsq**0.5)+1):
        if sieve[i]: sieve[i*i::i]=False
    lkpp, wts_p = [], []
    for p in range(2, lsq+1):
        if not sieve[p]: continue
        pk=p
        while pk<=lsq:
            lkpp.append(math.log(pk)); wts_p.append(math.log(p)/math.sqrt(pk)); pk*=p
    lkpp=np.array(lkpp); wts_p=np.array(wts_p)
    K=len(lkpp)
    S=np.zeros(D); S_diag=np.zeros(D)
    for k0 in range(0,K,K_BATCH):
        k1=min(k0+K_BATCH,K)
        ph=TWO_PI*np.outer(ns,lkpp[k0:k1])/L
        S+=np.sin(ph)@wts_p[k0:k1]
        S_diag+=np.cos(ph)@(2*wts_p[k0:k1]*(1-lkpp[k0:k1]/L))
    Si=S[:,None]; Sj=S[None,:]
    dni=ns[:,None]-ns[None,:]
    with np.errstate(divide='ignore',invalid='ignore'):
        P=np.where(dni!=0,(Sj-Si)/(math.pi*dni),0.0)
    np.fill_diagonal(P,S_diag)
    tau = arch + B - P
    return (tau+tau.T)/2, arch, B, P, ns, L

def build_berry_keating(N, L):
    """
    Build H_BK = xp + px in the Fourier basis {e^{2πin·/L} / sqrt(L)} for n=-N..N.

    In this basis:
      x = position operator (diagonal in position space, off-diagonal in Fourier)
      p = -i ∂_x = (2πn/L) diagonal in Fourier basis

    Matrix elements in Fourier basis |n⟩ = e^{2πinx/L}/sqrt(L):

    ⟨n|p|m⟩ = (2πm/L) δ_{nm}         (p is diagonal in Fourier basis)
    ⟨n|x|m⟩ = ∫_0^L x e^{2πi(m-n)x/L} dx / L

    For m ≠ n: ⟨n|x|m⟩ = L/(2πi(m-n))  [standard Fourier coefficient of x]
    For m = n: ⟨n|x|n⟩ = L/2            [mean position]

    Then:
    H_BK = xp + px = xp + (xp)† + i[x,p] ... but [x,p] = i on L²
    Actually xp + px = 2 Re(xp) + i[x,p] = 2 Re(xp) - I (on L²)

    Standard: H = xp + px
    In 2nd quantization flavor, the BK conjecture uses H = -(1/2π)(xd/dx + d/dx·x)
    which = -(1/2π)(2x d/dx + 1)

    Let's compute H_BK = xp + px (with p = -i d/dx, ℏ=1) in Fourier basis.
    H_{nm} = ⟨n|xp + px|m⟩ = ⟨n|xp|m⟩ + ⟨n|px|m⟩

    ⟨n|px|m⟩ = Σ_k ⟨n|p|k⟩⟨k|x|m⟩ = (2πn/L)⟨n|x|m⟩
    ⟨n|xp|m⟩ = Σ_k ⟨n|x|k⟩⟨k|p|m⟩ = (2πm/L)⟨n|x|m⟩

    So H_{nm} = (2π/L)(m+n) x_{nm}
    where x_{nm} = ⟨n|x|m⟩

    x_{nm} = L/2 δ_{nm} + L/(2πi(n-m)) (1 - δ_{nm})
    """
    ns = np.arange(-N, N+1)
    D = len(ns)

    # Position operator matrix elements in Fourier basis
    # x_{nm} = L/2 if n==m, else L/(2πi(n-m))
    ni = ns[:,None]; nj = ns[None,:]
    k = ni - nj  # n-m

    x_mat = np.zeros((D,D), dtype=complex)
    np.fill_diagonal(x_mat, L/2)
    off = (k != 0)
    x_mat[off] = L / (2*math.pi*1j*k[off])

    # Momentum operator (diagonal in Fourier basis)
    p_diag = (TWO_PI/L) * ns  # eigenvalues of p

    # H_BK = xp + px
    # H_{nm} = (2π/L)(m+n) * x_{nm}
    nm_sum = ni + nj  # n+m
    H_BK = (TWO_PI/L) * nm_sum * x_mat  # complex matrix

    # Symmetrize: H should be Hermitian
    H_BK = (H_BK + H_BK.conj().T) / 2

    return H_BK.real  # should be nearly real for Hermitian observable

print("Berry-Keating operator construction + comparison with τ_N")
print("="*60)
print()

results = {}

for N, lsq in [(12, 50), (16, 50), (20, 50)]:
    label = f"N={N}, lsq={lsq}"
    print(f"\n── {label} ──")
    tau, arch, B, P, ns, L = build_tau_N(N, lsq)
    H_BK = build_berry_keating(N, L)
    D = 2*N+1

    evs_tau = scipy.linalg.eigvalsh(tau)
    evs_bk  = scipy.linalg.eigvalsh(H_BK)
    evs_bk_sorted = np.sort(evs_bk)

    print(f"  τ_N: min_ev={evs_tau[0]:.4e}, max_ev={evs_tau[-1]:.4e}")
    print(f"  H_BK: min_ev={evs_bk_sorted[0]:.4f}, max_ev={evs_bk_sorted[-1]:.4f}")
    print(f"  H_BK first 8 positive evs: {evs_bk_sorted[evs_bk_sorted>0][:8]}")
    print(f"  Known ζ zeros (first 8):   {KNOWN_ZEROS[:8]}")

    # Spacing statistics of H_BK vs ζ zeros
    bk_pos = evs_bk_sorted[evs_bk_sorted > 0]
    if len(bk_pos) > 3:
        bk_spacings = np.diff(bk_pos)
        mean_bk = bk_spacings.mean()
        print(f"  H_BK positive eigenvalue spacings: mean={mean_bk:.4f}")
        print(f"  ζ zero spacings (first 8): mean={np.diff(KNOWN_ZEROS[:8]).mean():.4f}")

    # Check if H_BK^2 ~ τ_N (up to rescaling)
    H_BK_sq = H_BK @ H_BK
    evs_bk_sq = scipy.linalg.eigvalsh(H_BK_sq)
    # Correlation between tau eigenvalues and H_BK^2 eigenvalues
    if len(evs_tau) == len(evs_bk_sq):
        corr = np.corrcoef(np.sort(evs_tau), np.sort(evs_bk_sq))[0,1]
        print(f"  Correlation τ_N evs vs H_BK² evs: r={corr:.4f}")

    # Check: is tau = a*H_BK + b*I for some scalars?
    # Flatten and regress
    tau_flat = tau.flatten()
    hbk_flat = H_BK.flatten()
    A = np.vstack([hbk_flat, np.ones_like(hbk_flat)]).T
    coeffs, res, _, _ = np.linalg.lstsq(A, tau_flat, rcond=None)
    a_lin, b_lin = coeffs
    tau_pred = a_lin * H_BK + b_lin * np.eye(D)
    residual_lin = np.linalg.norm(tau - tau_pred, 'fro') / np.linalg.norm(tau, 'fro')
    print(f"  Linear fit τ ≈ a*H_BK + b*I: a={a_lin:.4f}, b={b_lin:.4f}, rel_residual={residual_lin:.4f}")

    # Check: tau ≈ a*H_BK² + b*H_BK + c*I
    A2 = np.vstack([(H_BK@H_BK).flatten(), hbk_flat, np.ones_like(hbk_flat)]).T
    coeffs2, _, _, _ = np.linalg.lstsq(A2, tau_flat, rcond=None)
    tau_pred2 = coeffs2[0]*(H_BK@H_BK) + coeffs2[1]*H_BK + coeffs2[2]*np.eye(D)
    residual_quad = np.linalg.norm(tau - tau_pred2, 'fro') / np.linalg.norm(tau, 'fro')
    print(f"  Quadratic fit τ ≈ a*H²+b*H+c: a={coeffs2[0]:.4e}, b={coeffs2[1]:.4e}, c={coeffs2[2]:.4e}")
    print(f"    rel_residual={residual_quad:.4f}")

    results[label] = {
        "tau_min": float(evs_tau[0]), "tau_max": float(evs_tau[-1]),
        "bk_evs_pos": evs_bk_sorted[evs_bk_sorted>0].tolist()[:10],
        "corr_tau_bk2": float(corr) if len(evs_tau)==len(evs_bk_sq) else None,
        "linear_fit_residual": float(residual_lin),
        "quad_fit_residual": float(residual_quad),
        "linear_a": float(a_lin), "linear_b": float(b_lin),
    }

print()
print("=== SYNTHESIS ===")
print()
print("Berry-Keating conjecture: there exists H with eigenvalues = ζ imaginary parts.")
print("Our τ_N is NOT that H (it's Poisson, not GUE).")
print("Key question: is τ_N a polynomial function of H_BK?")
print("  - If τ_N = f(H_BK) with f ≥ 0 provable, we'd have τ_N ⪰ 0 from H_BK structure")
print("  - Linear fit residual ~ 0 would mean τ_N ≈ aH_BK (antisymmetric, not PSD)")
print("  - Quadratic residual ~ 0 would mean τ_N ≈ aH_BK² (always PSD if a>0!)")
print()
print("If τ_N ≈ a·H_BK² + corrections, then:")
print("  τ_N ⪰ 0  ↔  corrections are PSD")
print("  That might be easier to prove than the full statement")

with open("results/berry_keating.json","w") as fh:
    json.dump(results, fh, indent=2)
print("\nSaved → results/berry_keating.json")
