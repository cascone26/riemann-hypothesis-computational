#!/usr/bin/env python3
"""hp_scaling_compress.py -- §236: the AKM Hilbert-Polya test on our truncation.

AKM (2410.20758) mechanism, transported: Theta = scaling-flow generator (skew), acting on
Ker(Delta_F) = leafwise-harmonic forms; spectrum in iR = the zeros; reality automatic.
Our analog:
  - scaling flow on the log-line = translation u -> u+t; generator Theta0 = d/du = diag(2pi i n/L)
    in the Fourier basis e_n. SKEW-adjoint on L^2 (imaginary spectrum), but its raw spectrum is the
    Dirichlet lattice {2pi n/L}, NOT the zeros.
  - Ker(Delta_F) analog = the NEAR-NULL (radical) of the Weil form tau (connes_decisive_v3): the space
    where the arch+prime interaction has bent the lattice. This is the 'harmonic' subspace.
  - THE TEST: compress Theta0 to the near-null V:  Theta_V = V^H Theta0 V  (V = orthonormal near-null
    eigenvectors of tau). Compression of skew is skew => Theta_V AUTOMATICALLY skew-adjoint, spectrum
    imaginary -- exactly AKM's 'reality is automatic'. THE CONTENT: does Im(spec Theta_V) = the zeta
    zeros gamma_k? If yes, Theta_V IS the Hilbert-Polya operator on the truncation (skew + spectrum=zeros).
    If no, the zeros live in the eigenVECTOR structure (as connes_decisive finds them) but NOT as the
    compressed generator spectrum -- which itself localizes the obstruction.

This is DISTINCT from the room's generic-G feasibility gate: we don't solve for any skew-making G, we
use the natural translation-invariant (L^2) metric and the natural generator, per the AKM template.
Verify-before-ship: prints the raw lattice too, and sweeps near-null dimension.
"""
import numpy as np, mpmath as mp
from scipy.special import digamma as _digamma

mp.mp.dps = 30
LAMBDA_SQ = 14
L = float(2 * mp.log(mp.sqrt(LAMBDA_SQ)))
TWO_PI = 2.0 * np.pi
LOGPI = float(np.log(np.pi))
TRUE_ZEROS = [14.134725142, 21.022039639, 25.010857580, 30.424876126, 32.935061588,
              37.586178159, 40.918719012, 43.327073281, 48.005150881, 49.773832478,
              52.970321478, 56.446247697, 59.347044003, 60.831778525]

def vm_table(kmax):
    tab = []
    for k in range(2, kmax + 1):
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
        if ok and p > 0: tab.append((k, float(mp.log(p))))
    return tab
VM = vm_table(LAMBDA_SQ)

def omega_vec(r): return np.real(_digamma(0.25 + 0.5j * r)) - LOGPI

def e_hat_vec(r, n):
    s = TWO_PI * n / L; d = s - r
    out = np.empty_like(d, dtype=complex); small = np.abs(d) < 1e-12
    out[small] = np.sqrt(L) + 0j; dd = d[~small]
    out[~small] = (1.0/np.sqrt(L))*np.exp(0.5j*dd*L)*2.0*np.sin(0.5*dd*L)/dd
    return out

def q_kernel(n, m, y):
    if n == m: return 2.0*(1.0-abs(y)/L)*np.cos(TWO_PI*n*y/L)
    return (np.sin(TWO_PI*m*y/L)-np.sin(TWO_PI*n*y/L))/(np.pi*(n-m))

def prime_mat(idx):
    D=len(idx); P=np.zeros((D,D))
    for ii,n in enumerate(idx):
        for jj in range(ii,D):
            m=idx[jj]; s=0.0
            for k,Lam in VM: s+=Lam*k**(-0.5)*q_kernel(n,m,np.log(k))
            P[ii,jj]=P[jj,ii]=s
    return P

def w02_mat(idx):
    D=len(idx); W=np.zeros((D,D)); sh=np.sinh(L/4)**2
    for ii,n in enumerate(idx):
        for jj,m in enumerate(idx):
            num=(L**2-16*np.pi**2*m*n)
            den=(L**2+16*np.pi**2*m**2)*(L**2+16*np.pi**2*n**2)
            W[ii,jj]=32*L*sh*num/den
    return W

def build_tau(N,R,G):
    idx=list(range(-N,N+1)); D=len(idx)
    grid=np.linspace(-R,R,G); h=grid[1]-grid[0]
    w=np.full(G,h); w[0]=w[-1]=0.5*h
    Om=omega_vec(grid); E=np.empty((D,G),dtype=complex)
    for k,n in enumerate(idx): E[k]=e_hat_vec(grid,n)
    A=E*(w*Om); arch=(A@E.conj().T).real/(2.0*np.pi)
    tau=arch+w02_mat(idx)-prime_mat(idx); tau=(tau+tau.T)/2
    return tau, idx

if __name__=="__main__":
    N,R,G=32,60.0,24000
    print(f"§236 AKM Hilbert-Polya compression test: N={N} (dim {2*N+1}), L={L:.5f}, lattice 2pi/L={TWO_PI/L:.4f}")
    tau, idx = build_tau(N,R,G)
    idx=np.array(idx)
    ev, U = np.linalg.eigh(tau)         # ascending; U columns orthonormal in L2
    absev=np.abs(ev); order=np.argsort(absev)
    mx=absev.max()
    # scaling generator Theta0 = d/du = diag(2pi i n / L)
    theta0_diag = 1j*TWO_PI*idx/L
    print(f"\nRaw generator spectrum = Dirichlet lattice (NOT zeros): {sorted(set(np.round(np.abs(theta0_diag[idx>0]),3)))[:8]} ...")
    zpos=[z for z in TRUE_ZEROS if z<=abs(theta0_diag).max()]
    print(f"True zeros in-range: {zpos}\n")
    print(f"{'nn_dim':>7} {'compressed |Im(spec)| (sorted, positive half)':>50}   {'match to zeros'}")
    for K in [12,16,20,24,28,32]:
        sel=order[:K]                    # K smallest-|eig| eigenvectors = near-null V
        V=U[:,sel]                        # (D,K) orthonormal
        ThetaV = V.conj().T @ (theta0_diag[:,None]*V)   # V^H Theta0 V, (K,K) skew-Hermitian
        w2=np.linalg.eigvals(ThetaV)     # imaginary
        impos=sorted([abs(x.imag) for x in w2])[::-1] # largest first
        impos=[x for x in impos if x>1e-6][:14]
        impos=sorted(impos)
        # match: for each true zero, nearest compressed |Im|
        matches=[]
        for z in zpos[:8]:
            d=min(abs(z-x) for x in impos) if impos else 9.9
            matches.append(f"{z:.1f}->{d:.2f}")
        skew_err=np.linalg.norm(ThetaV+ThetaV.conj().T)
        print(f"{K:>7} {str([round(x,2) for x in impos]):>50}")
        print(f"        skew-err ||T+T^H||={skew_err:.1e}   nearest-zero gaps: {matches}")
    print("\nREAD: if some near-null dim gives compressed |Im(spec)| clustering at the gamma_k (small gaps),")
    print("Theta_V is the Hilbert-Polya operator on the truncation (skew auto + spectrum=zeros). Else the zeros")
    print("are vector-encoded (connes_decisive) but not generator-spectrum -> obstruction localized to the metric.")
