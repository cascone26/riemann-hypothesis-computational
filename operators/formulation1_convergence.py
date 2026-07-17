#!/usr/bin/env python3
"""
Formulation #1 convergence study with CORRECT G scaling.

G must scale as ~N^2 to sample the off-diagonal arch terms:
  E_n(r) × E_m*(r) oscillates with frequency ~ π(n-m)/L
  For n-m = 2N (extreme off-diagonal): freq = 2πN/L
  Need h << L/(2πN): G >> 4πN×R/L ≈ 4π×N×(2πN/L+50)/L ≈ 8π²N²/L²

Conservative: G = max(60000, ceil(500*N^2/L^2))

Tracking convergence of ev[0] at lsq=γ₁ (should → 0) vs lsq=13.8 (should → negative).

At N=64, G=60000 (from zero_crossing_scan.json): ev[0](γ₁) ≈ 1e-12 ≈ 0.
Goal: verify this transition and track the depth ratio.
"""
import math, numpy as np, scipy.linalg, json, time, sys, os
from scipy.special import digamma

TWO_PI = 2*math.pi
LOGPI  = math.log(math.pi)
GAMMA1 = 14.134725141734693
L_G1 = math.log(GAMMA1)  # ≈ 2.648

def G_needed(N, L):
    """Conservative G: sample off-diagonal E_n×E_m at least 20 pts/period."""
    R = TWO_PI*N/L + 50.0
    # Fastest oscillation: 2πN/L rad/unit, period = L/(N)
    h_needed = L / (N * 20)
    return max(60000, int(2*R / h_needed + 1))

def build_ev01_f1(N, lsq):
    """Build τ_N in formulation #1: lsq is small (14.13), all int primes ≤ lsq included."""
    L = math.log(lsq)  # = 2*log(sqrt(lsq))
    G = G_needed(N, L)
    ns = np.arange(-N, N+1); D = len(ns)
    R = TWO_PI*N/L + 50.0
    grid = np.linspace(-R, R, G)
    h = grid[1]-grid[0]; w=np.full(G,h); w[0]=w[-1]=0.5*h
    Om = np.real(digamma(0.25+0.5j*grid)) - LOGPI
    c = TWO_PI*ns[:,None]/L; d = c - grid[None,:]
    sm = np.abs(d)<1e-12
    E = np.where(sm, np.sqrt(L)+0j,
                 (1/np.sqrt(L))*np.exp(0.5j*d*L)*2*np.sin(0.5*d*L)/d)
    arch = (E*(w*Om)@E.conj().T).real/TWO_PI
    sh = np.sinh(L/4)**2; ni=ns[:,None]; nj=ns[None,:]
    B = 32*L*sh*(L**2-16*math.pi**2*ni*nj)/((L**2+16*math.pi**2*ni**2)*(L**2+16*math.pi**2*nj**2))
    lsq_int = int(lsq)
    sieve=np.ones(lsq_int+1,dtype=bool); sieve[0]=sieve[1]=False
    for i in range(2,int(lsq_int**0.5)+1):
        if sieve[i]: sieve[i*i::i]=False
    lkpp,wts_p=[],[]
    for p in range(2,lsq_int+1):
        if not sieve[p]: continue
        pk=p
        while pk<=lsq_int:
            lkpp.append(math.log(pk)); wts_p.append(math.log(p)/math.sqrt(pk)); pk*=p
    lkpp=np.array(lkpp); wts_p=np.array(wts_p)
    if len(lkpp)>0:
        ph=TWO_PI*np.outer(ns,lkpp)/L
        S=np.sin(ph)@wts_p; Sd=np.cos(ph)@(2*wts_p*(1-lkpp/L))
        Si=S[:,None]; Sj=S[None,:]
        dni=ns[:,None]-ns[None,:]
        with np.errstate(divide='ignore',invalid='ignore'):
            P=np.where(dni!=0,(Sj-Si)/(math.pi*dni),0.0)
        np.fill_diagonal(P,Sd)
    else:
        P=np.zeros((D,D))
    tau=(arch+B-P); tau=(tau+tau.T)/2
    ev=scipy.linalg.eigvalsh(tau, subset_by_index=[0,1])
    return float(ev[0]), float(ev[1]), G

print("=== Formulation #1 convergence with correct G ===")
print(f"G_correct(N, L=2.648) = max(60000, 500×N²/L²)")
print()
print(f"  {'N':>6}  {'G':>8}  {'ev0(γ₁)':>14}  {'ev0(13.8)':>14}  {'ratio':>8}  dt")
print("-"*75)
sys.stdout.flush()

data = []
for N in [10, 20, 30, 40, 50, 64, 80]:
    t0 = time.time()
    ev0_g1, ev1_g1, G_g1 = build_ev01_f1(N, GAMMA1)
    ev0_138, ev1_138, _ = build_ev01_f1(N, 13.8)
    dt = time.time()-t0
    ratio = abs(ev0_g1)/abs(ev0_138) if ev0_138!=0 else float('nan')
    print(f"  {N:>6}  {G_g1:>8}  {ev0_g1:>14.6e}  {ev0_138:>14.6e}  {ratio:>8.4f}  ({dt:.1f}s)")
    sys.stdout.flush()
    data.append({"N": N, "G": G_g1, "ev0_g1": ev0_g1, "ev0_138": ev0_138, "ratio": ratio})

print()
print("=== Convergence: ev0(γ₁) as N grows ===")
print("  If RH connection: ev0(γ₁) → 0 from below as N→∞")
print("  If no connection: ev0(γ₁) → negative constant")

os.makedirs("results", exist_ok=True)
with open("results/formulation1_convergence.json","w") as f:
    json.dump(data, f, indent=2)
print("\nSaved → results/formulation1_convergence.json")
