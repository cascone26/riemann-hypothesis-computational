#!/usr/bin/env python3
"""
Formulation #1 at large N: direct lsq = γ₁ = 14.1347...

Key question: as N grows, does λ_min(lsq=γ₁) → 0 from above
while λ_min(lsq=non-zero) → negative?

lsq=14.1347 means only primes {2,3,5,7,11,13} in the sum (6 primes + prime powers).
Matrix size D=2N+1. But the prime sum is trivial (≤ 50 terms total).
This is EXTREMELY fast even for N=5000!

The L = 2*log(sqrt(lsq)) = log(lsq) varies continuously.
At lsq=γ₁=14.134725: L = log(14.134725) = 2.6484...

Scan N = 10, 50, 100, 500, 1000, 5000, 10000, 50000
For each N, compute λ_min at:
  lsq = γ₁ = 14.134725 (Riemann zero — expect λ_min → 0 from above)
  lsq = 13.8 (non-zero — expect λ_min → negative)
  lsq = 13.0 (far from zero — more negative)
  lsq = 15.0 (between γ₁ and γ₂ — expect negative)
  lsq = 21.022 = γ₂ (second zero — expect λ_min → 0)
"""
import math, numpy as np, scipy.linalg, json, time, sys, os
from scipy.special import digamma

TWO_PI = 2*math.pi
LOGPI  = math.log(math.pi)
GAMMA1 = 14.134725141734693
GAMMA2 = 21.022039638771555

def build_ev0_formulation1(N, lsq, G_RATIO=50.0):
    """Formulation #1: lsq is small (e.g. 14.1347), N can be large.
    G is chosen relative to N to maintain integration quality."""
    G = max(2000, int(G_RATIO * N))  # scales with N for resolution
    L = 2*math.log(math.sqrt(lsq))
    ns = np.arange(-N, N+1)
    D = len(ns)
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

    # Prime sum: only primes ≤ lsq (at lsq=14.1347: just {2,3,5,7,11,13})
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

    K=len(lkpp)
    if K > 0:
        S=np.zeros(D); Sd=np.zeros(D)
        ph=TWO_PI*np.outer(ns,lkpp)/L
        S=np.sin(ph)@wts_p; Sd=np.cos(ph)@(2*wts_p*(1-lkpp/L))
        Si=S[:,None]; Sj=S[None,:]
        dni=ns[:,None]-ns[None,:]
        with np.errstate(divide='ignore',invalid='ignore'):
            P=np.where(dni!=0,(Sj-Si)/(math.pi*dni),0.0)
        np.fill_diagonal(P,Sd)
    else:
        P = np.zeros((D,D))

    tau=(arch+B-P); tau=(tau+tau.T)/2
    ev = scipy.linalg.eigvalsh(tau, subset_by_index=[0,1])
    return float(ev[0]), float(ev[1]), G

print("=== Formulation #1 at large N: lsq = γ₁ directly ===")
print(f"γ₁ = {GAMMA1}  (primes: {{2,3,5,7,11,13}})")
print(f"γ₂ = {GAMMA2}  (primes: {{2,3,5,7,11,13,17,19}})")
print()

LSQ_VALS = [
    ("γ₁", GAMMA1),
    ("13.8", 13.8),
    ("13.0", 13.0),
    ("15.0", 15.0),
    ("γ₂", GAMMA2),
]

N_VALS = [10, 50, 100, 500, 1000, 5000, 10000]

results = {name: [] for name, _ in LSQ_VALS}

for N in N_VALS:
    t0 = time.time()
    row = {}
    for name, lsq in LSQ_VALS:
        ev0, ev1, G = build_ev0_formulation1(N, lsq)
        row[name] = (ev0, ev1)
    dt = time.time()-t0

    print(f"\nN = {N}  (D={2*N+1})  (dt={dt:.1f}s)")
    print(f"  {'lsq':>8}  {'ev[0]':>14}  {'ev[1]':>14}")
    for name, lsq in LSQ_VALS:
        ev0, ev1 = row[name]
        results[name].append({"N": N, "ev0": ev0, "ev1": ev1})
        flag = "  ← should be ≥ 0" if name.startswith("γ") else ""
        print(f"  {name:>8}  {ev0:>14.6e}  {ev1:>14.6e}{flag}")
    sys.stdout.flush()

print()
print("=== Summary: λ_min trajectory with N ===")
print(f"  {'N':>8}  {'γ₁ ev[0]':>14}  {'13.8 ev[0]':>14}  {'ratio':>10}")
for i, N in enumerate(N_VALS):
    ev0_g1 = results["γ₁"][i]["ev0"]
    ev0_138 = results["13.8"][i]["ev0"]
    ratio = abs(ev0_g1)/abs(ev0_138) if ev0_138 != 0 else float('nan')
    print(f"  {N:>8}  {ev0_g1:>14.6e}  {ev0_138:>14.6e}  {ratio:>10.4f}")

os.makedirs("results", exist_ok=True)
with open("results/formulation1_large_N.json","w") as f:
    json.dump({"N_vals": N_VALS, "lsq_vals": {n: l for n,l in LSQ_VALS},
               "results": {k: v for k,v in results.items()}}, f, indent=2)
print("\nSaved → results/formulation1_large_N.json")
