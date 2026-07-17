#!/usr/bin/env python3
"""
Find the ACTUAL Dip1 for j=40 in formulation #2 (resonant lsq=52738111).

Prior claim: Dip1 at N=4200. But N_Dip1 ~ 138×log(lsq) scaling from j=20,25 data:
  j=20: N_Dip1=1220, log(7262)=8.89  → ratio=137.2
  j=25: N_Dip1≈1550, log(67039)=11.11 → ratio=139.5
  j=40: log(52738111)=17.78 → N_Dip1 ≈ 138×17.78 ≈ 2455

Scan N=1000..5000 step 500 to find the local minimum of ev[0](N).
Then fine-scan around the minimum.

Runs on Olympus (cobo@100.124.234.43). G=15000 (half of margin_j40 to run faster).
"""
import math, numpy as np, scipy.linalg, json, time, sys
from scipy.special import digamma

TWO_PI = 2*math.pi
LOGPI  = math.log(math.pi)
GAMMA1 = 14.134725141734693

J        = 40
G        = 15000   # Faster than margin_j40's G=30000
K_BATCH  = 5000
LSQ_RIEM = 52738111

def build_ev0(N, lsq=LSQ_RIEM, G=G):
    L = 2*math.log(math.sqrt(lsq))
    ns = np.arange(-N, N+1)
    D = len(ns)
    R = TWO_PI*N/L + 50.0
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
    tau=(arch+B-P); tau=(tau+tau.T)/2
    evs = scipy.linalg.eigvalsh(tau, subset_by_index=[0,1])
    return float(evs[0]), float(evs[1])

print(f"j=40 Dip1 scan  lsq={LSQ_RIEM}  G={G}")
print(f"log(lsq)={math.log(LSQ_RIEM):.4f}  predicted N_Dip1 ≈ {int(138*math.log(LSQ_RIEM))}")
print(f"{'N':>6}  {'ev[0]':>14}  {'ev[1]':>14}  dt")
print("-"*55)
sys.stdout.flush()

scan = []
best_ev0 = float('inf')
best_N = None

# Coarse scan N=1000..5000 step 500
for N in range(1000, 5500, 500):
    t0=time.time()
    ev0, ev1 = build_ev0(N)
    dt=time.time()-t0
    mark = " <-- MIN" if ev0 < best_ev0 else ""
    if ev0 < best_ev0: best_ev0=ev0; best_N=N
    print(f"{N:>6}  {ev0:>14.4e}  {ev1:>14.4e}  ({dt:.0f}s){mark}")
    sys.stdout.flush()
    scan.append({"N": N, "ev0": ev0, "ev1": ev1})

print()
print(f"Coarse Dip1: N={best_N}, ev[0]={best_ev0:.4e}")
print()
print("Fine scan ±500 around coarse Dip1:")
print(f"{'N':>6}  {'ev[0]':>14}  {'ev[1]':>14}  dt")
print("-"*55)
sys.stdout.flush()

fine_scan = []
N_lo = max(500, best_N - 500)
N_hi = best_N + 500
for N in range(N_lo, N_hi + 1, 100):
    t0=time.time()
    ev0, ev1 = build_ev0(N)
    dt=time.time()-t0
    mark = " <-- MIN" if ev0 < best_ev0 else ""
    if ev0 < best_ev0: best_ev0=ev0; best_N=N
    print(f"{N:>6}  {ev0:>14.4e}  {ev1:>14.4e}  ({dt:.0f}s){mark}")
    sys.stdout.flush()
    fine_scan.append({"N": N, "ev0": ev0, "ev1": ev1})

print()
print(f"CONFIRMED Dip1: N={best_N}, ev[0]={best_ev0:.4e}")
print(f"Predicted: N_Dip1 ≈ 138×log(lsq) = {int(138*math.log(LSQ_RIEM))}")
print(f"At N=4200 (prior claim): ev[0]=-14462")
print()
print("=== Recovery trajectory after Dip1 ===")
recovery_scan = []
for N in range(best_N, best_N + 3001, 250):
    t0=time.time()
    ev0, ev1 = build_ev0(N)
    dt=time.time()-t0
    print(f"{N:>6}  {ev0:>14.4e}  ({dt:.0f}s)")
    sys.stdout.flush()
    recovery_scan.append({"N": N, "ev0": ev0, "ev1": ev1})

import os; os.makedirs("results", exist_ok=True)
out = {"j": 40, "lsq": LSQ_RIEM, "G": G,
       "coarse_scan": scan, "fine_scan": fine_scan,
       "dip1_N": best_N, "dip1_ev0": best_ev0,
       "recovery_scan": recovery_scan}
with open("results/j40_dip1_scan.json","w") as f:
    json.dump(out, f, indent=2)
print("Saved → results/j40_dip1_scan.json")
