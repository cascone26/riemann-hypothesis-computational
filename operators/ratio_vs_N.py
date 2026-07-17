#!/usr/bin/env python3
"""
Does depth_ratio grow with N?

Key question: as N→∞, does depth_ratio(Δγ) = |ev0(γ₁)| / |ev0(γ̃)| grow or shrink?

If it grows → Riemann case becomes relatively shallower → suggestive of PSD at N=∞
If it shrinks → fake also deepens, ratio → 1
If it's flat → both converge to same asymptote (no selectivity at N→∞)

Compare j=25 at two Δγ values over N=500..4000:
  Riemann: lsq = round(exp(2π×25/γ₁)) ≈ 67039
  Fake +0.5: γ̃ = γ₁+0.5 = 14.635, lsq = round(exp(2π×25/14.635))
  Fake -0.5: γ̃ = γ₁-0.5 = 13.635, lsq = round(exp(2π×25/13.635))

Uses G=2500 and K_BATCH=2000 (faster) to scan N quickly.
"""
import math, numpy as np, scipy.linalg, json, time, sys, os
from scipy.special import digamma

TWO_PI = 2*math.pi
LOGPI  = math.log(math.pi)
GAMMA1 = 14.134725141734693
J = 25
G = 2500; K_BATCH = 2000

LSQ_RIEM = round(math.exp(TWO_PI*J/GAMMA1))
LSQ_P05  = round(math.exp(TWO_PI*J/(GAMMA1+0.5)))
LSQ_M05  = round(math.exp(TWO_PI*J/(GAMMA1-0.5)))
LSQ_P10  = round(math.exp(TWO_PI*J/(GAMMA1+1.0)))

print(f"j=25 depth_ratio vs N")
print(f"  Riemann lsq = {LSQ_RIEM}")
print(f"  fake +0.5   lsq = {LSQ_P05}")
print(f"  fake -0.5   lsq = {LSQ_M05}")
print(f"  fake +1.0   lsq = {LSQ_P10}")
print()

def build_ev0(N, lsq):
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
    return float(scipy.linalg.eigvalsh(tau, subset_by_index=[0,0])[0])

N_vals = list(range(200, 4001, 200))

print(f"{'N':>6}  {'ev0_riem':>12}  {'ev0_p05':>12}  {'ev0_m05':>12}  {'ev0_p10':>12}  {'R_p05':>8}  {'R_m05':>8}  {'R_p10':>8}  dt_total")
print("-"*110)
sys.stdout.flush()

results = []
for N in N_vals:
    t0 = time.time()
    e_riem = build_ev0(N, LSQ_RIEM)
    e_p05  = build_ev0(N, LSQ_P05)
    e_m05  = build_ev0(N, LSQ_M05)
    e_p10  = build_ev0(N, LSQ_P10)
    dt = time.time()-t0

    R_p05 = abs(e_riem)/abs(e_p05) if e_p05 != 0 else float('nan')
    R_m05 = abs(e_riem)/abs(e_m05) if e_m05 != 0 else float('nan')
    R_p10 = abs(e_riem)/abs(e_p10) if e_p10 != 0 else float('nan')

    print(f"{N:>6}  {e_riem:>12.4e}  {e_p05:>12.4e}  {e_m05:>12.4e}  {e_p10:>12.4e}  {R_p05:>8.4f}  {R_m05:>8.4f}  {R_p10:>8.4f}  ({dt:.1f}s)")
    sys.stdout.flush()
    results.append({"N": N, "ev0_riem": e_riem, "ev0_p05": e_p05, "ev0_m05": e_m05, "ev0_p10": e_p10,
                    "R_p05": R_p05, "R_m05": R_m05, "R_p10": R_p10})

os.makedirs("results", exist_ok=True)
with open("results/ratio_vs_N.json","w") as f:
    json.dump({"j": J, "lsqs": {"riem": LSQ_RIEM, "p05": LSQ_P05, "m05": LSQ_M05, "p10": LSQ_P10},
               "data": results}, f, indent=2)
print("\nSaved → results/ratio_vs_N.json")

# Summary: is ratio growing or shrinking?
print("\n=== Ratio trend: depth_ratio growing = Riemann case asymptotes to something less negative ===")
if len(results) >= 4:
    early = sum(r["R_p05"] for r in results[:4]) / 4
    late  = sum(r["R_p05"] for r in results[-4:]) / 4
    trend = "GROWING" if late > early else "SHRINKING" if late < early else "FLAT"
    print(f"  R_p05 early avg: {early:.4f}  late avg: {late:.4f}  → {trend}")
