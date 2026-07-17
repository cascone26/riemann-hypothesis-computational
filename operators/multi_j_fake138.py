#!/usr/bin/env python3
"""
Multi-j scan: can ANY j level discriminate γ=13.800 (fake) from γ₁=14.135?

Key problem: at j=40, lsq_fake(13.8) = 81.2M > lsq_riem(14.135) = 52.7M.
Larger lsq → more primes → also near-PSD. sel_m335 ≈ 2-35 at j=40.

Hypothesis: at j where lsq_fake(13.8) < lsq_riem(14.135), there might be
better discrimination. But since γ=13.8 < γ₁=14.135:
  lsq_fake(j) = exp(2πj/13.8) > exp(2πj/14.135) = lsq_riem(j)  for ALL j.
So fake_13.8 ALWAYS has larger lsq. No j gives lsq_fake < lsq_riem for γ=13.8.

This multi-j scan CONFIRMS this and shows: which j gives the BEST discrimination
(highest sel = |ev0_fake|/|ev0_riem|), even if small?

Also tests γ=21.022 = γ₂ (second Riemann zero) to verify it also shows near-PSD.
And tests γ=18.000 (between γ₁ and γ₂, NOT a Riemann zero) for comparison.

N=64, G=60000 (fast, accurate at small N).
"""
import math, numpy as np, scipy.linalg, time, sys, json, os
from scipy.special import digamma

TWO_PI=2*math.pi; LOGPI=math.log(math.pi)
GAMMA1=14.134725141734693
GAMMA2=21.022039638771555

G=60000; N=64; K_BATCH=10000

gammas = {
    "gamma1":   14.134725,   # First Riemann zero → expect near-PSD
    "fake13.8": 13.800000,   # NOT a Riemann zero, γ < γ₁ → lsq_fake > lsq_riem
    "fake18.0": 18.000000,   # NOT a Riemann zero, between γ₁ and γ₂
    "gamma2":   21.022040,   # Second Riemann zero → expect near-PSD
    "fake19.0": 19.000000,   # NOT a Riemann zero, near γ₂
}

J_vals = [5, 10, 15, 20, 25, 30, 35, 40, 50, 60]

def build_ev0_from_lsq(lsq):
    if lsq < 2: return float('nan')
    lsq_int = int(lsq)
    L=2*math.log(math.sqrt(lsq)); ns=np.arange(-N,N+1); D=len(ns)
    R=TWO_PI*N/L+50.0; grid=np.linspace(-R,R,G)
    h=grid[1]-grid[0]; w=np.full(G,h); w[0]=w[-1]=0.5*h
    Om=np.real(digamma(0.25+0.5j*grid))-LOGPI
    c=TWO_PI*ns[:,None]/L; d=c-grid[None,:]
    sm=np.abs(d)<1e-12
    E=np.where(sm,np.sqrt(L)+0j,(1/np.sqrt(L))*np.exp(0.5j*d*L)*2*np.sin(0.5*d*L)/d)
    arch=(E*(w*Om)@E.conj().T).real/TWO_PI
    sh=np.sinh(L/4)**2; ni=ns[:,None]; nj=ns[None,:]
    B=32*L*sh*(L**2-16*math.pi**2*ni*nj)/((L**2+16*math.pi**2*ni**2)*(L**2+16*math.pi**2*nj**2))
    sieve=np.ones(lsq_int+1,dtype=bool); sieve[0]=sieve[1]=False
    for i in range(2,int(lsq_int**0.5)+1):
        if sieve[i]: sieve[i*i::i]=False
    lkpp,wts_p=[],[]
    for p in range(2,lsq_int+1):
        if not sieve[p]: continue
        pk=p
        while pk<=lsq_int:
            lkpp.append(math.log(pk)); wts_p.append(math.log(p)/math.sqrt(pk)); pk*=p
    lkpp=np.array(lkpp); wts_p=np.array(wts_p); K=len(lkpp)
    S=np.zeros(D); Sd=np.zeros(D)
    for k0 in range(0,K,K_BATCH):
        k1=min(k0+K_BATCH,K)
        ph=TWO_PI*np.outer(ns,lkpp[k0:k1])/L
        S+=np.sin(ph)@wts_p[k0:k1]
        Sd+=np.cos(ph)@(2*wts_p[k0:k1]*(1-lkpp[k0:k1]/L))
    dn=ns[:,None]-ns[None,:]; Si=S[:,None]; Sj=S[None,:]
    with np.errstate(divide='ignore',invalid='ignore'):
        P=np.where(dn!=0,(Sj-Si)/(math.pi*dn),0.0)
    np.fill_diagonal(P,Sd)
    tau=(arch+B-P); tau=(tau+tau.T)/2
    return float(scipy.linalg.eigvalsh(tau,subset_by_index=[0,0])[0])

print(f"Multi-j scan: N={N}, G={G}")
print(f"  Testing discriminability of γ=13.800 (fake) vs γ₁=14.135 (Riemann)")
print()
print(f"  sel(γ) = |ev0(lsq_gamma)| / |ev0(lsq_riem)|")
print(f"  sel > 1 → γ is MORE negative (Riemann is closer to 0 → Weil discriminates)")
print(f"  sel < 1 → γ is LESS negative (Riemann is NOT closer to 0 → false positive)")
print()

gamma_list = list(gammas.items())
hdr = f"{'j':>5}  {'lsq_riem':>12}  {'ev0_g1':>12}"
for name, g in gamma_list[1:]:  # skip gamma1
    hdr += f"  {'ev0_'+name:>14}  {'sel_'+name:>10}"
print(hdr)
print("-"*120)
sys.stdout.flush()

results = []
for j in J_vals:
    t0=time.time()
    row = {"j": j}
    evs = {}
    for name, gval in gammas.items():
        lsq = round(math.exp(TWO_PI*j/gval))
        row[f"lsq_{name}"] = lsq
        ev = build_ev0_from_lsq(lsq)
        evs[name] = ev
        row[f"ev0_{name}"] = ev

    ev_riem = evs["gamma1"]
    row["lsq_riem"] = row["lsq_gamma1"]

    line = f"{j:>5}  {row['lsq_riem']:>12}  {ev_riem:>12.4e}"
    for name, gval in gamma_list[1:]:
        ev = evs[name]
        sel = abs(ev)/abs(ev_riem) if ev_riem!=0 and ev!=0 else float('nan')
        row[f"sel_{name}"] = sel
        line += f"  {ev:>14.4e}  {sel:>10.2f}"
    dt=time.time()-t0
    print(f"{line}  ({dt:.1f}s)")
    sys.stdout.flush()
    results.append(row)

print()
print("=== Best j for discriminating each fake ===")
for name, gval in gamma_list[1:]:
    sels = [(r["j"], r.get(f"sel_{name}", float('nan'))) for r in results]
    valid = [(j, s) for j, s in sels if not math.isnan(s)]
    if valid:
        best = max(valid, key=lambda x: x[1])
        print(f"  {name} (γ={gval}): best j={best[0]}, sel={best[1]:.2f}")

os.makedirs("results", exist_ok=True)
with open("results/multi_j_fake138.json", "w") as f:
    json.dump({"N": N, "G": G, "J_vals": J_vals, "gammas": gammas, "results": results}, f, indent=2)
print("\nSaved → results/multi_j_fake138.json")
