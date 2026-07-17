#!/usr/bin/env python3
"""
Margin analysis at j=20, N≈2200 — runs on Olympus (cobo@100.124.234.43).

The correct N scaling: N_Dip1(j) ≈ 105×j (linear in j, NOT sqrt(lsq)).
  j=40: N=4200 (confirmed). 4200/40 = 105.
  j=20: N≈2100 (confirmed by near_miss_tracker, actual Dip1 found at N=2200).

At j=20, lsq_riem = round(exp(2π×20/γ₁)) = 7262.
D = 2×2200+1 = 4401. Large but manageable on Olympus.

This script scans N=1900..2400 at lsq_riem to find exact Dip1,
then at Dip1 computes min_ev for multiple γ̃ offsets to measure the margin.
"""
import math, numpy as np, scipy.linalg, json, time, os
from scipy.special import digamma

TWO_PI = 2*math.pi
LOGPI  = math.log(math.pi)
GAMMA1 = 14.134725141734693
GAMMA2 = 21.022039638771555
GAMMA3 = 25.010857580145688

J = 20
G = 4000       # Integration grid (halved from 8000 to speed up)
K_BATCH = 3000
N_START = 1800  # Scan range for Dip1
N_END   = 2400
N_STEP  = 100   # Coarse step

def omega_vec(r):
    return np.real(digamma(0.25+0.5j*r)) - LOGPI

def build_tau_min2ev(N, lsq, G=G, K_BATCH=K_BATCH):
    """Build τ_N and return (ev[0], ev[1]) = bottom 2 eigenvalues."""
    L = 2*math.log(math.sqrt(lsq))
    ns = np.arange(-N, N+1)
    D = len(ns)
    R = TWO_PI*N/L + 50.0
    grid = np.linspace(-R, R, G)
    h = grid[1]-grid[0]; w=np.full(G,h); w[0]=w[-1]=0.5*h
    Om = omega_vec(grid)
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

lsq_riem = round(math.exp(TWO_PI*J/GAMMA1))
print(f"j={J}, lsq_riem={lsq_riem}, γ₁={GAMMA1}")
print(f"Dip1 expected near N≈{105*J}")
print(f"Scanning N={N_START}..{N_END} step {N_STEP}")
print()

os.makedirs("results", exist_ok=True)
scan = []
t_total = time.time()

# ═══════════════════════════════════════
# Phase 1: N-scan to find Dip1
# ═══════════════════════════════════════
print("Phase 1: N-scan at lsq_riem to find Dip1")
print(f"{'N':>8}  {'ev[0]':>14}  {'ev[1]':>14}  {'dt(s)':>8}")
print("-"*50)

min_ev0_seen = float('inf')
N_dip1 = None
for N in range(N_START, N_END+1, N_STEP):
    t0 = time.time()
    ev0, ev1 = build_tau_min2ev(N, lsq_riem)
    dt = time.time()-t0
    print(f"{N:>8}  {ev0:>14.4e}  {ev1:>14.4e}  {dt:>8.1f}")
    scan.append({"N": N, "ev0": ev0, "ev1": ev1})
    if ev0 < min_ev0_seen:
        min_ev0_seen = ev0
        N_dip1 = N

# Fine-scan around Dip1
print()
print(f"Coarse Dip1 at N={N_dip1}, ev[0]={min_ev0_seen:.4e}")
print("Fine-scanning ±100 around Dip1:")
N_fine_start = max(N_START, N_dip1 - 100)
N_fine_end   = min(N_END, N_dip1 + 100)
fine_scan = []
for N in range(N_fine_start, N_fine_end+1, 25):
    t0=time.time()
    ev0,ev1=build_tau_min2ev(N, lsq_riem)
    dt=time.time()-t0
    print(f"{N:>8}  {ev0:>14.6e}  {ev1:>14.6e}  {dt:>8.1f}")
    fine_scan.append({"N": N, "ev0": ev0, "ev1": ev1})
    if ev0 < min_ev0_seen:
        min_ev0_seen = ev0
        N_dip1 = N

print()
print(f"CONFIRMED Dip1 at N={N_dip1}, ev[0]_min={min_ev0_seen:.6e}")

# ═══════════════════════════════════════
# Phase 2: Margin scan at N_dip1
# ═══════════════════════════════════════
print()
print(f"Phase 2: Margin scan at N={N_dip1}")
print()

gamma_vals = [
    ("riem_baseline", GAMMA1),
    ("fake_-2.0",  GAMMA1-2.0),
    ("fake_-1.0",  GAMMA1-1.0),
    ("fake_-0.5",  GAMMA1-0.5),
    ("fake_-0.2",  GAMMA1-0.2),
    ("fake_-0.1",  GAMMA1-0.1),
    ("fake_+0.1",  GAMMA1+0.1),
    ("fake_+0.2",  GAMMA1+0.2),
    ("fake_+0.5",  GAMMA1+0.5),
    ("fake_+1.0",  GAMMA1+1.0),
    ("fake_+2.0",  GAMMA1+2.0),
    ("fake_13.8",  13.8),
    ("zero_g2",    GAMMA2),
    ("zero_g3",    GAMMA3),
]

margin_results = []
ev0_riem = None
print(f"{'γ̃':>14}  {'Δγ':>7}  {'lsq_γ̃':>10}  {'ev[0]':>14}  {'ev[1]':>14}  {'margin':>10}  label")
print("-"*90)
for label, gamma in gamma_vals:
    lsq_g = round(math.exp(TWO_PI*J/gamma))
    t0=time.time()
    ev0,ev1=build_tau_min2ev(N_dip1, lsq_g)
    dt=time.time()-t0
    if label=="riem_baseline": ev0_riem=ev0
    margin = ev0/abs(ev0_riem) if (ev0_riem and ev0_riem!=0) else float('nan')
    print(f"{gamma:>14.6f}  {gamma-GAMMA1:>+7.3f}  {lsq_g:>10}  {ev0:>14.4e}  {ev1:>14.4e}  {margin:>10.4f}x  {label}  ({dt:.0f}s)")
    margin_results.append({"label": label, "gamma": gamma, "delta": gamma-GAMMA1,
                           "lsq": lsq_g, "ev0": ev0, "ev1": ev1, "margin": margin})

# Power-law fit on fake points
fakes = [p for p in margin_results if "fake" in p["label"] and abs(p["delta"])>0 and p["margin"]>0]
print()
if len(fakes) >= 3:
    dg = np.array([abs(p["delta"]) for p in fakes])
    mg = np.array([p["margin"] for p in fakes])
    valid = (dg>0)&(mg>0)
    if valid.sum() >= 2:
        coeffs = np.polyfit(np.log(dg[valid]), np.log(mg[valid]), 1)
        alpha_fit = coeffs[0]; A_fit = math.exp(coeffs[1])
        print(f"Power-law fit M(Δγ) ≈ {A_fit:.4f} × |Δγ|^{alpha_fit:.4f}")
        if abs(alpha_fit - 2.0) < 0.5:
            print(f"*** α ≈ {alpha_fit:.2f} — QUADRATIC margin! M ≥ C(γ̃-γ_k)² ***")
        elif abs(alpha_fit - 1.0) < 0.3:
            print(f"*** α ≈ {alpha_fit:.2f} — LINEAR margin ***")
    else:
        alpha_fit, A_fit = float('nan'), float('nan')
else:
    alpha_fit, A_fit = float('nan'), float('nan')

out = {
    "j": J, "N_dip1": N_dip1, "lsq_riem": lsq_riem,
    "ev0_dip1": float(min_ev0_seen), "ev0_riem": float(ev0_riem or 0),
    "scan": scan, "fine_scan": fine_scan,
    "margin_results": margin_results,
    "power_law_alpha": float(alpha_fit), "power_law_A": float(A_fit),
    "machine": "olympus"
}
with open("results/margin_olympus_j20.json","w") as fh:
    json.dump(out, fh, indent=2)
print(f"\nTotal time: {time.time()-t_total:.0f}s")
print("Saved → results/margin_olympus_j20.json")
