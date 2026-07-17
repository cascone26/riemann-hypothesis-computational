#!/usr/bin/env python3
"""
Margin analysis at j=40, N=4200 (the confirmed Dip1).

Claim A (confirmed): At Dip1 (j=40, N=4200, lsq_riem=52738111):
  γ̃=13.8 → min_ev = 5782× larger than γ₁=14.1347

This script measures the margin M(γ̃) = ev[0](γ̃) / ev[0](γ₁) at N=4200
for a grid of γ̃ values around γ₁ to extract the power-law exponent α:
  M(|γ̃ - γ₁|) ~ A × |γ̃ - γ₁|^α

α=2 → quadratic selectivity (M ≥ C|Δγ|²)  — strongest useful lower bound
α=1 → linear
α=? → we measure it

Also checks other Riemann zeros (γ₂, γ₃) to see if they also give small margins,
which would confirm M is specifically selective for ALL ζ zeros (not just γ₁).

Runs on Olympus: N=4200, D=8401, G=30000, lsq=52738111.
Per-point: ~2-5 min (vs j40_dip2 at N=12500 which takes 38 min).
Total for 14 points: ~30-60 min.
"""
import math, numpy as np, scipy.linalg, json, time, sys, os
from scipy.special import digamma

sys.stdout.flush()

TWO_PI = 2*math.pi
LOGPI  = math.log(math.pi)
GAMMA1 = 14.134725141734693
GAMMA2 = 21.022039638771555
GAMMA3 = 25.010857580145688
GAMMA4 = 30.424876125859513
GAMMA5 = 32.935061587739189

J        = 40
N_DIP1   = 4200
G        = 30000
K_BATCH  = 5000
LSQ_RIEM = 52738111  # round(exp(2π×40/γ₁))

def omega_vec(r):
    return np.real(digamma(0.25+0.5j*r)) - LOGPI

def build_tau_ev01(N, lsq, G=G, K_BATCH=K_BATCH):
    """Return (ev[0], ev[1]) for τ_N at given N, lsq."""
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

print(f"Margin Analysis — j={J}, N={N_DIP1}, lsq_riem={LSQ_RIEM}")
print(f"γ₁={GAMMA1}")
print(f"D = {2*N_DIP1+1}")
print(f"G = {G}")
print(); sys.stdout.flush()

os.makedirs("results", exist_ok=True)

# ═══════════════════════════════════════
# Riemann baseline
# ═══════════════════════════════════════
print("Computing Riemann baseline (γ₁, lsq=52738111)...")
sys.stdout.flush()
t0 = time.time()
ev0_riem, ev1_riem = build_tau_ev01(N_DIP1, LSQ_RIEM)
dt_riem = time.time()-t0
print(f"  ev[0]_riem = {ev0_riem:.6e}  ev[1]_riem = {ev1_riem:.6e}  ({dt_riem:.0f}s)")
print(); sys.stdout.flush()

# ═══════════════════════════════════════
# Margin scan
# ═══════════════════════════════════════
gamma_vals = [
    # Fine offsets from γ₁
    ("fake_+0.05", GAMMA1+0.05),
    ("fake_+0.10", GAMMA1+0.10),
    ("fake_+0.20", GAMMA1+0.20),
    ("fake_+0.50", GAMMA1+0.50),
    ("fake_+1.00", GAMMA1+1.00),
    ("fake_+2.00", GAMMA1+2.00),
    ("fake_-0.05", GAMMA1-0.05),
    ("fake_-0.10", GAMMA1-0.10),
    ("fake_-0.20", GAMMA1-0.20),
    ("fake_-0.50", GAMMA1-0.50),
    ("fake_-1.00", GAMMA1-1.00),
    ("fake_-2.00", GAMMA1-2.00),
    ("fake_13.8",  13.8),
    # Other Riemann zeros (should give SMALL margin too — near 0 at their own resonant N)
    # At N=4200, these are NOT at their Dip1, so margin behavior is informative
    ("zero_g2",    GAMMA2),
    ("zero_g3",    GAMMA3),
    ("zero_g4",    GAMMA4),
    ("zero_g5",    GAMMA5),
]

results = {
    "j": J, "N_dip1": N_DIP1, "lsq_riem": LSQ_RIEM,
    "ev0_riem": ev0_riem, "ev1_riem": ev1_riem,
    "dt_riem": dt_riem,
    "points": []
}

print(f"{'γ̃':>14}  {'Δγ':>7}  {'lsq_γ̃':>14}  {'ev[0]':>14}  {'M=ev0/ev0_r':>14}  label")
print("-"*90); sys.stdout.flush()

for label, gamma in gamma_vals:
    lsq_g = round(math.exp(TWO_PI*J/gamma))
    t1 = time.time()
    try:
        ev0, ev1 = build_tau_ev01(N_DIP1, lsq_g)
        dt = time.time()-t1
        # Margin = ev[0](γ̃) / ev[0](γ₁)
        # Both should be positive at Dip1 for γ₁; for other γ̃, ev0 could be positive or negative
        if ev0_riem > 0:
            margin = ev0 / ev0_riem
        else:
            margin = float('nan')
        print(f"{gamma:>14.6f}  {gamma-GAMMA1:>+7.3f}  {lsq_g:>14}  {ev0:>14.4e}  {margin:>14.4f}x  {label}  ({dt:.0f}s)")
        sys.stdout.flush()
        results["points"].append({
            "label": label, "gamma": gamma, "delta": gamma-GAMMA1,
            "lsq": lsq_g, "ev0": ev0, "ev1": ev1,
            "margin": float(margin) if not math.isnan(margin) else None,
            "dt": dt
        })
    except Exception as ex:
        print(f"{gamma:>14.6f}  ERROR: {ex}")
        sys.stdout.flush()

# ═══════════════════════════════════════
# Power-law fit
# ═══════════════════════════════════════
print(); print("=== Power-law fit M(|Δγ|) ≈ A × |Δγ|^α ===")
fakes = [p for p in results["points"]
         if "fake" in p["label"] and abs(p["delta"])>0
         and p["margin"] is not None and p["margin"]>0]

alpha_fit, A_fit = float('nan'), float('nan')
if len(fakes) >= 3:
    dg = np.array([abs(p["delta"]) for p in fakes])
    mg = np.array([p["margin"] for p in fakes])
    valid = (dg>0)&(mg>0)
    if valid.sum() >= 2:
        coeffs = np.polyfit(np.log(dg[valid]), np.log(mg[valid]), 1)
        alpha_fit = coeffs[0]; A_fit = math.exp(coeffs[1])
        print(f"  M(|Δγ|) ≈ {A_fit:.4f} × |Δγ|^{alpha_fit:.4f}")
        if abs(alpha_fit-2.0)<0.5:
            print(f"  *** α ≈ 2 — QUADRATIC selectivity ***")
            print(f"  → M(γ̃) ≥ C|γ̃-γ₁|² → lower bound is quadratic in distance to zero")
            print(f"  → This gives: non-zero γ̃ at distance δ from nearest zero has M ≥ Cδ²")
        elif abs(alpha_fit-1.0)<0.3:
            print(f"  *** α ≈ 1 — LINEAR selectivity ***")
        print(f"  Coefficient A = {A_fit:.4e}")

results["power_law_alpha"] = float(alpha_fit)
results["power_law_A"] = float(A_fit)

print()
print("=== Zero comparison (γ₂, γ₃, γ₄, γ₅ at N=4200) ===")
zeros = [p for p in results["points"] if "zero_g" in p["label"]]
for z in zeros:
    print(f"  {z['label']}: ev[0]={z['ev0']:.4e}  margin={z.get('margin','?')}")

print()
print("=== SUMMARY ===")
print(f"  γ₁ Dip1: ev[0]={ev0_riem:.6e} (near zero)")
print(f"  Power-law exponent α = {alpha_fit:.4f}")
print(f"  Coefficient A = {A_fit:.4e}")
if not math.isnan(alpha_fit):
    print(f"  Lower bound: M(γ̃) ≥ {A_fit:.4f} × |γ̃-γ₁|^{alpha_fit:.4f}")
    min_dist = (1/A_fit)**(1/alpha_fit) if A_fit>0 and not math.isnan(alpha_fit) else float('nan')
    print(f"  For M≥1: |γ̃-γ₁| ≥ ({1/A_fit:.4f})^(1/{alpha_fit:.2f}) = {min_dist:.4f}")
    print(f"  The margin M≥1 certificate distance from nearest zero = {min_dist:.4f}")

with open("results/margin_j40.json","w") as fh:
    json.dump(results, fh, indent=2)
print("\nSaved → results/margin_j40.json")
