#!/usr/bin/env python3
"""
Plan #2: Claim A margin function M(γ̃, N).

Claim A (confirmed): At Riemann-resonant N=4200 for j=40 (lsq=52M),
non-Riemann γ̃=13.8 gives min_ev 5782.7× larger than Riemann's γ₁=14.1347.

This script probes a SMALLER version (lower j) of the same effect:
for each fake zero γ̃, find the Riemann-resonant N_j(γ̃) = round(exp(2πj/γ̃))
and the Riemann-resonant N_j(γ₁), then compare min_ev(τ_{N_Riem}, lsq_γ̃) vs
min_ev(τ_{N_Riem}, lsq_γ₁).

Goal: fit M(γ̃) = min_ev(τ_{N_Riem}, lsq_γ̃) / min_ev(τ_{N_Riem}, lsq_γ₁)
and find:
  (a) Does M → 1 only when γ̃ ∈ {ζ zeros}?
  (b) What is the functional form of M(γ̃)? Power law? |γ̃ - γ₁|^α?
  (c) Is there an analytic lower bound M(γ̃) ≥ C(γ̃-γ₁)² for γ̃ not a zero?

Use j=20 (small enough to run on Mac in a few minutes per point).
lsq_20(γ₁) = round(exp(2π*20/14.1347)) = round(exp(8.886)) = round(7208) = 7208
"""
import math, numpy as np, scipy.linalg, json
from scipy.special import digamma

TWO_PI = 2*math.pi
LOGPI  = math.log(math.pi)
GAMMA1 = 14.134725141734693
GAMMA2 = 21.022039638771555
GAMMA3 = 25.010857580145688
J = 20   # j value to use

def lsq_resonant(gamma, j=J):
    """Riemann-resonant lsq for zero gamma at level j."""
    return round(math.exp(TWO_PI * j / gamma))

N_DIP1_RATIO = 4200 / math.sqrt(52738111)  # 0.5783 from j=40 confirmed Dip1

def N_for_lsq(lsq):
    return max(10, round(N_DIP1_RATIO * math.sqrt(lsq)))

def omega_vec(r):
    return np.real(digamma(0.25+0.5j*r)) - LOGPI

def build_tau_min_ev(N, lsq, G=6000, K_BATCH=3000):
    """Compute min_ev of τ_N for given N and lsq."""
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
    return float(scipy.linalg.eigvalsh(tau)[0])

print(f"Claim A Margin Function — j={J}")
print("="*60)

# Riemann-resonant N for γ₁ at j=J
lsq_riem = lsq_resonant(GAMMA1, J)
N_riem   = N_for_lsq(lsq_riem)
print(f"γ₁={GAMMA1:.6f}, j={J}: lsq_riem={lsq_riem}, using N={N_riem}")

import time
t0 = time.time()
min_ev_riem = build_tau_min_ev(N_riem, lsq_riem)
print(f"min_ev(τ_N, lsq_riem) = {min_ev_riem:.6e}  ({time.time()-t0:.1f}s)")

# Test γ̃ values: actual ζ zeros (should match riem) + non-zeros
gamma_vals = [
    # Non-Riemann (should give larger margin)
    ("fake_-1.0",  GAMMA1 - 1.0),
    ("fake_-0.5",  GAMMA1 - 0.5),
    ("fake_-0.2",  GAMMA1 - 0.2),
    ("fake_-0.1",  GAMMA1 - 0.1),
    ("fake_+0.1",  GAMMA1 + 0.1),
    ("fake_+0.2",  GAMMA1 + 0.2),
    ("fake_+0.5",  GAMMA1 + 0.5),
    ("fake_+1.0",  GAMMA1 + 1.0),
    ("fake_13.8",  13.8),
    ("fake_15.0",  15.0),
    # Actual ζ zeros (should give SMALL margin — close to Riemann baseline)
    ("zero_gamma2", GAMMA2),
    ("zero_gamma3", GAMMA3),
]

results = {"j": J, "N_riem": N_riem, "lsq_riem": lsq_riem,
           "min_ev_riem": min_ev_riem, "points": []}

print()
print(f"  {'γ̃':>12}  {'Δγ':>8}  {'lsq_γ̃':>12}  {'min_ev_γ̃':>14}  {'margin':>10}  label")
print(f"  {'-'*75}")
for label, gtilde in gamma_vals:
    lsq_g = lsq_resonant(gtilde, J)
    t1 = time.time()
    try:
        min_ev_g = build_tau_min_ev(N_riem, lsq_g)
        margin = min_ev_g / abs(min_ev_riem) if min_ev_riem != 0 else float('inf')
        dt = time.time()-t1
        print(f"  {gtilde:>12.6f}  {gtilde-GAMMA1:>+8.4f}  {lsq_g:>12}  {min_ev_g:>14.6e}  {margin:>10.2f}x  {label}  ({dt:.1f}s)")
        results["points"].append({
            "label": label, "gamma": gtilde, "delta_gamma": float(gtilde-GAMMA1),
            "lsq": lsq_g, "min_ev": float(min_ev_g), "margin": float(margin)
        })
    except Exception as ex:
        print(f"  {gtilde:>12.6f}  ERROR: {ex}")

# Fit margin(Δγ) = A * |Δγ|^α for non-zero points
pts = [p for p in results["points"] if "fake" in p["label"] and p["margin"] > 0]
if len(pts) >= 3:
    print()
    dg = np.array([abs(p["delta_gamma"]) for p in pts])
    mg = np.array([p["margin"] for p in pts])
    valid = (dg > 0) & (mg > 0)
    if valid.sum() >= 2:
        log_dg = np.log(dg[valid])
        log_mg = np.log(mg[valid])
        coeffs = np.polyfit(log_dg, log_mg, 1)
        alpha_fit = coeffs[0]
        A_fit = math.exp(coeffs[1])
        print(f"Power-law fit M(Δγ) ≈ {A_fit:.4f} × |Δγ|^{alpha_fit:.4f}")
        print(f"Exponent α={alpha_fit:.4f}")
        if abs(alpha_fit - 2.0) < 0.3:
            print(f"  *** α ≈ 2 — margin grows QUADRATICALLY with distance from zero! ***")
            print(f"  This would mean: M(γ̃) ≥ C(γ̃-γ_k)² for nearest zero γ_k")
            print(f"  → The lower bound M ≥ 1 gives |γ̃-γ_k| ≤ 1/√C")
            print(f"  → This is a quantitative uniqueness statement for ζ zeros!")
        results["power_law_alpha"] = float(alpha_fit)
        results["power_law_A"]     = float(A_fit)

with open("results/claim_a_margin.json","w") as fh:
    json.dump(results, fh, indent=2)
print("\nSaved → results/claim_a_margin.json")
