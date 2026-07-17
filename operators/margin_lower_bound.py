#!/usr/bin/env python3
"""
Plan #5: Margin lower bound — functional form and exponent extraction.

Extend Claim A to multiple j levels and multiple γ̃ offsets.
For each j: lsq_j = round(exp(2πj/γ₁)), N_j chosen to capture the Dip1.

Measure M(j, δ) = min_ev(τ_{N_j}, lsq(γ₁+δ)) / |min_ev(τ_{N_j}, lsq(γ₁))|
for small offsets δ = γ̃ - γ₁ ∈ {-0.5, -0.2, -0.1, +0.1, +0.2, +0.5}.

Goals:
  1. Confirm margin ~ |δ|^α — what is α?
  2. Does α depend on j (scale dependence)?
  3. Lower bound: M(j, δ) ≥ f(j) × |δ|^α — what is f(j)?
  4. If α=2: lower bound is M ≥ C|δ|² → uniqueness certificate for zeros

Use small j (j=15,17,20) to stay on Mac.
"""
import math, numpy as np, scipy.linalg, json, time
from scipy.special import digamma

TWO_PI = 2*math.pi
LOGPI  = math.log(math.pi)
GAMMA1 = 14.134725141734693

def lsq_for_gamma(gamma, j):
    return round(math.exp(TWO_PI * j / gamma))

N_DIP1_RATIO = 4200 / math.sqrt(52738111)  # 0.5783, from j=40 confirmed Dip1

def N_for_lsq(lsq):
    return max(10, round(N_DIP1_RATIO * math.sqrt(lsq)))

def min_ev_tau(N, lsq, G=5000, K_BATCH=3000):
    L = 2*math.log(math.sqrt(lsq))
    ns = np.arange(-N, N+1); D = len(ns)
    R = TWO_PI*N/L + 40.0
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
    return float(scipy.linalg.eigvalsh(tau)[0])

print("Margin Lower Bound — Cross-j Analysis")
print("="*60)
print()

J_VALS  = [15, 17, 20]   # j levels to test
DELTAS  = [-0.5, -0.2, -0.1, +0.1, +0.2, +0.5]   # offsets from γ₁

all_results = {}

for j in J_VALS:
    lsq_riem = lsq_for_gamma(GAMMA1, j)
    N_j      = N_for_lsq(lsq_riem)
    print(f"\n── j={j}: lsq_riem={lsq_riem}, N_j={N_j} ──")

    t0 = time.time()
    mev_riem = min_ev_tau(N_j, lsq_riem)
    print(f"  Riemann baseline: min_ev={mev_riem:.4e} ({time.time()-t0:.1f}s)")

    pts = []
    for delta in DELTAS:
        gtilde = GAMMA1 + delta
        lsq_g  = lsq_for_gamma(gtilde, j)
        t1 = time.time()
        try:
            mev_g = min_ev_tau(N_j, lsq_g)
            margin = mev_g / abs(mev_riem) if mev_riem != 0 else float('inf')
            print(f"  δ={delta:+.2f}: lsq={lsq_g}, min_ev={mev_g:.4e}, margin={margin:.2f}x ({time.time()-t1:.1f}s)")
            pts.append({"delta": float(delta), "margin": float(margin), "min_ev": float(mev_g)})
        except Exception as ex:
            print(f"  δ={delta:+.2f}: ERROR {ex}")

    # Fit M ~ A|δ|^α
    valid = [(abs(p["delta"]), p["margin"]) for p in pts if p["margin"] > 0 and abs(p["delta"]) > 0]
    if len(valid) >= 3:
        dabs = np.array([v[0] for v in valid])
        mabs = np.array([v[1] for v in valid])
        coeffs = np.polyfit(np.log(dabs), np.log(mabs), 1)
        alpha = coeffs[0]; A_coeff = math.exp(coeffs[1])
        print(f"  Power law fit: M ≈ {A_coeff:.4f} × |δ|^{alpha:.4f}")
        all_results[f"j{j}"] = {
            "j": j, "lsq_riem": lsq_riem, "N_j": N_j,
            "mev_riem": float(mev_riem),
            "alpha": float(alpha), "A": float(A_coeff),
            "points": pts
        }

# Cross-j comparison: does α vary with j?
if len(all_results) >= 2:
    print()
    print("=== Cross-j Exponent Analysis ===")
    for k, v in all_results.items():
        print(f"  j={v['j']}: α={v['alpha']:.4f}, A={v['A']:.4e}, lsq={v['lsq_riem']}")
    alphas = [v["alpha"] for v in all_results.values()]
    print(f"  α values: {[f'{a:.4f}' for a in alphas]}")
    print(f"  Mean α = {np.mean(alphas):.4f}, std = {np.std(alphas):.4f}")
    if np.std(alphas) < 0.2:
        print(f"  α is STABLE across j — universal exponent")
        mean_a = np.mean(alphas)
        if abs(mean_a - 2.0) < 0.3:
            print(f"  *** α ≈ 2 (quadratic) — M(δ) ≥ C|δ|² ***")
            print(f"  This gives: non-zero γ̃ always has margin ≥ C(γ̃-nearest_zero)²")
            print(f"  Lower bound is TIGHT near zeros and grows away from them.")
            print(f"  Proof target: show this quadratic lower bound analytically.")
        elif abs(mean_a - 1.0) < 0.2:
            print(f"  *** α ≈ 1 (linear) — M(δ) ≥ C|δ| ***")
        elif abs(mean_a - 0.5) < 0.2:
            print(f"  *** α ≈ 1/2 (square root) — M(δ) ≥ C√|δ| ***")

# Theoretical implication summary
print()
print("=== THEORETICAL IMPLICATIONS ===")
print()
print("If M(γ̃) ≥ C|γ̃ - nearest_zero|^α with α confirmed:")
print()
print("  α = 2 (quadratic): Taylor expansion M ≈ (∂²M/∂γ²)|_{γ=zero} × δ²")
print("    This means M has a SIMPLE ZERO at each ζ zero.")
print("    Analytic proof: show ∂M/∂γ = 0 and ∂²M/∂γ² > 0 at all ζ zeros.")
print("    M = 0 at ζ zeros means min_ev(τ) = 0 there — which is the Limit Theorem!")
print()
print("  The Limit Theorem (already confirmed): min_ev(τ_N) → 0 as N → lsq_j(γ₁)")
print("  This means M = min_ev_off / |min_ev_on| and M → ∞ as τ_on → 0.")
print("  But the RATIO M(γ̃) = off_ev / |on_ev|: both go to 0, ratio measures selectivity.")
print()
print("  KEY INSIGHT: The margin M is not a simple distance — it's a RATIO of infinitesimals.")
print("  The exponent α characterizes how much faster the non-zero γ̃ drops vs γ₁.")
print("  α > 0 means the minimum is SHARPER at the actual ζ zeros.")

with open("results/margin_lower_bound.json","w") as fh:
    json.dump(all_results, fh, indent=2)
print("\nSaved → results/margin_lower_bound.json")
