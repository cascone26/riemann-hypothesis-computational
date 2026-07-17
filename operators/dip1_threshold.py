#!/usr/bin/env python3
"""
Dip1 threshold: at what j does min_ev(τ_N) first approach zero?

j=20: scan showed ev[0]≈-145 everywhere (no Dip1 in accessible N range).
j=40: confirmed Dip1 at N=4200, ev[0]=1.021e-8.

This script scans j=10,15,20,25,30,35 and for each j:
  - Uses N=10,...,N_MAX in steps
  - Finds the N that minimizes |ev[0]|
  - Reports how close |ev[0]_min| gets to zero

The threshold j* where Dip1 first appears tells us:
  - How large lsq must be for the Weil criterion to become accessible
  - Whether the Dip1 is a finite-j artifact or a generic large-j phenomenon

Also studies the minimum eigenvector: which modes does it live in?
"""
import math, numpy as np, scipy.linalg, json, time, sys
from scipy.special import digamma

TWO_PI = 2*math.pi
LOGPI  = math.log(math.pi)
GAMMA1 = 14.134725141734693

def omega_vec(r):
    return np.real(digamma(0.25+0.5j*r)) - LOGPI

def build_tau_minev(N, lsq, G=3000, K_BATCH=2000):
    """Return min eigenvalue and its eigenvector for τ_N."""
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
    ev_all, vecs = scipy.linalg.eigh(tau, subset_by_index=[0,1])
    return float(ev_all[0]), float(ev_all[1]), vecs[:,0], ns

print("=== Dip1 Threshold Analysis ===")
print()

results = {}

# ═══════════════════════════════════════════════════════════
# For small j: scan N over a wide range
# ═══════════════════════════════════════════════════════════
J_CASES = [
    # (j, N_vals_to_scan, G)
    (5,  list(range(5,101,5)), 500),       # lsq=9, tiny
    (10, list(range(5,201,10)), 800),      # lsq=86, small
    (15, list(range(10,501,25)), 1500),    # lsq=795, small
    (20, list(range(20,2001,100)), 2000),  # lsq=7262, medium
    (25, list(range(50,5001,250)), 2500),  # lsq≈68K, larger
]

for j, N_vals, G in J_CASES:
    lsq = round(math.exp(TWO_PI*j/GAMMA1))
    print(f"\n── j={j}: lsq={lsq}, scanning N={N_vals[0]}..{N_vals[-1]} ──")
    print(f"  {'N':>6}  {'ev[0]':>14}  {'ev[1]':>14}  {'|ev[0]|':>14}  dt")

    best_N = None; best_abslam = float('inf')
    row_data = []
    for N in N_vals:
        t0 = time.time()
        try:
            ev0, ev1, vec, ns = build_tau_minev(N, lsq, G=G)
            dt = time.time()-t0
            abslam = abs(ev0)
            mark = " <-- MIN" if abslam < best_abslam else ""
            if abslam < best_abslam:
                best_abslam = abslam; best_N = N
                best_vec = vec.copy(); best_ns = ns.copy()
            print(f"  {N:>6}  {ev0:>14.4e}  {ev1:>14.4e}  {abslam:>14.4e}  ({dt:.1f}s){mark}")
            sys.stdout.flush()
            row_data.append({"N": N, "ev0": ev0, "ev1": ev1})
        except Exception as ex:
            print(f"  {N:>6}  ERROR: {ex}")
            sys.stdout.flush()

    print()
    print(f"  Best |ev[0]| = {best_abslam:.4e} at N={best_N}")

    # Mode decomposition of minimum eigenvector
    if best_N is not None:
        ev0b, ev1b, vec_best, ns_best = build_tau_minev(best_N, lsq, G=G)
        # Which modes dominate the minimum eigenvector?
        top5_idx = np.argsort(np.abs(vec_best))[::-1][:5]
        top5_modes = ns_best[top5_idx]
        top5_weights = np.abs(vec_best[top5_idx])**2
        print(f"  Minimum eigenvector top modes (at N={best_N}):")
        for k in range(5):
            print(f"    n={top5_modes[k]:+4d}: weight={top5_weights[k]:.4f}")

    results[f"j{j}"] = {
        "j": j, "lsq": lsq, "best_N": best_N,
        "best_abs_ev0": float(best_abslam),
        "scan": row_data
    }

# ═══════════════════════════════════════════════════════════
# Summary: how does min |ev[0]| scale with j?
# ═══════════════════════════════════════════════════════════
print()
print("=== SUMMARY: How does min|ev[0]| scale with j? ===")
print()
print(f"  {'j':>4}  {'lsq':>12}  {'N_best':>8}  {'min|ev[0]|':>14}  {'log(min|ev[0]|)':>18}")
for k, v in results.items():
    j = v['j']
    lam = v['best_abs_ev0']
    N = v['best_N']
    lsq = v['lsq']
    log_lam = math.log(lam) if lam > 0 else float('nan')
    print(f"  {j:>4}  {lsq:>12}  {N:>8}  {lam:>14.4e}  {log_lam:>18.4f}")

print()
print("j=40 (from prior confirmation): min|ev[0]| = 1.021e-8 at N=4200")
print()

# Check: does min|ev[0]| decay like exp(-c×j) or j^{-α}?
print("Interpretation:")
print("  If Dip1 exists and min|ev[0]| → 0 as j→∞ (Limit Theorem):")
print("  The approach rate determines how feasible it is to verify RH numerically.")
print()
print("  The Limit Theorem predicts: for Riemann-resonant lsq, min_ev(τ_N) → 0.")
print("  The question is WHICH j values show this numerically (accessible Dip1).")

with open("results/dip1_threshold.json","w") as fh:
    json.dump(results, fh, indent=2)
print("\nSaved → results/dip1_threshold.json")
