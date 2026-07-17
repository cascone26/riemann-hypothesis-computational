#!/usr/bin/env python3
"""
May's 1972 stability criterion applied to tau_N.

PRISM III marquee finding: May's criterion sigma*sqrt(S*C) < 1 for ecosystem stability.
Real food webs survive via non-random structure: weak interactions dominate,
modular structure, few strong channels.

Application to tau_N: is positivity surviving because the prime-comb coupling
is "food-web-structured" (weak-interaction-dominant, modular) rather than generic?

Metrics computed:
1. sigma = std dev of off-diagonal entries of P (prime-sum coupling matrix)
2. S = matrix sparsity (fraction of "significant" off-diagonal entries)
3. C = connectivity (fraction of entry pairs that are correlated)
4. May_criterion = sigma * sqrt(S * C) — if < 1, structurally stable
5. Off-diagonal decay profile: |P[n, n+k]| vs k — food-web predicts rapid decay
6. Comparison: random matrix with same sigma → would May's criterion hold?
7. Block structure / modularity detection on P

Compare: Riemann lsq vs fake lsq — is Riemann MORE food-web-structured?
"""
import math, numpy as np, scipy.linalg, time, sys, json, os
from scipy.special import digamma

TWO_PI=2*math.pi; LOGPI=math.log(math.pi)
GAMMA1=14.134725141734693
G=30000; N=128; K_BATCH=3000

lsq_list = [
    ("riem_j20", round(math.exp(TWO_PI*20/GAMMA1))),   # 7262
    ("fake_p05", round(math.exp(TWO_PI*20/(GAMMA1+0.5)))),  # 5360
    ("fake_m05", round(math.exp(TWO_PI*20/(GAMMA1-0.5)))),  # 10061
]

def build_matrices(N, lsq):
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
    sieve=np.ones(lsq+1,dtype=bool); sieve[0]=sieve[1]=False
    for i in range(2,int(lsq**0.5)+1):
        if sieve[i]: sieve[i*i::i]=False
    lkpp,wts_p=[],[]
    for p in range(2,lsq+1):
        if not sieve[p]: continue
        pk=p
        while pk<=lsq:
            lkpp.append(math.log(pk)); wts_p.append(math.log(p)/math.sqrt(pk)); pk*=p
    lkpp=np.array(lkpp); wts_p=np.array(wts_p); K=len(lkpp)
    S=np.zeros(D); Sd=np.zeros(D)
    for k0 in range(0,K,K_BATCH):
        k1=min(k0+K_BATCH,K)
        ph=TWO_PI*np.outer(ns,lkpp[k0:k1])/L
        S+=np.sin(ph)@wts_p[k0:k1]
        Sd+=np.cos(ph)@(2*wts_p[k0:k1]*(1-lkpp[k0:k1]/L))
    dn=ns[:,None]-ns[None,:]; Si_mat=S[:,None]; Sj_mat=S[None,:]
    with np.errstate(divide='ignore',invalid='ignore'):
        P=np.where(dn!=0,(Sj_mat-Si_mat)/(math.pi*dn),0.0)
    np.fill_diagonal(P,Sd)
    tau=(arch+B-P); tau=(tau+tau.T)/2
    return arch, B, P, tau, ns

print(f"May's Stability Criterion Diagnostic (N={N}, G={G})")
print(f"Testing: is prime-comb coupling matrix food-web-structured?")
print()

results = {}
for label, lsq in lsq_list:
    print(f"--- {label} (lsq={lsq}) ---")
    t0=time.time()
    arch, B, P, tau, ns = build_matrices(N, lsq)
    D = tau.shape[0]

    # Extract off-diagonal P
    P_od = P.copy(); np.fill_diagonal(P_od, 0.0)
    arch_od = arch.copy(); np.fill_diagonal(arch_od, 0.0)

    # May criterion: sigma * sqrt(S * C)
    # sigma = std dev of ALL off-diagonal entries
    # S = "sparsity" (May used C = connectance = fraction nonzero, we use fraction > threshold)
    # C = matrix size normalization (May uses matrix size N^2)
    # May: stable if sigma * sqrt(N * C) < 1 where C=connectance, N=number of species

    all_od = P_od.flatten()
    sigma = float(np.std(all_od))
    connectance = float(np.sum(np.abs(all_od) > 1e-10) / (D*(D-1)))  # fraction "significant"
    may_param = sigma * math.sqrt(D * connectance)

    # Off-diagonal decay profile: |P[0, k]| vs k (row n=0, columns n=k)
    n0_idx = N  # index of n=0 mode
    decay = [float(abs(P_od[n0_idx, n0_idx+k])) for k in range(1, min(20, D-n0_idx))]

    # Random matrix comparison with same sigma
    rng = np.random.default_rng(42)
    P_rand_od = rng.normal(0, sigma, (D, D)); np.fill_diagonal(P_rand_od, 0)
    P_rand = P_rand_od + np.diag(np.diag(P))  # same diagonal as real P
    tau_rand = (arch + B - P_rand); tau_rand = (tau_rand + tau_rand.T)/2

    ev0_real = float(scipy.linalg.eigvalsh(tau, subset_by_index=[0,0])[0])
    ev0_rand = float(scipy.linalg.eigvalsh(tau_rand, subset_by_index=[0,0])[0])

    # Power-law decay fit: |P[0, n+k]| ~ A * k^(-beta)
    ks = np.arange(1, len(decay)+1)
    valid = np.array(decay) > 1e-15
    if valid.sum() > 3:
        log_k = np.log(ks[valid]); log_d = np.log(np.array(decay)[valid])
        beta_coef = np.polyfit(log_k, log_d, 1)
        beta = -beta_coef[0]  # negative exponent → positive beta means decay
    else:
        beta = float('nan')

    # Modularity: block structure — compare |P[n,m]| for |n-m|=1,5,10,50
    block_diag = float(np.mean([abs(P_od[n0_idx, n0_idx+k]) for k in [1,2,3]]))
    block_far  = float(np.mean([abs(P_od[n0_idx, n0_idx+k]) for k in range(20, 30) if n0_idx+k < D]))

    print(f"  ev0(tau)         = {ev0_real:.4e}  (real prime-comb P)")
    print(f"  ev0(tau_rand)    = {ev0_rand:.4e}  (same sigma, random P off-diag)")
    print(f"  sigma_P_od       = {sigma:.4e}  (std dev of P off-diagonal)")
    print(f"  connectance C    = {connectance:.4f}  (fraction |P_od| > 1e-10)")
    print(f"  May param        = sigma*sqrt(D*C) = {may_param:.4f}")
    print(f"  May stable?      = {'YES' if may_param < 1 else 'NO'} (< 1 → stable)")
    print(f"  Off-diag decay β = {beta:.3f}  (|P[0,n+k]| ~ k^(-β))")
    print(f"  Near-diag |P|    = {block_diag:.4e}  (|k|=1,2,3)")
    print(f"  Far-diag  |P|    = {block_far:.4e}   (|k|=20..30)")
    print(f"  Locality ratio   = {block_diag/block_far if block_far>0 else float('nan'):.1f}  (near/far — food-web has high ratio)")
    print(f"  Duration: {time.time()-t0:.1f}s")
    print()
    sys.stdout.flush()

    results[label] = {
        "lsq": lsq, "D": D,
        "ev0_tau": ev0_real, "ev0_tau_rand": ev0_rand,
        "sigma_P_od": sigma, "connectance": connectance, "may_param": may_param,
        "decay_beta": float(beta) if not math.isnan(beta) else None,
        "near_diag_P": block_diag, "far_diag_P": block_far,
        "locality_ratio": float(block_diag/block_far) if block_far>0 else None,
        "decay_profile": decay,
    }

print("=== Summary: May's criterion comparison ===")
print(f"{'label':>12}  {'ev0(tau)':>12}  {'ev0(rand)':>12}  {'may_param':>10}  {'decay_β':>8}  {'locality':>10}")
for label, r in results.items():
    print(f"{label:>12}  {r['ev0_tau']:>12.4e}  {r['ev0_tau_rand']:>12.4e}  "
          f"{r['may_param']:>10.4f}  {r['decay_beta'] if r['decay_beta'] else float('nan'):>8.3f}  "
          f"{r['locality_ratio'] if r['locality_ratio'] else float('nan'):>10.1f}")

print()
print("PRISM III prediction: Riemann lsq should have HIGHER locality ratio,")
print("  MORE food-web-structured coupling (weak-interaction-dominant, modular)")
print("  than fake lsq. May param should be further below 1 for Riemann.")

os.makedirs("results", exist_ok=True)
with open("results/may_stability_diagnostic.json", "w") as f:
    json.dump(results, f, indent=2)
print("\nSaved → results/may_stability_diagnostic.json")
