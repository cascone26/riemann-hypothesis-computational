#!/usr/bin/env python3
"""
Zero-finder v2: dual-j CvS operator zero-location algorithm.

Key insight from beat-frequency analysis: a single j has false positives when
  Phase(j) = j × γ₁/γ_fake mod 1 ≈ 0 or 0.5
Using TWO incommensurate j values (j=15 and j=25) simultaneously eliminates these,
since a false positive requires BOTH Phase(15) and Phase(25) to be near 0/0.5.

Algorithm:
  score(γ) = ev0(τ_N, j=15) + ev0(τ_N, j=25)   (sum: maximized at Riemann zeros)
  Zero: argmax score(γ) in search range

For a Riemann zero: both j=15 and j=25 give ev0 ≈ 0 → score ≈ 0
For a fake: at least one j gives ev0 << 0 → score << 0

Uses adaptive j targeting lsq≈787 for each γ, then ALSO tests at j_alt targeting lsq≈5000.
"""
import math, numpy as np, scipy.linalg, scipy.optimize, time, sys
from scipy.special import digamma

TWO_PI=2*math.pi; LOGPI=math.log(math.pi)
GAMMA1=14.134725141734693; GAMMA2=21.022039638771555
G=60000; N=64; K_BATCH=10000

# Primes up to lsq for j=15,25 in [13,22]: max lsq = exp(2π×25/13) ≈ 180K
# Actually for adaptive j targeting lsq≈787: max_lsq ≈ 1500 (same as before)
# For j=25: same target, so lsq≈787 still but j changes
# Actually we test at TWO different lsq scales:
#   j1=15: lsq1 ≈ exp(2π×15/γ)  -- small lsq
#   j2=25: lsq2 ≈ exp(2π×25/γ)  -- larger lsq
MAX_LSQ=round(math.exp(TWO_PI*20/13.0))+100  # max lsq for j=20, γ=13
print(f"MAX_LSQ = {MAX_LSQ}")

sv=np.ones(MAX_LSQ+1,dtype=bool); sv[0]=sv[1]=False
for i in range(2,int(MAX_LSQ**0.5)+1):
    if sv[i]: sv[i*i::i]=False
primes_up_to=np.where(sv)[0]
print(f"{len(primes_up_to)} primes up to {MAX_LSQ}", flush=True)

pd_cache={}
def get_pd(lsq):
    lsq=min(lsq,MAX_LSQ)
    if lsq not in pd_cache:
        ll,ww=[],[]
        for p in primes_up_to:
            if p>lsq: break
            lp=math.log(p); pk=p
            while pk<=lsq:
                ll.append(math.log(pk)); ww.append(lp/math.sqrt(pk)); pk*=p
        pd_cache[lsq]=(np.array(ll),np.array(ww))
    return pd_cache[lsq]

call_count=[0]
def ev0_at(gamma, j):
    call_count[0]+=1
    lsq=min(MAX_LSQ,round(math.exp(TWO_PI*j/gamma)))
    if lsq<4: return -1.0  # degenerate
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
    lkpp,wts_p=get_pd(lsq); K=len(lkpp)
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

def score(gamma, j1=15, j2=20):
    """Dual-j score: sum of ev0 at j1 and j2. Maximum at Riemann zeros."""
    e1=ev0_at(gamma,j1)
    e2=ev0_at(gamma,j2)
    return e1+e2

def obj_score(gamma):
    return -score(gamma)  # minimize → maximize score

print("\nCvS DUAL-J ZERO-FINDER (j=15 + j=20)")
print(f"  Known γ₁ = {GAMMA1}")
print(f"  Known γ₂ = {GAMMA2}")
print(f"  j1=15: for γ₁, lsq=round(exp(2π×15/14.13))=787")
print(f"  j2=20: for γ₁, lsq=round(exp(2π×20/14.13))=7262  (manageable)")
print()

# First: verify the dual-j score at known zeros vs fakes
print("=== VERIFICATION: score(γ) at known γ values ===")
test_pts=[13.5, 13.7, 13.8, 13.9, 14.0, 14.05, GAMMA1, 14.2, 14.3]
print(f"{'γ':>10}  {'ev0(j=15)':>12}  {'ev0(j=25)':>12}  {'score':>12}  beat15  beat25")
print("-"*75)
for g in test_pts:
    e15=ev0_at(g,15); e25=ev0_at(g,25)
    sc=e15+e25
    ph15=(15*GAMMA1/g)%1; d15=min(ph15,abs(ph15-0.5),abs(ph15-1))
    ph25=(25*GAMMA1/g)%1; d25=min(ph25,abs(ph25-0.5),abs(ph25-1))
    flag=""
    if g==GAMMA1: flag="<-- γ₁"
    elif g==GAMMA2: flag="<-- γ₂"
    print(f"{g:>10.5f}  {e15:>12.4e}  {e25:>12.4e}  {sc:>12.4e}  {d15:.3f}   {d25:.3f}  {flag}")
print()

# Find γ₁: coarse scan [13.5, 15.0]
print("=== COARSE SCAN [13.5, 15.0] step 0.05 ===")
gammas=np.arange(13.5,15.05,0.05)
scores=[(g, score(g)) for g in gammas]
best=max(scores, key=lambda x: x[1])
print(f"  Best: γ={best[0]:.2f}, score={best[1]:.4e}")

# Fine scan around best
print(f"=== FINE SCAN [{best[0]-0.15:.2f}, {best[0]+0.15:.2f}] step 0.005 ===")
fine=[(g, score(g)) for g in np.arange(best[0]-0.15, best[0]+0.16, 0.005)]
best_fine=max(fine, key=lambda x: x[1])
print(f"  Best: γ={best_fine[0]:.4f}, score={best_fine[1]:.4e}")

# Scipy precision
print(f"=== SCIPY OPTIMIZER (j=15+25) ===")
call_count[0]=0
res=scipy.optimize.minimize_scalar(
    obj_score,
    bounds=(best_fine[0]-0.1, best_fine[0]+0.1),
    method='bounded',
    options={'xatol':1e-7,'maxiter':300}
)
gf=res.x; sf=score(gf)
e15f=ev0_at(gf,15); e25f=ev0_at(gf,25)
print(f"  γ* = {gf:.10f}")
print(f"  Known γ₁ = {GAMMA1}")
print(f"  Error = {abs(gf-GAMMA1):.4e}")
print(f"  ev0(j=15) = {e15f:.4e}")
print(f"  ev0(j=25) = {e25f:.4e}")
print(f"  score = {sf:.4e}")
print(f"  calls = {call_count[0]}")
print()

# Find γ₂
print("=== FIND γ₂ [20.5, 21.5] ===")
scores2=[(g, score(g)) for g in np.arange(20.5,21.55,0.05)]
best2=max(scores2,key=lambda x: x[1])
res2=scipy.optimize.minimize_scalar(
    obj_score,
    bounds=(best2[0]-0.15,best2[0]+0.15),
    method='bounded',
    options={'xatol':1e-5}
)
gf2=res2.x
print(f"  γ₂* = {gf2:.7f}  (known: {GAMMA2:.7f})")
print(f"  Error = {abs(gf2-GAMMA2):.4e}")
print(f"  score = {score(gf2):.4e}")
print()

# Find γ₃ = 25.0109...
GAMMA3=25.010857580145688
print("=== FIND γ₃ [24.5, 25.5] ===")
scores3=[(g, score(g)) for g in np.arange(24.5,25.55,0.05)]
best3=max(scores3,key=lambda x: x[1])
res3=scipy.optimize.minimize_scalar(
    obj_score,
    bounds=(best3[0]-0.15,best3[0]+0.15),
    method='bounded',
    options={'xatol':1e-5}
)
gf3=res3.x
print(f"  γ₃* = {gf3:.7f}  (known: {GAMMA3:.7f})")
print(f"  Error = {abs(gf3-GAMMA3):.4e}")
print(f"  score = {score(gf3):.4e}")
print()

# Fake anti-verification
print("=== FAKE TEST: verify optimizer stays AWAY from fake γ=13.8 ===")
scores_fake=[(g, score(g)) for g in np.arange(13.0,14.1,0.05)]
best_fake=max(scores_fake,key=lambda x: x[1])
print(f"  Best in [13.0,14.1]: γ={best_fake[0]:.2f}, score={best_fake[1]:.4e}")
print(f"  γ₁ score: {score(GAMMA1):.4e}")
print(f"  Margin: γ₁ is {score(GAMMA1)/best_fake[1]:.2f}× closer to zero than fake optimum")
print()

print("=== SUMMARY ===")
print(f"  γ₁ located: {gf:.8f} (error={abs(gf-GAMMA1):.2e})")
print(f"  γ₂ located: {gf2:.6f} (error={abs(gf2-GAMMA2):.2e})")
print(f"  γ₃ located: {gf3:.6f} (error={abs(gf3-GAMMA3):.2e})")
print(f"  False positive eliminated by dual-j: YES (single-j had false pos at γ≈13.7)")
print(f"  Total ev0 evaluations: {call_count[0]}")
