#!/usr/bin/env python3
"""
Mechanistic analysis: why does j=15 give 1.79M× selectivity for fake_13.8?

At j=15:
  lsq_riem = round(exp(2π×15/14.135)) = 787  (138 primes, 149 prime powers)
  lsq_fake = round(exp(2π×15/13.800)) = 923  (160 primes, 172 prime powers)

At j=20:
  lsq_riem = 7262  (~906 primes)
  lsq_fake = round(exp(2π×20/13.800)) = 9988 (~1228 primes)
  → fake MORE negative at j=15 but LESS negative at j=20??

Analysis:
1. Compute P[n,m] off-diagonal norm for riem vs fake at j=15 and j=20 (N=32, small)
2. Compute arch+B+P vs arch+B alone → cancellation metric
3. Compute S[n] = prime sum function — show the PATTERN of contributions
4. Show which primes (787→923) disrupt the cancellation vs which don't
5. Compare |P_od| distribution (Frobenius norm per off-diagonal band) riem vs fake
"""
import math, numpy as np, scipy.linalg, time, sys
from scipy.special import digamma

TWO_PI=2*math.pi; LOGPI=math.log(math.pi)
GAMMA1=14.134725141734693
GAMMA2=13.800

G=30000; N=32  # small N for visualization, large G for accuracy

def build_all(j, gamma, N=32):
    lsq = round(math.exp(TWO_PI*j/gamma))
    L = 2*math.log(math.sqrt(lsq))
    ns = np.arange(-N,N+1); D=len(ns)
    R = TWO_PI*N/L+50.0; grid=np.linspace(-R,R,G)
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
        lp=math.log(p); pk=p
        while pk<=lsq:
            lkpp.append(math.log(pk)); wts_p.append(lp/math.sqrt(pk)); pk*=p
    lkpp=np.array(lkpp); wts_p=np.array(wts_p); K=len(lkpp)
    S=np.zeros(D); Sd=np.zeros(D)
    for k0 in range(0,K,10000):
        k1=min(k0+10000,K)
        ph=TWO_PI*np.outer(ns,lkpp[k0:k1])/L
        S+=np.sin(ph)@wts_p[k0:k1]
        Sd+=np.cos(ph)@(2*wts_p[k0:k1]*(1-lkpp[k0:k1]/L))
    dn=ns[:,None]-ns[None,:]; Si=S[:,None]; Sj=S[None,:]
    with np.errstate(divide='ignore',invalid='ignore'):
        P=np.where(dn!=0,(Sj-Si)/(math.pi*dn),0.0)
    np.fill_diagonal(P,Sd)
    tau=(arch+B-P); tau=(tau+tau.T)/2
    ev0=float(scipy.linalg.eigvalsh(tau,subset_by_index=[0,0])[0])
    return {
        'j':j,'gamma':gamma,'lsq':lsq,'L':L,'D':D,'K':K,
        'arch':arch,'B':B,'P':P,'tau':tau,'S':S,'Sd':Sd,
        'lkpp':lkpp,'wts_p':wts_p,'ev0':ev0,
        'arch00':arch[N,N],'B00':B[N,N],'P00':P[N,N],'tau00':tau[N,N],
    }

print("Building all 4 cases (j=15,20 × gamma1,fake13.8)...\n")
cases = {}
for j in [15, 20]:
    for label, gamma in [('riem', GAMMA1), ('fake13.8', GAMMA2)]:
        key = f"j{j}_{label}"
        t0=time.time()
        r = build_all(j, gamma, N)
        print(f"{key}: lsq={r['lsq']}, K={r['K']} prime powers, ev0={r['ev0']:.4e} ({time.time()-t0:.1f}s)")
        cases[key] = r
print()

# Analysis 1: diagonal comparison (arch+B-P cancellation at n=0)
print("=== 1. Diagonal cancellation at n=0 ===")
for key, r in cases.items():
    diag = r['arch00'] + r['B00'] - r['P00']
    print(f"  {key}: arch={r['arch00']:.4f}, B={r['B00']:.4f}, P={r['P00']:.4f}, tau[0,0]={diag:.6e}")
print()

# Analysis 2: off-diagonal Frobenius norms per diagonal band
print("=== 2. Off-diagonal P norm per band (|k| = offset) ===")
for j in [15, 20]:
    print(f"  j={j}:")
    for offset in [1, 2, 5, 10, 20]:
        riem_key = f"j{j}_riem"
        fake_key = f"j{j}_fake13.8"
        P_r = cases[riem_key]['P']; N_r = cases[riem_key]['P'].shape[0]//2
        P_f = cases[fake_key]['P']

        def band_norm(P, k, n_center):
            idxs = [(n_center+i, n_center+i+k) for i in range(-n_center, n_center-k+1)]
            vals = [P[a,b] for a,b in idxs if 0<=a<len(P) and 0<=b<len(P)]
            return np.sqrt(sum(v**2 for v in vals)) if vals else 0.0

        rn = band_norm(P_r, offset, N)
        fn = band_norm(P_f, offset, N)
        print(f"    |k|={offset}: riem={rn:.4e}, fake={fn:.4e}, ratio={fn/rn:.2f}x")
print()

# Analysis 3: S[n] vector comparison
print("=== 3. S[n] vector: riem vs fake at j=15, n=0..5 ===")
for j in [15, 20]:
    r_r = cases[f"j{j}_riem"]; r_f = cases[f"j{j}_fake13.8"]
    N_c = N  # center index
    print(f"  j={j}:")
    for i in range(6):
        sr = r_r['S'][N_c+i]; sf = r_f['S'][N_c+i]
        print(f"    S[{i}]: riem={sr:.6f}, fake={sf:.6f}, diff={sf-sr:.6f}")
print()

# Analysis 4: P off-diagonal distribution (histogram shape)
print("=== 4. P off-diagonal statistics ===")
for key, r in cases.items():
    P = r['P'].copy(); np.fill_diagonal(P, 0)
    vals = P[P!=0]
    print(f"  {key}: std={np.std(vals):.4e}, max|P|={np.max(np.abs(vals)):.4e}, "
          f"skew(3rd/std^3)={(np.mean(vals**3)/np.std(vals)**3):.3f}")
print()

# Analysis 5: which primes at j=15 are RESPONSIBLE for the discrimination
# Compute |P[0,n+1]| = |S[1]-S[0]|/π incrementally as we add primes
print("=== 5. P[0,1] buildup as primes added (j=15, riem lsq=787) ===")
r_r15 = cases['j15_riem']
lkpp = r_r15['lkpp']; wts = r_r15['wts_p']
L = r_r15['L']; N_c = N
S0 = 0.0; S1 = 0.0
K = len(lkpp)
print(f"  L={L:.4f}, total prime powers: {K}")
print(f"  {'n_pp':>6}  {'p^k':>8}  {'S[0]':>10}  {'S[1]':>10}  {'P[0,1]':>12}  {'change':>10}")

# Process in small batches to show buildup
milestones_lkpp = sorted(set([0, K//5, 2*K//5, 3*K//5, 4*K//5, K-1]))
ns_01 = np.array([0, 1])
ph = TWO_PI*np.outer(ns_01, lkpp)/L  # shape (2, K)
S_cumul = np.cumsum(np.sin(ph)*wts[None,:], axis=1)  # shape (2, K)
for idx in milestones_lkpp:
    s0 = S_cumul[0, idx]; s1 = S_cumul[1, idx]
    p01 = (s1-s0)/math.pi
    pk_val = round(math.exp(lkpp[idx]))
    print(f"  {idx+1:>6}  {pk_val:>8}  {s0:>10.4f}  {s1:>10.4f}  {p01:>12.6e}  ")

# Last few primes around 787 (the boundary)
print(f"\n  Last 10 prime powers (approaching lsq=787):")
for idx in range(max(0,K-10), K):
    s0 = S_cumul[0, idx]; s1 = S_cumul[1, idx]
    p01 = (s1-s0)/math.pi
    pk_val = round(math.exp(lkpp[idx]))
    ds0 = math.sin(TWO_PI*0*lkpp[idx]/L)*wts[idx]
    ds1 = math.sin(TWO_PI*1*lkpp[idx]/L)*wts[idx]
    print(f"  {idx+1:>6}  {pk_val:>8}  {s0:>10.4f}  {s1:>10.4f}  {p01:>12.6e}  dS1={ds1:.4e}")

print("\nDone.")
