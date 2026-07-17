#!/usr/bin/env python3
"""
Operator component analysis: arch vs B vs P contributions to ev[0].

Answers the key question: WHY does ev[0] → -∞ as N grows for ALL lsq?
Which term (arch, B, or P) drives the growing negativity?

For each N: compute ev[0] of arch alone, B alone, arch+B, and tau=arch+B-P.
Also report arch[0,0], B[0,0], P[0,0] (diagonal at n=0 mode).
"""
import math, numpy as np, scipy.linalg, time, sys
from scipy.special import digamma

TWO_PI=2*math.pi; LOGPI=math.log(math.pi)
LSQ=7262; L=math.log(LSQ); K_BATCH=3000; G=30000

print(f"Component analysis: tau = arch + B - P  (lsq={LSQ}, L={L:.4f})")
print(f"  G={G} (fixed for all N)")
print()

def analyze(N):
    ns=np.arange(-N,N+1); D=len(ns)
    R=TWO_PI*N/L+50.0
    grid=np.linspace(-R,R,G)
    h=grid[1]-grid[0]; w=np.full(G,h); w[0]=w[-1]=0.5*h
    Om=np.real(digamma(0.25+0.5j*grid))-LOGPI
    c=TWO_PI*ns[:,None]/L; d=c-grid[None,:]
    sm=np.abs(d)<1e-12
    E=np.where(sm,np.sqrt(L)+0j,(1/np.sqrt(L))*np.exp(0.5j*d*L)*2*np.sin(0.5*d*L)/d)
    arch=(E*(w*Om)@E.conj().T).real/TWO_PI
    sh=np.sinh(L/4)**2; ni=ns[:,None]; nj=ns[None,:]
    B=32*L*sh*(L**2-16*math.pi**2*ni*nj)/((L**2+16*math.pi**2*ni**2)*(L**2+16*math.pi**2*nj**2))
    sieve=np.ones(LSQ+1,dtype=bool); sieve[0]=sieve[1]=False
    for i in range(2,int(LSQ**0.5)+1):
        if sieve[i]: sieve[i*i::i]=False
    lkpp,wts_p=[],[]
    for p in range(2,LSQ+1):
        if not sieve[p]: continue
        pk=p
        while pk<=LSQ:
            lkpp.append(math.log(pk)); wts_p.append(math.log(p)/math.sqrt(pk)); pk*=p
    lkpp=np.array(lkpp); wts_p=np.array(wts_p)
    K=len(lkpp); S=np.zeros(D); Sd=np.zeros(D)
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
    arch_s=(arch+arch.T)/2; B_s=(B+B.T)/2

    idx0 = N  # index of n=0 mode
    return {
        "N": N, "D": D,
        "arch00": arch[idx0,idx0], "B00": B[idx0,idx0], "P00": P[idx0,idx0],
        "arch_diag_min": float(np.min(np.diag(arch))),
        "arch_diag_max": float(np.max(np.diag(arch))),
        "B_diag_min": float(np.min(np.diag(B))),
        "B_diag_max": float(np.max(np.diag(B))),
        "P_diag_min": float(np.min(np.diag(P))),
        "P_diag_max": float(np.max(np.diag(P))),
        "ev0_arch": float(scipy.linalg.eigvalsh(arch_s, subset_by_index=[0,0])[0]),
        "ev0_archB": float(scipy.linalg.eigvalsh((arch_s+B_s)/2*2, subset_by_index=[0,0])[0]),
        "ev0_tau": float(scipy.linalg.eigvalsh(tau, subset_by_index=[0,0])[0]),
        "arch_offdiag_norm": float(np.linalg.norm(arch - np.diag(np.diag(arch)), 'fro')),
        "P_offdiag_norm": float(np.linalg.norm(P - np.diag(np.diag(P)), 'fro')),
    }

N_vals=[64,128,256,512,768,1024]
print(f"{'N':>6}  {'arch00':>10}  {'B00':>10}  {'P00':>10}  {'tau00':>10}  "
      f"{'ev0(arch)':>12}  {'ev0(a+B)':>12}  {'ev0(tau)':>12}  "
      f"{'|arch_od|':>12}  {'|P_od|':>12}")
print("-"*130)

for N in N_vals:
    t0=time.time()
    r=analyze(N)
    print(f"{N:>6}  {r['arch00']:>10.3f}  {r['B00']:>10.3f}  {r['P00']:>10.3f}  "
          f"{r['arch00']+r['B00']-r['P00']:>10.3f}  "
          f"{r['ev0_arch']:>12.4e}  {r['ev0_archB']:>12.4e}  {r['ev0_tau']:>12.4e}  "
          f"{r['arch_offdiag_norm']:>12.4f}  {r['P_offdiag_norm']:>12.4f}  ({time.time()-t0:.1f}s)")
    sys.stdout.flush()

print()
print("Legend: tau00 = arch[0,0]+B[0,0]-P[0,0] (diagonal at n=0 mode)")
print("  arch_od = Frobenius norm of off-diagonal arch (grows with N → drives negative ev)")
print("  P_od = Frobenius norm of off-diagonal P (prime sum contribution)")
