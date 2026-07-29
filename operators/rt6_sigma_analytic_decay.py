# Riemann: analytic decay of sigma_off(N). Sieve the OFF-DIAGONAL numerator E[u_t u_t']
# and the DIAGONAL E[u_t^2] SEPARATELY vs N, for a fixed worst-cell shape-mismatched pair.
# If numerator ~ N-independent (const) and diagonal ~ (log N)^3  =>  sigma_off ~ log-family (fork II).
# If numerator grows with a compensating power  =>  power law (fork I).  Let the object decide.
import numpy as np
def sieve_Lambda(M):
    spf=np.zeros(M+1,dtype=np.int64)
    for i in range(2,M+1):
        if spf[i]==0:
            spf[i::i]=np.where(spf[i::i]==0,i,spf[i::i])
    Lam=np.zeros(M+1)
    for n in range(2,M+1):
        p=spf[n]; m=n; k=0
        while m%p==0: m//=p; k+=1
        if m==1: Lam[n]=np.log(p)   # prime power
    return Lam

# worst-cell pair (offsets): t=(1,7,40), t'=(5,11,32)
tA=[1,7,40]; tB=[5,11,32]
allh=tA+tB
M=10_000_000
print("sieving Lambda up to",M,"...")
Lam=sieve_Lambda(M)
dL=Lam-1.0     # delta-Lambda = Lambda - 1
hmax=max(allh)

Ns=[100_000,300_000,1_000_000,3_000_000,10_000_000-hmax]
print(f"{'N':>11} {'numer E[u_t u_t7]':>18} {'diag E[u_t^2]':>14} {'m2':>7} {'R=num/diag':>11}")
rows=[]
for N in Ns:
    idx=np.arange(1,N)   # n from 1
    # numerator: product of dL at the 6 offsets
    num=np.ones(N-1)
    for h in allh: num=num*dL[idx+h]
    numer=num.mean()
    # diagonal for t=(1,7,40): product of dL^2 at the 3 offsets
    dg=np.ones(N-1)
    for h in tA: dg=dg*dL[idx+h]**2
    diag=dg.mean()
    m2=(dL[idx]**2).mean()
    R=numer/diag if diag!=0 else float('nan')
    rows.append((N,numer,diag,m2,R))
    print(f"{N:>11} {numer:>18.5f} {diag:>14.3f} {m2:>7.3f} {R:>11.6f}")

# fit scalings: numerator vs N (power), diagonal vs log N (power in logN), R vs N
lnN=np.log([r[0] for r in rows])
lnlogN=np.log(np.log([r[0] for r in rows]))
numer=np.array([abs(r[1]) for r in rows])
diag=np.array([r[2] for r in rows])
Rabs=np.array([abs(r[4]) for r in rows])
def slope(x,y): return np.polyfit(x,np.log(y),1)[0]
print(f"\nnumerator  |E[u_t u_t']| : power-in-N exponent = {slope(lnN,numer):+.3f}  (0 => N-independent HL const)")
print(f"diagonal   E[u_t^2]      : power-in-logN exponent = {slope(lnlogN,diag):+.2f}  (3 => (logN)^3)")
print(f"R=num/diag              : power-in-N exponent = {slope(lnN,Rabs):+.3f}  ,  power-in-logN = {slope(lnlogN,Rabs):+.2f}")
print("\nREAD: if numerator exponent ~0 and R tracks (logN)^-k => LOG-FAMILY => fork II.")
print("      if numerator carries a compensating +N power => power law => fork I lives.")
import json
json.dump({'pair':[tA,tB],'M':M,'rows':[{'N':r[0],'numer':r[1],'diag':r[2],'m2':r[3],'R':r[4]} for r in rows],
           'numer_powN':float(slope(lnN,numer)),'diag_powLogN':float(slope(lnlogN,diag)),
           'R_powN':float(slope(lnN,Rabs)),'R_powLogN':float(slope(lnlogN,Rabs))},
          open('results/rt6_sigma_analytic_decay.json','w'),indent=2)
print("\nCHECKPOINT -> results/rt6_sigma_analytic_decay.json")
