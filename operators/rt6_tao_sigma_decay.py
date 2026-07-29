"""Tao: discriminate sigma_off(N) POWER vs LOG, with sampling-floor null-subtraction.
Fixed sector s=42 (n=36 accessible). Chunked Gram + per-entry variance -> analytic
sampling floor -> sigma_arith = sqrt(sigma_off^2 - sigma_floor^2). Pre-registered
discriminator: local log-log slope of sigma_arith(N). CONSTANT => power (fork I open);
DRIFTING toward 0 => log (fork II, per-s only)."""
import json, numpy as np
from math import log
def sieve_Lambda(M):
    spf=np.zeros(M+1,dtype=np.int64)
    for i in range(2,M+1):
        if spf[i]==0:
            for j in range(i,M+1,i):
                if spf[j]==0: spf[j]=i
    L=np.zeros(M+1)
    for n in range(2,M+1):
        p=int(spf[n]); m=n
        while m%p==0: m//=p
        if m==1: L[n]=log(p)
    return L
def triples_at_sum(s,lo=1):
    return [(a,b,s-a-b) for a in range(lo,s) for b in range(a+1,s) if s-a-b>b]
def admissible(t):
    for p in (3,5,7):
        if len(set(x%p for x in t))==p: return False
    return True

s=42
T=[t for t in triples_at_sum(s) if admissible(t)]; n=len(T)
NMAX=8_000_000; MAXO=60
print(f"sieving to {NMAX+MAXO} ...")
Lam=sieve_Lambda(NMAX+MAXO+5); dL=Lam-1.0
off=np.array([[1 if i!=j else 0 for j in range(n)] for i in range(n)],dtype=bool)

def measure(N,chunk=250_000):
    G=np.zeros((n,n)); H=np.zeros((n,n))
    done=0
    while done<N:
        c=min(chunk,N-done); base=np.arange(done,done+c)
        Vc=np.array([dL[base+t[0]]*dL[base+t[1]]*dL[base+t[2]] for t in T])  # n x c
        G+=Vc@Vc.T; H+=(Vc**2)@(Vc**2).T; done+=c
    G/=N; H/=N
    d=np.sqrt(np.diag(G)); D=np.outer(d,d)
    R=G/D
    var=(H-G**2)                       # per-entry Var(u_t u_t')
    floor=np.sqrt(np.maximum(var,0)/N)/D  # standard error of each R entry
    so=np.sqrt((R[off]**2).mean())
    sf=np.sqrt((floor[off]**2).mean())
    sa=np.sqrt(max(so**2-sf**2,0.0))
    return so,sf,sa

Ns=[100_000,200_000,400_000,800_000,1_600_000,3_200_000,6_400_000]
rows=[]
print(f"{'N':>10} {'sig_off':>9} {'sig_floor':>10} {'sig_arith':>10}")
for N in Ns:
    so,sf,sa=measure(N); rows.append((N,so,sf,sa))
    print(f"{N:>10} {so:>9.5f} {sf:>10.5f} {sa:>10.5f}")
# local log-log slopes of sigma_arith
print("\nlocal log-log slope d(log sig_arith)/d(log N):")
sl_arith=[]; sl_off=[]
for k in range(1,len(rows)):
    N0,so0,sf0,sa0=rows[k-1]; N1,so1,sf1,sa1=rows[k]
    sla=(log(sa1)-log(sa0))/(log(N1)-log(N0))
    slo=(log(so1)-log(so0))/(log(N1)-log(N0))
    sl_arith.append(sla); sl_off.append(slo)
    print(f"  {N0:>8}->{N1:<8}  arith={sla:+.3f}   raw_off={slo:+.3f}")
# discriminator: is arith-slope constant (power) or trending to 0 (log)?
drift=sl_arith[-1]-sl_arith[0]
print(f"\narith-slope first={sl_arith[0]:+.3f} last={sl_arith[-1]:+.3f} drift={drift:+.3f}")
print("POWER if |drift|~0 & slope stable; LOG if slope drifts toward 0 as N grows")
json.dump({"object":"Tao sigma_off decay discriminator, null-subtracted",
 "version":"rt6-tao-sigma-decay-2026-07-28","sector_s":s,"n":n,
 "ladder":[{"N":r[0],"sig_off":r[1],"sig_floor":r[2],"sig_arith":r[3]} for r in rows],
 "loglog_slope_arith":sl_arith,"loglog_slope_rawoff":sl_off,"arith_slope_drift":drift},
 open("results/rt6_tao_sigma_decay.json","w"),indent=2)
print("wrote results/rt6_tao_sigma_decay.json")
