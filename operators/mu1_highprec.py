import numpy as np, mpmath as mp
from suzuki_screw_lambda_a import screw_g_no_primes, lambda_a
import suzuki_screw_lambda_a as S
mp.mp.dps=25
C0=float(mp.euler); A=0.5*(float(mp.log(2*mp.pi))+C0-1)
const=float(-mp.log(2*mp.pi)+mp.digamma(2)-1)
S.screw_g=screw_g_no_primes
print("high-precision mu1 (arch-only, small a, extrapolate a->0):",flush=True)
print("A=%.6f  const=-2A-1=%.6f"%(A,const),flush=True)
aa=[0.05,0.07,0.09,0.12,0.15,0.18,0.22,0.27]
mus=[]
for a in aa:
    M=min(7000,max(4000,int(300.0/a)))
    la=lambda_a(a,M=M,K=90)
    mu=(la-np.log(1.0/a))-const
    mus.append((a,mu)); print("  a=%.3f  lam=%.6f  mu1_est=%.6f"%(a,la,mu),flush=True)
a_arr=np.array([x[0] for x in mus]); m_arr=np.array([x[1] for x in mus])
for deg,lab in [(1,"linear"),(2,"quad")]:
    coef=np.polyfit(a_arr,m_arr,deg)
    print("  %s extrap mu1(a->0)=%.5f (O(a) coef=%.3f)"%(lab,coef[-1],coef[-2]),flush=True)
mu1=np.polyfit(a_arr,m_arr,2)[-1]
print("\nEXTRACTED mu1 ~= %.5f"%mu1,flush=True)
L2=float(mp.log(2)); Lpi=float(mp.log(mp.pi))
cands={"gamma/2":C0/2,"1-gamma":1-C0,"gamma":C0,"log2-1/2":L2-0.5,"1-log2":1-L2,
       "log2-gamma/2":L2-C0/2,"gamma^2":C0**2,"gamma-1/3":C0-1/3,"1/2-log2+gamma/2... ":0,
       "3/2-log(pi)":1.5-Lpi,"1-gamma/2-log2/2... ":0,"gamma-log2... ":C0-L2+0.4,
       "log(4/pi)-1/2":float(mp.log(4/mp.pi))-0.5,"(gamma+log2)/... ":0,"1-gamma+... ":0}
print("closed-form candidates near mu1:",flush=True)
for k,v in sorted(cands.items(),key=lambda kv:abs(kv[1]-mu1)):
    if v!=0: print("  %-20s=%.5f  |diff|=%.5f"%(k,v,abs(v-mu1)),flush=True)
