import numpy as np, time
from suzuki_screw_lambda_a import lambda_a
t0=time.time()
print("R7 cross-check: is λ_a>0, decreasing, smooth through 1st prime threshold a=½log2=%.4f?"%(0.5*np.log(2)))
print(f"{'a':>7} {'λ_a':>14} {'log λ_a':>10}")
res=[]
for a in [0.20,0.28,0.347,0.40,0.50,0.60,0.70,0.80,0.90]:
    M=max(8000,int(600/a)); la=lambda_a(a,M=M,K=80)
    res.append((a,la)); print(f"{a:7.3f} {la:14.6e} {np.log(la) if la>0 else float('nan'):10.4f}",flush=True)
# fit log λ vs a and vs a^2 to test exponential vs superexponential
a=np.array([r[0] for r in res]); lam=np.array([r[1] for r in res])
ok=lam>0; a,lam=a[ok],lam[ok]; L=np.log(lam)
# exp fit L=c0+c1*a ; super fit L=c0+c1*a+c2*a^2
c_e=np.polyfit(a,L,1); c_s=np.polyfit(a,L,2)
re=L-np.polyval(c_e,a); rs=L-np.polyval(c_s,a)
print(f"\nexp fit  logλ={c_e[1]:.3f}+({c_e[0]:.3f})a   resid RMS={np.sqrt(np.mean(re**2)):.4f}")
print(f"super fit logλ={c_s[2]:.3f}+({c_s[1]:.3f})a+({c_s[0]:.3f})a²  resid RMS={np.sqrt(np.mean(rs**2)):.4f}")
print(f"quadratic coeff c2={c_s[0]:.3f} (negative & curving down ⇒ faster-than-exponential = superexponential)")
print(f"all positive: {np.all(lam>0)}   monotone decreasing: {np.all(np.diff(lam)<0)}   [{time.time()-t0:.0f}s]")
