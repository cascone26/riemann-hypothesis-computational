import numpy as np
from suzuki_screw_lambda_a import screw_g, screw_g_no_primes, lambda_a
import suzuki_screw_lambda_a as S
thr=0.5*np.log(2)
print(f"FINE prime-onset scan, clean range (threshold={thr:.5f}, before λ_arch<0):",flush=True)
print(f"{'a':>7} {'delta=2a-log2':>13} {'Δλ':>13} {'Δλ/δ^2':>11} {'Δλ/δ^3':>11}",flush=True)
orig=S.screw_g; rows=[]
for a in [0.349,0.352,0.355,0.358,0.361,0.364,0.367]:
    M=max(9000,int(700/a)); K=100
    S.screw_g=orig;              lf=lambda_a(a,M=M,K=K)
    S.screw_g=screw_g_no_primes; la=lambda_a(a,M=M,K=K)
    S.screw_g=orig
    d=2*a-np.log(2); dl=lf-la; rows.append((d,dl))
    print(f"{a:7.4f} {d:13.5f} {dl:13.4e} {dl/d**2:11.4f} {dl/d**3:11.2f}",flush=True)
# log-log fit power
d=np.array([r[0] for r in rows]); dl=np.array([-r[1] for r in rows])  # -Δλ (positive)
ok=dl>0
p=np.polyfit(np.log(d[ok]),np.log(dl[ok]),1)
print(f"\nlog-log fit: -Δλ ~ δ^{p[0]:.3f} * exp({p[1]:.3f})   (power p; 2=quadratic, 3=cubic)",flush=True)
print(f"coefficient at δ^{round(p[0])}: {np.exp(p[1]):.4f}",flush=True)
