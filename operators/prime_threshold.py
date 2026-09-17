import numpy as np, mpmath as mp
from suzuki_screw_lambda_a import screw_g, screw_g_no_primes, lambda_a
import suzuki_screw_lambda_a as S
mp.mp.dps=20
thr=0.5*np.log(2)
print(f"λ_a across the FIRST prime threshold a=½log2={thr:.5f} (n=2 term turns on):",flush=True)
print(f"{'a':>7} {'λ_full':>12} {'λ_arch':>12} {'Δλ=prime corr':>14} {'(2a-log2)':>10}",flush=True)
orig=S.screw_g
for a in [0.30,0.33,0.345,0.3466,0.35,0.36,0.38,0.42,0.47]:
    M=max(6000,int(500/a)); K=80
    S.screw_g=orig;              lf=lambda_a(a,M=M,K=K)
    S.screw_g=screw_g_no_primes; la=lambda_a(a,M=M,K=K)
    S.screw_g=orig
    d=2*a-np.log(2)
    print(f"{a:7.4f} {lf:12.6f} {la:12.6f} {lf-la:14.3e} {d:10.4f}",flush=True)
print("\nΔλ should be 0 below threshold (2a<log2), turn on smoothly above (R7: 'passes smoothly').",flush=True)
print("Characterize onset: Δλ ~ C*(2a-log2)^p for a just above threshold.",flush=True)
