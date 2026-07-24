# Phase 2a: lambda_a on a FINE a-grid through prime thresholds a_n=(log n)/2. Continuous Suzuki object
# (domain grows with a). Does lambda_a DIP between thresholds + get RESCUED at each prime (relay), or stay pinned?
# Trust only a<=~0.9 (through prime 7); beyond = discretization leak (per suzuki_screw_lambda_a caveat).
import numpy as np, sys, os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),'..'))
import suzuki_screw_lambda_a as S
thr={2:0.3466,3:0.5493,4:0.6931,5:0.8047,7:0.9730}   # a_n=(log n)/2 for prime powers 2,3,4,5,7
print("a       lambda_a       nearest prime-threshold")
agrid=np.round(np.arange(0.20,0.95,0.02),3)
prev=None
for a in agrid:
    la=S.lambda_a(float(a), M=max(6000,int(500/a)), K=70)
    # tag if a just crossed a threshold
    tags=[f"a_{n}={t}" for n,t in thr.items() if abs(a-t)<0.011]
    d=f"  d(lam)/da~{(la-prev)/0.02:+.2f}" if prev is not None else ""
    print(f"{a:.2f}   {la:+.6f}   {' '.join(tags)}{d}", flush=True)
    prev=la
print("\nDips-between + rescue-at-threshold => RELAY (induction visible). Monotone/pinned => different structure.")
