# Claim: phi(0+) = (1/2)(gamma_E + log 2pi)  [Lerch asymptotic: S(t) ~ -(1/2)(log 2t + psi(1/4) + gamma);
#   Omega0 = psi(1/4) - log pi  =>  psi(1/4) cancels, gamma survives, then Term-B's -2*phi(0) kills Term-A's +gamma]
# =>  Q[k,k] = log(kb/(2pi)) - (wp*(2a-log2)/a)*cos(kb*log2) + O(1/k)
# with the leading 1/k piece containing -1/(2k) from Si(k*pi) = pi/2 + 1/(k*pi).
import numpy as np, math
from scipy.special import digamma

A=0.549; b=math.pi/(2*A); log2=math.log(2); wp=log2/math.sqrt(2)
gam=0.5772156649015329
phi0_pred=0.5*(gam+math.log(2*math.pi))
print(f"phi(0+) closed form = (1/2)(gamma+log 2pi) = {phi0_pred:.6f}")
print(f"  (measured H'-0.5log t at t=1e-2 was 1.1900; pred - O(t): consistent)")
amp=wp*(2*A-log2)/A
c0=-math.log(2*math.pi)
print(f"c0 = -log(2pi) = {c0:.6f}   (fit was -1.8527; diff 0.015 = absorbed 1/k terms)")
print(f"amp = wp(2a-log2)/a = {amp:.6f}  (fit 0.3608)")

Q=np.load('/Users/scones/tmp/Q51.npy'); ODD=list(range(1,52,2))
print("\n  k    Q[k,k]     closed-form    resid      resid*k")
for k in ODD:
    if k<15: continue
    q=Q[ODD.index(k),ODD.index(k)]
    pred=math.log(k*b/(2*math.pi))-amp*math.cos(k*b*log2)
    r=q-pred
    print(f" {k:3d}  {q:.5f}   {pred:.5f}   {r:+.5f}   {r*k:+.4f}")
