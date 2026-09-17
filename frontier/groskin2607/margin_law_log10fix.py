import sys, math
sys.path.insert(0,'.')
from arb_ldlt_certify import build_arb_tau
from flint import arb, ctx
import mpmath as mp, numpy as np

def spectrum(c,N,prec):
    ctx.prec=prec
    A,DIM=build_arb_tau(c,N,prec)
    mp.mp.dps=int(prec*0.29)
    Amp=mp.matrix(DIM,DIM)
    for i in range(DIM):
        for j in range(DIM): Amp[i,j]=mp.mpf(A[i,j].str(int(prec*0.28),radius=False))
    E=sorted(float(mp.log10(abs(e))) for e in mp.eigsy(Amp,eigvals_only=True))
    return E

# CORRECTED: fit log10(min-eig) vs log10(c), not vs ln(c) -- the original
# margin_law.py used L=math.log(c) (natural log) as the x-axis while the
# y-axis was already log10(min-eig). The fitted slope from that script is
# d(log10 minEig)/d(ln c), NOT the exponent alpha in min-eig ~ c^(-alpha).
# True alpha requires log10-log10 (or ln-ln) consistently on both axes.
print("CORRECTED B) true MIN-EIG vs log10(c) (N=8), consistent log10-log10 basis:")
L10s=[]; m2=[]
for c in [13,30,100,300,1000]:
    E=spectrum(c,8,1400); L10=math.log10(c); L10s.append(L10); m2.append(E[0])
    print(f"  c={c:5} log10(c)={L10:.4f} log10(min-eig)={E[0]:+.4f}")
lin_correct=np.polyfit(L10s,m2,1)
alpha=-lin_correct[0]
print(f"  log10(min-eig) = {lin_correct[0]:.4f}*log10(c) {lin_correct[1]:+.4f}")
print(f"  => min-eig ~ c^(-alpha) with alpha = {alpha:.4f}")
print(f"  resid: {np.std(np.polyval(lin_correct,L10s)-m2):.4f}")

# Sanity cross-check: convert the OLD natural-log slope to the same alpha
# convention, to confirm the two methods agree once units are fixed.
print("\nCross-check against natural-log fit converted to base-10 exponent:")
Ls=[]; m3=[]
for c in [13,30,100,300,1000]:
    E=spectrum(c,8,1400); Ls.append(math.log(c)); m3.append(E[0])
lin_natural=np.polyfit(Ls,m3,1)
alpha_from_natural = -lin_natural[0]*math.log(10)
print(f"  old-style slope wrt ln(c): {lin_natural[0]:.4f}")
print(f"  converted alpha = -slope*ln(10) = {alpha_from_natural:.4f}")
print(f"  (should match corrected alpha={alpha:.4f} above -- confirms conversion, not a new independent measurement)")
