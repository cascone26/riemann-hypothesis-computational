#!/usr/bin/env python3
"""§430: re-test the ρ_edge asymptotic trend with the EXTENDED instrument (150 zeros up to γ=318, from mpmath).
Now high bands have enough zeros above them for the explicit formula to close (valid residual). Question stands:
does ρ_edge=||prime||op/archmargin →1 (op-bound asymptotically tight) or level off >1 (irreducible deficit) as γ→∞?"""
import numpy as np, l2_fejer as LF
from scipy.special import digamma
Z=np.load('zeros150.npy')
LF.ZEROS=Z
# extend the integration grid to cover high bands (sinc^2 slow decay); rebuild OMEGA on it
LF.RGRID=np.arange(-320.0,320.0,0.006)
LF.OMEGA=np.real(digamma(0.25+0.5j*LF.RGRID))-np.log(np.pi)
def nz(band): return int(np.sum((Z>=band[0])&(Z<=band[1])))
print(f"§430 ρ_edge extended (150 zeros, γ<{Z.max():.0f}). W={LF.W}, prime cutoff n<{int(np.exp(4*LF.W))}.")
print(f"{'band':>14} {'γmid':>5} {'Nz':>3} {'Nabove':>6} {'resid':>9} {'archmarg':>9} {'||pr||op':>9} {'ρ_edge':>7} {'VALID':>6}")
for band in [(14,44),(40,70),(70,100),(100,130),(130,160),(160,190),(200,230)]:
    N=nz(band)
    if N<3: continue
    nabove=int(np.sum(Z>band[1]))
    s=(band[1]-band[0])/(N-1)   # ~N probes
    cs=LF.probes(band,s)
    f=LF.build(cs)
    resid=np.max(np.abs(f['full']-f['zero']))
    am=float(np.linalg.eigvalsh(f['archpole'])[0])
    ev=np.linalg.eigvalsh(f['prime']); opn=max(abs(ev[0]),abs(ev[-1]))
    print(f"{str(band):>14} {(band[0]+band[1])//2:5d} {len(cs):3d} {nabove:6d} {resid:9.1e} {am:+9.3e} "
          f"{opn:9.3e} {opn/am if am>1e-15 else np.inf:7.3f} {'YES' if resid<1e-5 else 'no':>6}",flush=True)
print("""
reading: over the now-VALID high bands, if ρ_edge keeps falling toward 1 => op-bound asymptotically tight
(archmargin~log γ catches the bounded prime norm). If it levels at ρ*>1 => irreducible cancellation deficit ρ*.
""")
