#!/usr/bin/env python3
"""§469 strengthening: PROVE the prime's bulk Jacobi imprint is the FIRST-ORDER Euler term.
|L_p(1/2-is)|^2 = 1/(1 - 2 p^{-1/2} cos(s ln p) + 1/p) = 1 + 2 p^{-1/2} cos(s ln p) + O(1/p).
So to leading order the semilocal measure is dm_0 * (1 + 2 eps cos(s ln p)) with eps = p^{-1/2}.
TEST: does the FULL |L_p|^2 imprint delta b_n equal the FIRST-ORDER-truncated-measure imprint to O(eps)?
If yes, the p^{-1/2} scaling is exactly first-order perturbation in the Euler factor's leading term, with a
FIXED kernel (independent of p): delta b_n(p) = eps * kappa_n + O(eps^2), eps = p^{-1/2}.
Also validates linearity: delta b at half eps should be ~half.
"""
import numpy as np
from scipy.special import loggamma

def w0(s):  return np.exp(2.0*np.real(loggamma(0.25+0.5j*s)))
def wfull(s,p): return w0(s)/(1.0 - 2.0*p**-0.5*np.cos(s*np.log(p)) + 1.0/p)      # exact |L_p|^2 * dm0
def wfo(s,p):   return w0(s)*(1.0 + 2.0*p**-0.5*np.cos(s*np.log(p)))              # first-order in p^-1/2
def weps(s,eps,lp): return w0(s)*(1.0 + 2.0*eps*np.cos(s*lp))                     # generic small-eps probe

def jac(weight,L=45.0,M=200000,K=16):
    s=np.linspace(-L,L,M); ds=s[1]-s[0]; w=weight(s)*ds
    Mt=w.sum(); qp=np.zeros_like(s); q=np.sqrt(w/Mt); b=np.zeros(K); bb=0.0
    for k in range(K):
        v=s*q; a=np.dot(q,v); v=v-a*q-bb*qp; v=v-q*np.dot(q,v)-qp*np.dot(qp,v)
        bb=np.linalg.norm(v)
        if k+1<K: b[k+1]=bb
        qp=q; q=v/bb
    return b

b0=jac(w0)
print("TEST 1 — full |L_p|^2 imprint vs first-order-truncated-measure imprint (should agree to O(1/p)):")
print(f"{'p':>4} {'n':>3} {'db_full':>10} {'db_1storder':>12} {'diff':>10}   (diff should be O(1/p), small)")
for p in [2,3,5,7]:
    bf=jac(lambda s:wfull(s,p)); b1=jac(lambda s:wfo(s,p))
    for n in [1,2,3,4]:
        print(f"{p:>4} {n:>3} {bf[n]-b0[n]:>10.4f} {b1[n]-b0[n]:>12.4f} {(bf[n]-b1[n]):>10.4f}")
    print()

print("TEST 2 — fixed kernel: delta b_n / p^{-1/2} should be ~p-independent (=kappa_n) at first order:")
lp=np.log(2.0)  # use a fixed frequency? No -- frequency IS s*ln p, so kernel depends on ln p too.")
# The kernel depends on the frequency ln p. So test linearity in eps at FIXED frequency instead:
print("  (kernel depends on frequency ln p, so instead verify LINEARITY in eps at fixed frequency ln 2)")
print(f"  {'eps':>8} {'db_1/eps':>10} {'db_2/eps':>10} {'db_3/eps':>10}")
kappa=None
for eps in [0.02,0.05,0.1,0.2]:
    be=jac(lambda s:weps(s,eps,np.log(2.0)))
    row=[ (be[n]-b0[n])/eps for n in [1,2,3] ]
    if kappa is None: kappa=row
    print(f"  {eps:>8.3f} {row[0]:>10.4f} {row[1]:>10.4f} {row[2]:>10.4f}")
print("  -> if db_n/eps is ~constant across eps, delta b is LINEAR in eps (first-order confirmed);")
print("     kappa_n (the fixed kernel at frequency ln2) = the small-eps limit above.")
print(f"\n  Cross-check: at eps=p^{{-1/2}} with p=2 (eps={2**-0.5:.4f}) the first-order prediction kappa*eps")
print(f"  db_1 pred = {kappa[0]*2**-0.5:.4f}  vs full |L_2|^2 db_1 = {jac(lambda s:wfull(s,2))[1]-b0[1]:.4f}")
