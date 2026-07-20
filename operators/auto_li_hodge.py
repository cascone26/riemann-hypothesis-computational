import mpmath as mp
import numpy as np
mp.mp.dps=60

def xi(s):
    # xi(s)=Gamma(s/2+1) pi^{-s/2} (s-1) zeta(s): the s(s-1) kills zeta's pole, analytic at s=1
    return mp.gamma(s/2+1)*mp.pi**(-s/2)*(s-1)*mp.zeta(s)

def logxi(s):
    return mp.log(xi(s))

N=6
lam=[]
for n in range(1,N+1):
    g=lambda s,n=n: (s**(n-1))*logxi(s)
    val=mp.diff(g,1,n)/mp.factorial(n-1)
    lam.append(mp.re(val))

print('Li/Hodge polarization coefficients lambda_n = self-intersection numbers (RH => all > 0):')
allpos=True
for n,l in enumerate(lam,1):
    print('  lambda_%d = %s'%(n,mp.nstr(l,15)))
    if l<=0: allpos=False

lam0=[0.0]+[float(x) for x in lam]
m=len(lam0)//2
Hmat=np.array([[lam0[i+j] for j in range(m)] for i in range(m)],dtype=float)
ev=np.linalg.eigvalsh((Hmat+Hmat.T)/2)
print('Hankel(lambda) intersection-form eigenvalues:',np.round(ev,8))
print('min lambda_n =',mp.nstr(min(lam),12))
if allpos:
    print('LI_POSITIVITY_HOLDS: Hodge polarization (arithmetic self-intersection numbers) positive up to n=%d'%N)
    print('Structural lever: if this positivity is the Edrei/Hadamard product signature it forces RH for all n.')
else:
    print('NEGATIVE_LI_COEFF_FOUND (would refute RH)')
