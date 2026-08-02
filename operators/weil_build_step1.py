import numpy as np
import mpmath as mp

# WEIL FORM BUILD, STEP 1: calibrate the Riemann-Weil explicit formula.
# Riemann-Weil: for even test Phi(t) with phi(u)=(1/2pi) int Phi(t) e^{-iut} dt,
#   sum_gamma Phi(gamma) = (1/2pi) int Phi(t) Omega(t) dt + [Phi(i/2)+Phi(-i/2)] - 2 sum_n Lambda(n)/sqrt(n) phi(log n)
# with Omega(t) = Re[psi(1/4 + i t/2)] - log(pi).
# CALIBRATE: LHS from 10^4 real zeros (both +-gamma) vs RHS; they must match if RH-side is consistent
# and the normalization is right. This is verify-before-ship on the formula before building the form.

mp.mp.dps = 15
Z = np.load("zeros10k.npy")   # positive gammas
gam = Z

def test_gaussian(c, w):
    # Phi(t) = exp(-w^2 t^2 / 2) * cos(c t)  (even);  phi(u) = (1/(2 sqrt(2pi) w)) [gauss(u-c)+gauss(u+c)]
    Phi = lambda t: np.exp(-0.5 * w**2 * t**2) * np.cos(c * t)
    phi = lambda u: (1.0/(2*w*np.sqrt(2*np.pi))) * (np.exp(-0.5*((u-c)/w)**2) + np.exp(-0.5*((u+c)/w)**2))
    return Phi, phi

def Omega(t):
    return float(mp.re(mp.digamma(mp.mpf(0.25) + 0.5j*t))) - np.log(np.pi)

def explicit_formula_check(c, w, NM=2000):
    Phi, phi = test_gaussian(c, w)
    # LHS: sum over +- gamma (Phi even so = 2 sum_{gamma>0})
    LHS = 2.0*np.sum(Phi(gam))
    # arch: (1/2pi) int Phi(t) Omega(t) dt  (Phi decays, integrate on [-R,R])
    R = 12.0/w
    tt = np.linspace(-R, R, 8000)
    Om = np.array([Omega(t) for t in tt])
    arch = np.trapezoid(Phi(tt)*Om, tt)/(2*np.pi)
    # pole terms Phi(i/2)+Phi(-i/2) = 2 Re Phi(i/2); Phi(i/2)=exp(-0.5 w^2 (i/2)^2) cos(c i/2)
    #   = exp(w^2/8) cosh(c/2)
    pole = 2.0*np.exp(w**2/8.0)*np.cosh(c/2.0)
    # prime: 2 sum_n Lambda(n)/sqrt(n) phi(log n)
    L = np.zeros(NM+1)
    for p in range(2, NM+1):
        if all(p % d for d in range(2, int(p**0.5)+1)):
            pk=p
            while pk<=NM:
                L[pk]=np.log(p); pk*=p
    ns = np.arange(2, NM+1)
    prime = 2.0*np.sum(L[2:NM+1]/np.sqrt(ns)*phi(np.log(ns)))
    RHS = arch + pole - prime
    return LHS, RHS, arch, pole, prime

print("Riemann-Weil calibration: LHS (zero sum) vs RHS (arch+pole-prime)")
print(f"{'c':>5} {'w':>5} {'LHS':>12} {'RHS':>12} {'diff':>10} {'arch':>9} {'pole':>9} {'prime':>9}")
for c, w in [(0.0,0.5),(1.0,0.4),(2.0,0.4),(0.5,0.3),(3.0,0.5)]:
    LHS,RHS,ar,po,pr = explicit_formula_check(c,w)
    print(f"{c:>5.1f} {w:>5.2f} {LHS:>12.4f} {RHS:>12.4f} {LHS-RHS:>10.4f} {ar:>9.3f} {po:>9.3f} {pr:>9.3f}", flush=True)
print("\nif diff ~ 0: normalization correct, ready to build the quadratic form arch+pole-prime.")
print("(residual = finite zero-count tail + prime cutoff; should shrink with more zeros/primes)")
