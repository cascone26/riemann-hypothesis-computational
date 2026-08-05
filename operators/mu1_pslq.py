"""High-precision mu1 (Suzuki's uncomputed constant) via archimedean-only Rayleigh quotient
of the fundamental mode (the near-null, established). For a<0.347 the kernel is purely
archimedean (no primes) -> arbitrary mpmath precision. Fit lambda_a = c log(1/a)+b+O(a),
mu1 = b - (log2pi - psi(2))/2, then PSLQ against archimedean constants."""
import mpmath as mp
mp.mp.dps = 40
PSI14_LOGPI = mp.digamma(mp.mpf(1)/4) - mp.log(mp.pi)   # Omega(0) = -5.372...
PHI = mp.lerchphi(1, 2, mp.mpf(1)/4)                      # Phi(1,2,1/4)
def phi_tail(t):
    return mp.lerchphi(mp.e**(-2*abs(t)), 2, mp.mpf(1)/4)
def g_arch(t):
    at = abs(t)
    boundary = -4*(mp.e**(t/2) + mp.e**(-t/2) - 2)
    core = -(at/2)*PSI14_LOGPI
    tail = -mp.mpf(1)/4*(PHI - mp.e**(-at/2)*phi_tail(t))
    return boundary + core + tail
def lambda_fundamental(a, M=120):
    # Rayleigh quotient Q(phi0)/||phi0'||^2-normalized for phi0 = sin(pi(x+a)/2a), phi0'=(pi/2a)cos(...)
    # Q = int int g(x-y) phi0'(x) phi0'(y); mass = a (sum c^2)=a for single mode. lambda ~ Q/a.
    # double integral via mp.quad
    k = mp.pi/(2*a)
    def inner(x):
        return mp.quad(lambda y: g_arch(x-y)*k*mp.cos(k*(y+a)), [-a, x, a])
    Q = mp.quad(lambda x: k*mp.cos(k*(x+a))*inner(x), [-a, 0, a])
    return Q/a
avals = [mp.mpf('0.02'), mp.mpf('0.04'), mp.mpf('0.06'), mp.mpf('0.09'), mp.mpf('0.13')]
las = []
for a in avals:
    la = lambda_fundamental(a)
    las.append(la)
    print(f"a={float(a):.3f}: lambda_a = {mp.nstr(la, 20)}", flush=True)
# fit c log(1/a) + b + e a  (3-param, high precision)
import numpy as np
A = mp.matrix([[mp.log(1/a), mp.mpf(1), a] for a in avals[:3]])
bb = mp.matrix([las[i] for i in range(3)])
sol = mp.lu_solve(A, bb)
c, b, e = sol[0], sol[1], sol[2]
const = (mp.log(2*mp.pi) - mp.digamma(2))/2
mu1 = b - const
print(f"\nc={mp.nstr(c,8)} (Suzuki=1) b={mp.nstr(b,12)} e={mp.nstr(e,8)}")
print(f"mu1 = b - (log2pi-psi(2))/2 = {mp.nstr(mu1,16)}")
# PSLQ hunt
print("\nPSLQ / closed-form candidates for mu1:")
cands = {'1/4':mp.mpf(1)/4,'gamma/2':mp.euler/2,'(1-gamma)/2':(1-mp.euler)/2,
         'gamma':mp.euler,'log2/2':mp.log(2)/2,'1-log2':1-mp.log(2),
         'zeta2/pi^2 stuff':mp.zeta(2)/mp.pi**2, 'gamma^2/2':mp.euler**2/2}
for nm,v in cands.items():
    print(f"  {nm}={mp.nstr(v,12)} diff={mp.nstr(abs(v-mu1),6)}")
try:
    rel = mp.pslq([mu1, mp.mpf(1), mp.euler, mp.log(2), mp.log(mp.pi)], maxcoeff=10**6)
    print(f"PSLQ [mu1,1,gamma,log2,logpi]: {rel}")
except Exception as ex:
    print("pslq:", ex)
print("DONE")
