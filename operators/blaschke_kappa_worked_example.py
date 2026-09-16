"""
riemann11 -- CORRECTED: the ORDINARY Blaschke factor B_a(p)=(1-pa)^{-1}(a-p) is a self-map of the
ball (zero at p=a, POLE at 1/a which is OUTSIDE |p|<1) -- an ordinary bounded/Pick function, kappa=0.
That's what test 1 measured (kappa=0 confirmed my formula is internally consistent, just applied to
the wrong object -- an honest, useful negative-control check).

The Krein-Langer "kappa=1 from a pole" object (matching the project's own cited fact, RH_CORE_
LOCALIZATION.md sec 8: "a Nevanlinna function with a single simple pole has kappa=1") is the
RECIPROCAL S(p) = 1/B_a(p) = (1-ap)/(a-p), which has an interior POLE at p=a (|a|<1) -- unbounded
inside the domain, the genuinely indefinite object.

S(p) = (1-ap)/(a-p) = sum_k s_k p^k for |p|<|a|:  s_0=1/a,  s_k = a^-(k+1) - a^-(k-1) for k>=1.
"""
import numpy as np

def recip_blaschke_coeffs(a, K):
    s = np.zeros(K+1)
    s[0] = 1.0/a
    for k in range(1, K+1):
        s[k] = a**(-(k+1)) - a**(-(k-1))
    return s

def verify_power_series(a, K):
    s = recip_blaschke_coeffs(a, K)
    ok = True
    for p in [0.2*a, -0.3*a, 0.4*a]:
        exact = (1 - a*p)/(a - p)
        series = sum(s[k]*p**k for k in range(K+1))
        if abs(exact-series) > 1e-6:
            print("  MISMATCH", a, p, exact, series); ok=False
    return s, ok

def kappa_of_series(s, MU_MAX):
    kappas = []
    for mu in range(0, MU_MAX+1):
        A = np.zeros((mu+1, mu+1))
        for n in range(mu+1):
            for m in range(mu+1):
                if n == m:
                    A[n,m] = (n+1) - sum(s[k]**2 for k in range(n+1))
                else:
                    A[n,m] = -sum(s[n-p]*s[m-p] for p in range(min(n,m)+1))
        ev = np.linalg.eigvalsh(A)
        kneg = int(np.sum(ev < -1e-9*max(1,np.abs(ev).max())))
        kappas.append(kneg)
    return kappas

for a in [0.3, 0.5, 0.7, -0.4, -0.6]:
    K = 30
    s, ok = verify_power_series(a, K)
    ks = kappa_of_series(s, MU_MAX=25)
    print(f"a={a:+.2f} series_ok={ok}  kappa(mu=0..25) = {ks}")
    print(f"          stabilized kappa = {ks[-1]}  (max = {max(ks)})")
