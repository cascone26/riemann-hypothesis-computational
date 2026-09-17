"""
cy460 — Arithmetic-Wronskian / Laguerre-inequality probe (opus roundtable residue C).

QUESTION: For the ground-state eigenfunction of the truncated Weil form QW_λ, form its
real trig-polynomial Fourier image F(x)=Σ v_n cos(2πnx/L). The Laguerre inequality
    L[F](x) := F'(x)^2 - F(x)F''(x) >= 0  on ℝ
is a NECESSARY condition for F to be in the Laguerre-Pólya class (all real zeros).
CvS proved the finite truncation IS real-rooted, so L[F]>=0 should hold.

The genuinely-informative test (non-circular question): is the MARGIN of L[F]>=0
  (a) O(1) / robust  -> hints a structural (arithmetic) real-rootedness certificate might exist (LEAD)
  (b) ~ ε_λ scale    -> the certificate tracks the delicate Weil-positivity balance (CIRCLES back to the wall)
We report min_x L[F], its scale vs the ground eigenvalue ε_λ (eps1) and the spectral gap, across λ.
Report-tier numerics (mpmath/numpy). NO proof claim.
"""
import sys, numpy as np, mpmath as mp
sys.path.insert(0, "/Users/scones/projects/riemann/operators")
from qw_lambda_faithful import build_QW

mp.mp.dps = 40

def probe(N, lam):
    Q, idx, L = build_QW(N, lam, verbose=False)
    n = len(idx)
    A = np.array([[float(Q[i,j]) for j in range(n)] for i in range(n)], dtype=float)
    A = 0.5*(A+A.T)
    ev, evec = np.linalg.eigh(A)
    eps1, eps2 = ev[0], ev[1]
    v = evec[:,0]                      # ground-state eigenvector, real (numpy eigh)
    ns = [int(k) for k in idx]         # integer frequencies n
    # even-fraction (parity) sanity
    vmap = {ns[i]: v[i] for i in range(n)}
    evenfrac = None
    if all((-k) in vmap for k in ns):
        num = sum((vmap[k]+vmap[-k])**2 for k in ns)
        den = sum((2*vmap[k])**2 for k in ns) or 1.0
        evenfrac = num/den
    # F(x)=Σ v_n cos(2π n x/L); F'=Σ v_n(-2πn/L)sin; F''=Σ v_n(-(2πn/L)^2)cos
    w = [2*np.pi*k/L for k in ns]
    def F(x):  return sum(v[i]*np.cos(w[i]*x)            for i in range(n))
    def dF(x): return sum(-v[i]*w[i]*np.sin(w[i]*x)      for i in range(n))
    def ddF(x):return sum(-v[i]*w[i]*w[i]*np.cos(w[i]*x) for i in range(n))
    xs = np.linspace(0, L, 4000)
    Lag = np.array([dF(x)**2 - F(x)*ddF(x) for x in xs])
    Fv  = np.array([F(x) for x in xs])
    # normalize the Laguerre margin by the natural scale of F'^2 (so it's a relative sign test)
    scaleF = np.max(np.abs(Fv))**2 * (max(abs(wi) for wi in w))**2 + 1e-300
    return dict(N=N, lam=lam, eps1=eps1, gap=eps2-eps1, evenfrac=evenfrac,
                minLag=float(Lag.min()), relminLag=float(Lag.min()/scaleF),
                fracNeg=float((Lag<0).mean()))

print(f"{'N':>4} {'lam':>6} {'eps1':>12} {'gap':>12} {'evenfrac':>9} {'min L[F]':>13} {'rel min L[F]':>13} {'%neg':>6}")
for (N,lam) in [(20, mp.e), (24, mp.sqrt(10)), (28, mp.sqrt(14)), (32, mp.sqrt(18))]:
    try:
        r = probe(N, float(lam))
        ef = f"{r['evenfrac']:.4f}" if r['evenfrac'] is not None else "  n/a "
        print(f"{r['N']:>4} {r['lam']:>6.3f} {r['eps1']:>12.3e} {r['gap']:>12.3e} {ef:>9} "
              f"{r['minLag']:>13.3e} {r['relminLag']:>13.3e} {100*r['fracNeg']:>5.1f}%")
    except Exception as e:
        print(f"{N:>4} {float(lam):>6.3f}  ERROR: {e}")
