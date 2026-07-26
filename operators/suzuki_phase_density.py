#!/usr/bin/env python3
"""suzuki_phase_density.py -- §247: the phase function Phi_a(z) and the phase-dynamics law.

Structure: W_theta = u + e^{i theta} w, u=(z-i)vhat_+, w=(z+i)vhat_-. On the real axis |u|=|w|
(de Branges reality), so e^{i Phi(z)} := -u/w defines the PHASE FUNCTION; every extension's spectrum
= level set Phi(z) = theta + pi(2n+1)... All theta-branches ride ONE Phi -> explains §241/§243
theta-universality immediately.
Phase-dynamics law to verify:  v_branch = -dPhi/da / dPhi/dz.
Then: parked v_min = |d_aPhi(gamma)| / Phi'(gamma), with Phi'(gamma) = pi/(zero step width).
MEASURE across a: (1) |u|/|w| on axis (structure check ~1), (2) Phi'(gamma_j) growth law (predict
~a^{2.4} from c(a)~a^{-1.44} if d_aPhi ~ gamma·a^0), (3) d_aPhi(gamma_j), (4) predicted vs measured v.
"""
import numpy as np, sys, os, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import suzuki_exact_Q as SQ
from scipy.linalg import solve
from scipy.interpolate import CubicSpline

_ts = np.logspace(-9, np.log10(9.0), 3000)
_ks = np.array([0.5*np.exp(-t/2)*SQ.lerch1_quarter_t(t) for t in _ts])
_spl = CubicSpline(np.log(_ts), _ks)
def tail_nodes_fast(a, NG=320):
    L = 2*a
    xs, ws = np.polynomial.legendre.leggauss(NG)
    s = 0.5*(xs+1.0); t = L*s**3
    return t, 0.5*ws*3*L*s**2*_spl(np.log(np.maximum(t, 1e-9)))
SQ.tail_nodes = tail_nodes_fast

GAMMAS = [14.134725, 21.022040, 25.010858, 30.424876, 32.935062]

def uw_funcs(a, K=None, NG=320):
    if K is None: K = int(np.ceil(50*2*a/np.pi)) + 40
    kk = np.arange(1, K+1)
    Q = SQ.build_Q(a, K, NG)
    lam_a = np.linalg.eigvalsh(Q)[0]/a
    Ta = Q - (lam_a-0.2)*a*np.eye(K)
    fp = SQ.E_p(1.0, a, kk); fm = SQ.E_p(-1.0, a, kk)
    cp = solve(Ta, fp, assume_a='sym'); cm = solve(Ta, fm, assume_a='sym')
    cm *= np.sqrt(abs(cp@fp))/np.sqrt(abs(cm@fm))
    om = kk*np.pi/(2*a)
    def uw(z):
        den = om*om - z*z; den = np.where(np.abs(den) < 1e-9, 1e-9, den)
        sh = np.exp(-1j*z*a)*om*(1.0-((-1.0)**kk)*np.exp(2j*z*a))/den
        vp = cp@sh; vm = cm@sh
        return (z-1j)*vp, (z+1j)*vm
    return uw

def Phi_at(uw, z):
    u, w = uw(z)
    return float(np.angle(-u/w))

def Phi_prime(uw, z, dz=2e-4):
    # unwrapped derivative
    p1, p2 = Phi_at(uw, z-dz), Phi_at(uw, z+dz)
    d = p2 - p1
    while d < -np.pi: d += 2*np.pi
    while d > np.pi: d -= 2*np.pi
    return d/(2*dz)

if __name__ == "__main__":
    print("=== structure check: |u|/|w| on axis (should be 1) ===")
    uw = uw_funcs(2.0)
    for z in [12.0, 14.1347, 20.0, 25.01, 33.0]:
        u, w = uw(z)
        print(f"  z={z:8.4f}  |u|/|w| = {abs(u)/abs(w):.6f}")
    print("\n=== Phi'(gamma_j) growth vs a + phase-dynamics law check ===")
    avals = [1.6, 2.0, 2.4, 2.8, 3.2, 3.6]
    da = 2e-3
    out = {}
    print(f"{'a':>5} {'gamma':>8} {'Phi_prime':>11} {'dPhi/da':>9} {'v_pred':>9}  (v_pred = -daPhi/Phi')")
    for a in avals:
        uw0 = uw_funcs(a-da); uw1 = uw_funcs(a); uw2 = uw_funcs(a+da)
        for j, gm in enumerate(GAMMAS):
            pp = Phi_prime(uw1, gm)
            d = Phi_at(uw2, gm) - Phi_at(uw0, gm)
            while d < -np.pi: d += 2*np.pi
            while d > np.pi: d -= 2*np.pi
            daPhi = d/(2*da)
            vpred = -daPhi/pp
            out.setdefault(str(a), {})[str(j)] = {"Phip": pp, "daPhi": daPhi, "vpred": vpred}
            print(f"{a:5.1f} {gm:8.3f} {pp:11.3f} {daPhi:9.3f} {vpred:+9.4f}")
    json.dump(out, open("../results/phase_density.json", "w"), indent=1)
    # growth law of Phi'(gamma)
    print("\nPhi'(gamma) growth: fit log Phi' vs log a per zero")
    for j, gm in enumerate(GAMMAS):
        A = np.array(avals); P = np.array([out[str(a)][str(j)]["Phip"] for a in avals])
        if np.all(P > 0):
            sl = np.polyfit(np.log(A), np.log(P), 1)[0]
            print(f"  gamma_{j+1}: Phi'({gm:.1f}) = {[round(float(p),2) for p in P]}  slope={sl:.2f}")
    print("\nCompare: c(a)~a^-1.44 predicts Phi'(gamma) ~ a^{2.44} IF daPhi ~ a^0*gamma; see daPhi column.")
