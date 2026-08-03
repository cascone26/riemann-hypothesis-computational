"""Universality of the plain-Fourier square-wave comb (LEDGER 610) for chi3:
J_chi(u) = (1/K) int sgn Z_chi(t) e^{i theta_chi(t)} e^{iut} dt  -> signs mu(n)chi3(n)?
Z_chi3 real (self-dual odd character mod 3). L via Hurwitz: L(s,chi)=q^{-s} sum chi(a) zeta(s,a/q)."""
import numpy as np, mpmath as mp
mp.mp.dps = 15
q = 3
chi = {1: 1, 2: -1}  # chi3
def theta(t):
    return float(mp.im(mp.loggamma((0.5+1j*t+1)/2)) - (t/2)*mp.log(mp.pi/q))
def Lline(t):
    s = mp.mpf(0.5) + 1j*t
    return q**(-s) * sum(chi[a]*mp.zeta(s, mp.mpf(a)/q) for a in chi)
T0, T, dt = 5.0, 130.0, 0.008
ts = np.arange(T0, T, dt)
th = np.array([theta(t) for t in ts])
Zv = np.array([complex(mp.e**(1j*mp.mpf(theta(t)))*Lline(t)).real if False else 0 for t in [0]])  # placeholder
Z = np.empty(len(ts))
for i, t in enumerate(ts):
    Z[i] = complex(mp.e**(1j*mp.mpf(th[i]))*Lline(t)).real
sq = np.sign(Z)
K = (th[-1]-th[0])/np.pi
print("grid %d pts, K~%.0f zeros' worth" % (len(ts), K))
mob = {2:-1,3:-1,4:0,5:-1,6:1,7:-1,10:1,11:-1,13:-1}
c = None
for n in sorted(mob):
    b = mob[n]*(chi[n % 3] if n % 3 else 0)
    J = np.trapezoid(sq*np.exp(1j*(th+np.log(n)*ts)), dx=dt)/K
    if c is None:
        anchor_n, anchor_b = n, b
        c = np.angle(J) - (0 if b > 0 else np.pi)
    v = (J*np.exp(-1j*c)).real
    ok = "OK" if (b == 0 and abs(v) < 0.04) or (b != 0 and np.sign(v) == np.sign(b)) else "MISS"
    print(f"  n={n:<3} b=mu*chi={b:+d}  Re J={v:+.4f}  {ok}")
Joff = np.trapezoid(sq*np.exp(1j*(th+1.35*ts)), dx=dt)/K
print("  off-support u=1.35: |J| = %.4f" % abs(Joff))
print("DONE")
