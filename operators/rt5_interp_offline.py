#!/usr/bin/env python3
"""§378: does the zero-interpolation corner mode (§375) DISTINGUISH on-line from off-line zeros?

§375: at criticality the Weil near-null mode phi = sum v_full F vanishes on the (on-line) zero
ordinates. RH-mechanism question: build the zero-side Gram for a FAKE world where one zero sits
off-line at 1/2+delta+i*gamma (a rank-2 real contribution, §364), find ITS near-null mode, and
ask: can the real even test class still interpolate (vanish on) the fake configuration as
cleanly? On-line zeros give real ordinates {gamma_k} (n points, rank-1 each). An off-line zero
replaces one rank-1 real-ordinate constraint by a rank-2 constraint (real+imag parts at the
complex ordinate) — so the null space DIMENSION drops by 1 per off-line zero at fixed probe
count. Measure: the smallest eigenvalue of the fake zero-side Gram at n=n_crit (on-line world
is rank-deficient ~0; off-line world should be LIFTED because the extra rank fills the deficiency)
and the interpolation quality of the fake near-null mode. If off-line consistently LIFTS the
floor and DEGRADES interpolability, that is the geometric signature of RH on this instrument
(conservative per Q5, but the mechanism made explicit). Run on BOTH instruments (Gaussian+Fejer).
"""
import numpy as np
import l2_edge_sweep as G
from scipy.special import digamma
import l2_fejer as F

def N_rvm(T):
    if T <= 2 * np.pi: return 0.0
    return T / (2 * np.pi) * np.log(T / (2 * np.pi)) - T / (2 * np.pi) + 7.0 / 8.0

def run(name, Frow, zeros, cs, w):
    """zero-side lmin & interp for on-line vs one-zero-off-line, at probe set cs."""
    def zside(zlist):
        # zlist entries: real g (on-line) or complex g+1j*delta (off-line, adds real+imag rank-2)
        n = len(cs); B = np.zeros((n, n))
        for z in zlist:
            if np.imag(z) == 0:
                v = np.array([Frow(c, np.array([np.real(z)]))[0] for c in cs])
                B += 2.0 * np.outer(v, v)
            else:
                # off-line: F evaluated at complex ordinate; even real test fn => use Re and Im parts
                vc = np.array([Frow(c, np.array([z]))[0] for c in cs])
                B += 2.0 * np.real(np.outer(vc, np.conj(vc)))
        return B
    online = [complex(g, 0) for g in zeros]
    lmin_on = float(np.linalg.eigvalsh(zside(online))[0])
    out = [f"  [{name}] n={len(cs)} on-line zero-side lmin = {lmin_on:+.3e}"]
    for delta in (0.05, 0.1, 0.2):
        # move the zero with largest projection on the on-line null mode
        Bon = zside(online); lo, Vo = np.linalg.eigh(Bon); vnull = Vo[:, 0]
        proj = [(abs(np.array([Frow(c, np.array([g]))[0] for c in cs]) @ vnull), g) for g in zeros]
        gsel = max(proj)[1]
        fake = [complex(g, 0) for g in zeros if g != gsel] + [complex(gsel, delta)]
        lmin_off = float(np.linalg.eigvalsh(zside(fake))[0])
        out.append(f"    delta={delta}: off-line lmin = {lmin_off:+.3e}  lift = {lmin_off - lmin_on:+.3e}")
    return "\n".join(out)

# Gaussian instrument
w = 3.0
ncrit_g = 9
cs_g = np.linspace(14, 44, ncrit_g)
def Grow(c, r): return np.exp(-((r - c) ** 2) / (2 * w * w)) + np.exp(-((r + c) ** 2) / (2 * w * w))
zg = [g for g in G.ZEROS if 8 <= g <= 50]
print("§378 off-line signature on the zero-side Gram at criticality (on-line rank-deficient ~0):")
print(run("Gaussian w=3", Grow, zg, cs_g, w))
# Fejer instrument
ncrit_f = 7
cs_f = np.linspace(14, 44, ncrit_f)
def Frow_f(c, r): return F.fej(r - c) + F.fej(r + c)
zf = [g for g in F.ZEROS if 8 <= g <= 50]
print(run("Fejer W=1.2", Frow_f, zf, cs_f, w))
print("\nreading: if off-line consistently LIFTS the zero-side lmin above the on-line ~0 floor (a")
print("POSITIVE lift growing with delta) on BOTH instruments, the rank-2 off-line contribution")
print("provably cannot be interpolated away by the real even class — the geometric RH signature")
print("(§364 rank-2 + §375 interpolation unified). Lift ~ c*delta^2 expected (§364 c(delta) law).")
