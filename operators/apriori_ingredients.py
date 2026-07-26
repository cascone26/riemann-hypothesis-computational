#!/usr/bin/env python3
"""§340b: A-PRIORI INGREDIENT VERIFICATION for the adaptive chain (band [14,44], w=3).
Each chain ingredient gets an explicit, provable-shape bound; verify each against the measured value.

I1 (overlap Gram floor, CLASSICAL): the overlap Gram is (near-)Toeplitz with symbol
    f(theta) = w*sqrt(pi) * sum_k e^{-k^2 s^2/4w^2} e^{ik theta},  min f = w*sqrt(pi)*theta4(q), q=e^{-s^2/4w^2}
    => lambda_min(G_n) >= w*sqrt(pi)*theta4(q) - (e_far correction).   [finite Toeplitz sections >= ess-inf f]
I2 (arch floor): A >= (Omega_0/2pi)*G - eps_tail*I, Omega_0 = min_{|r|>=r0} Omega(r), with the |r|<r0 tail
    bounded explicitly by Gaussian mass:  eps_tail = (max|Omega|_{<r0}/pi) * n * max_j Int_{|r|<r0} F_j^2.
I3 (tail norm): lambda_max(P_off) <= Gershgorin: max_j sum_k sum_{n notin core} 2*Lam(n)/sqrt(n)*|g_pair| with
    |cos|<=1 (triangle) — the PB-style bound, explicit.
Open (flagged, not verified here): a-priori lower bound on lambda_{k+1}(A) for the gap corner —
    needs finite-n Toeplitz eigenvalue bounds (Bogoya-Boettcher-Grudsky) on the arch symbol.
"""
import numpy as np
import l2_edge_sweep as L
from resonant_cert import prime_terms

w = 3.0
r0 = 7.0   # below r0, Omega can be negative/small

def theta4(q, terms=60):
    return 1 + 2 * sum((-1) ** k * q ** (k * k) for k in range(1, terms))

def theta3(q, terms=60):
    return 1 + 2 * sum(q ** (k * k) for k in range(1, terms))

# Omega_0 = min Omega on |r|>=r0 (Omega increasing for r>~7; min at r0)
from scipy.special import digamma
Omega = lambda r: np.real(digamma(0.25 + 0.5j * r)) - np.log(np.pi)
Omega0 = Omega(r0)
OmegaMaxLow = np.max(np.abs(Omega(np.linspace(0, r0, 500))))
print(f"Omega_0 = Omega({r0}) = {Omega0:.4f};  max|Omega| on [0,{r0}] = {OmegaMaxLow:.4f}")

print(f"\n{'s':>5} {'lmin(G) meas':>12} {'I1 theta4 bnd':>13} {'lmin(A) meas':>12} {'I2 bnd':>10} "
      f"{'lmax(Poff) meas':>15} {'I3 bnd':>10}")
for s in (6.0, 5.5, 5.0, 4.5, 4.0, 3.5):
    cs = L.probes(L.BAND, s)
    n = len(cs)
    f = L.build_forms(cs, w)
    G = f['overlap']; A = f['archpole']; P = f['prime']
    q = np.exp(-s * s / (4 * w * w))
    I1 = w * np.sqrt(np.pi) * theta4(q)
    # I2: A >= (Omega0/2pi) G - eps_tail
    # gaussian mass below r0 for probe at c: 2*int_{-r0}^{r0} F^2 dr <= 4*w*sqrt(pi)*exp(-(c-r0)^2/w^2)... explicit:
    from scipy.special import erf
    def lowmass(c):
        # int_{-r0}^{r0} F_c(r)^2 dr, F=exp(-(r-c)^2/2w^2)+exp(-(r+c)^2/2w^2); bound by 4*int of nearest bump sq
        z = np.linspace(-r0, r0, 2000)
        Fv = np.exp(-(z - c) ** 2 / (2 * w * w)) + np.exp(-(z + c) ** 2 / (2 * w * w))
        return np.trapezoid(Fv * Fv, z)
    eps_tail = (Omega0 / (2 * np.pi) + OmegaMaxLow / (2 * np.pi)) * n * max(lowmass(c) for c in cs)
    I2 = (Omega0 / (2 * np.pi)) * w * np.sqrt(np.pi) * theta4(q) - eps_tail
    # I3: Gershgorin triangle bound on the tail beyond core {2,3,4,5} (R=4 reference core)
    terms = prime_terms(cs)
    core = [2, 3, 4, 5]
    Poff = P - sum(terms[nn] for nn in core if nn in terms)
    meas_off = float(np.linalg.eigvalsh(Poff)[-1])
    pref = w * np.sqrt(np.pi) / (2 * np.pi)
    I3 = 0.0
    for i, nn in enumerate(L.NS):
        if L.LAMW[i] == 0.0 or int(nn) in core:
            continue
        u = L.LOGN[i]
        # row-sum bound: 2*LAMW * pref * e^{-u^2w^2/4} * 2 * (theta3 row-overlap sum + e_far tiny)
        I3 += 2 * L.LAMW[i] * pref * np.exp(-u * u * w * w / 4) * 2 * theta3(np.exp(-s * s / (4 * w * w)))
    print(f"{s:5.1f} {np.linalg.eigvalsh(G)[0]:12.4e} {I1:13.4e} {np.linalg.eigvalsh(A)[0]:12.4e} "
          f"{I2:+10.3e} {meas_off:15.4e} {I3:10.4e}")

print("""
Reading: I1 valid iff bound <= measured (Toeplitz floor). I2 valid iff bound <= lmin(A) — and useful
iff positive. I3 valid iff bound >= lmax(Poff). The open piece (flagged): a-priori lambda_{k+1}(A)
for the gap corner — Bogoya-Boettcher-Grudsky finite-n Toeplitz eigenvalue theory is the named tool.""")
