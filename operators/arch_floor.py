#!/usr/bin/env python3
"""§341: localized a-priori arch floor (ingredient I2-fixed).

Pointwise-minorant principle: for the quadratic form <v, arch v> = (1/2pi) Int phi^2 Omega dr
(phi = sum v_j F_j >= 0 pointwise-squared), any pointwise split Omega >= Psi - Err charges
form-wise:  arch >= M_Psi - (1/2pi)*(max Err on zone)*(zone mass form).

Three-zone construction (all explicit / provable-shape):
  Z1 |r| < r1(=7):    Omega can be negative; charge max|Omega - Omega_floor| x Gaussian zone mass
                      (tiny for probes at c>=14: erfc tails).
  Z2 r1 <= |r| < c0:  Omega in [Omega(r1), Omega(c0)]; charge (Omega(c0)-Omega(r1)) x zone mass
                      (only EDGE probes have real mass here — the identified refinement point).
  Z3 |r| >= c0:       Omega >= Omega(c0) =: Om0 (monotonicity, verified);
                      floor (Om0/2pi)*G with the theta4 Toeplitz bound.
Variants measured:
  V1: flat floor Om0 = Omega(c0) + zone charges           [min over probe centers]
  V2: D-weighted: arch >= lambda_min(G)*(1/2pi)*D - charges, D_j = Omega(c_j)  [matrix floor kept
      diagonal — better on interior probes; enters the chain as a matrix, evaluated here by its
      lambda_min for standalone reach]
Verify Omega monotone increasing on [r1, 200] numerically (provable via psi asymptotics).
"""
import numpy as np
from scipy.special import digamma, erf
import l2_edge_sweep as L

w = 3.0
r1 = 7.0
Omega = lambda r: np.real(digamma(0.25 + 0.5j * np.asarray(r, dtype=float))) - np.log(np.pi)

# monotonicity check on [r1, 200]
rr = np.linspace(r1, 200, 4000)
dOm = np.diff(Omega(rr))
print(f"Omega monotone increasing on [{r1},200]: {bool(np.all(dOm > 0))} (min step {dOm.min():.2e})")

def theta4(q, terms=80):
    return 1 + 2 * sum((-1) ** k * q ** (k * k) for k in range(1, terms))

def zone_mass_matrix(cs, a, b):
    """Z[j,k] = int_{a<=|r|<b} F_j F_k dr, exact via erf (Gaussian products)."""
    n = len(cs)
    Z = np.zeros((n, n))
    s2 = np.sqrt(2.0)
    def bump_pair(cj, ck, lo, hi):
        # int_lo^hi e^{-(r-cj)^2/2w^2} e^{-(r-ck)^2/2w^2} dr
        m = 0.5 * (cj + ck)
        pref = np.exp(-(cj - ck) ** 2 / (4 * w * w)) * w * np.sqrt(np.pi) / s2 * np.sqrt(2)
        # int e^{-(r-m)^2/w^2} dr = (w sqrt(pi)/ ... ) use erf: = (sqrt(pi) w / 2)*[erf((hi-m)/w)-erf((lo-m)/w)]
        return np.exp(-(cj - ck) ** 2 / (4 * w * w)) * (np.sqrt(np.pi) * w / 2) * (
            erf((hi - m) / w) - erf((lo - m) / w))
    for j in range(n):
        for k in range(j, n):
            v = 0.0
            for (a1, b1) in [(a, b), (-b, -a)]:
                for (sj, sk) in [(1, 1), (1, -1), (-1, 1), (-1, -1)]:
                    v += bump_pair(sj * cs[j], sk * cs[k], a1, b1)
            Z[j, k] = Z[k, j] = v
    return Z

print(f"\n{'s':>5} {'lminA meas':>11} {'V1 flat':>10} {'V2 wtd':>10} {'Z1 chg':>8} {'Z2 chg':>8} "
      f"{'V2/meas':>8}")
c0 = 14.0
Om_r1, Om_c0 = Omega(r1), Omega(c0)
maxAbsLow = float(np.max(np.abs(Omega(np.linspace(0, r1, 800)) - Om_r1)))
for s in (6.0, 5.5, 5.0, 4.5, 4.0, 3.5, 3.0):
    cs = L.probes(L.BAND, s)
    f = L.build_forms(cs, w)
    A = f['archpole']; G = f['overlap']
    lminA = float(np.linalg.eigvalsh(A)[0])
    q = np.exp(-s * s / (4 * w * w))
    lminG_bnd = w * np.sqrt(np.pi) * theta4(q)          # I1, verified valid
    Z1 = zone_mass_matrix(cs, 0.0, r1)
    Z2 = zone_mass_matrix(cs, r1, c0)
    chg1 = maxAbsLow / (2 * np.pi) * float(np.linalg.eigvalsh(Z1)[-1])
    chg2 = (Om_c0 - Om_r1) / (2 * np.pi) * float(np.linalg.eigvalsh(Z2)[-1])
    V1 = (Om_c0 / (2 * np.pi)) * lminG_bnd - chg1 - chg2
    # V2: arch >= (1/2pi)[ lminG_bnd * D_min-per-direction ] - charges; standalone reach via
    # lambda_min of the matrix bound (1/2pi)*lminG_bnd*D - chg1 - chg2  (D diagonal Omega(c_j))
    D = np.diag(Omega(cs))
    M = (lminG_bnd / (2 * np.pi)) * D
    V2 = float(np.linalg.eigvalsh(M)[0]) - chg1 - chg2   # = same min for diagonal, but kept for chain use
    print(f"{s:5.1f} {lminA:11.4e} {V1:+10.3e} {V2:+10.3e} {chg1:8.4f} {chg2:8.4f} {max(V1,V2)/lminA:8.2f}")

print("""
Reading: valid iff bound <= measured lminA; useful iff positive. Z2 charge = the edge-probe cost
(only probes near c0 have mass in [r1,c0]) — if it dominates, the refinement is an edge-block
treatment (2-block chain with the edge probes as the block), not a better global floor.""")
