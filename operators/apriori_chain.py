#!/usr/bin/env python3
"""§341b: FULLY A-PRIORI finite-band chain. Replace the arch integral by an explicit
piecewise-constant Riemann MINORANT (rigorous by Omega-monotonicity on [r1,inf), explicit floor
on [0,r1]) so every chain ingredient is an elementary closed formula (erf/exp/cos/Lambda) —
no numeric integration, no zeros anywhere.

  M_Psi[j,k] = (1/2pi) * [ Omega_neg_floor * Z(0,r1) + sum_i Omega(a_i) * Z(a_i, a_{i+1}) ]
  with Z(a,b)[j,k] = int_{a<=|r|<b} F_j F_k dr  (exact erf),  a_i partition of [r1, Rmax], step h.
  Omega evaluated at LEFT endpoints => M_Psi <= arch as forms (pointwise minorant, phi^2 >= 0).
Then run the §339 adaptive chain on  B_lower = pole + M_Psi - P  (B_lower <= B_full as forms, so
certifying B_lower PD certifies B_full PD).  Measure reach vs partition step h.
"""
import numpy as np
from scipy.special import digamma, erf
import l2_edge_sweep as L
from resonant_cert import prime_terms
from band_sweep_cert import gap_bound

w = 3.0
r1 = 7.0
Rmax = 120.0
Omega = lambda r: np.real(digamma(0.25 + 0.5j * np.asarray(r, dtype=float))) - np.log(np.pi)
Om_neg_floor = float(np.min(Omega(np.linspace(1e-6, r1, 4000))))   # most negative value on [0,r1]

def zone_masses_vec(cs, edges):
    """Z_i[j,k] = int over |r| in [edges[i], edges[i+1]) of F_j F_k dr, all i at once.
    F_j F_k = sum over 4 sign combos of Gaussian bumps; each integral is erf-explicit."""
    n = len(cs)
    E = np.asarray(edges)
    combos = []
    for sj in (1, -1):
        for sk in (1, -1):
            cj = sj * cs[:, None]      # (n,1)
            ck = sk * cs[None, :]      # (1,n)
            m = 0.5 * (cj + ck)        # (n,n)
            amp = np.exp(-(cj - ck) ** 2 / (4 * w * w))
            combos.append((m, amp))
    nz = len(E) - 1
    Z = np.zeros((nz, n, n))
    for (m, amp) in combos:
        # int_a^b e^{-(r-m)^2/w^2} dr = (sqrt(pi) w/2)(erf((b-m)/w)-erf((a-m)/w)); plus mirror [-b,-a]
        for (lo, hi) in [(E[:-1], E[1:]), (-E[1:], -E[:-1])]:
            t1 = erf((hi[:, None, None] - m[None]) / w)
            t0 = erf((lo[:, None, None] - m[None]) / w)
            Z += amp[None] * (np.sqrt(np.pi) * w / 2) * (t1 - t0)
    return Z

def build_MPsi(cs, h):
    edges = np.arange(r1, Rmax + h, h)
    Z = zone_masses_vec(cs, edges)
    om = Omega(edges[:-1])                      # LEFT endpoints (minorant)
    M = np.tensordot(om, Z, axes=(0, 0))
    Zlow = zone_masses_vec(cs, np.array([0.0, r1]))[0]
    M += Om_neg_floor * Zlow                    # rigorous floor on the low zone
    return M / (2 * np.pi)

print("A-PRIORI CHAIN (pole + Riemann-minorant arch - primes; all entries elementary formulas)")
print("first: minorant quality — ||arch - M_Psi|| vs h at s=4.0:")
cs = L.probes(L.BAND, 4.0)
f = L.build_forms(cs, w)
for h in (1.0, 0.3, 0.1, 0.03, 0.01):
    MP = build_MPsi(cs, h)
    gap = np.linalg.eigvalsh(f['archpole'] - (MP + f['archpole'] - f['archpole']))  # placeholder
    diff = f['archpole'] - (np.real(MP) + (f['archpole'] - f['archpole']))
    arch_only = f['archpole'] - f['overlap'] * 0  # archpole includes pole; compare arch part
    D = f['archpole'] - MP  # pole included in archpole, not in MP -> D >= pole >= 0 expected
    print(f"  h={h:5.2f}:  min eig(archpole - MPsi) = {np.linalg.eigvalsh(D)[0]:+.3e}  "
          f"(must be >= 0: minorant + pole)  maxdiag gap {np.max(np.diag(D)):.3e}")

print(f"\n{'s':>5} {'n':>3} {'true':>10} " + " ".join(f"h={h:<11}" for h in (0.1, 0.01)) + "  R_used")
for s in (5.0, 4.5, 4.2, 4.0, 3.9, 3.8, 3.7, 3.6, 3.5, 3.4):
    cs = L.probes(L.BAND, s)
    f = L.build_forms(cs, w)
    P = f['prime']
    # pole matrix (exact closed form, a-priori)
    z = 0.5j
    fp = np.exp(-((z - cs) ** 2) / (2 * w * w)) + np.exp(-((z + cs) ** 2) / (2 * w * w))
    pole = 2.0 * np.real(np.outer(fp, fp))
    true_m = float(np.linalg.eigvalsh(f['full'])[0])
    terms = prime_terms(cs)
    order = sorted(terms, key=lambda n: -np.linalg.norm(terms[n], 2))
    line = f"{s:5.1f} {len(cs):3d} {true_m:+10.3e} "
    Rrec = []
    for h in (0.1, 0.01):
        A_lower = pole + build_MPsi(cs, h)
        best = None
        for R in range(1, min(13, len(order) + 1)):
            Pres = sum(terms[n] for n in order[:R])
            b = gap_bound(A_lower - Pres, P - Pres, k=3)
            if best is None or b > best[1]:
                best = (R, b)
            if b > 0:
                break
        line += f"{best[1]:+12.3e} "
        Rrec.append(best[0])
    print(line + f"  R={Rrec}", flush=True)
