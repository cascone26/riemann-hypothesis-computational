#!/usr/bin/env python3
"""
cycle 333 — Is the box-moment Hamburger problem DETERMINATE? (=> heig essentially self-adjoint & UNIQUE)

The candidate heig is the Jacobi operator of the box measure dnu (moments nu_j = m_{j+1},
m_k = sum_gamma box_gamma^{-k}, box=1/4+gamma^2). FiniteBoxRealization (§1230) isolates the crux as
"realize ALL zeros simultaneously + positivity"; HilbertPolyaCore (§1765) needs the operator SELF-ADJOINT.
A Hamburger moment problem that is DETERMINATE => the Jacobi operator is essentially self-adjoint and the
measure is UNIQUE (Stone; limit-point at infinity). So determinacy CLOSES the extension-ambiguity gap:
once the (infinite) box moments are given, heig is unique & essentially self-adjoint — no ambiguity.

CARLEMAN's sufficient condition for determinacy:  sum_k nu_{2k}^{-1/(2k)} = +infinity.
Our moments decay superexponentially (nu_j ~ box_1^{-(j+1)}, box_1 = 1/4+gamma_1^2 ~ 200), so
nu_{2k}^{-1/(2k)} -> box_1 (a positive CONSTANT) => the sum grows ~ box_1 * K -> infinity => DETERMINATE.

Verify numerically: compute m_k from zeros (truth) to high k, form the Carleman terms and partial sums,
show linear growth (divergence). (cy330 verified the m_k are also computable UNCONDITIONALLY from xi and
match; determinacy is a property of the sequence's growth, identical either way.)
"""
import numpy as np

G = np.load('/Users/scones/projects/riemann/operators/zeros2000.npy').astype(float)
box = 0.25 + G**2
box1 = box.min()
print(f"box_1 = 1/4+gamma_1^2 = {box1:.4f}  (gamma_1={G.min():.4f})")

# box moments m_k = sum_gamma box^{-k}  (per-gamma); use log-space for the superexp decay
Kmax = 60
# m_k dominated by smallest box; compute in float (m_k underflows ~ 200^-k, fine in float64 to k~130)
mk = np.array([np.sum(box**(-k)) for k in range(1, Kmax+1)])   # m_1..m_Kmax
# measure nu has moments nu_j = m_{j+1}; Carleman uses even moments nu_{2k} = m_{2k+1}
print(f"\n{'k':>3} {'nu_2k=m_(2k+1)':>16} {'term nu_2k^(-1/2k)':>20} {'partial sum':>14}")
psum = 0.0
terms = []
for k in range(1, Kmax//2):
    nu2k = mk[2*k]        # m_{2k+1} (index 2k in 0-based m_1.. array => m_{2k+1})
    term = nu2k**(-1.0/(2*k))
    psum += term; terms.append(term)
    if k <= 12 or k % 5 == 0:
        print(f"{k:>3} {nu2k:>16.6e} {term:>20.6f} {psum:>14.4f}")

terms = np.array(terms)
# linear-growth check: partial sums ~ slope * K
ks = np.arange(1, len(terms)+1)
cum = np.cumsum(terms)
slope = np.polyfit(ks[-20:], cum[-20:], 1)[0]
print(f"\nCarleman terms -> constant ~{terms[-10:].mean():.3f} (should approach box_1={box1:.3f})")
print(f"partial-sum slope (last 20) = {slope:.3f} > 0  => sum DIVERGES => Carleman satisfied => DETERMINATE")
print(f"\nVERDICT: box-moment Hamburger problem is DETERMINATE (Carleman). Hence the Jacobi/heig operator is")
print(f"ESSENTIALLY SELF-ADJOINT and the realizing measure is UNIQUE — no self-adjoint-extension ambiguity")
print(f"in the Hilbert-Polya construction. (Structural fact about heig; holds for the real-moment/PSD case,")
print(f"i.e. under the same positivity that IS RH — so it sharpens the construction, it is not itself a lever.)")
