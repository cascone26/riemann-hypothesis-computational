#!/usr/bin/env python3
"""RIGOROUS (Arb) single-Gaussian Weil functional D(a)=arch(a,a)+pole(a,a)−prime(a,a).
Float sweep (weil_diagonal_curve.py) showed D(a)≈0 to machine-eps for a below the first zero (~14.13)
then D>0 above — the explicit formula D(a)=Σ_γ h_a(γ), a positive zero-sum, grazing 0 where φ_a
overlaps no zeros. The SIGN in the grazing zone is the RH-content; float can't resolve ±1e-15.
Reuse organ2_certify's Arb entries to CERTIFY the sign/enclosure of D(a) at test points."""
from flint import arb
import organ2_certify as O

print(f"{'a':>6} {'arch':>14} {'pole':>14} {'prime':>14} {'D(a) enclosure':>26} {'sign':>10}")
for av in [3, 5, 8, 10, 12, 13, 14, 15, 20]:
    a = arb(av)
    ar = O.arch_entry(a, a); po = O.pole_entry(a, a); pr = O.prime_entry(a, a)
    D = ar + po - pr
    mid = float(D.mid()); rad = float(D.rad())
    lo, hi = mid - rad, mid + rad
    sign = "≥0 CERT" if lo > 0 else ("≤0 CERT" if hi < 0 else "ambiguous")
    print(f"{av:>6} {float(ar.mid()):>+14.8f} {float(po.mid()):>+14.8f} {float(pr.mid()):>+14.8f} "
          f"[{lo:+.3e},{hi:+.3e}] {sign:>10}")
