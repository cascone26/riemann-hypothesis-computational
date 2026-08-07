"""
G2 (BARRIER_MAP): explicit averaging-wall witness family — continuous-level demo.

Family of K continuous Beurling-type zero configurations (zero sets are ours to
prescribe; discretization to genuine g-prime systems with these zeta zeros is
Broucke-Debruyne-Revesz arXiv:2309.01567 Thm 1.2, cited not executed here).

  - K-1 "good" members: all zeros on the critical line, low-lying (Beurling
    zetas can and do have low zeros), ordinates from ~1 upward.
  - 1 "bad" member: an off-line pair rho = sigma0 +/- i*t0 with LOW ordinate t0,
    its own on-line zeros starting higher (t0 below the first on-line ordinate
    is what lets the off-line term dominate its own form). sigma0 near 1 is
    exactly Diamond-Montgomery-Vorhauer territory (zeros approaching sigma=1).

Weil zero-side form for Gaussian h(r) = exp(-r^2/w^2):
  W_j(h) = sum over zeros of h(gamma), gamma = (rho-1/2)/i; pairs +/-.
  On-line zero: gamma real -> positive term.
  Off-line pair: gamma = t0 -/+ i*delta, delta = sigma0-1/2, contributes
    2*exp((delta^2 - t0^2)/w^2) * cos(2*t0*delta/w^2)  (sign-indefinite).

Demonstrated witness statement:
  mean_j W_j(h) > 0, W_bad(h) < 0, and the outlier's deviation is only a few
  sigma — i.e. a Chebyshev/variance (large-sieve-strength) bound CANNOT exclude
  one violator in a family this size. The averaging wall (LEDGER section 62)
  realized as an explicit family, not just a logic-ceiling.
"""
import numpy as np

rng = np.random.default_rng(8)

K = 12                       # family size
NZ = 40                      # on-line zeros per member

sigma0, t0 = 0.95, 2.0       # off-line pair (DMV territory: sigma near 1, low t)
delta = sigma0 - 0.5

members = []
for j in range(K - 1):       # good members: on-line zeros from ~1 upward
    ords = np.cumsum(0.7 + rng.exponential(0.5, NZ)) + 0.3 + rng.normal(0, 0.03)
    members.append({"ordinates": ords, "off": None})
# bad member: off-line pair at low ordinate; own on-line zeros start at 10
bad_ords = np.cumsum(0.7 + rng.exponential(0.5, NZ)) + 10.0
members.append({"ordinates": bad_ords, "off": (t0, delta)})

def W(member, w):
    s = 2.0 * np.sum(np.exp(-(member["ordinates"] / w) ** 2))
    if member["off"] is not None:
        t, d = member["off"]
        s += 2.0 * np.exp((d * d - t * t) / (w * w)) * np.cos(2.0 * t * d / (w * w))
    return s

print(f"family K={K}; bad member: off-line pair rho = {sigma0} +/- {t0}i, "
      f"on-line zeros from {bad_ords[0]:.1f}")
print(f"{'w':>6} {'W_bad':>12} {'mean_W':>12} {'std_W':>12} {'outlier_sig':>12}")
found = None
for w in np.arange(0.4, 3.0, 0.02):
    vals = np.array([W(m, w) for m in members])
    bad, mean, std = vals[-1], vals.mean(), vals.std()
    if bad < 0 and mean > 0:
        nsig = abs(bad - mean) / std if std > 0 else np.inf
        if found is None:
            found = (w, vals)
        print(f"{w:6.2f} {bad:12.4e} {mean:12.4e} {std:12.4e} {nsig:12.2f}")

if found is None:
    print("no width found -- adjust geometry")
else:
    w, vals = found
    bad, mean, std = vals[-1], vals.mean(), vals.std()
    nsig = abs(bad - mean) / std
    cheb = 1.0 / nsig**2
    print(f"\nWITNESS at first admissible width w = {w:.2f}:")
    print(f"  W_bad      = {bad:.6e}  < 0   (individual member VIOLATES positivity)")
    print(f"  mean_j W_j = {mean:.6e}  > 0   (family average positivity HOLDS)")
    print(f"  outlier at {nsig:.2f} sigma; Chebyshev bound on outlier prob = "
          f"{cheb:.3f} vs family rate 1/K = {1/K:.3f}")
    if cheb > 1.0 / K:
        print("  => a variance/Chebyshev (large-sieve-strength) bound CANNOT exclude")
        print("     this single violator: mean + variance information never reaches")
        print("     the individual. Averaging wall realized as explicit family.")
    else:
        print("  NOTE: Chebyshev would exclude this outlier at this width -- "
              "geometry needs adjusting for the full witness statement.")
