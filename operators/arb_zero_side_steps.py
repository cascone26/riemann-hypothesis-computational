"""S3023: RIGOROUS zero side of the detector-step identity (lane 2).
For the Arb-certified prime-side steps s_j (S3015/S3021 dumps, lower bound + mid), compute the zero-tail sum
P_j(T) = sum_{V+1 < gamma <= T} (4/L) sin^2(gamma L/2) [ 1/(nu_j-gamma) - 1/(nu_{j+1}-gamma) + 1/(nu_j+gamma) - 1/(nu_{j+1}+gamma) ]
in ball arithmetic with CERTIFIED zero ordinates (acb.zeta_zeros) and certified completeness (arb.zeta_nzeros(T)),
and a RIGOROUS upper bound U_j for the gamma > T remainder:
  each term <= (4/L) * 2*dnu/(gamma - V)^2 ; sum_{gamma>T} f(gamma) <= int_T^inf N(t)|f'(t)| dt with
  N(t) <= (t/2pi) log(t/2pi) (valid for t >= 12000 by Trudgian's explicit N(T) bound) and t/(t-V) <= k := T/(T-V):
  U = (8/pi) dnu k^3 (log(T/2pi)+1) / (L T).
Checks, for every j:  0 < D_j := s_j - P_j(T)  and  D_j <= U_j  (identity A = G_inf with all zeros <= T on the line).
usage: python3 arb_zero_side_steps.py lam2 N dumpfile T"""
import sys, math
from flint import arb, acb, ctx
ctx.prec = 200
lam2, N, dump, T = sys.argv[1], int(sys.argv[2]), sys.argv[3], int(sys.argv[4])
from fractions import Fraction as F
q = F(lam2); L = (arb(q.numerator) / arb(q.denominator)).log(); PI = arb.pi()
rows = [l.split() for l in open(dump) if not l.startswith('#')]
nu = [2 * PI * n / L for n in range(-N, N + 1)]
assert len(rows) == 2 * N and all(abs(float(nu[j].mid()) - float(rows[j][1])) < 1e-6 for j in range(2 * N))
# certified steps: [lower, mid] widened by 1e-11 relative for the 13-digit dump rounding
S = []
for r in rows:
    lo, mid = F(r[4]), F(r[3])
    a = arb(lo.numerator) / arb(lo.denominator); b = arb(mid.numerator) / arb(mid.denominator)
    S.append((a * (1 - arb('1e-11'))).union(b * (1 + arb('1e-11'))))
V = nu[-1]; Vp = V + 1
n0 = int(round(float(arb.zeta_nzeros(Vp).mid()))); n1b = arb.zeta_nzeros(arb(T)); n1 = int(round(float(n1b.mid())))
assert float(n1b.rad()) < 0.5
zs = [z.imag for z in acb.zeta_zeros(n0 + 1, n1 - n0)]
assert all(z.overlaps(z) and (z > Vp) and (z < T) for z in zs)
P = [arb(0)] * (2 * N)
for z in zs:
    w = (4 / L) * (z * L / 2).sin() ** 2
    for j in range(2 * N):
        P[j] += w * (1 / (nu[j] - z) - 1 / (nu[j + 1] - z) + 1 / (nu[j] + z) - 1 / (nu[j + 1] + z))
dnu = nu[1] - nu[0]; k = arb(T) / (arb(T) - V)
U = (8 / PI) * dnu * k ** 3 * ((arb(T) / (2 * PI)).log() + 1) / (L * T)
D = [S[j] - P[j] for j in range(2 * N)]
pos = all(d > 0 for d in D); le = all(d < U for d in D)
dl = [float(d.lower()) for d in D]; du = [float(d.upper()) for d in D]
print(f"lam2={lam2} N={N} T={T}: certified zeros #{n0+1}..#{n1} (all on the line, gamma in ({float(Vp.mid()):.2f},{T})); "
      f"D_j in [{min(dl):.6e}, {max(du):.6e}], spread<={max(du)-min(dl):.2e}; U={float(U.upper()):.6e}; "
      f"CERTIFIED 0<D_j: {pos}; CERTIFIED D_j<=U: {le}", flush=True)
