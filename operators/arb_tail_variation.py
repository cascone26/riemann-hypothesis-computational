"""S3027: RIGOROUS bound theta on the variation ACROSS THE GRID of the gamma > T zero tail in the detector steps (lane 2).
Per zero, step profile e_j(gamma) = dnu*h(nu_j), h(nu) = g(nu) + g(-nu-dnu), g(x) = 1/((gamma-x)(gamma-x-dnu)).
h(nu) = H(nu + dnu/2) with H even, so |h(nu) - h(nu')| <= sup|H''| * U^2 (U = V + dnu), and with p = gamma - x >= gamma - V'':
|g''| <= 6/(p - dnu)^4 ⇒ sup|H''| <= 12/(gamma - V'')^4, V'' = V + 2 dnu.  Weight w <= 4/L.  Partial summation with
N(t) <= (t/2pi) log(t/2pi) (t >= 12000, Trudgian) and t/(t-V'') <= k:
 theta <= (4/L) dnu 12 U^2 * (4 k^5 / 2pi) * (3 log(T/2pi) + 1)/(9 T^3)."""
import sys, math
from fractions import Fraction as F
lam2, N, T = sys.argv[1], int(sys.argv[2]), float(sys.argv[3])
q = F(lam2); L = math.log(q.numerator / q.denominator); dnu = 2 * math.pi / L; V = dnu * N
U = V + dnu; V2 = V + 2 * dnu; k = T / (T - V2)
theta = (4 / L) * dnu * 12 * U * U * (4 * k ** 5 / (2 * math.pi)) * (3 * math.log(T / (2 * math.pi)) + 1) / (9 * T ** 3)
theta *= 1.0000001  # float safety margin
print(f"lam2={lam2} N={N} T={T:.0f}: V={V:.3f} dnu={dnu:.4f}  rigorous tail-variation bound theta <= {theta:.4e}")
