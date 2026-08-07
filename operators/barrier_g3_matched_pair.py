"""
G3 (BARRIER_MAP): the matched-local-data pair, realized at finite truncation.

Two Beurling systems discretized by the SAME deterministic greedy rule:
  System A (RH-true target):   Psi_A(x) = x - x0                    on [x0, X]
  System B (off-line zero):    Psi_B(x) = x - x0                    on [x0, Xs]
                               Psi_B(x) = (x-x0) - 2Re[(x^rho - Xs^rho)/rho]
                                                                    on (Xs, X]
  rho = 0.75 + 2i, x0 = 16, split Xs = 1000, X = 1e6.

The zero-at-rho term is restricted to x > Xs; the log-singularity of
int_Xs^inf t^{rho-1-s} dt/... at s = rho comes from the tail at infinity
(verified analytically, BARRIER_G_RESULTS.md G3), so zeta_B still has its zero
pair at rho while zeta_A is zero-free off-line. Greedy discretization is
deterministic, so identical targets on [x0, Xs] give IDENTICAL g-primes there.

Claims verified below:
  (1) g-primes of A and B agree EXACTLY below Xs (matched local data);
  (2) psi_B - x carries the x^{0.75} oscillation above Xs, psi_A - x does not;
  (3) both systems' discretization residual stays O(log x).

Consequence (the barrier): any RH-certificate reading g-prime data up to Xs
returns the same verdict for A and B — and their RH-analogue verdicts differ.
Xs is arbitrary, so no finite-local-data certificate decides RH-analogues.
"""
import numpy as np

sigma0, t0 = 0.75, 2.0
rho = complex(sigma0, t0)
x0, Xs, X = 16.0, 1000.0, 1.0e6

def PsiA(x):
    return x - x0
def dPsiA(x):
    return 1.0
def PsiB(x):
    x = np.asarray(x, dtype=float)
    extra = np.where(x > Xs, -2.0 * ((np.power(x + 0j, rho) - Xs**rho) / rho).real, 0.0)
    return (x - x0) + extra
def dPsiB(x):
    if np.isscalar(x):
        x = float(x)
        return 1.0 if x <= Xs else 1.0 - 2.0 * x**(sigma0 - 1.0) * np.cos(t0 * np.log(x))
    raise ValueError

def discretize(target, dtarget, x_start, x_end):
    primes, S = [], 0.0
    p = x_start * 1.0001
    while True:
        y = max(p * 1.0001, x_start * 1.0001)
        for _ in range(60):
            f = float(target(y)) - S - np.log(y)
            fp = dtarget(y) - 1.0 / y
            step = f / fp if fp != 0 else f
            y_new = y - step
            if y_new <= x_start:
                y_new = 0.5 * (y + x_start)
            if abs(y_new - y) < 1e-9 * y:
                y = y_new
                break
            y = y_new
        if y > x_end:
            break
        primes.append(y)
        S += np.log(y)
        p = y
    return np.array(primes)

print("discretizing system A (no zero) and system B (zero at rho above Xs=%.0f)..." % Xs)
gpA = discretize(PsiA, dPsiA, x0, X)
gpB = discretize(PsiB, dPsiB, x0, X)
print("  A: %d g-primes, B: %d g-primes" % (len(gpA), len(gpB)))

# (1) matched local data below the split
nA = np.searchsorted(gpA, Xs)
nB = np.searchsorted(gpB, Xs)
match = (nA == nB) and np.allclose(gpA[:nA], gpB[:nB], rtol=0, atol=1e-9)
print("\n(1) g-primes below Xs: A has %d, B has %d, identical: %s"
      % (nA, nB, match))
if match and nA > 0:
    print("    max |p_A - p_B| below split: %.2e" % np.abs(gpA[:nA] - gpB[:nB]).max())

# (2)+(3) divergent zero content above the split
logsA, logsB = np.log(gpA), np.log(gpB)
def psi(logs, gp, x):
    return logs[: np.searchsorted(gp, x)].sum()

print("\n(2) psi(x) - (x - x0) above the split:")
print("%10s %14s %14s %14s %12s" % ("x", "A", "B", "B target osc", "B resid"))
bad = False
for x in [2e3, 1e4, 5e4, 2e5, 1e6 * 0.99]:
    dA = psi(logsA, gpA, x) - (x - x0)
    dB = psi(logsB, gpB, x) - (x - x0)
    osc = -2.0 * ((x**rho - Xs**rho) / rho).real
    residB = dB - osc
    print("%10.0f %14.2f %14.2f %14.2f %12.3f" % (x, dA, dB, osc, residB))
    if abs(residB) > 3 * np.log(X) or abs(dA) > 3 * np.log(X):
        bad = True

print("\nRESULT:")
if match and not bad:
    print("  MATCHED PAIR REALIZED: identical g-primes below Xs=%.0f;" % Xs)
    print("  above it, B's psi carries the prescribed off-line-zero oscillation")
    print("  (~x^0.75) while A's stays within O(log x). Any certificate reading")
    print("  local data below Xs cannot distinguish A (RH-true) from B (RH-false).")
else:
    print("  check failed -- do not claim; investigate.")
