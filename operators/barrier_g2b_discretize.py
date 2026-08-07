"""
G2b/G3 (BARRIER_MAP): discretize a continuous Beurling system with a PRESCRIBED
off-line zero into an actual g-prime multiset, and verify the witness signal.

Continuous target (bad member of the G2 family, now at sigma0=0.75 so the prime
density is positive for x >= x0 = 2^{1/(1-sigma0)} = 16):
  d Psi(x) = (1 - 2 x^{sigma0-1} cos(t0 log x)) dx   on [x0, X],
  Psi(x) = x - 2 Re[(x^rho - x0^rho)/rho] - (x0 - ...)   (exact integral below),
  rho = sigma0 + i t0.  Its zeta has a zero pair at rho, rho_bar; the explicit-
  formula shadow is the x^{sigma0} oscillation in Psi(x) - x.

Greedy discretization (deterministic, error <= log x — sharper than BDR's
random Thm 1.2, which we cite for the N(x) side):
  place g-primes p_1 < p_2 < ... so that after each placement
  sum_{i<=n} log p_i tracks Psi(p_n):  solve Psi(p_{n+1}) = S_n + log p_{n+1}
  by Newton (Psi' = density > 0 on [x0, X]).

Verification (the point):
  psi_P(x) - x for the DISCRETE multiset must reproduce -2Re[x^rho/rho] —
  an x^{0.75}-size oscillation with period 2*pi/t0 in log x — far above the
  O(log x) discretization error. A control system (density dx, no zero) must
  show error O(log x) with NO such oscillation.

G3 tie-in: a second discretization sharing the target on [x0, X_split] and
diverging above realizes the matched-local-data pair at finite truncation.
"""
import numpy as np

sigma0, t0 = 0.75, 2.0
rho = complex(sigma0, t0)
x0 = 2.0 ** (1.0 / (1.0 - sigma0))     # = 16: density nonnegative from here
X = 1.0e6

def Psi(x):
    # integral of (1 - 2 u^{sigma0-1} cos(t0 log u)) du from x0 to x
    # = (x - x0) - 2*Re[ (x^rho - x0^rho)/rho ]
    return (x - x0) - 2.0 * ((x**rho - x0**rho) / rho).real

def dPsi(x):
    return 1.0 - 2.0 * x**(sigma0 - 1.0) * np.cos(t0 * np.log(x))

def discretize(target, dtarget, x_start, x_end):
    """Greedy: Psi(p_next) = S + log(p_next), Newton from previous prime."""
    primes = []
    S = 0.0
    p = x_start * 1.0001
    while True:
        # solve f(y) = target(y) - S - log(y) = 0, f' = dtarget - 1/y
        y = max(p * 1.0001, x_start * 1.0001)
        for _ in range(60):
            f = target(y) - S - np.log(y)
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

print("discretizing witness system (prescribed zero rho = %.2f + %.1fi) ..." % (sigma0, t0))
gp = discretize(Psi, dPsi, x0, X)
print("  g-primes generated: %d (up to X=%.0e)" % (len(gp), X))

print("discretizing control system (no zero, density du) ...")
gp_ctrl = discretize(lambda x: x - x0, lambda x: 1.0, x0, X)
print("  control g-primes: %d" % len(gp_ctrl))

# --- verify the witness signal ------------------------------------------------
logs = np.log(gp)
logs_ctrl = np.log(gp_ctrl)
def psi_discrete(logs_sorted, gp_sorted, x):
    return logs_sorted[: np.searchsorted(gp_sorted, x)].sum()

print("\n%10s %14s %14s %12s %12s" % ("x", "psi_P-x", "target osc", "resid", "ctrl err"))
sample = np.exp(np.linspace(np.log(50), np.log(X * 0.99), 25))
resids, oscs = [], []
for x in sample:
    d = psi_discrete(logs, gp, x) - (x - x0)          # discrete psi minus linear part
    osc = -2.0 * ((x**rho - x0**rho) / rho).real       # prescribed-zero oscillation
    resid = d - osc
    dc = psi_discrete(logs_ctrl, gp_ctrl, x) - (x - x0)
    resids.append(resid); oscs.append(osc)
    print("%10.0f %14.2f %14.2f %12.3f %12.3f" % (x, d, osc, resid, dc))

resids = np.array(resids); oscs = np.array(oscs)
amp = np.abs(oscs)[-5:].mean()
print("\nRESULT:")
print("  oscillation amplitude at top: %.1f  (~x^0.75 as prescribed)" % amp)
print("  max |residual| after removing prescribed-zero term: %.3f (vs log X = %.1f)"
      % (np.abs(resids).max(), np.log(X)))
if np.abs(resids).max() < 2 * np.log(X) and amp > 100 * np.abs(resids).max():
    print("  => DISCRETE g-prime system CARRIES the off-line-zero signal: psi_P(x)-x")
    print("     equals -2Re[x^rho/rho] up to discretization error O(log x).")
    print("     The G2/G3 witness machinery is now handled at finite truncation,")
    print("     not just cited.")
else:
    print("  => signal/error separation not clean; investigate before claiming.")
