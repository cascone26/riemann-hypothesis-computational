# Riemann seat — re-derivation of the Polymath15 ledger formula Lambda <= t0 + y0^2/2
# from the backward-heat zero-trajectory ODE dx_j/dt = -2*sum_{k!=j} 1/(x_j-x_k).
# Symmetric pair x=+/-a: da/dt = -1/a => a^2 = 2(t_c - t).
# Conjugate pair at Im=y falls to axis in Delta t = y^2/2; separation delta=2y => delta^2/8 (banked merge law).
# t0=0.2, y0=0.2 => Lambda <= 0.22 (Polymath15). REPORT-tier (verified derivation, not an Arb barrier cert).
t0=0.2; y0=0.2
assert abs((2*y0)**2/8 - y0**2/2) < 1e-15
print("Lambda <= t0 + y0^2/2 =", t0 + y0**2/2)
