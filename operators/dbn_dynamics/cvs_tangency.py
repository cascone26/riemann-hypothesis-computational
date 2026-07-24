#!/usr/bin/env python3
"""
cvs_tangency.py — LEDGER §166: the §154 TANGENCY test on the CORRECT (connes-cvs) operator.

Our exploratory finding (§154, on the sine-Galerkin Weil operator): with g_t = g_arch + t*g_prime,
lambda_min(t) is concave and peaks at the TRUE prime coupling t=1 — a cusp/eigenvalue-crossing (§155).
Question now: does that survive on the EXACT CvS truncated Weil operator (the frontier object whose
ground state reproduces the Riemann zeros to 1e-17)?

Structural claim under test: lambda_min(t) <= 0 for all t, tangent to 0 EXACTLY at t=1 (the arithmetic).
If true, it is a variational characterization: the true von Mangoldt coupling is the unique value where
the Weil form is marginally positive.

EFFICIENCY: only psi_prime scales with t. psi_arch (expensive digamma quadrature) + psi_pole are
t-independent, so we compute the three components ONCE per basis index, then the whole t-sweep is
cheap reassembly + even-sector eigsy. Cost ~= one matrix build.

Usage: cvs_tangency.py [N] [dps] [T] [c]
"""
import sys
import mpmath as mp
sys.path.insert(0, ".venv_cvs/lib/python3.14/site-packages")
import connes_cvs.operator as O

N   = int(sys.argv[1]) if len(sys.argv) > 1 else 30
dps = int(sys.argv[2]) if len(sys.argv) > 2 else 60
T   = int(sys.argv[3]) if len(sys.argv) > 3 else 200
c   = int(sys.argv[4]) if len(sys.argv) > 4 else 13

mp.mp.dps = dps
L = mp.log(c)
prime_data, _ = O.prime_powers_up_to(int(c))

print(f"§166 CvS tangency test — c={c} N={N} dps={dps} T={T}  (L=log c={float(L):.5f})")
print(f"decompose psi = t*prime + pole + arch; sweep t; watch lambda_min(t) for a peak/cusp at t=1")
print("precomputing t-independent components (arch is the slow part)...", flush=True)

# --- precompute the three psi-components (value + derivative) per basis index, arch paid ONCE ---
idx = list(range(-N, N + 1))
pri  = {}; prid = {}
pol  = {}; pold = {}
arc  = {}; arcd = {}
for n in idx:
    x = mp.mpf(n)
    pri[n]  = O.psi_prime(x, L, prime_data)
    prid[n] = O.psi_prime_deriv(x, L, prime_data)
    pol[n]  = O.psi_pole(x, L)
    pold[n] = O.psi_pole_deriv(x, L)
    O._hplus_cache_clear()
    arc[n]  = O.psi_arch(x, L, T, dps)
    arcd[n] = O.psi_arch_deriv(x, L, T, dps)
    O._hplus_cache_clear()
print("  ...components ready. sweeping t.\n", flush=True)

DIM = 2 * N + 1

def lam_min(t):
    """lambda_min of the even-sector CvS operator with prime piece scaled by t."""
    psi  = {n: t * pri[n]  + pol[n]  + arc[n]  for n in idx}
    psid = {n: t * prid[n] + pold[n] + arcd[n] for n in idx}
    Q = mp.matrix(DIM, DIM)
    for i in range(DIM):
        m = i - N
        for j in range(DIM):
            nn = j - N
            Q[i, j] = psid[nn] if m == nn else (psi[m] - psi[nn]) / (m - nn)
    for i in range(DIM):
        for j in range(i + 1, DIM):
            a = (Q[i, j] + Q[j, i]) / 2
            Q[i, j] = a; Q[j, i] = a
    lam, _ = O.compute_ground_state(Q)
    return lam

# sweep: dense near t=1 (the claimed tangency point), plus wider context
ts = [mp.mpf(s) for s in
      ["0", "0.5", "0.9", "0.99", "0.9999", "0.99999", "1", "1.00001", "1.0001", "1.01", "1.1", "1.5", "2"]]
print(f"{'t':>8}   {'lambda_min':>16}   sign")
rows = []
for t in ts:
    lm = lam_min(t)
    sign = "+" if lm > 0 else ("0" if lm == 0 else "-")
    rows.append((t, lm))
    print(f"{float(t):>8.4f}   {mp.nstr(lm, 8):>16}   {sign}", flush=True)

print("\nREAD:")
print(" - peak of lambda_min(t) exactly at t=1 (and only there) => tangency / arithmetic is the critical coupling")
print(" - lambda_min<0 for t!=1, ~0 at t=1 => variational characterization of the primes")
print(" - sharp corner (not smooth max) at t=1 => cusp / eigenvalue crossing (matches §155)")
