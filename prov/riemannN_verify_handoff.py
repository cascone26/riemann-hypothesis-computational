#!/usr/bin/env python3
"""Independent verification of §426/§427 handoff (Riemann-N seat).
Consumes /tmp/ops_X{60,100,140}_exact.npy (op-eigenvalues at three truncations).
Checks: (a) delta_n = gamma_op - gamma_true, stats; (b) X-stability of delta_n
(is the O(1) bound a genuine operator property or an extraction artifact?);
(c) the algebraic identity F_op - S_true == N_op - N_true (smooth terms cancel),
i.e. ask(2) reduces to ask(1)."""
import numpy as np

Z = np.load('operators/zeros500.npy')

def Nsmooth(T):
    return (T/(2*np.pi))*np.log(T/(2*np.pi*np.e)) + 7.0/8.0

def match_delta(roots, gmax=95.0):
    """Match each op-root to nearest true zero, return (gamma_true, delta)."""
    zt = Z[Z < gmax]
    out = []
    for r in roots:
        if r < 13.0 or r > gmax:  # ignore below grid start / above range
            continue
        j = np.argmin(np.abs(zt - r))
        out.append((zt[j], r - zt[j]))
    return np.array(out)

Xs = [60, 100, 140]
data = {}
for X in Xs:
    try:
        data[X] = np.load(f'/tmp/ops_X{X}_exact.npy')
    except FileNotFoundError:
        print(f"MISSING /tmp/ops_X{X}_exact.npy");

print("=== op-root counts ===")
for X in Xs:
    if X in data:
        print(f"  X={X}: {len(data[X])} roots in [13.5,96];  true zeros<95 = {int((Z<95).sum())}")

# (a) delta stats per X
print("\n=== delta_n = gamma_op - gamma_true  (per truncation X) ===")
for X in Xs:
    if X not in data: continue
    md = match_delta(data[X])
    g, d = md[:,0], md[:,1]
    # slope vs index and vs gamma
    idx = np.arange(len(d))
    s_idx = np.polyfit(idx, d, 1)[0]
    s_g = np.polyfit(g, d, 1)[0]
    imax = np.argmax(np.abs(d))
    print(f"  X={X}: N={len(d)}  max|d|={np.abs(d).max():.4f} at g={g[imax]:.3f}"
          f"  slope/idx={s_idx:+.4f}  slope/gamma={s_g:+.5f}  mean|d|={np.abs(d).mean():.4f}")

# (b) X-stability: match op-roots ACROSS truncations (same true zero) and see if delta drifts.
print("\n=== X-STABILITY of delta_n (the real question): delta at each true zero, X=60 vs 100 vs 140 ===")
if all(X in data for X in Xs):
    md = {X: dict((round(gt,4), dd) for gt, dd in match_delta(data[X])) for X in Xs}
    common = sorted(set(md[60]) & set(md[100]) & set(md[140]))
    print(f"  {'gamma_true':>10} {'d(X=60)':>9} {'d(X=100)':>9} {'d(X=140)':>9} {'drift 60->140':>13}")
    drifts=[]
    for gt in common:
        d60,d100,d140 = md[60][gt], md[100][gt], md[140][gt]
        drift = d140-d60
        drifts.append(drift)
        print(f"  {gt:10.4f} {d60:+9.4f} {d100:+9.4f} {d140:+9.4f} {drift:+13.4f}")
    drifts=np.array(drifts)
    print(f"\n  |drift 60->140|: max={np.abs(drifts).max():.4f} mean={np.abs(drifts).mean():.4f}")
    print(f"  --> if drift ~ max|delta|, the 'O(1) bound' is X-truncation noise, NOT a fixed operator property.")
    print(f"  --> if drift << max|delta|, delta_n is X-stable and the O(1) bound is real.")

# (c) algebraic reduction: F_op - S_true = N_op - N_true  (Nsmooth cancels identically)
print("\n=== ask(2) reduces to ask(1): F_op(T)-S_true(T) == N_op(T)-N_true(T) (smooth term cancels) ===")
if 100 in data:
    roots = data[100]
    grid = np.linspace(14.0, 95.0, 400)
    Ns = Nsmooth(grid)
    Nop = np.searchsorted(np.sort(roots), grid, side='right')
    Ntr = np.searchsorted(Z, grid, side='right')
    Fop = Nop - Ns; Str = Ntr - Ns
    lhs = Fop - Str; rhs = Nop - Ntr
    print(f"  max|(F_op-S_true) - (N_op-N_true)| = {np.max(np.abs(lhs-rhs)):.2e}  (should be ~0: identity)")
    print(f"  So comparing F_op to S_true is IDENTICALLY comparing the two counting functions;")
    print(f"  the smooth <N> term is common and cancels -> NOT an independent physical 'decider'.")
    print(f"  RMS(N_op - N_true) = {np.sqrt(np.mean(rhs**2)):.4f} integer-step counting discrepancy")
    print(f"  max|N_op - N_true| = {int(np.max(np.abs(rhs)))} (# of order-1 count displacements)")
