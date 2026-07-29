#!/usr/bin/env python3
"""§406 (handoff iii): the never-built transversal c(delta). Measure the frame lower bound
inf_f [sum_gamma |f(gamma+i*delta)|^2] / ||f||^2 for f in the probe span, vs delta, on the ZERO set.
Fit c(delta): is it ~ (1 - k*delta^2), exp(-c*delta), or worse? This is the empirical shape of the
open harmonic-analysis problem (§396). If c(delta) >= 1 - k*delta^2 cleanly, the transversal bound has
a simple form worth trying to prove; if it collapses, the difficulty is confirmed structural."""
import numpy as np
import l2_edge_sweep as L
w = 3.0
def frame_lb(cs, delta):
    # sum_gamma F(gamma+i*delta) F(gamma+i*delta)^H over probe basis; lower frame bound = lmin of that
    # Gram divided by the overlap Gram (generalized eig) = inf over f of ratio.
    zin = L.ZEROS[(L.ZEROS>=cs.min()-2*w)&(L.ZEROS<=cs.max()+2*w)]
    n = len(cs); G = np.zeros((n,n), dtype=complex)
    for g in zin:
        Fv = np.exp(-((g+1j*delta-cs)**2)/(2*w*w)) + np.exp(-((g+1j*delta+cs)**2)/(2*w*w))
        G += np.outer(Fv, np.conj(Fv))
    G = np.real(G)
    f = L.build_forms(cs, w); Ov = f['overlap']
    ev = np.linalg.eigvalsh(np.linalg.solve(Ov, G))  # generalized: G v = lam Ov v
    return float(ev.real.min())
print("transversal frame lower bound c(delta) = inf_f sum|f(gamma+i delta)|^2 / ||f||^2 (SUB-critical fill).")
print(f"{'s':>5} {'delta=0':>9} {'0.05':>9} {'0.1':>9} {'0.2':>9} {'c/delta^2 @0.1':>13}")
for s in (5.0, 4.5, 4.0):
    cs = L.probes(L.BAND, s)
    c0 = frame_lb(cs, 0.0)
    vals = [frame_lb(cs, d) for d in (0.05,0.1,0.2)]
    print(f"{s:5.1f} {c0:9.3e} {vals[0]:9.3e} {vals[1]:9.3e} {vals[2]:9.3e} {(c0-vals[1])/0.01:13.3e}", flush=True)
print("\nreading: (c0 - c(delta))/delta^2 constant => c(delta)=c0(1-k delta^2), a clean quadratic transversal")
print("bound (provable target). If it blows up/collapses => the open problem is genuinely hard (structural).")
