import numpy as np, os
import mpmath as mp

# Landau: sum_gamma x^{i gamma} spikes at x = p^m. The string's phase field twists
# this: G(u) = sum_k e^{-2i arg zeta'(rho_k)} e^{i u gamma_k}. Where does the
# TWISTED sum spike? Compare |G_twisted| vs |G_plain| on u in [0.5, 4.5]
# (log2=0.693 ... log90~4.5), all 500 zeros, with prime-power marks.

mp.mp.dps = 16
DIR = os.path.dirname(os.path.abspath(__file__))
Z = np.load(os.path.join(DIR, "zeros500.npy"))
fn = os.path.join(DIR, "argzp500.npy")
if os.path.exists(fn):
    ph = np.load(fn)
else:
    ph = np.empty(len(Z))
    for i, g in enumerate(Z):
        ph[i] = float(mp.arg(mp.zeta(mp.mpc(0.5, float(g)), derivative=1)))
        if i % 50 == 0:
            print(f"  zeta' {i}/{len(Z)}", flush=True)
    np.save(fn, ph)

us = np.linspace(0.5, 4.5, 4001)
E = np.exp(1j * np.outer(us, Z))
Gp = np.abs(E.sum(axis=1)) / len(Z)
Gt = np.abs((E * np.exp(-2j * ph)[None, :]).sum(axis=1)) / len(Z)

def isp(q): return q > 1 and all(q % d for d in range(2, int(q ** 0.5) + 1))
pps = []
for p in [q for q in range(2, 100) if isp(q)]:
    pm = p
    while np.log(pm) < 4.5:
        if np.log(pm) > 0.5: pps.append(pm)
        pm *= p
pps = sorted(set(pps))
marks = np.log(pps)

def peak_report(name, G):
    from scipy.signal import argrelmax
    ix = argrelmax(G, order=25)[0]
    top = ix[np.argsort(G[ix])[::-1][:15]]
    top = np.sort(top)
    print(f"\n{name}: top-15 peaks (u, height, nearest p^m):")
    hits = 0
    for i in top:
        u = us[i]
        j = int(np.argmin(np.abs(marks - u)))
        d = u - marks[j]
        hit = abs(d) < 0.02
        hits += hit
        print(f"  u={u:.3f} h={G[i]:.3f}  nearest log({pps[j]})={marks[j]:.3f} d={d:+.3f} {'HIT' if hit else ''}")
    print(f"  {name} prime-power hits: {hits}/15")
    return hits

h1 = peak_report("PLAIN Landau", Gp)
h2 = peak_report("TWISTED (e^-2i arg zeta')", Gt)
print(f"\nmean level: plain {Gp.mean():.3f}, twisted {Gt.mean():.3f}")
print("if twisted hits vanish: the zeta'-phase field is exactly what hides the primes in the string.")
print("if twisted peaks land elsewhere systematically: NEW spectral comb — chase what those u are.", flush=True)
