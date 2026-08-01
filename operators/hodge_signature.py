import numpy as np
from l2_edge_sweep import probes, build_forms

# Hodge-index signature hunt: in FF-land RH follows from the intersection form on
# Num(CxC) being Lorentzian (1, rho-1) plus one known positive class. Number-field
# analogue: inertia (n+, n0, n-) of the explicit-formula form's pieces on prime-
# localized probe spaces. Question: does the GEOMETRIC side (archpole + prime,
# i.e. full minus the zero term) show the Lorentzian "exactly one positive
# direction" pattern, stably across band and probe count?

def inertia(M, tol=1e-9):
    ev = np.linalg.eigvalsh(M)
    s = tol * max(1.0, np.abs(ev).max())
    return int(np.sum(ev > s)), int(np.sum(np.abs(ev) <= s)), int(np.sum(ev < -s)), ev

print(f"{'band':>12} {'w':>3} {'n':>3} | {'full':>10} | {'geom=arch+pole+prime':>20} | {'archpole':>10} | {'-prime':>10} | {'zero':>10}")
for band, w, s in [((14, 44), 3.0, 2.0), ((14, 44), 3.0, 1.0), ((10, 60), 3.0, 2.0),
                   ((20, 80), 4.0, 2.5), ((14, 44), 2.0, 1.5), ((30, 90), 3.0, 2.0)]:
    cs = probes(band, s)
    F = build_forms(cs, w)
    geom = F["full"] - F["zero"]
    rows = {}
    for k, M in [("full", F["full"]), ("geom", geom), ("archpole", F["archpole"]),
                 ("-prime", -F["prime"]), ("zero", F["zero"])]:
        p, z, m, ev = inertia(M)
        rows[k] = (p, z, m, ev)
    n = len(cs)
    fmt = lambda k: f"({rows[k][0]},{rows[k][1]},{rows[k][2]})"
    print(f"{str(band):>12} {w:>3.0f} {n:>3} | {fmt('full'):>10} | {fmt('geom'):>20} | {fmt('archpole'):>10} | {fmt('-prime'):>10} | {fmt('zero'):>10}")
    if band == (14, 44) and s == 2.0:
        ev = rows["geom"][3]
        print(f"    geom eigenvalues: top3={np.round(ev[-3:], 4)}  bottom3={np.round(ev[:3], 4)}")
        v = np.linalg.eigh(geom)[1][:, -1]
        print(f"    top geom eigenvector (positive direction) coeffs vs probe centers:")
        print("    " + "  ".join(f"{c:.0f}:{x:+.2f}" for c, x in zip(cs, v)))
print("\nFF prediction if the 'missing surface' shadows through: geom side Lorentzian-like —")
print("a FIXED small number of positive directions (ideally 1 = the pole class) as n grows.")
