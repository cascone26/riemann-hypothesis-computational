# Probe #2b (decisive): is |CROSS_K| < |LOCAL_K| an RH-DETECTOR?  I.e. does the definite-term dominance found on
# real zeros (screw_vs_crossterm.py: ratio 0.15-0.93 in ALL windows) BREAK when a zero goes off the real axis?
# Mechanism = dBN Lehmer-pair pinch-off: a close real pair {g_k,g_{k+1}} merges to c +/- i*delta (complex conjugate).
# If off-line => ratio > 1 (cross overtakes local => monotonicity/positivity fails), then |CROSS|<|LOCAL| <=> reality.
# That would be the deterministic Lambda<=0 handle RT's statistical method lacks. Complex arithmetic throughout;
# the RT energy identity (iii) holds formally over C, we track |.| of each (complex) piece and the REAL parts.
import numpy as np
G = np.loadtxt("/Users/scones/projects/riemann/data/zeros_100k.txt", max_rows=400)

def build_X(gammas, pinch=None, delta=0.0):
    """x_j = 2 gamma_j (+ mirrors). pinch=(i) merges gammas[i],gammas[i+1] to their mean +/- i*delta/2 (dBN pinch)."""
    g = gammas.astype(complex).copy()
    if pinch is not None and delta>0:
        i = pinch; c = 0.5*(g[i]+g[i+1])
        g[i]   = c + 1j*delta/2
        g[i+1] = c - 1j*delta/2
    xp = 2.0*g
    return np.concatenate([-xp[::-1], xp]), len(xp)

def cross_and_local(X, Kidx):
    K = np.array(Kidx); others=[i for i in range(len(X)) if i not in set(K.tolist())]
    xK = X[K]; loc=0j; cross=0j
    for a in range(len(K)):
        for b in range(len(K)):
            if a==b: continue
            d2=(xK[a]-xK[b])**2; s=2.0/d2
            for c2 in range(len(K)):
                if c2==a or c2==b: continue
                s-=1.0/((xK[c2]-xK[a])*(xK[c2]-xK[b]))
            loc+=s*s
            for j in others:
                cross+=4.0/(d2*(xK[a]-X[j])*(xK[b]-X[j]))
    return cross, -2.0*loc

# pick the tightest (Lehmer-like) low pair to pinch: smallest gap among first ~60
gaps = np.diff(G[:60]); ip = int(np.argmin(gaps))
print(f"Tightest low pair: g{ip+1},g{ip+2} = {G[ip]:.4f},{G[ip+1]:.4f}  gap={gaps[ip]:.4f} (mean gap~{gaps.mean():.3f})")
print("Window K = 5 particles centered on the pinched pair. Track |CROSS|/|LOCAL| vs off-line displacement delta.\n")
print(f"{'delta':>7} {'|CROSS|':>12} {'|LOCAL|':>12} {'|CROSS|/|LOCAL|':>15} {'Re(CROSS+LOCAL)':>16} {'monotone? (<0)':>14}")
start=max(0,ip-2)
Kidx=[ (len(G)+j) for j in range(start,start+5) ]   # positive-side indices into X (len X = 2*len(G))
for delta in [0.0, 0.05, 0.1, 0.2, 0.4, 0.8, 1.6, 3.0]:
    X,_ = build_X(G, pinch=ip, delta=delta)
    cross, loc = cross_and_local(X, Kidx)
    ratio = abs(cross)/abs(loc) if abs(loc)>0 else np.inf
    tot = (cross+loc).real
    print(f"{delta:>7.2f} {abs(cross):>12.4e} {abs(loc):>12.4e} {ratio:>15.4f} {tot:>16.4e} {'YES' if tot<0 else 'NO ':>14}", flush=True)
print("\nREAD: delta=0 is the real config (baseline ratio<1). If ratio climbs PAST 1 as delta grows (zero goes")
print("off-line), then |CROSS|<|LOCAL| is an RH-DETECTOR: reality <=> definite-term dominance <=> monotonicity.")
print("If ratio stays <1 regardless, the dominance is generic (not RH-specific) and this isn't the handle.")
