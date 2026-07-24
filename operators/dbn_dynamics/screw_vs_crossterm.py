# Probe #2 (BRIDGE.md §4): does screw/Weil POSITIVITY control RT's indefinite CROSS-TERM deterministically,
# replacing RT's statistical (pair-correlation) control? RT energy identity (iii), arXiv:1801.05914 line 758:
#   d/dt Sum_{k!=k' in K} E_{kk'} = CROSS_K  -  2*Sum_{k!=k' in K}( 2/(x_k-x_k')^2 - Sum_{k'' in K}1/((x_k''-x_k)(x_k''-x_k')) )^2
#   CROSS_K = Sum_{j notin K, k!=k' in K} 4/((x_k-x_k')^2 (x_k-x_j)(x_k'-x_j))     [INDEFINITE - RT control statistically]
# x_j = 2*gamma_j (RT scaling), and the FULL symmetric config includes negatives x_{-j} = -2*gamma_j.
# Question: is |CROSS_K| bounded by a screw-positive quantity for EVERY window K? If yes -> deterministic handle (crack).
# If no / needs the tail -> confirms the infinite-system frontier (89/91) is genuine. Either way = real result.
import numpy as np

# --- real zeta zeros ---
G = np.loadtxt("/Users/scones/projects/riemann/data/zeros_100k.txt", max_rows=400)   # gamma_j > 0
# RT particle positions: x_j = 2 gamma_j, plus mirror negatives. Index by signed j.
xs_pos = 2.0*G
X = np.concatenate([-xs_pos[::-1], xs_pos])          # full symmetric configuration, sorted
Npos = len(xs_pos)
def idx_pos(j): return Npos + j                      # j=0..Npos-1 -> positive side index into X

def cross_and_local(Kidx, envcap=None):
    """RT identity (iii) RHS pieces for finite window K (indices into X). env = all other particles (capped)."""
    K = np.array(Kidx)
    others = np.array([i for i in range(len(X)) if i not in set(K.tolist())])
    if envcap is not None:
        # nearest envcap environment particles to the window centroid (RT truncate; we test tail-sensitivity)
        c = X[K].mean(); others = others[np.argsort(np.abs(X[others]-c))[:envcap]]
    xK = X[K]
    # local definite term  -2 * sum_{k!=k'} ( 2/(x_k-x_k')^2 - sum_{k'' in K} 1/((x_k''-x_k)(x_k''-x_k')) )^2
    loc = 0.0
    for a in range(len(K)):
        for b in range(len(K)):
            if a==b: continue
            d2 = (xK[a]-xK[b])**2
            s = 2.0/d2
            for c2 in range(len(K)):
                if c2==a or c2==b: continue
                s -= 1.0/((xK[c2]-xK[a])*(xK[c2]-xK[b]))
            loc += s*s
    loc *= -2.0
    # cross term  sum_{j notin K, k!=k'} 4/((x_k-x_k')^2 (x_k-x_j)(x_k'-x_j))
    cross = 0.0
    for a in range(len(K)):
        for b in range(len(K)):
            if a==b: continue
            d2 = (xK[a]-xK[b])**2
            for j in others:
                cross += 4.0/(d2*(xK[a]-X[j])*(xK[b]-X[j]))
    return cross, loc

# --- screw-positivity margin of the window: the Weil/screw form Sum_gamma (2/gamma^2)|vhat(gamma)|^2 restricted ---
# Under RH the screw g(t)=c+Sum_gamma (2/gamma^2)cos(gamma t) (LEDGER 82). The window's positive screw "budget"
# is P_K = Sum_{gamma in window} 2/gamma^2  (the total positive atomic mass the window carries). This is the
# deterministic positivity resource. Test whether |CROSS_K| <= const * P_K * (local scale) uniformly.
def screw_budget(Kpos_gammas):  # gammas (positive) in the window
    return np.sum(2.0/np.asarray(Kpos_gammas)**2)

print("Probe #2: RT cross-term vs screw-positivity budget, sliding windows on real zeta zeros.")
print("x_j=2*gamma_j. K = m consecutive POSITIVE-side particles. env truncation tests tail-sensitivity.\n")
print(f"{'window(gammas)':>16} {'m':>2} {'CROSS_K':>12} {'LOCAL':>12} {'|CROSS|/|LOCAL|':>15} {'screwP_K':>10} {'|CROSS|/screwP':>14}")
for start in [0, 5, 20, 50, 100]:
    for m in [3, 5, 8]:
        Kidx = [idx_pos(j) for j in range(start, start+m)]
        gammas = G[start:start+m]
        # tail sensitivity: full env vs capped env
        cross_full, loc = cross_and_local(Kidx, envcap=None)
        cross_cap,  _   = cross_and_local(Kidx, envcap=40)
        Pk = screw_budget(gammas)
        ratio_loc = abs(cross_full)/abs(loc) if loc!=0 else np.inf
        ratio_scr = abs(cross_full)/Pk
        tail = abs(cross_full-cross_cap)/max(abs(cross_full),1e-12)
        print(f"g{start+1}-g{start+m:>3} {m:>2} {cross_full:>12.4e} {loc:>12.4e} {ratio_loc:>15.4f} "
              f"{Pk:>10.4e} {ratio_scr:>14.4e}  tailΔ={tail:.2%}", flush=True)

print("\nREAD: (a) is |CROSS|<|LOCAL| (definite part dominates -> monotone WITHOUT tail control)?")
print("     (b) does |CROSS|/screwP stay BOUNDED across windows/heights (deterministic positivity handle)?")
print("     (c) tailΔ small => cross-term is LOCAL (crack plausible); large => needs infinite tail (RT frontier real).")
