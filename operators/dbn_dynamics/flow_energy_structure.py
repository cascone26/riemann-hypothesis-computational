# Flow energy dE/ds = sum_j v_j^2, v_j=2 sum_{k!=j}1/(g_j-g_k) [= Xi''/Xi'].  Structure:
#  - <v^2> vs height (typical flow speed of a zero, prime-determined)
#  - decompose sum v_j^2 = 4[ two-point Tr: sum_{j!=k}1/(g_j-g_k)^2 ] + 4[ three-point: sum_{k!=l,both!=j} 1/((g_j-g_k)(g_j-g_l)) ]
import numpy as np
G=[]
with open('/Users/scones/projects/riemann/data/zeros_100k.txt') as f:
    for line in f:
        s=line.strip()
        if s:
            try: G.append(float(s.split()[0]))
            except: pass
        if len(G)>=6000: break
G=np.array(sorted(set(G)))
allz=np.concatenate([-G[::-1],G])   # symmetric
def vel_j(j_idx_in_allz):
    gj=allz[j_idx_in_allz]; d=gj-allz; d[j_idx_in_allz]=np.inf; return 2*np.sum(1.0/d)
# <v^2> in windows at different heights
print("height-window   mean|v|    <v^2>    (typical flow speed of a zero vs height)")
for lo,hi in [(0,200),(1000,1200),(3000,3200),(5000,5200)]:
    # indices in allz of the positive zeros G[lo:hi]
    base=len(G)  # allz[base:] = positive G
    vs=[vel_j(base+i) for i in range(lo,hi)]
    vs=np.array(vs); ht=G[(lo+hi)//2]
    print(f"  ~height {ht:7.0f}: mean|v|={np.mean(np.abs(vs)):.4f}   <v^2>={np.mean(vs**2):.4f}   (log-density (logT/2pi)^? )")
# two-point vs three-point for a mid window
base=len(G); idxs=[base+i for i in range(1000,1120)]
tot=0.0; two=0.0
for jj in idxs:
    gj=allz[jj]; d=gj-allz; d[jj]=np.inf; inv=1.0/d
    vj=2*np.sum(inv); tot+=vj**2
    two+=4*np.sum(inv**2)
three=tot-two
print(f"\nmid window (120 zeros): sum v^2={tot:.2f}  = two-point 4*sum 1/gap^2={two:.2f} + three-point={three:.2f}")
print(f"  two-point fraction={two/tot:.3f}, three-point fraction={three/tot:.3f}")
print("\ntwo-point = Tr(L)/... = pair-correlation (close-pair dominated); three-point = triple correlation.")
print("<v^2> scaling with height reveals the prime-determined flow-speed law (density = logT/2pi).")
