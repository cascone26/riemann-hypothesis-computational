"""Montgomery pair correlation F(alpha) from 100k zeros + the ZERO-PRIME cross-correlation
(the sign theorem's crux, 698b). Reliable numerics. Does F(alpha) match Montgomery's
conjecture 1-(sin pi alpha/pi alpha)^2 for |alpha|<1, and what does the prime cross-term show?"""
import numpy as np
z = np.load("zeros100k.npy")
# unfold: w_n = theta(gamma_n)/pi ~ (gamma/2pi)(log(gamma/2pi)-1); use smooth mean density
g = z
# unfold to unit mean spacing
dens = np.log(g/(2*np.pi))/(2*np.pi)
w = np.cumsum(1/np.gradient(np.arange(len(g))))  # not needed; use direct unfolding
# proper unfolding: N_smooth(t) = (t/2pi)(log(t/2pi)-1) + 7/8
def Nsm(t): return (t/(2*np.pi))*(np.log(t/(2*np.pi))-1)+7/8
wn = Nsm(g)
# pair correlation: F(alpha) = (1/N*) sum_{j,k} w(gamma_j-gamma_k) e^{i alpha (wj-wk)} form
# use the standard: histogram of unfolded differences within a window
diffs=[]
W=40  # neighbor window
for i in range(len(wn)):
    for j in range(i+1, min(i+W, len(wn))):
        d=wn[j]-wn[i]
        if d<6: diffs.append(d)
diffs=np.array(diffs)
# pair correlation R2(s) = histogram normalized
hist,edges=np.histogram(diffs, bins=120, range=(0,6), density=False)
centers=0.5*(edges[:-1]+edges[1:])
R2 = hist/(len(wn))/(edges[1]-edges[0])
# Montgomery/GUE prediction: R2(s) = 1 - (sin(pi s)/(pi s))^2
gue = 1-(np.sin(np.pi*centers)/(np.pi*centers))**2
print("pair correlation vs GUE (s, R2_data, R2_GUE):")
for i in range(0,40,3):
    print(f"  s={centers[i]:.2f}: data {R2[i]:.3f}  GUE {gue[i]:.3f}")
# the small-s repulsion (R2->0) and first peak = GUE signature
print(f"small-s R2(0.1)={R2[2]:.3f} (GUE {gue[2]:.3f}); correlation with GUE over s in [0.2,4]:")
m=(centers>0.2)&(centers<4)
print(f"  corr(data,GUE) = {np.corrcoef(R2[m],gue[m])[0,1]:.4f}")
print("DONE")
