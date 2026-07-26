#!/usr/bin/env python3
"""§318: DERIVE tau_d from first-order root dynamics. Flow c_k -> c_k e^{t g(k)}:
dp/dt|_{t=0} = G(x d/dx) p with G(k)=g(k) = (k-h)^3 - (3h^2/4)(k-h), h=d/2.
Root velocity: v_i = -[G(x∂)p](x_i)/p'(x_i). Predicted collapse: t* = min over adjacent pairs of
(x_{i+1}-x_i)/max(v_i - v_{i+1},0)  (pair closing). Compare t* to measured tau_d (0.6/d^3 law) and
report WHICH pair binds (edge or bulk). First-order theory - expect O(1) agreement, not exact."""
import mpmath as mp
mp.mp.dps=40
def hermite_coeffs(d):
    P=[mp.mpf(1)]; Q=[mp.mpf(0),mp.mpf(1)]
    for k in range(1,d):
        R=[mp.mpf(0)]+Q
        for i in range(len(P)): R[i]-=k*P[i]
        P,Q=Q,R
    return Q
print(" d :  t*_pred     tau_meas    ratio   binding pair (i from edge)",flush=True)
meas={8:0.0023437,10:0.00087473,12:0.00042793,14:0.00024488,16:0.00015518,18:0.00010561,20:7.5769e-5}
for d in [8,10,12,14,16,18,20]:
    co=hermite_coeffs(d); h=mp.mpf(d)/2
    g=lambda k:(k-h)**3-(3*h*h/4)*(k-h)
    Gco=[c*g(k) for k,c in enumerate(co)]            # G(x∂)p coefficients
    dco=[co[k]*k for k in range(1,len(co))]          # p' coefficients (low->high)
    rts=sorted(mp.re(r) for r in mp.polyroots(list(reversed(co)),maxsteps=300,extraprec=150))
    def ev(cs,x): return sum(c*x**j for j,c in enumerate(cs))
    v=[-ev(Gco,x)/ev(dco,x) for x in rts]
    best=None
    for i in range(len(rts)-1):
        close=v[i]-v[i+1]     # positive => pair closing
        if close>0:
            t=(rts[i+1]-rts[i])/close
            if best is None or t<best[0]: best=(t,i)
    t_pred,i0=best
    edge_idx=min(i0,d-2-i0)   # distance of binding pair from nearest edge
    tm=mp.mpf(str(meas[d]))
    print(f"{d:2d}: {mp.nstr(t_pred,5):>10} {mp.nstr(tm,5):>10}  {mp.nstr(t_pred/tm,4):>6}   pair#{i0} ({edge_idx} from edge)",flush=True)
print("",flush=True)
print("if ratio ~const and binding pair is at fixed edge-distance: mechanism identified;",flush=True)
print("then tau_d ~ [edge gap]/[edge closing rate] -> derive 0.6/d^3 from Hermite edge asymptotics.",flush=True)
print("DONE",flush=True)
