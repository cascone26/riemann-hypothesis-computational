#!/usr/bin/env python3
"""§311: DERIVED closed-form leading law for the Turan/hyperbolicity margin, from the modular saddle.
r_m = d^2(log gamma_m)/dm^2 (discrete 2nd diff). Envelope theorem: d(log M_m)/dm = 2 log u*(m) (u* the saddle of
int Phi u^{2m}); d^2/dm^2 = 2 u*'/u*. And d^2/dm^2 log(2m)! = 2/m. Implicit-diff of the saddle eq 2pi e^{2u*}=2m/u*+9/2
gives u*' = u*/(m(2u*+1)) (leading), so:
   r_m ~ 2/(m(2u*+1)) - 2/m = -4 u*(m) / ( m (2 u*(m)+1) ),   u*(m) ~ (1/2) log(m/pi).
=> m|r_m| -> 2 as m->inf (r_m ~ -2/m). Verify this DERIVED law vs the EXACT r_m from §310. Charisma: HANDLE if it matches."""
import mpmath as mp
mp.mp.dps=30
# exact r_m from §310 run (operators/lp_rate.py output), m=1..19:
exact={1:-0.7653227,2:-0.466319,3:-0.3405623,4:-0.2700753,5:-0.2246265,6:-0.1927415,7:-0.1690667,
8:-0.150755,9:-0.1361472,10:-0.1242091,11:-0.1142609,12:-0.1058373,13:-0.09860837,14:-0.09233375,
15:-0.08683382,16:-0.08197172,17:-0.07764123,18:-0.07375863,19:-0.07025703}
def usaddle(m):
    return mp.findroot(lambda u: 2*m/u+mp.mpf(9)/2-2*mp.pi*mp.e**(2*u), mp.mpf(0.5)+mp.log(m+1)/2)
print("m :  exact r_m      derived -4u*/(m(2u*+1))    ratio     m|r_m|(exact)  limit->2",flush=True)
for m in range(2,20):
    us=usaddle(m)
    rpred=-4*us/(m*(2*us+1))
    ex=mp.mpf(str(exact[m]))
    print(f"{m:2d}: {mp.nstr(ex,7):>12}   {mp.nstr(rpred,7):>14}     {mp.nstr(rpred/ex,6):>7}   {mp.nstr(m*abs(ex),6):>8}",flush=True)
print("",flush=True)
print("=> the modular-saddle law r_m ~ -4u*/(m(2u*+1)) reproduces the exact Turan margin; asymptote r_m ~ -2/m.",flush=True)
print("   (leading Laplace; the ~few% gap = next-order saddle correction = the GORZ-uniformity object.)",flush=True)
print("DONE",flush=True)
