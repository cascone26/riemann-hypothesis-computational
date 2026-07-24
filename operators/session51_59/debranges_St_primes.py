# De Branges step 2: the arithmetic content of the phase = S(t) = N(t) - theta(t)/pi - 1.
# Explicit formula predicts S(t) ~ -(1/pi) sum_p sum_k sin(t k log p)/(k p^{k/2}). Extract the primes:
# correlate S(t) against sin(t log p) over a t-window; recovered amplitude should ~ 1/(pi sqrt p).
import mpmath as mp
mp.mp.dps=20
theta=mp.siegeltheta
def N(t): return mp.nzeros(t)          # # zeros with 0<gamma<=t
def S(t): return N(t) - theta(t)/mp.pi - 1
# sample S(t) on a fine grid over a window well above small-t transients
a,b=20.0,80.0; M=1600
ts=[a+(b-a)*i/M for i in range(M+1)]
Sv=[float(S(mp.mpf(t))) for t in ts]   # S is a step-ish function; sample it
# subtract mean
mean=sum(Sv)/len(Sv); Sv=[v-mean for v in Sv]
dt=(b-a)/M
def corr(freq):  # projection of S onto sin(freq*t) and cos(freq*t), return amplitude
    cs=sum(Sv[i]*mp.cos(freq*ts[i]) for i in range(M+1))*dt
    sn=sum(Sv[i]*mp.sin(freq*ts[i]) for i in range(M+1))*dt
    return float(mp.sqrt(cs*cs+sn*sn))*2/(b-a)
primes=[2,3,5,7,11,13]
print("prime  log p    recovered-amp   predicted 1/(pi*sqrt p)   ratio")
for p in primes:
    amp=corr(mp.log(p)); pred=1/(mp.pi*mp.sqrt(p))
    print(f"  {p:2d}  {float(mp.log(p)):.4f}    {amp:.5f}          {float(pred):.5f}            {amp/float(pred):.2f}")
# control: non-log-prime frequencies should give small amplitude
print("\ncontrols (non-log-prime freqs):")
for f in [0.5,1.5,2.5]:
    print(f"  freq={f}: amp={corr(mp.mpf(f)):.5f}")
