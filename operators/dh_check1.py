import mpmath as mp
mp.mp.dps=15
chiv={1:1,2:1j,3:-1j,4:-1,0:0}
def chi(n): return chiv[n%5]
def Lc(s,conj=False):
    return 5**(-s)*sum((mp.conj(chi(a)) if conj else chi(a))*mp.zeta(s,mp.mpf(a)/5) for a in range(1,5))
kap=(mp.sqrt(10-2*mp.sqrt(5))-2)/(mp.sqrt(5)-1)
def F(s): return (1-1j*kap)/2*Lc(s)+(1+1j*kap)/2*Lc(s,conj=True)
print("kappa=",mp.nstr(kap,6),flush=True)
print("|F(0.8085+85.699i)| =",mp.nstr(abs(F(mp.mpc('0.8085','85.699'))),4),flush=True)
print("|F(0.5+85.699i)|    =",mp.nstr(abs(F(mp.mpc('0.5','85.699'))),4),flush=True)
z=mp.findroot(F,mp.mpc('0.808','85.70'))
print("findroot ->",mp.nstr(z,8)," |F|=",mp.nstr(abs(F(z)),3),flush=True)
