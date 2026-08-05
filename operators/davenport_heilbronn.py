"""Davenport-Heilbronn: functional equation but NO Euler product, and PROVEN off-line zeros.
Test: does the sign-theorem readout FAIL on DH's zeros (esp. off-line ones)? If off-line zeros
break the mu-reading, the readout is RH-SENSITIVE (a criterion), not just consequence-side."""
import mpmath as mp
mp.mp.dps = 25
# character chi mod 5 with chi(2)=i (order-4). chi: 1->1,2->i,3->-i,4->-1,0->0
chiv = {1:1, 2:1j, 3:-1j, 4:-1, 0:0}
def chi(n): return chiv[n % 5]
def L(s, conj=False):
    # L(s,chi) via Hurwitz zeta: L(s,chi)=5^{-s} sum_{a=1}^4 chi(a) zeta(s, a/5)
    tot = mp.mpc(0)
    for a in range(1,5):
        c = mp.conj(chi(a)) if conj else chi(a)
        tot += c*mp.zeta(s, mp.mpf(a)/5)
    return 5**(-s)*tot
kappa = (mp.sqrt(10-2*mp.sqrt(5)) - 2)/(mp.sqrt(5)-1)
def F(s):
    return (1-1j*kappa)/2*L(s) + (1+1j*kappa)/2*L(s, conj=True)
# verify DH has off-line zeros: known off-line zero near s = 0.808517 + 85.699348 i (Titchmarsh)
print("kappa =", mp.nstr(kappa, 10))
print("F(0.5+i*85.7) =", mp.nstr(F(0.5+85.7j), 6))
for test in [mp.mpc('0.808517','85.699348'), mp.mpc('0.5','85.7')]:
    print(f"  |F({mp.nstr(test,8)})| = {mp.nstr(abs(F(test)),4)}")
# find a few zeros on the line (Z-like) at low height + check the known off-line one
print("\nsearching DH zeros on the line (low height):")
zs=[]
prev=F(mp.mpc('0.5','1'))
t=mp.mpf('1')
while t < 40 and len(zs)<15:
    t2=t+mp.mpf('0.1')
    cur=F(mp.mpc('0.5',str(t2)))
    if mp.re(prev)*mp.re(cur)<0 or mp.im(prev)*mp.im(cur)<0:
        try:
            z=mp.findroot(lambda x: F(mp.mpc('0.5',str(x))), float(t))
            if abs(F(mp.mpc('0.5',str(z))))<1e-8 and (not zs or abs(float(z)-zs[-1])>0.01):
                zs.append(float(z))
        except: pass
    prev=cur; t=t2
print("on-line zeros found:", [round(z,3) for z in zs[:12]])
print("DONE")
