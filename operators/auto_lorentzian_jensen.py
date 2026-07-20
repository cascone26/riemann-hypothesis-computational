import mpmath as mp, numpy as np, math
mp.mp.dps=60
def xi(s): return 0.5*s*(s-1)*mp.power(mp.pi,-s/2)*mp.gamma(s/2)*mp.zeta(s)
f=lambda t: xi(mp.mpf('0.5')+1j*t)
ORD=60
try:
    coef=mp.taylor(f,0,ORD,method='quad',radius=5)
except Exception:
    coef=mp.taylor(f,0,ORD)
Nmax=ORD//2
a=[]
for n in range(Nmax+1):
    cn=mp.re(coef[2*n]); an=((-1)**n)*mp.factorial(2*n)*cn; a.append(an)
af=np.array([float(x) for x in a])
print('a_0..a_5:',[float(x) for x in a[:6]])
print('all a_n>0:',bool(np.all(af>0)))
# Turan / Newton (Hodge-Riemann rank-2) margins
tur=[float(a[n]*a[n]-a[n-1]*a[n+1]) for n in range(1,Nmax)]
print('min Turan margin a_n^2-a_{n-1}a_{n+1} =',min(tur))
# normalized Hankel positive-definiteness
m=min(Nmax//2,14)
H=np.zeros((m+1,m+1))
for i in range(m+1):
    for j in range(m+1): H[i,j]=af[i+j]
D=np.sqrt(np.abs(np.diag(H))); Hn=H/np.outer(D,D)
print('normalized Hankel(a) min eig =',np.linalg.eigvalsh(Hn).min())
# Jensen hyperbolicity across shifts/degrees (Lorentzian => real-rooted)
maxim=0.0; checked=0
for d in range(1,11):
    for n in range(0,min(Nmax-d,28)):
        cc=[mp.binomial(d,k)*a[n+k] for k in range(d+1)]
        hi=[cc[d-k] for k in range(d+1)]
        try:
            rts=mp.polyroots(hi,maxsteps=120,extraprec=60)
        except Exception:
            continue
        im=max(abs(float(mp.im(z))/(1+abs(float(mp.re(z))))) for z in rts)
        maxim=max(maxim,im); checked+=1
print('Jensen polys checked:',checked,' max rel |Im root| (0 => all hyperbolic):',maxim)
if np.all(af>0) and min(tur)>0 and np.linalg.eigvalsh(Hn).min()>-1e-9 and maxim<1e-6:
    print('LORENTZIAN_SEQUENCE_CONFIRMED_ON_RANGE: Turan>0, Hankel PSD, Jensen hyperbolic -> Hodge-Riemann signature holds on computable range')
