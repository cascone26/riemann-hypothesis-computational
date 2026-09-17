import numpy as np
from scipy.special import digamma
from scipy.linalg import eigh
def vonmangoldt(N):
    lam=np.zeros(N+1);isc=np.zeros(N+1,dtype=bool)
    for p in range(2,N+1):
        if not isc[p]:
            isc[2*p::p]=True;lp=np.log(p);pk=p
            while pk<=N: lam[pk]=lp;pk*=p
    return lam
NMAX=3000;LAM=vonmangoldt(NMAX);NS=np.arange(2,NMAX+1);LOGN=np.log(NS);LAMW0=LAM[2:]/np.sqrt(NS)
a,b,W=14.0,44.0,1.0;fourW=4*W;dr=0.004
RGRID=np.arange(-(b+120),(b+120),dr);OMEGA=np.real(digamma(0.25+0.5j*RGRID))-np.log(np.pi)
def fej(x):
    y=W*x;out=np.ones_like(y);nz=np.abs(y)>1e-12;out[nz]=(np.sin(y[nz])/y[nz])**2;return out
def fc(c):
    y=W*(0.5j-c);y2=W*(0.5j+c);return (np.sin(y)/y)**2+(np.sin(y2)/y2)**2
mask=LOGN<fourW;ln=LOGN[mask];nn=NS[mask];w0=LAMW0[mask];cosu=np.cos(np.outer(ln,RGRID))
print("Stability of RH-critical near-null mode vs N (centers). small-prime dominance persists?")
for N in [13,20,30,40]:
    CEN=np.linspace(a,b,N)
    P=np.array([fej(RGRID-c)+fej(RGRID+c) for c in CEN]);Gm=(P@P.T)*dr
    arch=(P*OMEGA[None,:])@P.T*dr/(2*np.pi)
    fp=np.array([fc(c) for c in CEN]);pole=2.0*np.real(np.outer(fp,fp))
    Mn=[]
    for jj in range(len(ln)):
        col=cosu[jj];Mk=np.zeros((N,N))
        for j in range(N): Mk[j,:]=2.0*(((P[j][None,:]*P)@col)*dr/(2*np.pi))
        Mn.append(0.5*(Mk+Mk.T))
    prime=np.tensordot(w0,np.array(Mn),axes=1)
    try:
        w,V=eigh(pole+arch-prime,Gm);v=V[:,0];Qv=w[0]
    except Exception:
        print(f"  N={N:3d}: SINGULAR G (overcomplete basis)",flush=True);continue
    av=float(v@arch@v);pv=float(v@prime@v)
    contribs=[(int(nn[jj]),float(v@Mn[jj]@v)*w0[jj]) for jj in range(len(ln))]
    contribs.sort(key=lambda t:-abs(t[1]));tot=pv
    top=", ".join(f"n{n}:{100*c/tot:.0f}%" for n,c in contribs[:4])
    print(f"  N={N:3d}: Q(v)={Qv:+.3e} arch={av:.4f} prime={pv:.4f} arch-prime={av-pv:+.3e} | top: {top}",flush=True)
