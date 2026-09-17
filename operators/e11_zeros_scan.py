import numpy as np, mpmath as mp
mp.mp.dps=30
atil=np.load('E11_an.npy')
N=11;K=2;eps=1;sN=mp.sqrt(N);TP=2*mp.pi;NM=32
araw=[mp.mpf(0)]+[mp.mpf(atil[n])*mp.sqrt(n) for n in range(1,NM+1)]
def ReLam(t):
    s=mp.mpc(1,t);tot=mp.mpc(0)
    for n in range(1,NM+1):
        x=TP*n/sN
        tot+=araw[n]*((x)**(-s)*mp.gammainc(s,x)+eps*(x)**(s-K)*mp.gammainc(K-s,x))
    return float(mp.re(tot))
ts=np.arange(0.3,32,0.1);vals=[ReLam(t) for t in ts]
zeros=[]
for i in range(len(vals)-1):
    if vals[i]*vals[i+1]<0:
        a,b=ts[i],ts[i+1];fa=vals[i]
        for _ in range(30):
            m=(a+b)/2;fm=ReLam(m)
            if fa*fm<0:b=m
            else:a,fa=m,fm
        zeros.append((a+b)/2)
zeros=np.array(zeros)
print(f"E11a on-line zeros in [0,32]: {len(zeros)}",flush=True)
print(f"  {np.round(zeros,4).tolist()}",flush=True)
np.save('/tmp/E11_zeros.npy', zeros)
