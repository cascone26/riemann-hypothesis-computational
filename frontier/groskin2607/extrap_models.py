import numpy as np
# data: log10(min-eig) vs N per c
data={
13:[(8,-22.11),(12,-28.67),(16,-34.07),(20,-38.81),(24,-42.51),(28,-45.96),(32,-48.65),(36,-51.08),(40,-53.02),(44,-54.68)],
30:[(8,-26.73),(12,-35.52),(16,-43.72),(20,-50.31),(24,-56.78),(28,-63.15),(32,-68.36),(36,-73.26),(40,-78.03),(44,-82.50)],
100:[(8,-31.36),(12,-42.66),(16,-52.31),(20,-61.52),(24,-70.45),(28,-79.11),(32,-86.80),(36,-94.50),(40,-101.88),(44,-108.84)],
}
for c,pts in data.items():
    N=np.array([p[0] for p in pts]); y=np.array([p[1] for p in pts])
    # model 1: diverges linearly y=a+b N
    m1=np.polyfit(N,y,1); r1=np.std(np.polyval(m1,N)-y)
    # model 2: diverges as sqrt: y=a+b sqrt(N)
    X2=np.vstack([np.ones_like(N),np.sqrt(N)]).T; c2,_,_,_=np.linalg.lstsq(X2,y,rcond=None); r2=np.std(X2@c2-y)
    # model 3: CONVERGES to floor A: y=A - B N^{-q}, fit via grid on q
    best=None
    for q in np.linspace(0.2,3,80):
        X=np.vstack([np.ones_like(N),-N.astype(float)**(-q)]).T
        cc,_,_,_=np.linalg.lstsq(X,y,rcond=None); r=np.std(X@cc-y)
        if best is None or r<best[0]: best=(r,q,cc[0],cc[1])
    r3,q,A,B=best
    print(f"c={c:3}: linear(diverge) resid={r1:.3f} | sqrt(diverge) resid={r2:.3f} | converge-to-floor A={A:.1f} q={q:.2f} resid={r3:.3f}")
    print(f"        -> best model: {'CONVERGE floor~10^%.0f'%A if r3<min(r1,r2)*0.9 else 'DIVERGE (min-eig->0 at fixed c)' if min(r1,r2)<r3*0.9 else 'AMBIGUOUS'}")
