import sys, mpmath as mp
sys.path.insert(0,'operators'); from ccm_parts import parts
mp.mp.dps=40
out=[]
for lam2,N in [(10,10),(20,12),(28,12),(36,12)]:
    W02,WR,WP,L,pw=parts(lam2,N,dps=40,as_mp=True); mp.mp.dps=40
    A=W02-WR-WP; M=2*N+1; nu=[2*mp.pi*n/L for n in range(-N,N+1)]; c=N  # center row n=0
    g=[mp.mpf(0) if j==c else mp.pi*(nu[c]-nu[j])*A[c,j] for j in range(M)]
    # consistency: residual of Loewner relation on all off-diagonal pairs
    res=max(abs((g[j]-g[i]) - mp.pi*(nu[i]-nu[j])*A[i,j]) for i in range(M) for j in range(M) if i!=j)
    V=max(abs(x) for x in nu)
    zs=[]; k=1
    while True:
        z=mp.im(mp.zetazero(k))
        if z>V+1: break
        zs.append(z); k+=1
    gt=[g[j]-mp.pi*sum((4/L)*mp.sin(z*L/2)**2*(1/(nu[j]-z)+1/(nu[j]+z)) for z in zs) for j in range(M)]
    steps=[-(gt[j+1]-gt[j])/mp.pi for j in range(M-1)]
    mn=min(steps)
    out.append((lam2,N,len(zs),mp.nstr(res,3),mp.nstr(mn,6),all(s>0 for s in steps)))
    print(f'lam2={lam2} N={N} removed zeros={len(zs)} loewner-residual={mp.nstr(res,3)} min step(-g_tail/pi)={mp.nstr(mn,8)} monotone={all(s>0 for s in steps)}',flush=True)
