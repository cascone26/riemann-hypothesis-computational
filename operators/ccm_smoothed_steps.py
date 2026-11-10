"""S3033: Fejer-smoothed detector-step decomposition (lane 2). Smooth steps with triangular (Fejer) weights of half-width h
grid steps; report, per h, the smoothed components s02, -sR, -sP, -sZ and the residual total, at centre and edges.
Question: does smoothing shrink the O(1) prime oscillation below the explicit (W02 + W_R) part so that positivity of the
smoothed tail could be certified from the explicit pieces + a bound on the smoothed prime sum?"""
import sys, mpmath as mp
sys.path.insert(0,'operators'); from ccm_parts import parts
lam2,N=int(sys.argv[1]),int(sys.argv[2])
W02,WR,WP,L,pw=parts(lam2,N,dps=25,as_mp=True); mp.mp.dps=25
M=2*N+1; nu=[2*mp.pi*n/L for n in range(-N,N+1)]; c=N
gen=lambda X:[mp.mpf(0) if j==c else mp.pi*(nu[c]-nu[j])*X[c,j] for j in range(M)]
g02,gR,gP=gen(W02),gen(WR),gen(WP)
V=nu[-1]; zs=[]; k=1
while True:
    z=mp.im(mp.zetazero(k))
    if z>V+1: break
    zs.append(z); k+=1
gz=[mp.pi*sum((4/L)*mp.sin(z*L/2)**2*(1/(nu[j]-z)+1/(nu[j]+z)) for z in zs) for j in range(M)]
st=lambda g:[float(-(g[j+1]-g[j])/mp.pi) for j in range(M-1)]
S02,SR,SP,SZ=st(g02),[-x for x in st(gR)],[-x for x in st(gP)],[-x for x in st(gz)]
tot=[a+b+cc+d for a,b,cc,d in zip(S02,SR,SP,SZ)]
def sm(a,h,j):
    w=[(h+1-abs(i)) for i in range(-h,h+1)]; num=sum(wi*a[j+i] for wi,i in zip(w,range(-h,h+1))); return num/sum(w)
print(f"lam2={lam2} N={N}: components smoothed with Fejer half-width h (steps); j0=centre, jq=quarter, edge excluded")
for h in [int(x) for x in (sys.argv[3] if len(sys.argv)>3 else "0,1,2,4,8").split(",")]:
    for name,j in [("centre",N),("quarter",N+N//2),("8th",N+N//4)]:
        if j+h>=M-1 or j-h<0: continue
        r=[sm(x,h,j) for x in (S02,SR,SP,SZ,tot)]
        print(f"  h={h:2d} {name:7s} s02={r[0]:+.4e} -sR={r[1]:+.4e} -sP={r[2]:+.4e} -sZ={r[3]:+.4e} total={r[4]:+.4e}  |sP|/total={abs(r[2])/r[4]:.1f}")
