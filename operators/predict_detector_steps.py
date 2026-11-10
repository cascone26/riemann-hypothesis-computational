# S3021: predict detector steps from the zeros ABOVE the window (explicit-formula tail), compare with Arb-certified steps.
import sys
from flint import arb, acb, ctx
ctx.prec=120
lam2, N, dump, T = sys.argv[1], int(sys.argv[2]), sys.argv[3], float(sys.argv[4])
from fractions import Fraction as F
q=F(lam2); L=(arb(q.numerator)/arb(q.denominator)).log(); PI=arb.pi()
rows=[l.split() for l in open(dump) if not l.startswith('#')]
nu=[arb(r[1]) for r in rows]+[arb(rows[-1][2])]
cert=[float(r[3]) for r in rows]
V=2*PI*N/L; Vp=V+1
n0=int(round(float(arb.zeta_nzeros(Vp).mid()))); n1=int(round(float(arb.zeta_nzeros(arb(T)).mid())))
zs=[z.imag for z in acb.zeta_zeros(n0+1, n1-n0)]
def gt(x):
    s=arb(0)
    for z in zs: s+=(4/L)*(z*L/2).sin()**2*(1/(x-z)+1/(x+z))
    return PI*s
G=[gt(x) for x in nu]
pred=[float((-(G[j+1]-G[j])/PI).mid()) for j in range(len(rows))]
# tail beyond T: d/dnu of sum ~ -sum w*2/gamma^2 (|nu|<<T); mean w = 2/L, density log(g/2pi)/(2pi)
dnu=float((nu[1]-nu[0]).mid()); Tl=__import__('math').log(T/(2*3.141592653589793))
tail=dnu*(2/float(L.mid()))*2*(Tl+1)/(2*3.141592653589793*T)
print(f"lam2={lam2} N={N}: zeros used #{n0+1}..#{n1} (gamma in ({float(Vp.mid()):.2f},{T}]), tail est per step ~{tail:.2e}")
mx=0
for j in range(0,len(rows),max(1,len(rows)//8)):
    print(f"  j={j:3d} nu={float(nu[j].mid()):9.3f} cert={cert[j]:.8e} pred={pred[j]:.8e} cert-pred={cert[j]-pred[j]:+.3e}")
d=[cert[j]-pred[j] for j in range(len(rows))]
print(f"  max|cert-pred|={max(abs(x) for x in d):.4e}  spread(max-min)={max(d)-min(d):.3e}  mean={sum(d)/len(d):.6e}  tail_est={tail:.6e}  mean/tail={sum(d)/len(d)/tail:.5f}")
