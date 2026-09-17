import sys,os,math
sys.path.insert(0, os.path.join('..','external','groskin_2607.02828','code'))
from flint import arb, ctx
import arb_ldlt_certify as G
ctx.prec=200
Cat=0.915965594177219015   # Catalan
# numeric (validated) (kappa+J)_DH targets, per c
data={14:0.5957, 25:0.61048, 50:0.6174, 75:0.61916, 100:0.61984}
print("Test closed form: (k+J)_DH ?= [arb_kappa+arb_J](L) + pi - 8*Catalan/L - log5")
print(f"{'c':>4} {'L':>7} {'G14=k+J(1/4)':>13} {'predicted':>11} {'numeric':>9} {'diff':>10}")
for c,num in data.items():
    L=arb(c).log(); Lf=float(L.mid())
    g14=float((G.arb_kappa(L)+G.arb_J(L)).mid())
    pred=g14 + math.pi - 8*Cat/Lf - math.log(5)
    print(f"{c:>4} {Lf:>7.4f} {g14:>13.5f} {pred:>11.5f} {num:>9.5f} {pred-num:>+10.5f}")
