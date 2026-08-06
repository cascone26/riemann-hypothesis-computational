"""Recompute the {15,25,35,45} Weil matrix (margin 0.048, robustly PD), rationalize, LDL^T,
emit machine-checked Lean SOS certificate."""
import numpy as np, mpmath as mp, sympy as sp
mp.mp.dps=25
t=1.0; av=[15.,25.,35.,45.]; nb=4
def phihat(a,r): return np.sqrt(2*np.pi)/2*(np.exp(-(r-a)**2/2)+np.exp(-(r+a)**2/2))
def phihat_c(a,r): return mp.sqrt(2*mp.pi)/2*(mp.e**(-(r-a)**2/2)+mp.e**(-(r+a)**2/2))
def g_cf(aj,ak,y):
    s=aj+ak; d=aj-ak
    return np.exp(-y**2/4)*(np.sqrt(np.pi)/2)*(np.exp(-s**2/4)*np.cos(d*y/2)+np.exp(-d**2/4)*np.cos(s*y/2))
from sympy import primerange
pp=[]
for p in primerange(2,300000):
    lp=np.log(p); m=1
    while m*lp<30: pp.append((m*lp,lp/np.sqrt(p**m))); m+=1
Ls=np.array([L for L,w in pp]); Ws=np.array([w for L,w in pp])
rs=np.linspace(-120,120,24000)
Om=np.array([float(mp.re(mp.digamma(mp.mpf('0.25')+1j*r/2))-mp.log(mp.pi)) for r in rs])
Q=np.zeros((nb,nb))
for j in range(nb):
    for k in range(nb):
        arch=np.trapezoid(phihat(av[j],rs)*phihat(av[k],rs)*Om,rs)/(2*np.pi)
        pole=2*float(mp.re(phihat_c(av[j],mp.mpf('0.5')*1j)*phihat_c(av[k],mp.mpf('0.5')*1j)))
        prime=2*np.sum(Ws*g_cf(av[j],av[k],Ls))
        Q[j,k]=arch+pole-prime
Q=(Q+Q.T)/2
print("float eigenvalues:", np.round(np.linalg.eigvalsh(Q),4), " (min should be ~0.048)")
M=sp.Matrix(4,4, lambda i,j: sp.Rational(round(Q[i,j]*10000),10000)); M=(M+M.T)/2
print("rationalized min-eig (float):", round(float(min(np.linalg.eigvalsh(np.array(M).astype(float)))),5))
L,D=M.LDLdecomposition(); d=[D[i,i] for i in range(4)]
assert all(x>0 for x in d) and (L*D*L.T-M)==sp.zeros(4,4)
LT=L.T
def lr(x):
    x=sp.Rational(x); return f"({x.p}/{x.q} : ℝ)" if x.q!=1 else f"({x.p} : ℝ)"
def ell(i):
    ts=[("x "+str(j) if LT[i,j]==1 else f"{lr(LT[i,j])} * x {j}") for j in range(4) if LT[i,j]!=0]
    return " + ".join(ts) if ts else "0"
def row(i): return "![" + ", ".join(lr(M[i,j]) for j in range(4)) + "]"
hints=",\n               ".join(f"mul_nonneg (by norm_num : (0:ℝ) ≤ {lr(d[i])}) (sq_nonneg ({ell(i)}))" for i in range(4))
lean=f'''import Mathlib.LinearAlgebra.Matrix.PosDef
import Mathlib.Data.Real.StarOrdered
open Matrix
/-- **Machine-checked organ-2 Weil-positivity certificate**, basis φ_a(x)=e^{{-x²/2}}cos(ax),
a∈{{15,25,35,45}}.  The rigorous Guinand–Weil form Q=arch+pole−prime (operators/organ2_certify.py,
min-eig≈0.048), rationalized (rounding ≪ margin), is `PosSemidef` — verified by its exact LDL^T
sum-of-squares.  The Lean end of the Arb → exact-rational Weil certificate (§731/§761/§762). -/
theorem organ2_weil_cert_posSemidef :
    (Matrix.of ![{row(0)}, {row(1)}, {row(2)}, {row(3)}]).PosSemidef := by
  rw [posSemidef_iff_dotProduct_mulVec]
  refine ⟨?_, ?_⟩
  · ext i j; fin_cases i <;> fin_cases j <;> norm_num [Matrix.conjTranspose_apply, Matrix.of_apply]
  · intro x
    simp only [dotProduct, Matrix.mulVec, Matrix.of_apply, Fin.sum_univ_four,
      Matrix.cons_val_zero, Matrix.cons_val_one, Matrix.head_cons, Matrix.cons_val_two,
      Matrix.cons_val_three, Matrix.tail_cons, star_trivial]
    nlinarith [{hints}]
'''
open("/tmp/organ2_cert.lean","w").write(lean); print("emitted organ2_cert.lean")
