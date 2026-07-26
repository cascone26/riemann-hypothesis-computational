#!/usr/bin/env python3
"""STORM #2: Castelnuovo-Severi arithmetic-surface construction on the Weil space.
Realize Δ, Γ (Frobenius/scaling), fiber as divisor classes; compute intersections via Weil form;
test Hodge-index (signature) and Castelnuovo-Severi inequality."""
import mpmath as mp
mp.mp.dps=22
def Omega(r): return mp.re(mp.digamma(mp.mpf(1)/4+1j*r/2))-mp.log(mp.pi)
w=mp.mpf(3)
def F(c,r): return mp.e**(-(r-c)**2/(2*w**2))+mp.e**(-(r+c)**2/(2*w**2))
def Fc(c,z): return mp.e**(-(z-c)**2/(2*w**2))+mp.e**(-(z+c)**2/(2*w**2))
def bpts(a,b):
    p=[mp.mpf(0)]
    for c in (a,b):
        for x in (c-4*w,c,c+4*w):
            if x>0:p.append(x)
    return sorted(set(p))+[mp.inf]
Lam={}
for pp in range(2,50):
    if all(pp%d for d in range(2,int(pp**.5)+1)):
        q=pp
        while q<50: Lam[q]=mp.log(pp);q*=pp
centers=[mp.mpf(14+3*i) for i in range(13)]  # [14,50], M=13
M=len(centers)
print(f"STORM#2 Castelnuovo-Severi: M={M}, band [14,50]",flush=True)

# Build Weil form B
B=mp.zeros(M)
for j in range(M):
    for k in range(M):
        pts=bpts(centers[j],centers[k])
        arch=(1/(2*mp.pi))*2*mp.quad(lambda r:F(centers[j],r)*F(centers[k],r)*Omega(r),pts)
        pole=mp.re(Fc(centers[j],1j/2)*Fc(centers[k],1j/2))*2
        pr=-2*sum(Lam[n]/mp.sqrt(n)*(1/(2*mp.pi))*2*mp.quad(lambda r:F(centers[j],r)*F(centers[k],r)*mp.cos(mp.log(n)*r),pts) for n in sorted(Lam))
        B[j,k]=pole+arch+pr
    print(f"  row {j}",flush=True)

# Divisor classes (vectors in the Gaussian basis):
# FIBER: constant function => (1,1,...,1)
f_fiber=mp.matrix([mp.mpf(1) for _ in range(M)])
f_fiber=f_fiber/mp.sqrt(mp.mpf(M))

# SCALING (Frobenius/height): weighted by center => (c_1, c_2, ..., c_M)
f_scal=mp.matrix([centers[j] for j in range(M)])
f_scal=f_scal/mp.sqrt(sum(centers[j]**2 for j in range(M)))

# DIAGONAL: the min-eigenvector (near-saturating function from storm#3, which encodes the zeros)
ev,V=mp.eigsy(B)
f_diag=mp.matrix([V[i,0] for i in range(M)])  # min-eig direction

# Intersection numbers via Weil pairing
def inner(f,g): return float((f.T * B * g)[0,0])  # f^T B g

II={}
for (name_a,f_a) in [("D",f_diag),("G",f_scal),("f",f_fiber)]:
    for (name_b,f_b) in [("D",f_diag),("G",f_scal),("f",f_fiber)]:
        II[name_a+"."+name_b]=inner(f_a,f_b)

print("\n--- Intersection numbers (classes: D=diagonal, G=Frobenius/scaling, f=fiber) ---",flush=True)
print(f"  D·D = {mp.nstr(II['D.D'],5)}")
print(f"  G·G = {mp.nstr(II['G.G'],5)}")
print(f"  f·f = {mp.nstr(II['f.f'],5)}")
print(f"  D·G = {mp.nstr(II['D.G'],5)}")
print(f"  D·f = {mp.nstr(II['D.f'],5)}")
print(f"  G·f = {mp.nstr(II['G.f'],5)}")

# Signature of the 3x3 intersection matrix
Mx=mp.matrix([[II['D.D'],II['D.G'],II['D.f']],
              [II['G.D'],II['G.G'],II['G.f']],
              [II['f.D'],II['f.G'],II['f.f']]])
evals=mp.eigsy(Mx,eigvals_only=True)
pos=sum(1 for e in evals if e>mp.mpf(10)**(-10))
neg=sum(1 for e in evals if e< -mp.mpf(10)**(-10))
print(f"\n  Intersection matrix eigenvalues: {[mp.nstr(e,4) for e in evals]}")
print(f"  Signature: {pos} positive, {neg} negative => {'HYPERBOLIC (1 pos, rest neg)' if pos==1 and neg==2 else 'NOT hyperbolic'}")

# Castelnuovo-Severi: primitive class D = Frobenius - (Frobenius·fiber)/(fiber·fiber) fiber - (Frobenius·height)/(height·height) height
# Actually, in the classical form: D = Gamma_Frobenius - (Gamma·v)h - (Gamma·h)v, where v,h are fiber directions.
# Here, simplify: let D_prim = Gamma - (projection onto fiber) = Gamma - (Gamma·f/f·f) f
coeff_fib=II['G.f']/II['f.f']
f_prim=f_scal - coeff_fib*f_fiber
D_prim_D_prim=inner(f_prim,f_prim)
print(f"\n--- Castelnuovo-Severi inequality ---",flush=True)
print(f"  Primitive class D_prim = Γ - {mp.nstr(coeff_fib,4)} f")
print(f"  D_prim · D_prim = {mp.nstr(D_prim_D_prim,5)}  => {'≤0 (Castelnuovo holds)' if D_prim_D_prim <= mp.mpf(10)**(-10) else '>0 (Castelnuovo FAILS)'}")

# rank=N(band) connection
rank_B=sum(1 for e in ev if e>mp.mpf(10)**(-11))
print(f"\n--- Rank vs #zeros ---",flush=True)
print(f"  Weil-form rank (M={M}) = {rank_B} (eigenvalues: min={mp.nstr(ev[0],4)}, max={mp.nstr(ev[-1],4)})")
print(f"  In-band zero count = 10 (from storm#3)")
print(f"  rank ≈ M - #zeros = {M} - 10 = 3  => actual rank {rank_B} {'consistent' if rank_B<=4 else 'INCONSISTENT'}")

print("\nDONE",flush=True)
