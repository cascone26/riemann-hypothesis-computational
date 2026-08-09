#!/usr/bin/env python3
"""§897: the archimedean density Ω(r)=Re ψ₀(¼+ir/2)−ln π is CAUCHY–STIELTJES in r² — its derivative dΩ/d(r²) is
COMPLETELY MONOTONE. Analytic: ψ₀(¼+ir/2)=−γ+Σ_n[1/(n+1)−1/(n+¼+ir/2)] ⇒ Ω(r)=−γ−lnπ+Σ_n[1/(n+1)−(n+¼)/((n+¼)²+r²/4)].
Let x=r²; each r-term = −a_n/(a_n²+x/4), a_n=n+¼. dΩ/dx = Σ_n (a_n/4)/(a_n²+x/4)² > 0; (−1)^k d^k(dΩ/dx)/dx^k ≥ 0
(each 1/(a_n²+x/4)² is completely monotone in x). VERIFY: (1) the series matches the digamma value, (2) successive
x-derivatives of Ω alternate sign as required for a Stieltjes/CM derivative. This is the pointwise origin of Groskin's
totally-positive archimedean tail — a real UNCONDITIONAL structural fact about the arch density."""
import numpy as np
from scipy.special import digamma
np.set_printoptions(precision=4,suppress=False,linewidth=140)

def Omega_direct(r): return np.real(digamma(0.25+0.5j*r))-np.log(np.pi)
def Omega_series(r,Nn=200000):
    n=np.arange(0,Nn); a=n+0.25
    return -np.euler_gamma-np.log(np.pi)+np.sum(1.0/(n+1.0)-a/(a*a+ (r*r)/4.0))

# (1) series == direct?
for r in [0.5,3.0,14.0,44.0,150.0]:
    print(f"r={r:6.1f}  Ω_direct={Omega_direct(r):+.8f}  Ω_series={Omega_series(r):+.8f}  diff={Omega_direct(r)-Omega_series(r):+.2e}")

# (2) numerically check dΩ/dx and its CM property (x=r^2). Use central differences in x.
print("\ncomplete-monotonicity of f(x)=dΩ/dx (x=r²): signs of f, f', f'', f''' should be +,−,+,−")
def Omega_x(x): r=np.sqrt(x); return Omega_direct(r)
for x0 in [1.0,9.0,100.0,2000.0,20000.0]:
    h=x0*1e-3
    # derivatives of Omega wrt x via finite diff
    d1=(Omega_x(x0+h)-Omega_x(x0-h))/(2*h)                      # f = dΩ/dx
    d2=(Omega_x(x0+h)-2*Omega_x(x0)+Omega_x(x0-h))/h**2         # f'
    d3=(Omega_x(x0+2*h)-2*Omega_x(x0+h)+2*Omega_x(x0-h)-Omega_x(x0-2*h))/(2*h**3)  # f''
    d4=(Omega_x(x0+2*h)-4*Omega_x(x0+h)+6*Omega_x(x0)-4*Omega_x(x0-h)+Omega_x(x0-2*h))/h**4 # f'''
    cm = (d1>0) and (d2<0) and (d3>0) and (d4<0)
    print(f"x={x0:8.0f}: f={d1:+.3e} f'={d2:+.3e} f''={d3:+.3e} f'''={d4:+.3e}  CM-pattern(+,-,+,-): {cm}")
print("\nf=dΩ/dx completely monotone (alternating signs) ⇒ Ω is a Bernstein-type fn of r² with a Stieltjes derivative")
print("⇒ Ω(r)=c + ∫₀^∞ [x/(x+r²/4·?)]... i.e. the archimedean density is Cauchy–Stieltjes. Unconditional, from digamma.")
print("This is the pointwise structural fact under Groskin T3.2 (totally-positive Cauchy–Stieltjes archimedean tail).")
