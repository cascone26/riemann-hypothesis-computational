"""Verify the Connes-Moscovici operator FORM via its semiclassical eigenvalue density.
Symbol H(x,p)=(x²-λ²)p²-(2π)²x² on |x|>λ; density of states (μ=-E>0):
  ρ(μ) = (1/π)∫_λ^∞ dx / sqrt((μ+(2π)²x²)(x²-λ²)).
Zeros have density in μ: ρ_zeros(μ) = log(√μ/2π)/(4π√μ)  [from N(γ)~(γ/2π)(log(γ/2π)-1), μ=γ²].
If ρ(μ) ~ ρ_zeros(μ), the operator FORM reproduces the zero density (validates the approach)."""
import numpy as np
from scipy import integrate
def rho_op(mu, lam):
    f=lambda x: 1.0/np.sqrt((mu+(2*np.pi)**2*x**2)*(x**2-lam**2))
    val,_=integrate.quad(f, lam, np.inf, limit=200)
    return val/np.pi
def rho_zeros(mu):
    g=np.sqrt(mu)
    return np.log(g/(2*np.pi))/(4*np.pi*g)
print(f"{'μ':>10} {'√μ=γ':>8} {'ρ_op(μ)':>13} {'ρ_zeros(μ)':>13} {'ratio':>8}")
for lam in [1.0, np.sqrt(2)]:
    print(f"--- λ={lam:.4f} ---")
    for mu in [200,500,1000,4000,1e4,1e5,1e6,1e8]:
        ro=rho_op(mu,lam); rz=rho_zeros(mu)
        print(f"{mu:>10.0f} {np.sqrt(mu):>8.1f} {ro:>13.4e} {rz:>13.4e} {ro/rz:>8.4f}")
