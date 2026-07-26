#!/usr/bin/env python3
"""§313A: third-order saddle law (degree-d control). From the saddle eq 2*pi*e^{2u}=2m/u+9/2:
  u' = 2u/D, D = 2m(2u+1)+9u^2;  u'' = [2u'D - 2u*D']/D^2, D' = 2(2u+1)+u'(4m+18u).
F(m)=log M_m, F'=2 log u*, F''=2u'/u, F'''=2(u''u-u'^2)/u^2.
log gamma_m = F(m) - log Gamma(2m+1):  r_m (2nd) = F'' - 4*psi'(2m+1),  c_m (3rd) = F''' - 8*psi''(2m+1).
Verify BOTH against exact central differences of log gamma_m (gamma_m to m=24, dps 45).
c_m is the object whose smallness vs the Hermite hyperbolicity margin makes GORZ effective. Charisma: HANDLE if matches."""
import mpmath as mp
mp.mp.dps=45
def Phi(u):
    s=mp.mpf(0); e2u=mp.e**(2*u)
    for n in range(1,90):
        term=(2*mp.pi**2*n**4*mp.e**(mp.mpf(9)*u/2)-3*mp.pi*n*n*mp.e**(mp.mpf(5)*u/2))*mp.e**(-mp.pi*n*n*e2u)
        s+=term
        if abs(term)<mp.mpf(10)**(-50) and n>3: break
    return s
Mmax=25
gam=[]
for m in range(Mmax):
    M=mp.quad(lambda u: Phi(u)*u**(2*m),[0,0.4,0.8,1.2,2,4])
    gam.append(M/mp.factorial(2*m))
lg=[mp.log(g) for g in gam]
def usaddle(m):
    return mp.findroot(lambda u: 2*m/u+mp.mpf(9)/2-2*mp.pi*mp.e**(2*u), mp.mpf(0.5)+mp.log(m+1)/2)
def saddle_derivs(m):
    u=usaddle(m)
    D=2*m*(2*u+1)+9*u*u
    up=2*u/D
    Dp=2*(2*u+1)+up*(4*m+18*u)
    upp=(2*up*D-2*u*Dp)/D**2
    F2=2*up/u
    F3=2*(upp*u-up*up)/(u*u)
    r_pred=F2-4*mp.polygamma(1,2*m+1)
    c_pred=F3-8*mp.polygamma(2,2*m+1)
    return u,r_pred,c_pred
print("m :  r_m exact(2nd cdiff)   r_m saddle    ratio  |  c_m exact(3rd cdiff)    c_m saddle    ratio",flush=True)
for m in range(3,Mmax-2):
    r_ex=lg[m+1]-2*lg[m]+lg[m-1]
    c_ex=(lg[m+2]-2*lg[m+1]+2*lg[m-1]-lg[m-2])/2
    u,r_pr,c_pr=saddle_derivs(m)
    print(f"{m:2d}: {mp.nstr(r_ex,7):>13} {mp.nstr(r_pr,7):>13}  {mp.nstr(r_pr/r_ex,5):>7} | {mp.nstr(c_ex,6):>12} {mp.nstr(c_pr,6):>12}  {mp.nstr(c_pr/c_ex,5):>7}",flush=True)
print("",flush=True)
print("asymptotics: r_m ~ -2/m ; c_m ~ +2/m^2 (from -8psi''(2m+1)~2/m^2 dominating F'''):",flush=True)
for m in [10,15,20]:
    u,r_pr,c_pr=saddle_derivs(m)
    print(f"  m={m}: m*|r|={mp.nstr(m*abs(r_pr),5)}  m^2*c={mp.nstr(m*m*c_pr,5)}",flush=True)
print("DONE",flush=True)
