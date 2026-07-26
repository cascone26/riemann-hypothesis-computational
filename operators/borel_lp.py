#!/usr/bin/env python3
"""§323: is B(s)=Σ_k γ_k (−s²)^k/k! UNCONDITIONALLY in Laguerre-Pólya? Structural test.
Write c_k = γ_k/k! (Taylor coeffs of B in x=−s²: B=Σ c_k x^k). LP for an even entire fn of genus 0 with
positive-x zeros needs the Turán/Hankel positivity of {c_k}: LAGUERRE inequalities c_k² − c_{k-1}c_{k+1} ≥ 0
(necessary for LP) AND the stronger log-concavity of k! c_k etc. Key structural fact to test:
  c_k = γ_k/k! = M_k/[(2k)! k!],  M_k=∫Φ(u)u^{2k}du.  Since Φ≥0, {M_k} is a STIELTJES moment sequence
  ⇒ its Hankel matrices are PSD (unconditional, from Φ≥0 alone!). Question: does dividing by (2k)!k! (a
  double-factorial / Bessel weight) PRESERVE the total positivity that forces real zeros?
Tests (all unconditional, from Φ≥0):
 (1) Hankel PSD of {M_k} (Stieltjes) — must hold (sanity).
 (2) Laguerre ineq for c_k (LP-necessary).
 (3) The DECISIVE one: is {c_k} itself a moment/Hankel-PSD sequence? If {γ_k/k!} is Hankel-PSD it is Stieltjes
     ⇒ B(x)=Σc_k x^k has only NEGATIVE real zeros in x ⇒ B(s) only REAL zeros in s — UNCONDITIONALLY.
 (4) double-Hankel / Hankel of {k! c_k} = {γ_k}.
"""
import mpmath as mp, pickle
mp.mp.dps=60
gam=[mp.mpf(g) for g in pickle.load(open("operators/_gam130.pkl","rb"))]
K=40
c=[gam[k]/mp.factorial(k) for k in range(K)]       # Taylor coeffs of B in x
M=[gam[k]*mp.factorial(2*k) for k in range(K)]      # moments M_k (=γ_k (2k)!)
def hankel_min_eig(seq,size):
    H=mp.matrix(size,size)
    for i in range(size):
        for j in range(size): H[i,j]=seq[i+j]
    ev,_=mp.eigsy((H+H.T)/2)
    return ev[0]
def hankel_shifted_min_eig(seq,size):  # shifted Hankel [a_{i+j+1}] (Stieltjes second condition)
    H=mp.matrix(size,size)
    for i in range(size):
        for j in range(size): H[i,j]=seq[i+j+1]
    ev,_=mp.eigsy((H+H.T)/2)
    return ev[0]
print("(1) Stieltjes Hankel of M_k (must be PSD, Φ≥0):",flush=True)
for s in (5,10,15):
    print(f"   size {s}: min eig H[M]={mp.nstr(hankel_min_eig(M,s),4)}  H1[M]={mp.nstr(hankel_shifted_min_eig(M,s),4)}",flush=True)
print("(2) Laguerre inequalities c_k²−c_{k-1}c_{k+1} ≥ 0 (LP-necessary):",flush=True)
lag=[c[k]**2-c[k-1]*c[k+1] for k in range(1,20)]
print("   "+", ".join(mp.nstr(x,3) for x in lag[:12]),flush=True)
print(f"   all ≥0 (k=1..19)? {all(x>=0 for x in lag)}",flush=True)
print("(3) DECISIVE: Hankel PSD of c_k={γ_k/k!} (PSD ⇒ Stieltjes ⇒ B has only real zeros, UNCONDITIONAL):",flush=True)
for s in (5,8,12,16):
    e0=hankel_min_eig(c,s); e1=hankel_shifted_min_eig(c,s)
    print(f"   size {s}: min eig H[c]={mp.nstr(e0,4)}  H1[c]={mp.nstr(e1,4)}  {'PSD' if e0>-mp.mpf('1e-40') and e1>-mp.mpf('1e-40') else 'NOT PSD'}",flush=True)
print("(4) Hankel PSD of γ_k (=k! c_k):",flush=True)
for s in (5,10,15):
    print(f"   size {s}: min eig H[γ]={mp.nstr(hankel_min_eig(gam,s),4)}  H1[γ]={mp.nstr(hankel_shifted_min_eig(gam,s),4)}",flush=True)
print("READ: (3) PSD at all sizes ⇒ {γ_k/k!} Stieltjes ⇒ B∈LP unconditionally (the free/RH separation).",flush=True)
print("DONE",flush=True)
