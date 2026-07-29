#!/usr/bin/env python3
"""§434 (C-C route, correct grounding): the classical Arakelov/van der Geer-Schoof Riemann-Roch on Spec ℤ.
h⁰(D)=log θ(e^{-2t}) [t=deg D, θ(y)=Σ_{n∈ℤ}e^{-πn²y} = effectivity/theta], h¹(D)=h⁰(−D). CLAIM (RR): 
h⁰(D)−h¹(D)=deg(D), and this IS the theta functional equation θ(1/y)=√y·θ(y) = the ζ functional equation.
Verify numerically. Then note: the self-square Spec ℤ×_{F₁}Spec ℤ needs the 2D analog + diagonal self-
intersection (the OPEN Hodge piece). This grounds the C-C route in the (correct, classical) functional equation."""
import numpy as np
def theta(y):  # Σ_{n∈ℤ} e^{-π n² y}, converges fast for y>0
    n=np.arange(-2000,2001); return np.sum(np.exp(-np.pi*n**2*y))
def h0(t): return np.log(theta(np.exp(-2*t)))    # h⁰(D), deg D = t
print("§434 Arakelov RR on Spec ℤ:  h⁰(D) − h⁰(−D) =?= deg(D) = t   (= theta/ζ functional equation)")
print(f"{'t=degD':>8} {'h0(D)':>10} {'h1=h0(-D)':>10} {'h0-h1':>10} {'deg':>6} {'match':>7}")
for t in [-1.5,-1.0,-0.5,0.0,0.5,1.0,1.5,2.0,3.0]:
    a=h0(t); b=h0(-t); diff=a-b
    print(f"{t:8.2f} {a:10.5f} {b:10.5f} {diff:10.5f} {t:6.2f} {'YES' if abs(diff-t)<1e-9 else 'no':>7}")
# the functional equation origin: θ(1/y)=√y θ(y)
y=0.3
print(f"\ntheta functional eq check: θ(1/y)={theta(1/y):.6f}, √y·θ(y)={np.sqrt(y)*theta(y):.6f} (equal => RR holds)")
print("""
reading: h⁰(D)−h¹(D)=deg(D) holds EXACTLY (=theta functional equation θ(1/y)=√y θ(y) = ζ's functional equation
ξ(s)=ξ(1-s)). So RR on Spec ℤ IS the functional equation — the SYMMETRY half of RH (banked, classical). What RR
does NOT give is the POSITIVITY/Hodge-index half (zeros ON the line) — that needs the SELF-SQUARE's diagonal
self-intersection (C-C open step). This correctly grounds the C-C route: functional equation = 1D RR (have it);
RH positivity = 2D Hodge index on the self-square (open). Matches §326 (explicit formula = Lefschetz) + §422.
""")
