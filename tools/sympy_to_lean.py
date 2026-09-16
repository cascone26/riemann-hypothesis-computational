import re
import sympy as sp

def sin_lean(i):
    return f"(((Real.sqrt 2 : ℝ) : ℂ) * Complex.sin (L * z{i} / 2))"

def csin_lean(i):
    return f"(starRingEnd ℂ (((Real.sqrt 2 : ℝ) : ℂ) * Complex.sin (L * z{i} / 2)))"

def sympy_to_lean(expr):
    s = str(sp.expand(expr))
    s = s.replace('**', '^')
    for i in range(1, 5):
        s = re.sub(rf'\bw{i}c\b', csin_lean(i), s)
        s = re.sub(rf'\bw{i}\b', sin_lean(i), s)
        s = re.sub(rf'\bz{i}c\b', f'(starRingEnd ℂ z{i})', s)
    return s
