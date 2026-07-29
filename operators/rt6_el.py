#!/usr/bin/env python3
"""§394 (angle 7): Euler-Lagrange criticality — is the margin-gradient orthogonal to the pair-preserving
tangent at the TRUE Lambda? Gradient (Hellmann-Feynman): g(n) = -(1/sqrt n)<v,G_n v>, v = bottom eig of
B_full. Constraint gradients: const 1 (mean) + f_h(n)=Lambda(n+h)+Lambda(n-h) (pair lags h). Regress g on
span{1, f_h}: R^2 = criticality degree; coefficients lambda(h) = the LAGRANGE MULTIPLIER — the d