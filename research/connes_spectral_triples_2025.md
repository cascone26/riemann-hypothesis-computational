# Connes, Consani, Moscovici -- "Zeta Spectral Triples" (arXiv:2511.22755, Nov 2025)

34 pages, 4 figures. To appear in EMS Lecture Notes in Mathematics.
MSC: 58B34, 11M06, 11M55, 33D60, 34B20

---

## 1. The Big Picture

The paper constructs explicit self-adjoint operators D_log^(lambda,N) whose eigenvalues
converge numerically to the nontrivial zeros of zeta(1/2 + is). The operators are
rank-one perturbations of the scaling operator on a logarithmic circle. The construction
uses only finitely many primes (those <= lambda^2). A rigorous proof that the eigenvalues
converge to zeta zeros as N, lambda -> infinity would prove the Riemann Hypothesis.

---

## 2. The Hilbert Space

**H = L^2([lambda^{-1}, lambda], d*u)**

where d*u = du/u is the multiplicative Haar measure (logarithmic measure).

Set L = 2 log(lambda). Under the change of variable x = log(lambda * u), the interval
[lambda^{-1}, lambda] maps to [0, L], and d*u becomes dx. The space becomes L^2([0,L], dx)
which is a circle of circumference L (periodic boundary conditions).

**Orthonormal basis:** V_n(u) = U_n(log(lambda * u)) where U_n(x) = (1/sqrt(L)) * exp(2*pi*i*n*x / L)

These are eigenfunctions of the scaling operator with eigenvalue 2*pi*n/L.

**Finite-dimensional truncation:** E_N = span{V_n : |n| <= N}, dimension 2N+1.

---

## 3. The Base Operator: Scaling Operator

**D_log^(lambda) = -i * u * d/du = -i * d/d(log u)**

with periodic boundary conditions on [lambda^{-1}, lambda].

Eigenvalues: s_n = 2*pi*n / L = pi*n / log(lambda), for n in Z.
Eigenfunctions: V_n as above.

This is the "unperturbed" Dirac operator of the spectral triple.

---

## 4. The Weil Quadratic Form Q_W

The Weil quadratic form Q_W_lambda(f,f) encodes arithmetic information. It has three
contributions:

### 4.1 Archimedean contribution (from the real place)

W_R(V_n, V_m) involves integrals of the form:

    integral_0^L g(x) * rho(x) dx

where rho(x) = exp(x/2) / (exp(x) - exp(-x)) and g involves sin, cos, x*sin, x*cos
terms with frequencies 2*pi*n/L and 2*pi*m/L.

The explicit formulas (Proposition 4.2) use:
- Hypergeometric functions 2F1
- Hurwitz-Lerch transcendent Phi(e^{-2L}, 2, ...)
- Digamma function psi(z) and its derivative psi^(1)(z)

The fast convergence factor is e^{-2L} which for lambda = sqrt(14) gives
e^{-2L} = e^{-4 log(sqrt(14))} = 1/14^2 ~ 0.005, ensuring rapid series convergence.

### 4.2 Non-archimedean contribution (from primes p <= lambda^2)

    sum_{p primes} W_p(V_n, V_m) = sum_{1 < k <= exp(L)} Lambda(k) * k^{-1/2} * q(U_n, U_m)(log k)

where Lambda(k) is the von Mangoldt function and:

**For n != m:**
    q(U_n, U_m)(y) = [sin(2*pi*m*y/L) - sin(2*pi*n*y/L)] / [pi*(n-m)]

**For n = m:**
    q(U_n, U_n)(y) = 2*(1 - |y|/L) * cos(2*pi*n*y/L)

Since L = 2*log(lambda), the sum runs over k <= lambda^2, meaning only primes
p <= x = lambda^2 contribute. This is the Euler product truncation.

### 4.3 Boundary evaluation terms

    W_{0,2}(V_n, V_m) = 32*L * sinh^2(L/4) * (L^2 - 16*pi^2*m*n) /
                         [(L^2 + 16*pi^2*m^2) * (L^2 + 16*pi^2*n^2)]

---

## 5. The Toeplitz-like Matrix Structure

The truncated Weil form Q_W_lambda^N restricted to E_N gives a (2N+1) x (2N+1) matrix
tau with special structure (Lemma 5.1):

    tau_{i,i} = a_i           (diagonal)
    tau_{i,j} = (b_i - b_j) / (i - j)   for i != j

with symmetries: a_{-j} = a_j (even) and b_{-j} = -b_j (odd).

The off-diagonal structure (b_i - b_j)/(i-j) is the key: this is a generalized Toeplitz
form that connects to the Caratheodory-Fejer theorem.

**The b_j values:**
    b_j = -(1/pi) * integral_0^L sin(2*pi*j*y/L) * D(y) dy

where D is the distribution log_*(Psi^#) encoding the arithmetic data.

**Commutator relation:**
    [D_log, tau] = |beta><eta| - |eta><beta|

where beta = sum_j b_j V_j and eta = sum_j V_j. This rank-two commutator structure
is what makes the rank-one perturbation work.

---

## 6. The Rank-One Perturbation

### 6.1 Finding the minimal eigenvector xi

1. Compute the (2N+1) x (2N+1) matrix tau (the truncated Weil form).
2. Find its smallest eigenvalue epsilon_N and corresponding eigenvector xi.
3. Require xi to be even: xi_{-j} = xi_j.
4. Normalize so that delta_N(xi) = 1.

### 6.2 The Dirichlet kernel

    delta_N = (1/sqrt(L)) * sum_{n=-N}^{N} V_n

This approximates a delta function at the boundary (point evaluation at u = lambda).
As N -> infinity, <delta_N | f> -> f(lambda).

### 6.3 The perturbed operator

    D_log^{lambda,N} = D_log^{lambda} - |D_log^{lambda} xi><delta_N|

In other words: for any f in E_N,

    D_log^{lambda,N}(f) = D_log^{lambda}(f) - <delta_N | f> * D_log^{lambda}(xi)

This is a rank-one perturbation. It ensures D_log^{lambda,N}(xi) = 0 (since <delta_N|xi> = 1).

### 6.4 Self-adjointness

D_log^{lambda,N} is NOT self-adjoint with respect to the standard inner product.
It IS self-adjoint with respect to the modified inner product on E_N' = E_N / C*xi:

    <f | g>_T = <(tau - epsilon_N * I) f | g>

This is where the Caratheodory-Fejer extension (from Connes-van Suijlekom, CMP 2025,
arXiv:2511.23257) is essential: the generalized Toeplitz structure of tau guarantees
that this modified inner product makes the operator self-adjoint and the spectrum real.

---

## 7. Eigenvalues = Zeros of an Entire Function

**Theorem 5.10:**

The Fourier transform of xi,

    xi_hat(z) = (2/sqrt(L)) * sin(z*L/2) * sum_{j=-N}^{N} xi_j / (z - 2*pi*j/L)

is an entire function. ALL its zeros are real, and they coincide exactly with the
spectrum of D_log^{lambda,N}.

**Regularized determinant:**

    det_reg(D_log^{lambda,N} - z) = -i * lambda^{-iz} * xi_hat(z)

So finding eigenvalues = finding real zeros of xi_hat(z).

This is the function to implement numerically: compute xi_j (the eigenvector components),
then evaluate xi_hat(z) and find its zeros.

---

## 8. Connection to Zeta Zeros

The zeros of xi_hat(z) approximate the imaginary parts of the nontrivial zeros of
zeta(1/2 + is). As N, lambda -> infinity, the regularized determinants (suitably
normalized) converge toward the Riemann Xi function:

    Xi(s) = (1/2) * s * (s-1) * pi^{-s/2} * Gamma(s/2) * zeta(s)

Since Xi(1/2 + it) is real for real t and its zeros are exactly the nontrivial zeta zeros
(on RH), convergence of xi_hat toward Xi would prove RH.

---

## 9. Numerical Results

### 9.1 Parameters used

- N = 120 (so 241 x 241 matrices)
- lambda = sqrt(12), sqrt(13), sqrt(14)
- This means x = lambda^2 = 12, 13, 14 (primes used: {2, 3, 5, 7, 11, 13})
- Computations done in 200-digit precision arithmetic

### 9.2 Error table (absolute difference: |eigenvalue - zeta zero|)

| Zero # | lambda=sqrt(12) | lambda=sqrt(13) | lambda=sqrt(14) |
|--------|-----------------|-----------------|-----------------|
| 1      | ~3.4e-50        | ~mid            | ~1.1e-60        |
| 5      | ~7.8e-41        | ~mid            | ~3.8e-51        |
| 25     | ~1.9e-15        | ~mid            | ~3.9e-24        |
| 50     | ~9.0e-2         | ~mid            | ~4.8e-6         |

The accuracy is extraordinary for low zeros and degrades for higher zeros (as expected
since higher zeros need larger lambda to resolve).

### 9.3 Statistical significance

The probability of this agreement occurring by chance is estimated at ~10^{-1235}.

### 9.4 Case study: lambda = 3, N = 120

Figure 1 in the paper shows differences between the first 20 zeta zeros and eigenvalues,
with errors in the range 10^{-20} to 10^{-10}.

---

## 10. Recipe for Numerical Implementation

### Step 1: Choose parameters
- Pick lambda (e.g., sqrt(14), so L = 2*log(sqrt(14)) = log(14) ~ 2.6391)
- Pick N (e.g., 120)
- Identify primes: all p <= lambda^2 = 14, so {2, 3, 5, 7, 11, 13}

### Step 2: Build the (2N+1) x (2N+1) Weil matrix tau
For each pair (n, m) with -N <= n, m <= N, compute:
- tau_{n,m} = W_R(V_n, V_m) + W_{0,2}(V_n, V_m) - sum_p W_p(V_n, V_m)

The archimedean part W_R requires evaluating hypergeometric/digamma functions.
The prime part uses von Mangoldt weights and the q(U_n,U_m) formulas above.
The boundary part W_{0,2} has a closed form.

### Step 3: Find the minimal eigenvector
- Diagonalize tau (symmetric matrix).
- Take the smallest eigenvalue epsilon_N and its eigenvector xi.
- Ensure xi is even (xi_{-j} = xi_j). If the eigenspace is > 1-dimensional, pick the even one.
- Normalize: sum_{j=-N}^{N} xi_j / sqrt(L) = 1, i.e., delta_N(xi) = 1.

### Step 4: Compute xi_hat(z) and find its zeros
    xi_hat(z) = (2/sqrt(L)) * sin(z*L/2) * sum_{j=-N}^{N} xi_j / (z - 2*pi*j/L)

Find real zeros of this function (e.g., Newton's method starting near known zeta zeros,
or a root-finding sweep).

### Step 5: Compare to known zeta zeros
The real zeros of xi_hat should approximate: 14.134725..., 21.022040..., 25.010858..., etc.

### Key numerical notes
- Use high-precision arithmetic (200+ digits) since the matrix entries involve cancellations
- The factor e^{-2L} ensures fast convergence of the special function series
- For L ~ 2.6 (lambda = sqrt(14)), e^{-2L} ~ 0.005 which is already small
- Larger lambda gives more primes and better accuracy for higher zeros
- Larger N gives more basis functions and better resolution

---

## 11. Follow-up Paper

**"Spectral Analysis of the D_log^{(lambda,N)} Operators"** by Dominik Sliwinski
(arXiv:2601.12133, January 2026, 7 pages)

Key finding: The mean absolute error satisfies a lower bound

    epsilon(lambda, N) >= 1 / (4 * ln(lambda))

This inverse-logarithmic bound connects spectral deviation to prime distribution.
Numerical tests with kappa from 50 to 7500 against first 1000 zeta zeros confirm
consistent inverse-logarithmic convergence behavior.

Conjecture: lim_{k->inf} E(k) * ln(k) exists, which would imply RH.

---

## 12. Predecessor Papers

1. **"Spectral triples and zeta-cycles"** -- Connes & Consani, Enseign. Math. 69 (2023),
   93-148. arXiv:2106.01715. Introduces the framework of spectral triples on logarithmic
   circles and the connection to zeta-cycles.

2. **"Quadratic Forms, Real Zeros and Echoes of the Spectral Action"** -- Connes &
   van Suijlekom, Commun. Math. Phys. 406 (2025). arXiv:2511.23257. Extends
   Caratheodory-Fejer theorem: if a Hermitian positive semidefinite matrix has the
   generalized Toeplitz structure tau_{i,j} = (b_i - b_j)/(i-j), and xi is in its kernel,
   then the associated polynomial has all zeros on the unit circle. This is the theoretical
   engine guaranteeing self-adjointness and real spectrum.

3. **"Zeta zeros and prolate wave operators"** -- Connes & Moscovici,
   arXiv:2310.18423 (2023). Earlier approach using prolate spheroidal wave operators.

4. **"The Riemann Hypothesis: Past, Present and a Letter Through Time"** -- Connes,
   arXiv:2602.04022 (Feb 2026), 42 pages. Survey covering 165 years of approaches,
   includes a "Letter to Riemann" demonstrating the optimization procedure.

---

## 13. Community Discussion

As of March 2026, I found NO MathOverflow threads, blog posts, or accessible reviews
specifically discussing the D_log^{lambda,N} construction from this paper. The paper is
very recent (Nov 2025) and highly technical. The Sliwinski follow-up (Jan 2026) is the
only independent analysis I found.

Connes gave a talk at RISM (Riemann International School of Mathematics) on June 4, 2025,
titled "Topos and noncommutative geometry: two views on space and numbers," which likely
covered related material. A 2025 AOFA Best Paper Award webinar by Connes was hosted by
Springer Nature but content details are not publicly available.

---

## 14. Summary for Implementation

To reproduce the numerical results, you need to:

1. **Build the Weil matrix** (hardest part): requires implementing hypergeometric functions,
   Hurwitz-Lerch transcendents, digamma functions, and the von Mangoldt-weighted prime sums.
   All at high precision (200+ digits). Use mpmath or Arb.

2. **Eigendecompose** the (2N+1) x (2N+1) symmetric matrix at high precision.

3. **Evaluate xi_hat(z)** -- this is a rational function times sin(zL/2), straightforward
   once you have the eigenvector components xi_j.

4. **Root-find** xi_hat(z) = 0 on the real line.

The non-archimedean (prime) part of the matrix is actually the simplest -- it is a finite
sum over prime powers k <= lambda^2 with explicit trigonometric weights. The archimedean
part involving special functions is the implementation bottleneck.

---

Sources:
- https://arxiv.org/abs/2511.22755
- https://arxiv.org/html/2511.22755v1
- https://arxiv.org/abs/2511.23257
- https://arxiv.org/abs/2601.12133
- https://arxiv.org/abs/2602.04022
- https://arxiv.org/abs/2106.01715
