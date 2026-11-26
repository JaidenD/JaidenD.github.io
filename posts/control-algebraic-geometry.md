---
title: Linear Control Theory as Algebraic Geometry
date: 2026-11-26
summary: Rephrasing the main results of linear control theory in the language of algebraic geometry, following Hermann and Martin's observation that a transfer function is a rational curve in a Grassmannian.
toc: true
draft: true
---

In 1977 Robert Hermann and Clyde Martin pointed out that the transfer function of a linear system with $m$ inputs and $p$ outputs is the same thing as a rational curve in the Grassmannian of $m$-planes in $\mathbb C^{m+p}$. Once you see this, a lot of linear control theory turns into classical algebraic geometry: controllability becomes surjectivity of a map of modules, the McMillan degree becomes the degree of a curve, and pole placement by output feedback becomes Schubert calculus.

In this post we'll go through the main results of linear control theory one at a time and rephrase each of them in algebro-geometric terms. Nothing here needs more than the definitions of a module, a Grassmannian and a vector bundle on $\mathbb P^1$.

## Linear systems in one page

A linear time-invariant system is
$$
\dot x = Ax + Bu,
\qquad
y = Cx,
$$
with state $x \in k^n$, input $u \in k^m$, output $y \in k^p$, and $k = \mathbb R$ or $\mathbb C$. Taking Laplace transforms, the output is $y(s) = G(s)u(s)$ with transfer function
$$
G(s) = C(sI - A)^{-1}B,
$$
a $p \times m$ matrix of strictly proper rational functions. The classical results we want to rephrase are:

1. **Controllability.** $(A,B)$ is controllable (every state can be reached from $0$) if and only if $\operatorname{rank}\,[B \;\, AB \;\, \cdots \;\, A^{n-1}B] = n$, if and only if $\operatorname{rank}\,[\lambda I - A \;\; B] = n$ for every $\lambda \in \mathbb C$ (the Popov–Belevitch–Hautus test). Observability of $(C,A)$ is the dual statement.
2. **Realization.** Every strictly proper rational $G$ is $C(sI-A)^{-1}B$ for some $(A,B,C)$. A realization has the smallest possible $n$ exactly when it is controllable and observable, and any two minimal realizations differ by a change of state coordinates $(A,B,C) \mapsto (TAT^{-1}, TB, CT^{-1})$ (Kalman). That smallest $n$ is the *McMillan degree* $\delta(G)$.
3. **Coprime factorizations.** $G = N(s)D(s)^{-1}$ with $N$, $D$ polynomial matrices which are right coprime: there are polynomial $X$, $Y$ with $XN + YD = I$ (a Bezout identity). Then $\deg\det D = \delta(G)$.
4. **Feedback.** State feedback $u = Fx$ can put the eigenvalues of $A + BF$ anywhere if and only if $(A,B)$ is controllable (Wonham). Output feedback $u = Ky$ is much harder, and that's where the geometry really pays off.

## Controllability is surjectivity

Make $k^n$ into a module over the polynomial ring $k[s]$ by letting $s$ act as $A$. The input matrix $B$ then extends to a module map
$$
\beta \colon k[s]^m \to k^n,
\qquad
(f_1, \dots, f_m) \mapsto f_1(A)\,b_1 + \cdots + f_m(A)\,b_m,
$$
where $b_i$ are the columns of $B$. Its image is spanned by the columns of $B, AB, A^2B, \dots$, which is the reachable subspace. So:

::: proposition
$(A,B)$ is controllable if and only if $\beta$ is surjective.
:::

The two classical tests are two ways of checking surjectivity. The Kalman rank condition checks it directly (by Cayley–Hamilton, powers of $A$ beyond $A^{n-1}$ add nothing new). The PBH test is Nakayama's lemma: the cokernel $M = k^n/\operatorname{im}\beta$ is a finite-dimensional $k[s]$-module, so it is zero if and only if $M/(s - \lambda)M = 0$ for every point $\lambda$ of the affine line, and
$$
M/(s-\lambda)M = k^n \big/ \big(\operatorname{im}B + (A - \lambda I)k^n\big),
$$
which vanishes exactly when $[\lambda I - A \;\; B]$ has rank $n$. So the PBH test is "a map of modules is surjective if and only if it is surjective at every point".

Now look at the kernel. Since $\beta$ is surjective, $K = \ker\beta$ is a submodule of $k[s]^m$ of codimension $n$. Submodules of free modules over a principal ideal domain are free, so $K$ is the column span of an $m \times m$ polynomial matrix $D(s)$, and
$$
k^n \cong k[s]^m / D(s)\,k[s]^m
$$
as $k[s]$-modules. In particular $\det D(s)$ is (up to a constant) the characteristic polynomial of $A$. This $D$ is the denominator of the coprime factorization in item 3 above.

::: example The double integrator
For $A = \begin{pmatrix} 0 & 1 \\ 0 & 0 \end{pmatrix}$ and $b = \begin{pmatrix} 0 \\ 1 \end{pmatrix}$, the map $\beta \colon k[s] \to k^2$ sends $1 \mapsto b$ and $s \mapsto Ab = \begin{pmatrix} 1 \\ 0 \end{pmatrix}$, so it is surjective, and its kernel is the ideal $(s^2)$. In general, a controllable system with one input is the module $k[s]/(\chi_A)$, and writing $A$ in the basis $1, s, \dots, s^{n-1}$ gives the companion matrix: the controllable canonical form.
:::

This also tells us what the space of all controllable systems looks like. Changing state coordinates doesn't change the module $K$, and $K$ determines $(A,B)$ up to a change of coordinates, so
$$
\{\text{controllable } (A,B)\}/GL_n \;\cong\; \{\text{submodules } K \subset k[s]^m \text{ of codimension } n\}.
$$
The right-hand side is the Quot scheme of length-$n$ quotients of $\mathcal O^m$ on the affine line, a smooth variety of dimension $nm$. When $m = 1$, a submodule is an ideal $(d(s))$ with $d$ monic of degree $n$, and the Quot scheme is just $k^n$, the space of coefficients of $d$. That's why single-input systems have a canonical form and multi-input systems are harder.

## Realization and the space of all systems

Kalman's realization theorem says that $GL_n$ acts freely on the triples $(A,B,C)$ that are controllable and observable, and that two such triples have the same transfer function exactly when they're in the same orbit. So the transfer function identifies the orbit space with
$$
\mathrm{Rat}_n^{p \times m} = \{\, p \times m \text{ strictly proper rational } G \text{ with } \delta(G) = n \,\}.
$$
In algebraic geometry, this says the quotient of the controllable and observable triples by $GL_n$ is a smooth quasi-projective variety of dimension $(n^2 + nm + pn) - n^2 = n(m+p)$, and its points are transfer functions. (The real-analytic version of this is due to Clark, and the algebraic one to Hazewinkel and to Byrnes and Hurt.)

Two more classical results fit here:

- **Hankel rank.** Expanding $G(s) = \sum_{k \ge 1} H_k s^{-k}$ at infinity, the Markov parameters $H_k = CA^{k-1}B$ fill an infinite block Hankel matrix $\mathcal H = (H_{i+j-1})_{i,j \ge 1}$, and $\delta(G) = \operatorname{rank}\mathcal H$ (Kronecker, Ho–Kalman). So systems of order at most $n$ form a *determinantal variety*: the vanishing of all $(n+1) \times (n+1)$ minors of $\mathcal H$.
- **Canonical forms.** A canonical form is a way of picking one triple $(A,B,C)$ in each orbit, continuously, that is, a section of the principal $GL_n$-bundle of controllable and observable triples over $\mathrm{Rat}_n^{p \times m}$. Hazewinkel showed this bundle has a global continuous section if and only if $m = 1$ or $p = 1$. For multi-input, multi-output systems every parametrization used in system identification has to be local: you need several charts to cover all systems.

## The Hermann–Martin curve

Take a right coprime factorization $G = ND^{-1}$. For each $s \in \mathbb C$ let
$$
h_G(s) = \text{column span of } \begin{pmatrix} D(s) \\ N(s) \end{pmatrix} \subset \mathbb C^{m+p}.
$$
Three things make this a good definition:

- **It's always an $m$-plane.** The Bezout identity $XN + YD = I$ says $\begin{pmatrix} D(s) \\ N(s) \end{pmatrix}$ has a left inverse at every $s$, so it never drops rank. (Conversely, if a polynomial matrix has full rank at every point, it has a polynomial left inverse: this is the Nullstellensatz on the affine line.)
- **It doesn't depend on the factorization.** Two right coprime factorizations differ by a unimodular matrix $U(s)$, which doesn't change the column span.
- **It's the graph of $G$.** Away from the poles, the column span of $\begin{pmatrix} D \\ N \end{pmatrix}$ is the column span of $\begin{pmatrix} I \\ G(s) \end{pmatrix}$, i.e. the graph $\{(u, G(s)u)\}$. Since $G$ is strictly proper, as $s \to \infty$ this tends to $\mathbb C^m \oplus 0$.

So $h_G$ extends to a morphism
$$
h_G \colon \mathbb P^1 \to \mathrm{Gr}(m, m+p),
\qquad
h_G(\infty) = \mathbb C^m \oplus 0,
$$
the *Hermann–Martin curve* of $G$. At the poles, the graph of $G(s)$ blows up but the curve passes smoothly through planes that are no longer graphs.

::: theorem Hermann–Martin
The degree of $h_G$ in the Plücker embedding is the McMillan degree of $G$.
:::

::: proof
Away from the poles, the Plücker coordinates of $h_G(s)$ are the $m \times m$ minors of $\begin{pmatrix} I \\ G(s) \end{pmatrix}$, which are $1$ and the minors of $G(s)$ of every size. Clearing denominators, we multiply all of them by the least common denominator of the minors of $G$, and get polynomials without a common root. The degree of the curve is the largest of their degrees, which is the degree of that common denominator (the minors of $G$ vanish at infinity, so the $1$ wins). One of the classical characterizations of the McMillan degree is exactly this: the degree of the least common denominator of all minors of $G$.
:::

Conversely, a degree-$n$ curve $h \colon \mathbb P^1 \to \mathrm{Gr}(m, m+p)$ with $h(\infty) = \mathbb C^m \oplus 0$ is a graph at all but finitely many points (being a graph is an open condition and holds at $\infty$), so it is the Hermann–Martin curve of a unique $G$. Therefore:
$$
\mathrm{Rat}_n^{p \times m} \;\cong\; \{\text{based curves } \mathbb P^1 \to \mathrm{Gr}(m, m+p) \text{ of degree } n\}.
$$
Transfer functions *are* rational curves in a Grassmannian.

::: example One input, one output
For $m = p = 1$, $G = n(s)/d(s)$ and $\mathrm{Gr}(1,2) = \mathbb P^1$, so $h_G(s) = [\,d(s) : n(s)\,]$ is a degree-$n$ rational map $\mathbb P^1 \to \mathbb P^1$ sending $\infty$ to $[1:0]$. Take $G(s) = 1/(s^2 + s)$, so $h_G(s) = [\,s^2 + s : 1\,]$. With output feedback $u = ky$ the closed-loop characteristic polynomial is $s^2 + s - k$, whose roots are exactly the $s$ with $h_G(s) = [\,k : 1\,]$. So the closed-loop poles at gain $k$ are the fibre $h_G^{-1}([k:1])$, and the root locus is the family of fibres of the Hermann–Martin map. A degree-$n$ map has $n$ points in each fibre: $n$ closed-loop poles.
:::

### Kronecker indices are a splitting type

Pull the tautological subbundle $S \subset \mathbb C^{m+p} \times \mathrm{Gr}(m, m+p)$ back along $h_G$. This gives a rank-$m$ vector bundle on $\mathbb P^1$ of degree $-n$, and by Grothendieck's theorem every such bundle splits as
$$
h_G^*S \cong \mathcal O(-\kappa_1) \oplus \cdots \oplus \mathcal O(-\kappa_m),
\qquad
\kappa_1 + \cdots + \kappa_m = n.
$$
Martin and Hermann showed that the $\kappa_i$ are the controllability (Kronecker) indices of a minimal realization: the column degrees of $D$ when it is column reduced. Dually, the quotient bundle $h_G^*(\mathbb C^{m+p}/S)$ splits according to the observability indices.

This puts two classical facts side by side. Brunovsky's theorem says that controllable pairs are classified up to feedback equivalence by their controllability indices, and Grothendieck's theorem says that vector bundles on $\mathbb P^1$ are classified by their splitting type. And "a generic system has controllability indices as equal as possible" is the fact that a generic bundle on $\mathbb P^1$ is balanced: its $\kappa_i$ differ by at most one.

### Topology of real systems

Over $\mathbb R$, a scalar transfer function of degree $n$ gives a map from the real line $\mathbb{RP}^1$ to itself, and that map has a degree. Brockett showed that the space of real scalar transfer functions of McMillan degree $n$ has exactly $n + 1$ connected components, and they are distinguished by the Cauchy index of $G$, which is (up to sign) this degree. So we can't deform a system continuously into one with a different Cauchy index without its McMillan degree dropping on the way, which is a real obstruction for identification algorithms that search through systems of a fixed order. Over $\mathbb C$, the space of based rational maps is connected, and its topology was worked out by Segal.

## Feedback, poles and Schubert calculus

### State feedback

With $u = Fx$ the closed-loop matrix is $A + BF$. Wonham's theorem says the map
$$
k^{m \times n} \to k^n,
\qquad
F \mapsto \text{coefficients of } \det(sI - A - BF)
$$
is surjective if and only if $(A,B)$ is controllable. For one input it is an affine isomorphism: in the controllable canonical form, $F$ just adds to the coefficients of the characteristic polynomial.

### Output feedback is a Schubert problem

With output feedback $u = Ky$, where $K$ is an $m \times p$ matrix, the closed-loop system has characteristic polynomial
$$
\det(sI - A - BKC) = \det\big(D(s) - K N(s)\big)
= \det\begin{pmatrix} D(s) & K \\ N(s) & I_p \end{pmatrix},
$$
for a minimal realization, normalizing $\det D = \det(sI - A)$. (The first equality is Sylvester's determinant identity, and the second is a Schur complement.) The columns of $\begin{pmatrix} K \\ I_p \end{pmatrix}$ span a $p$-plane
$$
W_K = \{(Ky, y)\} \in \mathrm{Gr}(p, m+p),
$$
and the determinant vanishes exactly when the $m$-plane $h_G(s)$ and the $p$-plane $W_K$ intersect. So:

::: proposition
$s$ is a closed-loop pole under $u = Ky$ if and only if the Hermann–Martin curve at $s$ meets $W_K$.
:::

The $p$-planes meeting a fixed $m$-plane form a Schubert hypersurface in $\mathrm{Gr}(p, m+p)$, and in Plücker coordinates it is a hyperplane section: by Laplace expansion, the determinant above is linear in the Plücker coordinates of $W_K$. So placing the poles at $s_1, \dots, s_n$ means finding a point of $\mathrm{Gr}(p, m+p)$ on $n$ given hyperplanes, namely
$$
W_K \in \bigcap_{i=1}^{n} \{\, W : W \cap h_G(s_i) \ne 0 \,\}.
$$
The Grassmannian has dimension $mp$, so we should expect solutions when $n \le mp$, and finitely many when $n = mp$, as many as the degree of the Grassmannian. That's the theorem of Brockett and Byrnes:

::: theorem Brockett–Byrnes
For a generic system with $mp \ge n$, any closed-loop poles can be placed by complex output feedback. If $mp = n$, a generic choice of poles is achieved by exactly
$$
d(m,p) = \frac{(mp)!\;\, 1!\,2!\cdots(p-1)!}{m!\,(m+1)!\cdots(m+p-1)!}
$$
feedback matrices $K$, counted with multiplicity: the degree of $\mathrm{Gr}(p, m+p)$.
:::

A few values:

| $m$ | $p$ | $n = mp$ | $d(m,p)$ |
|----:|----:|---------:|---------:|
| 1 | any | $p$ | 1 |
| 2 | 2 | 4 | 2 |
| 2 | 3 | 6 | 5 |
| 2 | 4 | 8 | 14 |
| 3 | 3 | 9 | 42 |

Over the reals, complex solutions come in conjugate pairs, so if $d(m,p)$ is odd there is always a real one. When it's even there might not be. The first case is $m = p = 2$, $n = 4$: the Grassmannian $\mathrm{Gr}(2,4)$ is the Klein quadric in $\mathbb P^5$, four hyperplanes in $\mathbb P^5$ meet in a line, and a line meets a quadric in two points, which can be complex conjugates. Then no real feedback places those poles. More generally, Eremenko and Gabrielov showed that when $m$ and $p$ are both even and $n = mp$, real pole placement fails for an open set of systems. On the other hand, with one dimension to spare, $mp > n$, real output feedback can place the poles of a generic system (Wang).

The same picture works for dynamic compensators, with the Grassmannian replaced by a space of rational curves in it (Ravi, Rosenthal and Wang).

## Optimal control and the Lagrangian Grassmannian

The last classical result we'll rephrase is the linear-quadratic regulator. Minimizing $\int_0^\infty x^TQx + u^TRu\,dt$ leads to the algebraic Riccati equation
$$
A^TX + XA - XBR^{-1}B^TX + Q = 0
$$
for a symmetric $n \times n$ matrix $X$, and the optimal control is $u = -R^{-1}B^TXx$. Let
$$
H = \begin{pmatrix} A & -BR^{-1}B^T \\ -Q & -A^T \end{pmatrix}
$$
be the Hamiltonian matrix. Expanding $H\begin{pmatrix} I \\ X \end{pmatrix} = \begin{pmatrix} I \\ X \end{pmatrix}M$ shows:

::: proposition
$X$ solves the algebraic Riccati equation if and only if its graph $\{(x, Xx)\}$ is an $H$-invariant $n$-plane in $k^{2n}$.
:::

$X$ being symmetric means its graph is Lagrangian for the standard symplectic form, so the solutions of the Riccati equation are the $H$-invariant points of the Lagrangian Grassmannian that are graphs. The stabilizing solution is the one whose graph is the sum of the generalized eigenspaces of $H$ for eigenvalues in the left half-plane. The Riccati *differential* equation is the linear flow $e^{tH}$ acting on the Lagrangian Grassmannian: if $\frac{d}{dt}\begin{pmatrix} U \\ V \end{pmatrix} = H\begin{pmatrix} U \\ V \end{pmatrix}$, then $X = VU^{-1}$ satisfies $-\dot X = A^TX + XA - XBR^{-1}B^TX + Q$, and "finite escape time" is the graph leaving the chart where it's a graph. Hermann and Martin also studied the Riccati equation from exactly this point of view.

## A dictionary

| control theory | algebraic geometry |
|---|---|
| controllable pair $(A,B)$ | surjection $k[s]^m \to k^n$ of $k[s]$-modules |
| PBH test | Nakayama's lemma: surjective at every point |
| controllable systems up to coordinates | Quot scheme of length-$n$ quotients of $\mathcal O^m$ on $\mathbb A^1$ |
| minimal realization, unique up to similarity | free $GL_n$-orbit, geometric quotient |
| transfer function of McMillan degree $n$ | based rational curve of degree $n$ in $\mathrm{Gr}(m, m+p)$ |
| coprime factorization $G = ND^{-1}$ | homogeneous coordinates of the curve |
| Bezout identity | Nullstellensatz on the affine line |
| Kronecker indices | splitting type of the pulled-back bundle |
| Brunovsky classification | Grothendieck's theorem on $\mathbb P^1$ |
| Hankel rank | determinantal variety |
| no global canonical form for MIMO systems | nontrivial principal $GL_n$-bundle |
| Cauchy index | degree of the real map $\mathbb{RP}^1 \to \mathbb{RP}^1$ |
| root locus (SISO) | fibres of the Hermann–Martin map |
| static output feedback $K$ | $p$-plane $W_K \in \mathrm{Gr}(p, m+p)$ |
| closed-loop poles | where the curve meets a Schubert hypersurface |
| number of pole-placing feedbacks when $n = mp$ | degree of the Grassmannian |
| Riccati equation | invariant Lagrangian subspaces of $H$ |

## Further reading

- R. Hermann and C. F. Martin, "Applications of algebraic geometry to systems theory, Part I", *IEEE Transactions on Automatic Control* (1977).
- C. F. Martin and R. Hermann, "Applications of algebraic geometry to systems theory: the McMillan degree and Kronecker indices of transfer functions as topological and holomorphic system invariants", *SIAM Journal on Control and Optimization* (1978).
- R. W. Brockett, "Some geometric questions in the theory of linear systems", *IEEE Transactions on Automatic Control* (1976).
- R. W. Brockett and C. I. Byrnes, "Multivariable Nyquist criteria, root loci, and pole placement: a geometric viewpoint", *IEEE Transactions on Automatic Control* (1981).
- M. Hazewinkel, "Moduli and canonical forms for linear dynamical systems II: the topological case", *Mathematical Systems Theory* (1977).
- A. Eremenko and A. Gabrielov, "Pole placement by static output feedback for generic linear systems", *SIAM Journal on Control and Optimization* (2002).
