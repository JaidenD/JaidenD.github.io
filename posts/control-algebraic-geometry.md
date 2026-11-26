---
title: Linear Control Theory as Algebraic Geometry
date: 2026-11-26
summary: A way of formalizing linear control theory that I find really satisfying, built on Hermann and Martin's observation that a transfer function is a rational curve in a Grassmannian.
toc: true
draft: true
---

Most of us learn linear control theory as a collection of matrix tests. A system is controllable if a certain matrix has full rank, two realizations are the same if certain other matrices are similar, and the poles of a feedback loop are the roots of a certain determinant. Every one of these results has a proof, but it's easy to come away feeling that the subject is a list of tricks rather than a theory.

I find there's a much more satisfying way to organize it, and I want to present it today. In 1977 Robert Hermann and Clyde Martin noticed that the transfer function of a linear system with $m$ inputs and $p$ outputs is the same thing as a rational curve in a Grassmannian: the space of $m$-dimensional subspaces of $\mathbb C^{m+p}$. That one observation turns a surprising amount of linear control theory into classical algebraic geometry. Controllability becomes the statement that a certain map of modules is onto. The McMillan degree of a system, the number of states it really needs, becomes the degree of a curve. And the old question of where output feedback can put the poles of a system becomes a problem in Schubert calculus, the nineteenth-century art of counting how many lines meet a given collection of lines. (Hermann's work is a big part of why I got interested in applying geometry to engineering in the first place, so this post is also a bit of a tribute.)

The plan is to walk through the main results of linear control theory one at a time and see what each of them looks like from the algebro-geometric side. I'll try to say why each translation is more than a change of vocabulary. You won't need more than the definitions of a module, a Grassmannian and a vector bundle on the projective line, and I'll explain the rest as we go.

## The classical story

Let's start by recalling what we're going to translate. A linear time-invariant system is
$$
\dot x = Ax + Bu,
\qquad
y = Cx,
$$
with state $x \in k^n$, input $u \in k^m$ and output $y \in k^p$, where $k$ is $\mathbb R$ or $\mathbb C$. I like to think of $u$ as the knobs we can turn, $x$ as the internal state we usually can't see, and $y$ as what our sensors report. If we start at rest and take Laplace transforms, the inputs and outputs are related by $y(s) = G(s)u(s)$, where
$$
G(s) = C(sI - A)^{-1}B
$$
is the transfer function. All of the internal machinery has been packaged into a single $p \times m$ matrix of rational functions of $s$, and that matrix is strictly proper, meaning every entry goes to zero as $s \to \infty$.

Here are the four pieces of the classical theory I want to revisit.

1. **Controllability.** Can the inputs steer the state anywhere we like? The answer is yes, and $(A,B)$ is called controllable, exactly when $\operatorname{rank}\,[B \;\, AB \;\, \cdots \;\, A^{n-1}B] = n$. Equivalently, $\operatorname{rank}\,[\lambda I - A \;\; B] = n$ for every $\lambda \in \mathbb C$; this is the Popov–Belevitch–Hautus (PBH) test. Observability, whether the outputs determine the state, is the dual question about $(C,A)$.
2. **Realization.** Going the other way, every strictly proper rational $G$ can be written as $C(sI-A)^{-1}B$ for some $(A,B,C)$. Kalman's realization theorem says that a realization uses the fewest possible states exactly when it is controllable and observable, and that any two such minimal realizations differ only by a change of state coordinates, $(A,B,C) \mapsto (TAT^{-1}, TB, CT^{-1})$. The fewest possible number of states is called the *McMillan degree* $\delta(G)$.
3. **Coprime factorizations.** Just as a rational number can be written as a fraction in lowest terms, a transfer function can be written as $G = N(s)D(s)^{-1}$ with $N$ and $D$ polynomial matrices that have no common factor. The precise condition, right coprimeness, is a Bezout identity: there are polynomial matrices $X$ and $Y$ with $XN + YD = I$. When this holds, $\deg\det D = \delta(G)$.
4. **Feedback.** With state feedback $u = Fx$ we can put the eigenvalues of $A + BF$ wherever we like, as long as $(A,B)$ is controllable (Wonham's theorem). With output feedback $u = Ky$, which is what we can actually implement when we can't measure the whole state, things are much harder. This is where the geometry will really pay off.

## Controllability is surjectivity

Here's the first translation, and it's the one that sets up everything else. The idea is to stop thinking of $A$ as a matrix and start thinking of it as a way of making the state space into a module.

Let the polynomial ring $k[s]$ act on $k^n$ by letting $s$ act as $A$, so a polynomial $f(s)$ acts as the matrix $f(A)$. Once we do that, the input matrix $B$ extends to a map of $k[s]$-modules
$$
\beta \colon k[s]^m \to k^n,
\qquad
(f_1, \dots, f_m) \mapsto f_1(A)\,b_1 + \cdots + f_m(A)\,b_m,
$$
where $b_1, \dots, b_m$ are the columns of $B$. The image of $\beta$ is spanned by the columns of $B, AB, A^2B, \dots$, which is exactly the set of states we can reach from rest. So controllability, which sounds like a statement about trajectories, is really a statement about a single algebraic map:

::: proposition
$(A,B)$ is controllable if and only if $\beta$ is surjective.
:::

From this point of view, the two classical tests are just two ways of checking that a map is onto. The Kalman rank condition checks it directly: by Cayley–Hamilton, powers of $A$ beyond $A^{n-1}$ add nothing new, so the image is spanned by the columns of $[B \;\, AB \;\, \cdots \;\, A^{n-1}B]$.

The PBH test is more interesting, because it is Nakayama's lemma in disguise, which is the algebraist's version of "it's enough to check locally". The cokernel $M = k^n/\operatorname{im}\beta$ is a finite-dimensional $k[s]$-module, and such a module is zero if and only if $M/(s - \lambda)M = 0$ for every point $\lambda$ of the affine line. Unwinding the definitions,
$$
M/(s-\lambda)M = k^n \big/ \big(\operatorname{im}B + (A - \lambda I)k^n\big),
$$
and this vanishes exactly when $[\lambda I - A \;\; B]$ has rank $n$. So the PBH test is the statement that a map of modules is surjective if and only if it is surjective at every point.

The kernel of $\beta$ turns out to be just as interesting as its image. When $\beta$ is surjective, $K = \ker\beta$ is a submodule of $k[s]^m$ whose quotient has dimension $n$. Submodules of free modules over a principal ideal domain are free, so $K$ is the column span of a square $m \times m$ polynomial matrix $D(s)$, and
$$
k^n \cong k[s]^m / D(s)\,k[s]^m
$$
as $k[s]$-modules. Comparing characteristic polynomials on both sides, $\det D(s)$ is the characteristic polynomial of $A$, up to a constant. And this $D$ is precisely the denominator of the coprime factorization from item 3 above, so coprime factorizations aren't an extra trick: they come for free from the module structure.

::: example The double integrator
Take $A = \begin{pmatrix} 0 & 1 \\ 0 & 0 \end{pmatrix}$ and $b = \begin{pmatrix} 0 \\ 1 \end{pmatrix}$. The map $\beta \colon k[s] \to k^2$ sends $1 \mapsto b$ and $s \mapsto Ab = \begin{pmatrix} 1 \\ 0 \end{pmatrix}$, so it is surjective, and its kernel is the ideal $(s^2)$, the characteristic polynomial of $A$. The same thing happens for any controllable system with one input: as a module it is $k[s]/(\chi_A)$, and if we write $A$ in the basis $1, s, \dots, s^{n-1}$ we get the companion matrix. That's the controllable canonical form, and now we know where it comes from.
:::

This also answers a question that sounds like it ought to be hard: what does the set of all controllable systems look like, once we stop distinguishing between systems that differ by a change of coordinates? Changing coordinates doesn't change the submodule $K$, and $K$ determines $(A,B)$ up to a change of coordinates. So
$$
\{\text{controllable } (A,B)\}/GL_n \;\cong\; \{\text{submodules } K \subset k[s]^m \text{ of codimension } n\}.
$$
Algebraic geometers call the right-hand side the Quot scheme of length-$n$ quotients of $\mathcal O^m$ on the affine line, and it's a smooth variety of dimension $nm$. When $m = 1$, a submodule is just an ideal $(d(s))$ generated by a monic polynomial of degree $n$, and the Quot scheme is plain old $k^n$: the coefficients of $d$.

This is the first appearance of a theme that runs through the whole post. Single-input systems are easy because their spaces of systems are just affine spaces. Multi-input systems are hard because their spaces have interesting global geometry, and we'll keep running into consequences of that.

## Realization and the space of all systems

Kalman's realization theorem has a nice geometric reading too. It says two things: that $GL_n$ acts freely on the triples $(A,B,C)$ that are controllable and observable (a minimal realization has no nontrivial symmetries), and that two such triples have the same transfer function exactly when they lie in the same orbit. Put together, the transfer function identifies the space of orbits with
$$
\mathrm{Rat}_n^{p \times m} = \{\, p \times m \text{ strictly proper rational } G \text{ with } \delta(G) = n \,\}.
$$
In the language of algebraic geometry, the quotient of the controllable and observable triples by $GL_n$ is a smooth quasi-projective variety of dimension $(n^2 + nm + pn) - n^2 = n(m+p)$, and its points are transfer functions. (The real-analytic version of this is due to Clark, and the algebraic one to Hazewinkel and to Byrnes and Hurt.) This is the space that a system identification algorithm searches through when it fits an $n$-state model to data, so its geometry has real consequences.

Two more classical results fit naturally here.

The first is about recognizing the McMillan degree from data. Expanding the transfer function at infinity, $G(s) = \sum_{k \ge 1} H_k s^{-k}$, gives the Markov parameters $H_k = CA^{k-1}B$, which are the derivatives of the impulse response $Ce^{At}B$ at $t = 0$. Arrange them in an infinite block Hankel matrix $\mathcal H = (H_{i+j-1})_{i,j \ge 1}$. A theorem going back to Kronecker, and in this form to Ho and Kalman, says that $\delta(G) = \operatorname{rank}\mathcal H$. Geometrically, the systems with at most $n$ states form a *determinantal variety*: the common zeros of all the $(n+1) \times (n+1)$ minors of $\mathcal H$.

The second explains something that people doing system identification run into in practice. A *canonical form* is a rule that picks one triple $(A,B,C)$ out of each orbit, continuously, so that we can parametrize systems by the entries of that triple. Geometrically, that's a continuous section of the principal $GL_n$-bundle of controllable and observable triples over $\mathrm{Rat}_n^{p \times m}$. Hazewinkel showed this bundle has a global continuous section if and only if $m = 1$ or $p = 1$. So for genuinely multi-input, multi-output systems, any parametrization has to be local, and identification software needs several overlapping charts to cover every system. That isn't a failure of imagination on anyone's part: the bundle is topologically twisted.

## The Hermann–Martin curve

Now for the main construction. The idea fits in one sentence: for each frequency $s$, look at the graph of $G(s)$, the set of input–output pairs $(u, G(s)u)$. That's an $m$-dimensional subspace of $\mathbb C^{m+p}$, so as $s$ varies we get a curve in the Grassmannian $\mathrm{Gr}(m, m+p)$. The only problem is that the graph doesn't make sense at the poles, where $G(s)$ blows up, and the coprime factorization is exactly what fixes this.

Take a right coprime factorization $G = ND^{-1}$, and for each $s \in \mathbb C$ let
$$
h_G(s) = \text{column span of } \begin{pmatrix} D(s) \\ N(s) \end{pmatrix} \subset \mathbb C^{m+p}.
$$
Three things make this a good definition.

- **It's always an $m$-plane.** The Bezout identity $XN + YD = I$ says that $\begin{pmatrix} D(s) \\ N(s) \end{pmatrix}$ has a left inverse at every $s$, so it never drops rank. (The converse is also true: a polynomial matrix with full rank at every point has a polynomial left inverse. That's the Nullstellensatz on the affine line, and it's why "coprime" and "no common zeros" mean the same thing.)
- **It doesn't depend on the factorization.** Two right coprime factorizations differ by a unimodular matrix $U(s)$, and multiplying on the right by $U(s)$ doesn't change the column span.
- **It really is the graph of $G$.** Away from the poles, the column span of $\begin{pmatrix} D \\ N \end{pmatrix}$ is the column span of $\begin{pmatrix} I \\ G(s) \end{pmatrix}$, which is the graph $\{(u, G(s)u)\}$. Since $G$ is strictly proper, as $s \to \infty$ this tends to $\mathbb C^m \oplus 0$.

So $h_G$ extends to a morphism
$$
h_G \colon \mathbb P^1 \to \mathrm{Gr}(m, m+p),
\qquad
h_G(\infty) = \mathbb C^m \oplus 0,
$$
called the *Hermann–Martin curve* of $G$. I like to picture it this way: as $s$ moves along the Riemann sphere, the graph of $G(s)$ swings around inside $\mathbb C^{m+p}$, and at a pole, where the graph "blows up", the curve passes perfectly smoothly through subspaces that just don't happen to be graphs.

The first thing we can read off from the curve is its degree, and it's exactly the number we care about most.

::: theorem Hermann–Martin
The degree of $h_G$ in the Plücker embedding is the McMillan degree of $G$.
:::

::: proof
Away from the poles, the Plücker coordinates of $h_G(s)$ are the $m \times m$ minors of $\begin{pmatrix} I \\ G(s) \end{pmatrix}$, which are $1$ together with the minors of $G(s)$ of every size. To clear denominators, multiply all of them by the least common denominator of the minors of $G$. This gives polynomials without a common root, and the degree of the curve is the largest of their degrees. Since $G$ is strictly proper, its minors vanish at infinity, so the largest degree belongs to the coordinate that started out as $1$: it's the degree of the common denominator itself. One of the classical characterizations of the McMillan degree is exactly this: the degree of the least common denominator of all the minors of $G$.
:::

The construction also runs backwards. Suppose $h \colon \mathbb P^1 \to \mathrm{Gr}(m, m+p)$ is any curve of degree $n$ with $h(\infty) = \mathbb C^m \oplus 0$. Being a graph is an open condition, and it holds at $\infty$, so $h(s)$ is the graph of some $p \times m$ matrix at all but finitely many points, and that matrix depends rationally on $s$. So $h$ is the Hermann–Martin curve of a unique $G$, and
$$
\mathrm{Rat}_n^{p \times m} \;\cong\; \{\text{based curves } \mathbb P^1 \to \mathrm{Gr}(m, m+p) \text{ of degree } n\}.
$$
In other words, transfer functions *are* rational curves in a Grassmannian, and the McMillan degree is how curved they are.

::: example One input, one output
When $m = p = 1$ the Grassmannian $\mathrm{Gr}(1,2)$ is just $\mathbb P^1$, and if $G = n(s)/d(s)$ the curve is $h_G(s) = [\,d(s) : n(s)\,]$: a degree-$n$ rational map from the Riemann sphere to itself, sending $\infty$ to $[1:0]$.

Here's a nice way to see a root locus with this. Take $G(s) = 1/(s^2 + s)$, so $h_G(s) = [\,s^2 + s : 1\,]$. If we close the loop with output feedback $u = ky$, the closed-loop characteristic polynomial is $s^2 + s - k$, and its roots are exactly the values of $s$ where $h_G(s) = [\,k : 1\,]$. So the closed-loop poles at gain $k$ are the fibre $h_G^{-1}([k:1])$, and the root locus is nothing but the family of fibres of the Hermann–Martin map. A degree-$n$ map has $n$ points in each fibre, which is why there are always $n$ closed-loop poles.
:::

### Kronecker indices are a splitting type

The curve carries more refined information than its degree. Every point of the Grassmannian is a subspace, so pulling back the tautological subbundle $S \subset \mathbb C^{m+p} \times \mathrm{Gr}(m, m+p)$ along $h_G$ gives a rank-$m$ vector bundle on $\mathbb P^1$, of degree $-n$. Grothendieck's theorem says every vector bundle on $\mathbb P^1$ splits as a sum of line bundles, so
$$
h_G^*S \cong \mathcal O(-\kappa_1) \oplus \cdots \oplus \mathcal O(-\kappa_m),
\qquad
\kappa_1 + \cdots + \kappa_m = n.
$$
Martin and Hermann showed that the numbers $\kappa_i$ are the controllability (Kronecker) indices of a minimal realization, which in the polynomial picture are the column degrees of $D$ once it's put in column-reduced form. Dually, the quotient bundle $h_G^*(\mathbb C^{m+p}/S)$ splits according to the observability indices.

I really like how this puts two classical facts side by side. Brunovsky's theorem says that controllable pairs are classified, up to feedback equivalence, by their controllability indices. Grothendieck's theorem says that vector bundles on $\mathbb P^1$ are classified by their splitting type. And the engineering folklore that "a generic system has controllability indices as equal as possible" is the algebro-geometric fact that a generic bundle on $\mathbb P^1$ is balanced, meaning its $\kappa_i$ differ by at most one.

### Topology of real systems

So far we've mostly worked over $\mathbb C$, but real systems have some extra structure that's worth a look. A real scalar transfer function of degree $n$ gives a map from the real projective line $\mathbb{RP}^1$ to itself, and a map from a circle to a circle has a degree. Brockett showed that the space of real scalar transfer functions of McMillan degree $n$ has exactly $n + 1$ connected components, distinguished by the Cauchy index of $G$, which is (up to sign) this degree.

This has a very practical consequence. You can't continuously deform one real system into another with a different Cauchy index without the McMillan degree dropping somewhere along the way. So an identification algorithm that searches over $n$-state models by local improvements can get stuck in the wrong component and never find the best fit. Over $\mathbb C$, by contrast, the space of based rational maps is connected, and its topology was worked out in detail by Segal.

## Feedback, poles and Schubert calculus

For engineers, the real payoff comes when we close the loop.

### State feedback

The easy case first. With state feedback $u = Fx$ the closed-loop matrix is $A + BF$, and Wonham's theorem says the map
$$
k^{m \times n} \to k^n,
\qquad
F \mapsto \text{coefficients of } \det(sI - A - BF),
$$
is surjective if and only if $(A,B)$ is controllable. For a single input, the module picture makes this obvious: in the controllable canonical form, $F$ simply adds to the coefficients of the characteristic polynomial, so the map is an affine isomorphism.

### Output feedback is a Schubert problem

In practice we rarely get to measure the whole state, so the question engineers actually face is what we can do with output feedback $u = Ky$, where $K$ is an $m \times p$ matrix of gains. This turned out to be a hard problem, and the Hermann–Martin curve is what turned it into one we know how to think about.

For a minimal realization, normalized so that $\det D = \det(sI - A)$, the closed-loop characteristic polynomial is
$$
\det(sI - A - BKC) = \det\big(D(s) - K N(s)\big)
= \det\begin{pmatrix} D(s) & K \\ N(s) & I_p \end{pmatrix}.
$$
The first equality is Sylvester's determinant identity and the second is a Schur complement. Now look at the last matrix. Its left block spans the $m$-plane $h_G(s)$, and its right block spans the $p$-plane
$$
W_K = \{(Ky, y)\} \in \mathrm{Gr}(p, m+p),
$$
which is just the graph of the feedback law. The determinant of an $(m+p) \times (m+p)$ matrix vanishes exactly when its columns are dependent, that is, when these two subspaces intersect. So we've turned the pole placement question into a question about subspaces meeting:

::: proposition
$s$ is a closed-loop pole under $u = Ky$ if and only if the Hermann–Martin curve at $s$ meets $W_K$.
:::

The $p$-planes meeting a fixed $m$-plane form what's called a Schubert hypersurface in $\mathrm{Gr}(p, m+p)$, and in Plücker coordinates it's a hyperplane section: by Laplace expansion, the determinant above is linear in the Plücker coordinates of $W_K$. So placing the closed-loop poles at chosen points $s_1, \dots, s_n$ means finding a point of $\mathrm{Gr}(p, m+p)$ that lies on $n$ given hyperplanes:
$$
W_K \in \bigcap_{i=1}^{n} \{\, W : W \cap h_G(s_i) \ne 0 \,\}.
$$
Now we can simply count. The Grassmannian has dimension $mp$, and each hyperplane cuts the dimension down by one, so we should expect solutions when $n \le mp$, and finitely many of them when $n = mp$, as many as the degree of the Grassmannian. That's exactly the theorem of Brockett and Byrnes:

::: theorem Brockett–Byrnes
For a generic system with $mp \ge n$, any closed-loop poles can be placed by complex output feedback. If $mp = n$, a generic choice of poles is achieved by exactly
$$
d(m,p) = \frac{(mp)!\;\, 1!\,2!\cdots(p-1)!}{m!\,(m+1)!\cdots(m+p-1)!}
$$
feedback matrices $K$, counted with multiplicity: the degree of $\mathrm{Gr}(p, m+p)$.
:::

Here are a few values:

| $m$ | $p$ | $n = mp$ | $d(m,p)$ |
|----:|----:|---------:|---------:|
| 1 | any | $p$ | 1 |
| 2 | 2 | 4 | 2 |
| 2 | 3 | 6 | 5 |
| 2 | 4 | 8 | 14 |
| 3 | 3 | 9 | 42 |

The case $m = p = 2$ is worth stopping on, because it's a problem Schubert himself solved. A 2-plane in $\mathbb C^4$ is the same thing as a line in projective 3-space, and two 2-planes meet exactly when the corresponding lines meet. So placing the four poles of a system with two inputs, two outputs and four states is the same as finding the lines in space that meet four given lines, and the classical answer is that there are two of them.

That example also shows where things get subtle over the reals. Complex solutions come in conjugate pairs, so whenever $d(m,p)$ is odd, real data always has at least one real solution. When $d(m,p)$ is even there might not be one. In the $m = p = 2$ case, $\mathrm{Gr}(2,4)$ is the Klein quadric in $\mathbb P^5$, four hyperplanes in $\mathbb P^5$ meet in a line, and a line meets a quadric in two points, which can perfectly well be complex conjugates. When that happens, no real feedback places those poles. More generally, Eremenko and Gabrielov showed that when $m$ and $p$ are both even and $n = mp$, real pole placement fails for an open set of systems, so it isn't just a measure-zero accident. On the other hand, as soon as there's one dimension to spare, $mp > n$, real output feedback can place the poles of a generic system (Wang).

The same picture extends to dynamic compensators, where the feedback law is itself a linear system: the Grassmannian gets replaced by a space of rational curves in it (Ravi, Rosenthal and Wang).

## Optimal control and the Lagrangian Grassmannian

The last stop is optimal control, and it's a nice example of the Grassmannian showing up somewhere unexpected. The linear-quadratic regulator asks for the input that minimizes $\int_0^\infty x^TQx + u^TRu\,dt$. The answer is the feedback law $u = -R^{-1}B^TXx$, where $X$ is a symmetric $n \times n$ matrix solving the algebraic Riccati equation
$$
A^TX + XA - XBR^{-1}B^TX + Q = 0.
$$
The Riccati equation is quadratic in $X$, which makes it look harder than everything else in linear control. The trick is to linearize it by going up a dimension. Let
$$
H = \begin{pmatrix} A & -BR^{-1}B^T \\ -Q & -A^T \end{pmatrix}
$$
be the Hamiltonian matrix of the problem. Expanding $H\begin{pmatrix} I \\ X \end{pmatrix} = \begin{pmatrix} I \\ X \end{pmatrix}M$ and comparing blocks gives:

::: proposition
$X$ solves the algebraic Riccati equation if and only if its graph $\{(x, Xx)\}$ is an $H$-invariant $n$-plane in $k^{2n}$.
:::

So a quadratic matrix equation turns into a linear eigenvalue-type problem. The condition that $X$ is symmetric says its graph is Lagrangian for the standard symplectic form, so the solutions of the Riccati equation are the $H$-invariant points of the Lagrangian Grassmannian that happen to be graphs. The stabilizing solution, the one LQR actually uses, is the one whose graph is the sum of the generalized eigenspaces of $H$ for eigenvalues in the left half-plane.

The Riccati *differential* equation has an equally clean description: it's the linear flow $e^{tH}$, watched on the Lagrangian Grassmannian. If $\frac{d}{dt}\begin{pmatrix} U \\ V \end{pmatrix} = H\begin{pmatrix} U \\ V \end{pmatrix}$, then $X = VU^{-1}$ satisfies $-\dot X = A^TX + XA - XBR^{-1}B^TX + Q$. The notorious "finite escape time" of Riccati equations is just the moving subspace leaving the chart where it's a graph; nothing actually goes wrong upstairs. Hermann and Martin also studied the Riccati equation from exactly this point of view.

## A dictionary

To finish, here is everything we've translated, in one place:

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

What I like most about this point of view is that it explains *why* the matrix tests look the way they do, and it tells us where to expect trouble. Every place where multi-input systems are genuinely harder than single-input ones (no global canonical forms, real pole placement failing, identification getting stuck) is a place where the underlying space stops being an affine space and starts having interesting global geometry. That's a much better reason for a result to be hard than "the matrices got bigger".

## Further reading

- R. Hermann and C. F. Martin, "Applications of algebraic geometry to systems theory, Part I", *IEEE Transactions on Automatic Control* (1977).
- C. F. Martin and R. Hermann, "Applications of algebraic geometry to systems theory: the McMillan degree and Kronecker indices of transfer functions as topological and holomorphic system invariants", *SIAM Journal on Control and Optimization* (1978).
- R. W. Brockett, "Some geometric questions in the theory of linear systems", *IEEE Transactions on Automatic Control* (1976).
- R. W. Brockett and C. I. Byrnes, "Multivariable Nyquist criteria, root loci, and pole placement: a geometric viewpoint", *IEEE Transactions on Automatic Control* (1981).
- M. Hazewinkel, "Moduli and canonical forms for linear dynamical systems II: the topological case", *Mathematical Systems Theory* (1977).
- A. Eremenko and A. Gabrielov, "Pole placement by static output feedback for generic linear systems", *SIAM Journal on Control and Optimization* (2002).
