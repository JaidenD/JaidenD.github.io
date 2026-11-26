---
title: Symmetry of Differential Equations
date: 2023-06
summary: Finding the symmetries and conservation laws of a PDE by linearizing it, with worked examples for the heat and Burgers equations, and what they are good for in engineering.
---

## Introduction

During the current academic year I have been working on research pertaining to symmetries of PDEs. In this post I'd like to briefly explain the method, and then explain why the 'modern' method, utilizing the variational derivative, is computationally and conceptually clearer than the historical version (which is still often used).

## Algebraic equations

To start, consider a circle. We can rotate and reflect a circle and it will look the same. How can we express this mathematically?

In the plane we can express a circle with radius $r$ as the set of points satisfying the algebraic equation $x^2 + y^2 = r^2$, a symmetry of the circle should leave this algebraic equation unchanged. A rotation $R(\theta)$ in the plane is described by the linear transformation

$$
\begin{pmatrix}
x' \\ y'
\end{pmatrix}
=
\begin{pmatrix}
\cos\theta & -\sin\theta \\
\sin\theta & \cos\theta
\end{pmatrix}
\begin{pmatrix}
x \\ y
\end{pmatrix}
$$

and a reflection is described by

$$
\begin{pmatrix}
x' \\ y'
\end{pmatrix}
=
\begin{pmatrix}
\cos\theta & \sin\theta \\
\sin\theta & -\cos\theta
\end{pmatrix}
\begin{pmatrix}
x \\ y
\end{pmatrix}
$$

Let's check by substitution if rotation is indeed a symmetry:

$$
\begin{aligned}
x'^2+y'^2 &=(x\cos\theta -y\sin\theta)^2+ (x \sin\theta + y\cos\theta)^2\\ &=x^2(\cos^2\theta+\sin^2\theta)+y^2(\sin^2\theta+\cos^2\theta)+xy(-\sin\theta\cos\theta+\sin\theta\cos\theta)\\ &= x^2+y^2
\end{aligned}
$$

Next, we want to think about how a circle can be infinitesimally rotated, what does that even look like? The picture is we want each point to walk infinitesimally in the direction of its tangent vector. Mathematically this corresponds to the derivative of the rotation matrix evaluated at zero:

$$
\frac{d}{d\theta}\Bigr|_{\theta=0}
\begin{pmatrix}
\cos\theta & -\sin\theta \\
\sin\theta & \cos\theta
\end{pmatrix} =
\begin{pmatrix}
0 & -1 \\
1 & 0
\end{pmatrix}
$$

This transformation is called the infinitesimal generator of a rotation.

From this it is easy to see that $x^2+y^2$ remains invariant. An important thing to note here is that the expression for an infinitesimal rotation is simpler than that of a finite rotation.

Finally, we can represent an infinitesimal rotation as an action on functions, this will reduce the problem of finding symmetries to one of solving a PDE. First note that we can define a representation of $G$ to act on functions by $(\pi(g)f)(x) = f(g^{-1}x)$. So, an infinitesimal rotation acts on a function $f(x,y)$ by $f( \begin{pmatrix} 0 & 1 \\ -1 & 0 \end{pmatrix} \begin{pmatrix} x \\ y \end{pmatrix}) = f(y, -x)$. Since we want this transformation to be infinitesimal we compute the first order Taylor expansion of $f(y, -x)$, this gives $y\partial_x f - x\partial_y f$.

Therefore, we can conclude that an algebraic equation is rotationally symmetric if the rotational vector field vanishes on the equation.

![The vector field y∂x − x∂y drawn around a circle](../PDEsymm_blog/inf_rotation.png "The vector field $y\partial_x-x\partial_y$ describes an infinitesimal symmetry of the circle.")

## Symmetries of Differential equations

In the previous example, we saw how symmetries correspond to vector fields, which can be solved as a linear PDE to yield classes of algebraic equations that are invariant under that symmetry. We will now develop an analogous method for PDEs, such that if we know the symmetries of a PDE, we can generate classes of solutions invariant under those symmetries.

First order of business is to understand how symmetries act on PDEs, then develop a procedure which takes a PDE and gives its symmetries.

## Jet Space

We will denote the space of independent variables, like $x, y, z$ in $\mathbb{R}^3$, with $X$, and the space of dependent variables, that is variables $u$ which depend on $x \in X$, with $U$. We know how transformations act on $X$, so how do they act on $U$?

Consider the action of $SO(2)$ on $X \times U \simeq \mathbb{R} \times \mathbb{R}$ where $(x, u)$ is transformed into $(\tilde{x}, \tilde{u})$, then
$$
(\tilde{x}, \tilde{u}) = g(\theta)(x, u) = (x \cos \theta - u \sin \theta, x \sin \theta + u \cos \theta)
$$
Now for example, take $u = f(x) = ax + b$, under the same transformation this gives
$$
(\tilde{x}, \tilde{f}(\tilde{x})) = ((x \cos \theta - (ax + b) \sin \theta), x \sin \theta + (ax + b) \cos \theta)
$$
Isolating $x$ we have
$$
x = \frac{\tilde{x} + b \sin \theta}{\cos \theta - a \sin \theta},
$$
and subbing into the expression for $\tilde{f}(\tilde{x})$ we get that $f$ transformed under the action of $SO(2)$ is
$$
\tilde{f}(\tilde{x}) = \frac{\sin \theta + a \cos \theta}{\cos \theta - a \sin \theta} \tilde{x} + \frac{b}{\cos \theta - a \sin \theta}
$$

Next, we'd like to prolong group action into "jet space", that is, extend our transformations to the derivatives of $u$. This will allow us to apply group actions to PDEs and determine their invariance.

## Prolonging the action

Look again at what the rotation did to the line $u = ax + b$: it gave another line, with slope
$$
\tilde a = \frac{\sin\theta + a\cos\theta}{\cos\theta - a\sin\theta}.
$$
The slope of the graph is $u_x$, so the rotation also acts on first derivatives. Differentiating at $\theta = 0$ gives $\frac{d\tilde a}{d\theta}\big|_{\theta=0} = 1 + a^2$, so the infinitesimal rotation $X = -u\,\partial_x + x\,\partial_u$, extended to the coordinates $(x, u, u_x)$, is
$$
\operatorname{pr}X = -u\,\partial_x + x\,\partial_u + (1 + u_x^2)\,\partial_{u_x}.
$$
This is the *prolongation* of $X$. The space with coordinates $(x, u, u_x, u_{xx}, \dots)$ is jet space, and a differential equation is just an equation between these coordinates, like $u_{xx} = 0$.

For a general $X = \xi(x,u)\,\partial_x + \eta(x,u)\,\partial_u$ the prolongation is
$$
\operatorname{pr}X = \xi\,\partial_x + \eta\,\partial_u + \eta^{x}\,\partial_{u_x} + \eta^{xx}\,\partial_{u_{xx}} + \cdots,
\qquad
\eta^{x} = D_x\eta - u_x\,D_x\xi,
\qquad
\eta^{xx} = D_x\eta^{x} - u_{xx}\,D_x\xi,
$$
and so on, where $D_x = \partial_x + u_x\,\partial_u + u_{xx}\,\partial_{u_x} + \cdots$ is the total derivative: differentiate with respect to $x$, remembering that $u$ depends on $x$. For the rotation, $\eta^x = D_x(x) - u_x D_x(-u) = 1 + u_x^2$, as above.

An equation $\Delta = 0$ is invariant when $\operatorname{pr}X(\Delta) = 0$ whenever $\Delta = 0$. For $\Delta = u_{xx}$ and the rotation, $\operatorname{pr}X(\Delta) = \eta^{xx} = 3u_xu_{xx}$, which vanishes when $u_{xx} = 0$: rotations take lines to lines, as they should.

That's the historical method, and it works, but the formulas get out of hand quickly. With $\xi$ and $\eta$ left general,
$$
\eta^{xx} = \eta_{xx} + (2\eta_{xu} - \xi_{xx})\,u_x + (\eta_{uu} - 2\xi_{xu})\,u_x^2 - \xi_{uu}\,u_x^3 + (\eta_u - 2\xi_x)\,u_{xx} - 3\xi_u\,u_xu_{xx},
$$
and that's one independent variable and second derivatives. With $x$ and $t$ and third-order equations, the bookkeeping is most of the work.

## Symmetries as deformations of solutions

Here's the approach I learned from Stephen Anco. Instead of tracking how a symmetry moves $x$ and $u$ separately, track how it moves *solutions*. Under $X$, a point $(x, u)$ on the graph of a solution moves by $\varepsilon(\xi, \eta)$. Seen at a fixed $x$, the graph moves up by $\varepsilon\eta$, and sliding it sideways by $\varepsilon\xi$ looks like a vertical change of $-\varepsilon\xi u_x$. So to first order the solution changes by $\varepsilon Q$, where
$$
Q = \eta - \xi\,u_x
$$
is called the *characteristic* of $X$. For the rotation, $Q = x + uu_x$. On the line $u = ax + b$ this is $Q = (1 + a^2)x + ab$, so $u + \varepsilon Q$ has slope $a + \varepsilon(1+a^2)$ and intercept $b + \varepsilon ab$: exactly the derivatives at $\theta = 0$ of the slope and intercept of $\tilde f$ in the previous section.

The reason $Q$ is the right object is the identity
$$
\operatorname{pr}X = \operatorname{pr}X_Q + \xi\,D_x,
\qquad
\operatorname{pr}X_Q = Q\,\partial_u + (D_xQ)\,\partial_{u_x} + (D_x^2Q)\,\partial_{u_{xx}} + \cdots
$$
(compare coefficients to check it). The prolongation of $X_Q = Q\,\partial_u$ has no formula to remember: we just keep differentiating $Q$. And the extra piece $\xi D_x$ never matters, because $D_x\Delta = 0$ on solutions (if $\Delta$ vanishes along a solution, so does its derivative along the solution). So invariance only depends on $\operatorname{pr}X_Q(\Delta)$, which has a familiar name:
$$
\operatorname{pr}X_Q(\Delta) = \frac{d}{d\varepsilon}\Big|_{\varepsilon=0}\Delta[u + \varepsilon Q] = \Delta'(Q),
$$
the *linearization* (Fréchet derivative) of $\Delta$ in the direction $Q$.

::: definition
A *symmetry* of $\Delta = 0$ is a function $Q$ on jet space such that
$$
\Delta'(Q) = 0 \quad\text{whenever}\quad \Delta = 0.
$$
:::

In words, $u + \varepsilon Q[u]$ is a solution to first order in $\varepsilon$ whenever $u$ is: $Q$ is a tangent vector to the space of solutions. That's the coordinate-free picture. A symmetry is a direction in which we can move through the space of solutions, and the definition says nothing about how $x$, $t$ or $u$ move individually. Nothing forces $Q$ to have the form $\eta - \xi u_x$ either: letting $Q$ depend on higher derivatives gives *generalized* symmetries for free, which is how equations like KdV end up with infinitely many. Point symmetries are the special case
$$
Q = \eta(x,t,u) - \xi(x,t,u)\,u_x - \tau(x,t,u)\,u_t
$$
when there are two independent variables $x$ and $t$.

Finding symmetries is now a three-step recipe:

1. linearize: compute $\Delta'(Q)$;
2. use the equation to eliminate a derivative (say $u_t$) everywhere;
3. split: the result is a polynomial in the remaining derivatives of $u$, so every coefficient must vanish.

::: example Rotations of lines, again
For $\Delta = u_{xx}$ the linearization is $\Delta'(Q) = D_x^2 Q$. With $Q = x + uu_x$,
$$
D_x^2 Q = D_x\big(1 + u_x^2 + uu_{xx}\big) = 3u_xu_{xx} + uu_{xxx},
$$
which vanishes whenever $u_{xx} = 0$ (and hence $u_{xxx} = 0$). Same answer as before, with no prolongation formula.
:::

::: example A Galilean boost of Burgers' equation
Burgers' equation $\Delta = u_t + uu_x - u_{xx} = 0$ is a model for weak shock waves in gases and for traffic flow. Its linearization is
$$
\Delta'(Q) = D_tQ + D_x(uQ) - D_x^2Q.
$$
Switching to a frame moving at speed $\varepsilon$, $u(x,t) \mapsto u(x - \varepsilon t, t) + \varepsilon$, is generated by $X = t\,\partial_x + \partial_u$, with characteristic $Q = 1 - tu_x$. Then
$$
\Delta'(1 - tu_x) = -t\,\big(u_{xt} + u_x^2 + uu_{xx} - u_{xxx}\big) = -t\,D_x\Delta,
$$
which is zero on solutions but not zero identically. The words "whenever $\Delta = 0$" in the definition are doing real work.
:::

### Worked example: the heat equation

Let's find every point symmetry of the heat equation $u_t = u_{xx}$. (In physical units it's $u_t = \alpha u_{xx}$; rescaling $t$ sets $\alpha = 1$.) The equation is linear, so its linearization is the same operator:
$$
\Delta'(Q) = D_tQ - D_x^2Q.
$$
**Step 2.** On solutions we can replace $u_t$ by $u_{xx}$ everywhere, so we may take $Q = \eta - \xi u_x - \tau u_{xx}$, with $\xi, \tau, \eta$ functions of $(x, t, u)$, and $D_t$ acts on functions of $(x, t, u, u_x, u_{xx}, \dots)$ as
$$
D_t = \partial_t + u_{xx}\,\partial_u + u_{xxx}\,\partial_{u_x} + u_{xxxx}\,\partial_{u_{xx}} + \cdots.
$$
**Step 3.** Expanding $D_tQ - D_x^2Q$ gives a polynomial in $u_x, u_{xx}, u_{xxx}, u_{xxxx}$ whose coefficients are functions of $(x,t,u)$. The $u_{xx}^2$ and $u_{xxxx}$ terms cancel, and the remaining coefficients must vanish:

| term | coefficient |
|------|-------------|
| $u_xu_{xxx}$ | $2\tau_u$ |
| $u_{xxx}$ | $2\tau_x$ |
| $u_x^2u_{xx}$ | $\tau_{uu}$ |
| $u_xu_{xx}$ | $2\xi_u + 2\tau_{xu}$ |
| $u_{xx}$ | $\tau_{xx} + 2\xi_x - \tau_t$ |
| $u_x^3$ | $\xi_{uu}$ |
| $u_x^2$ | $2\xi_{xu} - \eta_{uu}$ |
| $u_x$ | $\xi_{xx} - \xi_t - 2\eta_{xu}$ |
| $1$ | $\eta_t - \eta_{xx}$ |

These are the *determining equations*, and they're solved from the top down. The first two give $\tau = \tau(t)$, then $\xi_u = 0$, and the $u_{xx}$ row gives $\xi = \frac12\tau'(t)\,x + \sigma(t)$. The $u_x^2$ row gives $\eta_{uu} = 0$, so $\eta = f(x,t)\,u + g(x,t)$, and the $u_x$ row gives
$$
f = -\tfrac18\tau''x^2 - \tfrac12\sigma'x + \rho(t).
$$
Finally the last row, split into the parts with and without $u$, says $g_t = g_{xx}$ and $f_t = f_{xx}$. The second of these forces $\tau''' = 0$, $\sigma'' = 0$ and $\rho' = -\frac14\tau''$. Writing $\tau = c_1 + 2c_2t + 4c_3t^2$, $\sigma = c_4 + 2c_5t$ and $\rho = c_6 - 2c_3t$, each constant gives one symmetry:

| | generator $X$ | characteristic $Q$ | meaning |
|---|---|---|---|
| $c_1$ | $\partial_t$ | $-u_t$ | time translation |
| $c_4$ | $\partial_x$ | $-u_x$ | space translation |
| $c_6$ | $u\,\partial_u$ | $u$ | scaling $u$ |
| $c_2$ | $x\,\partial_x + 2t\,\partial_t$ | $-xu_x - 2tu_t$ | scaling with $x \sim \sqrt t$ |
| $c_5$ | $2t\,\partial_x - xu\,\partial_u$ | $-2tu_x - xu$ | Galilean boost |
| $c_3$ | $4tx\,\partial_x + 4t^2\partial_t - (x^2 + 2t)u\,\partial_u$ | $-4txu_x - 4t^2u_t - (x^2+2t)u$ | projective |
| $g$ | $g(x,t)\,\partial_u$ | $g$ | superposition |

where $g$ is any solution of the heat equation. That's the full point symmetry algebra of the heat equation, from expanding one expression and reading off coefficients. Packages like GeM for Maple automate exactly these three steps.

## Adjoint symmetries and conservation laws

Symmetries are solutions of the linearized equation. *Adjoint symmetries* are solutions of its adjoint.

The formal adjoint of a linear differential operator $L$ is the operator $L^*$ with
$$
Q\,L(P) - P\,L^*(Q) = D_tA + D_xB
$$
for some $A$ and $B$: integrate by parts until every derivative is off $P$. For example $D_t^* = -D_t$, $(D_x^2)^* = D_x^2$, and $(D_x \circ f)^* = -f\,D_x$.

::: definition
An *adjoint symmetry* of $\Delta = 0$ is a function $Q$ on jet space such that
$$
\Delta'^*(Q) = 0 \quad\text{whenever}\quad \Delta = 0.
$$
:::

Why would we care? Because of conservation laws. A *conservation law* is a pair $(T, \Phi)$ with
$$
D_tT + D_x\Phi = 0 \quad\text{on solutions}.
$$
Integrating over $x$ shows that $\int T\,dx$ can only change through the flux $\Phi$ at the boundary, so if things decay at infinity, $\int T\,dx$ is constant. The way to find conservation laws by hand is to multiply the equation by something clever and try to write the result as a divergence:
$$
Q\,\Delta = D_tT + D_x\Phi \quad\text{identically}.
$$
Such a $Q$ is called a *multiplier*. Deciding whether an expression is a divergence is exactly what the Euler operator (the variational derivative)
$$
E_u = \partial_u - D_t\,\partial_{u_t} - D_x\,\partial_{u_x} + D_x^2\,\partial_{u_{xx}} - \cdots
$$
is for: $E_u(F) = 0$ identically if and only if $F$ is a total divergence. (It's the same operator that produces Euler–Lagrange equations, and a Lagrangian that is a divergence gives no equations of motion.) So multipliers are the solutions of
$$
E_u(Q\Delta) = 0,
$$
another linear equation that we can solve by expanding and splitting. This is Anco and Bluman's *direct method*, and unlike Noether's theorem it doesn't need a Lagrangian.

Here's where adjoint symmetries come in. A product rule for the Euler operator says
$$
E_u(Q\Delta) = \Delta'^*(Q) + Q'^*(\Delta),
$$
and the last term vanishes on solutions. So every multiplier is an adjoint symmetry: conservation laws live inside the solutions of the adjoint linearized equation, just as symmetries are solutions of the linearized one. When $\Delta$ is an Euler–Lagrange equation, $\Delta'$ is self-adjoint, adjoint symmetries and symmetries are the same thing, and we recover Noether's theorem. The heat and Burgers equations aren't Euler–Lagrange equations, and adjoint symmetries take over the job symmetries do in Noether's theorem.

::: example Conservation laws of the heat equation
Here $\Delta'^*(Q) = -D_tQ - D_x^2Q$. For a multiplier $Q(x,t)$ the condition $E_u\big(Q(u_t - u_{xx})\big) = -Q_t - Q_{xx} = 0$ says that $Q$ solves the backward heat equation. Its simplest solutions give:

- $Q = 1$: $\;\Delta = D_t(u) + D_x(-u_x)$, so the total heat $\int u\,dx$ is conserved.
- $Q = x$: $\;x\Delta = D_t(xu) + D_x(u - xu_x)$, so $\int xu\,dx$ is conserved: the centre of mass of the heat doesn't move.
- $Q = x^2 - 2t$: $\;(x^2 - 2t)\Delta = D_t\big((x^2 - 2t)u\big) + D_x\big(2xu - (x^2 - 2t)u_x\big)$, so $\int x^2u\,dx = 2t\int u\,dx + \text{const}$.

The last one says the variance of the temperature distribution grows like $2\alpha t$ in physical units: heat spreads a distance of order $\sqrt{\alpha t}$. We learned that without solving the heat equation.
:::

::: example Burgers' equation has one conservation law
For Burgers' equation $\Delta'^*(Q) = -D_tQ - uD_xQ - D_x^2Q$. Looking for multipliers $Q(x,t,u)$, the direct method gives
$$
E_u(Q\Delta) = -2Q_u\,u_{xx} - Q_{uu}\,u_x^2 - 2Q_{xu}\,u_x - Q_t - uQ_x - Q_{xx}.
$$
The $u_{xx}$ coefficient forces $Q_u = 0$, then the term proportional to $u$ forces $Q_x = 0$, and then $Q_t = 0$. So $Q = 1$ is the only multiplier of this kind, and it gives
$$
D_tu + D_x\big(\tfrac12u^2 - u_x\big) = 0:
$$
conservation of $\int u\,dx$.
:::

## What symmetries are good for

From an engineer's point of view, a symmetry group is a way of getting information out of a PDE without solving it. These are the uses I reach for most.

### Similarity solutions and scaling laws

A solution is invariant under a symmetry exactly when it doesn't move, that is when $Q[u] = 0$. That's an extra equation, and for two independent variables it turns the PDE into an ODE.

Take the scaling symmetry of the heat equation, $Q = -xu_x - 2tu_t$. Solutions of $Q = 0$ have the form $u = F(\eta)$ with $\eta = x/\sqrt t$, and substituting into $u_t = u_{xx}$ gives the ODE $F'' + \frac{\eta}{2}F' = 0$, so $F = A + B\operatorname{erf}(\eta/2)$. In physical units, this is a solid at temperature $T_i$ whose surface is suddenly held at $T_s$:
$$
T(x,t) = T_s + (T_i - T_s)\operatorname{erf}\!\left(\frac{x}{2\sqrt{\alpha t}}\right).
$$
Even before finding $F$, the symmetry tells us that the temperature depends on $x$ and $t$ only through $x/\sqrt{\alpha t}$. Heat penetrates a depth proportional to $\sqrt{\alpha t}$, so it takes four times as long to get twice as deep: for concrete ($\alpha \approx 10^{-6}\,\mathrm{m^2/s}$) that's roughly three hours for 10 cm and twelve hours for 20 cm. The heat flux through the surface falls off like $1/\sqrt t$:
$$
q(t) = \frac{k\,(T_s - T_i)}{\sqrt{\pi\alpha t}}.
$$

<figure>
<div class="collapse-fig"></div>
<figcaption>Temperature profiles $\theta = (T - T_i)/(T_s - T_i)$ in concrete at different times (left) all fall on the single curve $1 - \operatorname{erf}$ when plotted against $x/2\sqrt{\alpha t}$ (right). This "data collapse" is what a scaling symmetry looks like in an experiment.</figcaption>
</figure>

<script src="pde-symmetries/collapse.js"></script>

The same reduction, applied to the boundary layer equations, is how Blasius turned the flow over a flat plate into an ODE, and it's why a laminar boundary layer grows like $\sqrt{\nu x/U}$. One caveat: the boundary conditions have to be invariant too. A semi-infinite solid has no length scale, so the scaling symmetry survives. A slab of thickness $L$ doesn't, and the similarity solution is only good while $\sqrt{\alpha t} \ll L$.

Conservation laws tell us which scaling to use. Any combination $x\,\partial_x + 2t\,\partial_t + k\,u\,\partial_u$ is a symmetry, with invariant solutions $u = t^{k/2}F(x/\sqrt t)$. Total heat is conserved, and $\int t^{k/2}F(x/\sqrt t)\,dx = t^{(k+1)/2}\int F(\eta)\,d\eta$ is constant only for $k = -1$. That choice gives the ODE $F'' + \frac12(\eta F)' = 0$, so $F = Ae^{-\eta^2/4}$ and
$$
u = \frac{A}{\sqrt t}\,e^{-x^2/4t},
$$
the temperature after a burst of heat at a point. Symmetries give the family of candidate similarity solutions, and conservation laws pick out the physical one.

### Travelling waves and shocks

With viscosity $\nu$, Burgers' equation is $u_t + uu_x = \nu u_{xx}$. The combined translation $\partial_t + c\,\partial_x$ has $Q = -u_t - cu_x$, and its invariant solutions are travelling waves $u = F(x - ct)$. Solving the resulting ODE for a front going from $u_L$ down to $u_R$ gives
$$
u = \frac{u_L + u_R}{2} - \frac{u_L - u_R}{2}\tanh\!\left(\frac{(u_L - u_R)(x - ct)}{4\nu}\right),
\qquad
c = \frac{u_L + u_R}{2}.
$$
Two engineering facts fall out. The shock moves at the average of the two states. (That also follows from the conservation law alone: integrating $D_tu + D_x(\frac12u^2 - \nu u_x) = 0$ across the shock gives the Rankine–Hugoniot condition.) And the shock thickness is about $4\nu/(u_L - u_R)$, so stronger shocks are thinner. In traffic models, $u$ is an affine function of the density of cars and the shock is the back of a jam.

### New solutions from old

Symmetries map solutions to solutions. The Galilean boost turns a standing shock into a moving one: if $u(x,t)$ solves Burgers' equation, so does $u(x - \varepsilon t, t) + \varepsilon$. For the heat equation, the projective symmetry turns the boring solution $u = 1$ into
$$
u = \frac{1}{\sqrt{1 + 4\varepsilon t}}\exp\!\left(-\frac{\varepsilon x^2}{1 + 4\varepsilon t}\right),
$$
which, after shifting $t$, is the point-source solution again. The practical use is test problems: every symmetry turns one exact solution into a family of them, and a solver that passes one case but fails a boosted or rescaled version of it is mishandling its frame or its units.

### Conservation laws as checks

Conservation laws give integral checks on simulations: a diffusion code that loses $\int u\,dx$, or whose variance doesn't grow like $2\alpha t$, has a bug. They also decide which weak form to discretize. A finite volume scheme is built from the flux $\Phi$ of a conservation law, and using the wrong one gives shocks that move at the wrong speed.

### Reduction in mechanics and control

For ODEs the same machinery reduces dimension. A rotational symmetry of a mechanical system gives a conserved angular momentum, and fixing it lets us drop a coordinate (Routh reduction, or symplectic reduction in the geometric language). This is where symmetry analysis of PDEs meets the control problems I work on now, which I'll leave for another post.

## Further reading

- P. J. Olver, *Applications of Lie Groups to Differential Equations*. The standard reference for prolongation, characteristics, the Euler operator and Noether's theorem.
- G. W. Bluman, A. F. Cheviakov and S. C. Anco, *Applications of Symmetry Methods to Partial Differential Equations*. Symmetries, adjoint symmetries and conservation laws through multipliers, with many examples.
- S. C. Anco and G. Bluman, "Direct construction method for conservation laws of partial differential equations", Parts I and II, *European Journal of Applied Mathematics* 13 (2002).
