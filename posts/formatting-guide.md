---
title: How to write a post
date: 2026-11-26
summary: Everything a post can contain, with the Markdown that produces it.
tags: [meta]
toc: true
draft: true
---

A post is one Markdown file in the `posts/` folder. Math is ordinary LaTeX. This page is a
**draft**, so it only appears in the local preview, never on the live site. Read it next to its
source, `posts/formatting-guide.md`, to see how everything is written.

## The top of the file

Every post starts with a short block of settings, called front matter:

```yaml
---
title: Chow's Theorem, Physics, and Control
date: 2026-06-18
summary: One or two sentences for the list on the home page.
status: incomplete
---
```

- `title` and `date` are the important ones. Write the date as `2026-06-18`, or just `2026-06`
  for the month. Newest posts are listed first.
- `summary` is optional and shows under the title on the home page. Titles and summaries can
  contain math, e.g. `title: Representations of $SU(2)$`.
- `status: incomplete` (or any short word) puts a small label next to the title.
- `author: V. I. Arnold` is for reposts of someone else's writing.
- `tags: [geometry, control]` is optional and shows on the post's page.
- `toc: true` adds a table of contents built from your `##` and `###` headings, like the one above.
- `draft: true` keeps the post off the live site. You can still preview it locally.

The file name becomes the link: `posts/chows-theorem.md` lives at `post.html?p=chows-theorem`.

## Text

Paragraphs are separated by a blank line. Single line breaks inside a paragraph are ignored,
so you can wrap lines however you like.

```markdown
*italic*, **bold**, `code`, [a link](https://en.wikipedia.org/wiki/Lagrangian_mechanics),
a link to another post: [Chow's theorem](chows-theorem.md), and a footnote.[^lang]

[^lang]: Footnote text can go anywhere in the file. It is collected at the bottom.
```

*italic*, **bold**, `code`, [a link](https://en.wikipedia.org/wiki/Lagrangian_mechanics),
a link to another post: [Chow's theorem](chows-theorem.md), and a footnote.[^lang]

[^lang]: Footnote text can go anywhere in the file. It is collected at the bottom.

Use `##` for sections and `###` for subsections. The post's title is already the top
heading, so don't use a single `#`.

### A subsection

Lists, quotes and horizontal rules work as usual:

1. First, a numbered item with math: $G/\ker\varphi \cong \im \varphi$.
2. Second, a numbered item.
   - Nested bullets are indented by three spaces.

> The book of nature is written in the language of mathematics.

---

## Math

Inline math goes between single dollars: `$\dot q = \partial H / \partial p$` gives
$\dot q = \partial H / \partial p$. Display math goes between double dollars:

```latex
$$
\frac{d}{dt}\frac{\partial L}{\partial \dot q} = \frac{\partial L}{\partial q}
$$
```

$$
\frac{d}{dt}\frac{\partial L}{\partial \dot q} = \frac{\partial L}{\partial q}
$$

`\( ... \)` and `\[ ... \]` work as well. Underscores, asterisks and backslashes inside math
are left alone, so write LaTeX exactly as you would in a `.tex` file:
$a_1 * b_1 = \{x \in \R : \norm{x} < 1\}$. For a literal dollar sign, write `\$`:
the book cost \$40.

### Numbered equations and references

Environments like `equation` and `align` number themselves. Give an equation a `\label` and
refer to it with `\eqref` anywhere in the text:

```latex
\begin{equation}\label{eq:action}
  S[q] = \int_{t_1}^{t_2} L(q, \dot q, t)\, dt
\end{equation}

Requiring $\delta S = 0$ in \eqref{eq:action} gives
\begin{align}
  \delta S &= \int_{t_1}^{t_2} \left( \frac{\partial L}{\partial q}\delta q
              + \frac{\partial L}{\partial \dot q}\delta \dot q \right) dt \label{eq:var} \\
           &= \int_{t_1}^{t_2} \left( \frac{\partial L}{\partial q}
              - \frac{d}{dt}\frac{\partial L}{\partial \dot q} \right) \delta q \, dt. \nonumber
\end{align}
```

\begin{equation}\label{eq:action}
  S[q] = \int_{t_1}^{t_2} L(q, \dot q, t)\, dt
\end{equation}

Requiring $\delta S = 0$ in \eqref{eq:action} gives
\begin{align}
  \delta S &= \int_{t_1}^{t_2} \left( \frac{\partial L}{\partial q}\delta q
              + \frac{\partial L}{\partial \dot q}\delta \dot q \right) dt \label{eq:var} \\
           &= \int_{t_1}^{t_2} \left( \frac{\partial L}{\partial q}
              - \frac{d}{dt}\frac{\partial L}{\partial \dot q} \right) \delta q \, dt. \nonumber
\end{align}

Use the starred versions (`align*`, `equation*`) or `$$ ... $$` for unnumbered displays.
Matrices, cases and the other AMS environments all work:

$$
\rho(\sigma) = \begin{pmatrix} 0 & 1 \\ 1 & 0 \end{pmatrix},
\qquad
\chi(g) = \begin{cases} 2 & g = e, \\ -1 & g \text{ a 3-cycle}, \\ 0 & \text{otherwise}. \end{cases}
$$

So do commutative diagrams, with the `CD` environment:

$$
\begin{CD}
A @>{f}>> B \\
@V{g}VV @VV{h}V \\
C @>>{k}> D
\end{CD}
$$

### Shortcuts

These are defined for every post: $\R, \C, \Q, \Z, \N, \F$, $\eps$, $\abs{x}$, $\norm{v}$,
$\inner{u}{v}$, $\Hom$, $\End$, $\Aut$, $\Gal$, $\tr$, $\im$, $\id$. Add more in the
`MATH_MACROS` list at the top of `assets/blog.js`, or define them for one post with
`\newcommand` on a line of its own (it doesn't show up on the page):

```latex
$\newcommand{\Rep}{\operatorname{Rep}}$
```

$\newcommand{\Rep}{\operatorname{Rep}}$

Now $\Rep(G)$ works for the rest of this post.

## Theorems and proofs

Fence a block with three colons and name the environment:

```markdown
::: definition
A *representation* of a group $G$ on a vector space $V$ is a homomorphism
$\rho \colon G \to GL(V)$.
:::

::: theorem Maschke
If $G$ is finite and $\operatorname{char} \F \nmid |G|$, every representation of $G$
over $\F$ is a direct sum of irreducible ones.
:::

::: proof
Average any projection over $G$:
$$
\pi_0 = \frac{1}{|G|} \sum_{g \in G} \rho(g)\, \pi\, \rho(g)^{-1}.
$$
:::
```

::: definition
A *representation* of a group $G$ on a vector space $V$ is a homomorphism
$\rho \colon G \to GL(V)$.
:::

::: theorem Maschke
If $G$ is finite and $\operatorname{char} \F \nmid |G|$, every representation of $G$
over $\F$ is a direct sum of irreducible ones.
:::

::: proof
Average any projection over $G$:
$$
\pi_0 = \frac{1}{|G|} \sum_{g \in G} \rho(g)\, \pi\, \rho(g)^{-1}.
$$
:::

The available blocks are `theorem`, `lemma`, `proposition`, `corollary`, `conjecture`,
`definition`, `example`, `exercise` and `problem` (numbered, sharing one counter), plus
`claim`, `remark`, `note`, `proof` and `solution` (unnumbered). Short names work too:
`thm`, `lem`, `prop`, `cor`, `def`, `ex`, `rem`.

::: lemma Schur
A $G$-linear map between irreducible representations is either $0$ or an isomorphism.
:::

::: proof of Lemma 3
Its kernel and image are subrepresentations.
:::

::: remark
To put a block inside another (say a proof inside an example), give the outer one
more colons: open it with `::::` and close it with `::::`.
:::

## Code, tables, images

Code blocks get syntax highlighting when you name the language:

```python
import numpy as np

def euler_lagrange_residual(q, t, dt):
    """Discrete check of d/dt (dL/dq') - dL/dq = 0 for L = q'^2/2 - q^2/2."""
    qdot = np.gradient(q, dt)
    return np.gradient(qdot, dt) + q  # $q$ here is just a comment
```

| Group   | Order | Irreducible representations |
|---------|------:|-----------------------------|
| $S_3$   | 6     | 3                           |
| $Q_8$   | 8     | 5                           |
| $A_5$   | 60    | 5                           |

Images and videos: put the files in a folder next to the post, e.g. `posts/my-post/plot.png`,
and write `![description](my-post/plot.png "Caption")`. The caption is optional and can contain
math. For a video, use HTML:

```html
<figure>
<video controls src="my-post/animation.mp4"></video>
<figcaption>What the video shows</figcaption>
</figure>
```

<details>
<summary>Raw HTML works too, e.g. a fold-out section</summary>

Hidden until clicked, and it can hold math: $e^{i\pi} + 1 = 0$.

</details>
