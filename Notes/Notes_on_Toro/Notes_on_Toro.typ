#set page(margin:2.5cm)

#set text(font:"New Computer Modern", size:12pt)

#set par(justify: true)
#set text(hyphenate: false)

#show figure: set block(
  spacing: 2em,
)
// outline/table of content
#import "@preview/outrageous:0.4.0"
#show outline.entry: outrageous.show-entry.with(
  font: ("New Computer Modern",) 
)


#set heading(numbering: "1.1")

#show heading: it => {
  if it.level == 1 {v(1em)} else {v(0.5em)}
  it
  if it.level == 1 {v(1em)} else {v(0.5em)}
}

// #set math.equation(numbering: "(1.)")
// #set math.equation(supplement: none, numbering: it => {numbering("(1.1)", counter(heading).get().first(), it)})

//#set math.equation(numbering:none)
#show math.equation: set text(font: "New Computer Modern Math")



#page[#include "chapters/0-Cover.typ"]

#page[
#align(center)[#heading(numbering: none, outlined : false)[Abstract]]
These are my personal notes on the book from Toro during my internship at the Cagire Team.

#align(center+horizon)[#outline(
title: "Content Table",
depth: 2,
)
]]


#set page(numbering: "1")
#counter(page).update(1)
#counter(heading).update(0)
#show link: set text(fill: blue)

== Notations

- $u = u(x_1, x_2, ..., x_m, t)$

$u_t$ : Derivative of u by t

$nabla u = ((partial u)/(partial x_i))_i$ : Gradient of u

$nabla² u = ((partial^2 u)/(partial x_i partial x_j))_(i,j)$ :  Hessian of u

$Delta u = sum_i (partial^2 u)/(partial x_i^2) $ Laplacian of u

- U is a vector $U = (u_i)_i$ with $ u_i = u_i (x_1, x_2, ..., x_m, t)$

$nabla U = ((partial u_i)/(partial x_j))_(i,j) = (nabla u_i^T)_i$ : Jacobian of U is a matrix

In a conservation law $U$ is the conserved variable and $F$ is the flux.

#pagebreak()
#include "chapters/1-Equations_of_Fluid_Dynamics.typ"
#pagebreak()
#include "chapters/2-Notions_on_Hyperbolic_PDEs.typ"
#pagebreak()
#include "chapters/3-Some_Properties_of_the_Euler_Equations.typ"
#pagebreak()
#include "chapters/4-The_Riemann_Problem_for_the_Euler_Equations.typ"
