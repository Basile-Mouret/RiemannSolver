#set page(margin:2.5cm)

#set text(font:"New Computer Modern", size:12pt)

#set par(justify: true)
#set text(hyphenate: false)

#set heading(numbering: (..nums) => {
  let n = nums.pos()
  if n.len() == 1 {
    [Chapter #n.first() - ]
  } else {
    numbering("1.1", ..nums)
  }
})


#v(2fr)


#block[
#align(center)[
  #text("Notes on the Book", size: 2em, weight: "bold")\
  #text("Riemann Solvers and Numerical Methods for Fluid Dynamics by Toro", size: 2em, style: "italic")
]]

#v(3fr)



== Notations

- $ u = u(x_1, x_2, ..., x_m, t)$

$u_t$ : Derivative of u by t

$nabla u = ((partial u)/(partial x_i))_i$ : Gradient of u

$nabla² u = ((partial^2 u)/(partial x_i partial x_j))_(i,j)$ :  Hessian of u

$Delta u = sum_i (partial^2 u)/(partial x_i^2) $ Laplacian of u

- U is a vector $U = (u_i)_i$ with $ u_i = u_i (x_1, x_2, ..., x_m, t)$

$nabla U = ((partial u_i)/(partial x_j))_(i,j) = (nabla u_i^T)_i$ : Jacobian of U is a matrix

#include "chapters/1-Equations_of_Fluid_Dynamics.typ"
#include "chapters/2-Notions_on_Hyperbolic_PDEs.typ"
