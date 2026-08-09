#import "../packages/theorems/config.typ" : *
#show: thmrules.with(qed-symbol: $square$)

// ploting
#import "@preview/lilaq:0.6.0" as lq
#import "@preview/tiptoe:0.4.0"

// == Notations
//
// scalar : lower case $s$
//
// vector : bold lower case $bold(v)$
//
// bold upper case : matrix ($bold(M)$)
//
// $nabla =  vec(partial/(partial x), partial/(partial y), partial/(partial z))$
//
// $"grad"(u)=nabla u =  vec((partial u)/(partial x), (partial u)/(partial y), (partial u)/(partial z))$
//
// $"Jac"(bold(u)) = nabla bold(u) = vec(nabla u_1^T, dots, nabla u_n^T)$
//
// $"div"(u) = nabla dot u =  (partial u)/(partial x) + (partial u)/(partial y) + (partial u)/(partial z)$
//
// $"rot"(u) = nabla times u $
//
// total derivative in time

= The finite volume method for fluid dynamics

In this first part I will detail the theoretical part of this internship.
I will present the equation governing fluid motion, then I will look into their specificities and present numerical methods used for solving them.

==  The Euler equations <euler_equations>

Describing the dynamics of fluids by each particle composing the medium is too difficult.
Instead physicists use macroscopic quantities that result from taking the mean of microscopic values over many particles.

For compressible fluids these macroscopic quantities are :
- the *density* $rho (bold(x),t)$ in $"kg"dot "m"^(-3)$ resulting from taking the mean mass per volume.

- the *velocity vector* $bold(v) (bold(x),t) = mat(u, v, w)^T$ in $"m" dot "s"^(-1)$ resulting from taking the mean velocity of each particle per volume.

- the *pressure* $p(bold(x), t)$ in $"Pa" = "N" dot "m"^(-2)$ resulting from the mean magnitude of force per area.

These are called the _primitive_ variables, but in order to compute them we use the _conserved variables_ : 

- the *density* $rho (bold(x),t)$
- the *momentum vector* $bold(m)  (bold(x),t) = rho bold(v)$
- the *total energy per unit volume* $E(bold(x), t)$

Using the fundamental laws of conservation of mass, Newton second law and the conservation of Energy we get the Euler equations governing inviscid compressible fluids on the _conserved variables_ : 
$
cases(
  rho_t + bold(nabla) dot (rho bold(v)) = 0,
  (rho bold(v))_t + nabla dot (rho bold(v) times.o bold(v) + p I) = 0,
  E_t + nabla dot [(E+p) bold(v)] = 0
)
$ 

where $E = rho(1/2 bold(v) dot bold(v) + e)$, with $e$ the *specific internal energy* that is determined using thermodynamical laws (EOS). 
We will study ideal gases, for which we have $e(rho, p) = p/((gamma-1)rho)$.

The system can be rewritten in conservative form : 
$
bold(u)_t + bold(F)(bold(u))_x + bold(G)(bold(u))_y + bold(H)(bold(u))_z = 0
$<Euler_cons>
with $bold(u) = vec(rho, rho u, rho v, rho w, E), quad bold(F)(bold(u)) = vec(rho u, rho u^2 + p, rho u v, rho u w, (E + p) u), quad bold(G)(bold(u)) = vec(rho v, rho u v, rho v^2 + p, rho v w, (E + p) v) "and" bold(H)(bold(u)) = vec(rho w, rho u w, rho v w, rho w^2 + p, (E + p) w)$

$bold(u)$ is the column vector of conserved variables and $bold(F), bold(G), bold(H)$ are the _flux vectors_ in the $x, y, z$ directions respectively.

This is a system of differential equations and thus assumes smooth solutions (partial derivatives exist). In order to handle discontinuous solutions we will need to rewrite it in integral form.
=== Submodels <submodels>

There exists a lot of derived systems obtained from simplifications on the Euler equations.
These are useful to analyze each of aspect of the full system individually as it has a lot of specificities.

==== Linear Waves and Linear Advection

Let us consider a linearization of the euler equation around a state at rest $bold(v) = vec(rho_0, 0, p_0)$ in which we introduce a small disturbance : $bold(v') = vec(rho', u', p')$.

When simplifying derivatives of constants and higher order terms, the one dimensional Euler equation becomes 
$
cases(
rho'_t + rho_0 u'_x = 0,
u'_t + (1/rho_0) p'_x = 0,
p'_t + rho_0 a_0^2 u'_x = 0
)
$
with $a_0 = sqrt((gamma p_0)/rho_0)$ the speed of sound.

As $rho'$ only appear in the first equation, it is fully determined by integrating $rho_0 u'_x$ over time.

The remaining $2 times 2$ system is called the wave system. It's specificity is that it combines two underlying linear waves travelling in opposite direction. This drops the non linearity in the full euler system, as the coefficient $(1/rho_0)$ and $rho_0 a_0^2$ are constants independent from $bold(v')$. In order to solve it, it is useful to decompose it into two independent scalar equations.

==== Burgers' Equation

The burgers' equation on the other hand is a scalar non linear equation : 
$
u_t + u u_x = 0
$<burgers>
as the coefficient in front of $u_x$ depends on the state $u$.

This equation is the simplest non-linear conservation law, and as such is useful to study non-linearity.

==== Isentropic Euler Equations

Finally, a more complete system that is useful before studying the full Euler system are obtained when assuming a constant entropy, $dif S = 0$.
This leads to pressure being fully described by the density field, $p(rho) = kappa rho^gamma$, with $kappa$ and $gamma$ two constants. We can thus drop the last equation of the euler system.

For $gamma=2$ and $kappa=g/2$ we obtain the same equations as in the Shallow Water Equations.







== Hyperbolicity

A key aspect of all of these systems is hyperbolicity.
For smooth solutions, the chain rule lets us rewrite the one dimensional
conservative system @Euler_cons in _quasilinear form_
$
bold(u)_t + bold(A)(bold(u)) bold(u)_x = 0
$<quasilinear>
where $bold(A) = nabla_bold(u) F(bold(u))$ is the Jacobian matrix of the flux.

#definition[
  The system @quasilinear is *hyperbolic* if $bold(A)$ has real eigenvalues
  $lambda_1 <= dots.c <= lambda_m$ and a complete set of linearly independent
  right eigenvectors $bold(r)_1, dots, bold(r)_m$.
  It is *strictly hyperbolic* if in addition the eigenvalues are distinct.
]

Thus there exists an invertible matrix formed by the right eigenvectors $P(bold(u)) = (bold(r)_1, dots, bold(r)_m)$  and a diagonal matrix with the diagonal elements being the eigenvalues $D(bold(u)) = mat(lambda_1,,0 ;,dots.down,;0,, lambda_m)$ such that $A(bold(u)) = P(bold(u)) D(bold(u)) P(bold(u))^(-1)$.


This is especially useful for *linear* hyperbolic conservation laws as $P$ and $D$ do not depend on $bold(u)$.
In that case, we can use the change of variable
$
bold(v) = P^(-1) bold(u)
$
These are called _characteristic variables_ and using them, the quasilinear form can be rewritten as a system of $m$ independent advection equations : 
$
(partial v_i)/(partial t) + lambda_(i) (partial v_i)/(partial x) = 0
$<linear_decomp>

Where $lambda_(i)$ are the diagonal entries of $D$, i.e. the eigenvalues of $nabla_(bold(u)) F(bold(u))$ and are called the _characteristic speeds_.
In the more general non-linear case, this can't be done as the eigenvalues depend on the state so no such decoupling of the system exists.

=== Relation to the classification of second order equations

The _hyperbolic_ name comes from the classification of second order scalar equation into elliptic, parabolic and hyperbolic. 

A second order equation of the form $a u_(x x) + b u_(x y) + c u_(y y) = 0$  is classified using the discriminant $Delta = b^2 - 4 a c$.
If $Delta>0$ the equation is hyperbolic, if $Delta=0$ the equation is parabolic and if $Delta<0$ the equation is elliptic.

By setting $bold(w) = vec(u_x, u_y)$ and using the condition $(u_x)_t = (u_t)_x$, we get a first order system : 
$
bold(w)_x + A bold(w)_y = 0\ "with" A = mat(b/a, c/a; -1, 0)
$

The characteristic polynomial of $A$ is given by  
$lambda^2 - b/a lambda + c/a = 0$ and it has real roots when $(b^2 - 4 a c) / (a^2) > 0$, so the hyperbolicity condition is the same for both definitions.

The classification of the equations translate in their behaviour. Hyperbolic equations are wave like with a finite propagation speed. This isn't the case for parabolic equations whose diffusive behaviour affects the whole domain without any speed limit. Finally elliptic problems do not have any real characteristic direction and thus no timelike direction.

For example, steady state heat conductivity, modelized by Laplace's equation $- Delta u = 0$ is elliptic and adding a time dependency $u_t - Delta u = 0$ makes it parabolic.
Linear advection, $u_t + c u_x = 0$ on the other hand is Hyperbolic and so are the equations presented in @submodels.
The full Euler equation however are Hyperbolic only under certain condition and can become elliptic, for example in a steady subsonic compressible model.
Classification becomes even more complex when using the Navier Stokes equations that add a diffusive term.


=== Characteristic curves

#definition([characteristic curves])[
  Characteristics are curve in the $x-t$ plane where the PDE reduces to an ODE in time.
]

These curve are a very useful tool to visualize the wave like propagation of hyperbolic systems.

Let us first study the linear scalar advection $u_t + c u_x = 0$ and let us consider a curve $x = x(t)$. We can thus write $u$ as a function of t $u(t) = u(x(t), t)$. Differentiating $u$ over $t$ along the curve $x(t)$ gives :
$
(dif u)/(dif t) = (partial u)/(partial t) + (dif x)/(dif t) (partial u)/(partial x)
$
This is similar to the LHS of the PDE.
In the case where the curve satisfies $(dif x)/(dif t) = c$ we have $(dif u)/(dif t) = 0$.
This means that $u(t) = u_0$ is constant along the curve. Similarly to isolines for spatial visualization, the characteristics show us the path of constant values over time.
From $(dif x)/(dif t) = c$ we obtain $x(t) = c t + x_0$, on the $x-t$ plane these give lines with a slope $1/c$.



#let char-plot = it => {
  show: lq.set-diagram(
    ylabel: $t$,
    xlabel: $x$,
    xaxis: (subticks: none, ticks: none, tip: tiptoe.stealth, mirror: false, position: 0),
    yaxis: (subticks: none, ticks: none, tip: tiptoe.stealth, mirror: false, position: 0),
    xlim: (-1, 5),
    ylim: (-1, 3),
    aspect-ratio: 1,
    width: auto,
  )
  show: lq.set-label(angle: 0deg)
  show lq.selector(lq.label): set align(top + right)
  it
}

#figure(
  {
    show: char-plot
    lq.diagram(
      height: 5cm,
      lq.plot((0, 1), (2, 3), mark: none, color: black),
      lq.plot((0, 2), (1, 3), mark: none, color: black),
      lq.plot((0, 3), (0, 3), mark: none, color: black),
      lq.plot((1, 4), (0, 3), mark: none, color: black),
      lq.plot((2, 5), (0, 3), mark: none, color: black),
      lq.plot((3, 5), (0, 2), mark: none, color: black),
      lq.ellipse(2, 0, width: 0.2, height: 0.2, align: center + horizon),
      lq.place(2,-0.3)[$u(x_0,0)$],
      lq.ellipse(0, 1, width: 0.2, height: 0.2, align: center + horizon),
      lq.place(-0.6,1)[$u(0,t)$],
    )
  },
  caption: [Characteristic curves in the $(x,t)$ plane for linear advection.],
) <fig-characteristics>

The values of $u$ on the rays are given either by the initial condition or the boundary condition.

Linear hyperbolic systems can be decomposed into independent scalar advection equations on the characteristic variables as detailed in @linear_decomp. Each of these variables are then getting transported by their correspond speed given by the eigenvalues.


#figure(
  {
    show: char-plot
    lq.diagram(
      height: 5cm,
      lq.plot((0, 1), (2, 3), mark: none, color: red),
      lq.plot((0, 2), (1, 3), mark: none, color: red),
      lq.plot((0, 3), (0, 3), mark: none, color: red),
      lq.plot((1, 4), (0, 3), mark: none, color: red),
      lq.plot((2, 5), (0, 3), mark: none, color: red),
      lq.plot((3, 5), (0, 2), mark: none, color: red),

      lq.plot((1, 0), (0, 3), mark: none, color: blue),
      lq.plot((2, 1), (0, 3), mark: none, color: blue),
      lq.plot((3, 2), (0, 3), mark: none, color: blue),
      lq.plot((4, 3), (0, 3), mark: none, color: blue),
      lq.plot((5, 4), (0, 3), mark: none, color: blue),
    )
  },
  caption: [Characteristic curves for a $2 times 2$ hyperbolic system.],
)

This can be seen when simulating the wave system in 1D and watching the characteristic variables. Introducing a small difference in the initial condition, one can see one wave travelling to the left and the other to the right when viewing each characteristic variable. 

Let us now consider a scalar hyperbolic conservation law $u_t + f(u)_x = 0$. For $f$ differentiable, we can write it in quasilinear form $u_t + lambda(u) u_x = 0$ with $lambda(u) = f'(u)$. If we then consider characteristic curves satisfying 

$
(dif x)/(dif t) = lambda(u), quad x(0) = x_0
$

then the local derivative of $u(t) = u(x(t), t)$ along $x(t)$ gives
$
(dif u)/(dif t) = u_t + lambda(u)u_x = 0
$
i.e. $u$ is constant on this curve whose value is given by following the characteristic back in time, crossing either the initial value or a boundary condition. The slope can then be evaluated as $lambda(u_0)$ and the curve is given by $x(t) = x_0 + lambda(u_0(x_0)) t$. This means that two rays can have different slopes which leads to crossings, where a single point has multiple values and the model breaks at time $t_c$.


#figure(
  {
    show: char-plot
    lq.diagram(
      height: 5cm,
      lq.plot((1, 4), (0, 2), mark:none, color:black),
      lq.plot((2, 4), (0, 2), mark:none, color:black),
      lq.plot((3, 4), (0, 2), mark:none, color:black),
      lq.plot((0, 4), (2, 2), mark:none, color:gray, stroke:(dash:"dashed")),
      lq.place(-0.2,2)[$t_c$]

    )
  },
  caption: [Crossing of characteristic curves.],
)

In order to keep an inviscid model we have to allow discontinuities in the computed solutions.
This can be done by searching for weak solution using the integral form of the conservation law.

=== Weak solutions

Crossing characteristics force us to accept solutions that are not differentiable in the classical sense.
We only require $bold(u)$ to be locally integrable, this is possible thanks to distribution theory.

#definition([weak solution])[
  $bold(u)$ is a _weak solution_ of $bold(u)_t + bold(F)(bold(u))_x = 0$ if for every test function $phi in C_0^1 (RR times RR^+)$
  $
  integral_0^oo integral_RR (bold(u) phi_t + bold(F)(bold(u)) phi_x) dif x dif t + integral_RR bold(u)(x,0) phi(x,0) dif x = 0 .
  $
]

The weak formulation admits discontinuous data, which the classical theory
could not. The simplest such data is a single jump, and the resulting initial
value problem is central to everything that follows.

#definition([Riemann problem])[
  The _Riemann problem_ for a conservation law is the initial value problem
  $
  cases(
    bold(u)_t + bold(F)(bold(u))_x = 0,
    bold(u)(x,0) = cases(bold(u)_L "if" x<0, bold(u)_R "if" x>0)
  )
  $<riemann_problem>
  with two constant states $bold(u)_L$ and $bold(u)_R$.
]
This problem is the central building block of the Godunov method.

Now consider a discontinuity between two smooth states to travel at speed $S$ and evaluate the integral form on a control volume encompassing this discontinuity.
Taking its limit yields the Rankine-Hugoniot condition, also called the jump condition (see @toro2009Riemann, p. 70): 
$
bold(F)(bold(u)_R) - bold(F)(bold(u)_L) = S (bold(u)_R - bold(u)_L)
$<rankine-hugoniot>

=== Shocks and Rarefactions

When characteristics cross, we now approximate it by a discontinuity called a *shock*.
This happens when a characteristic on the left is faster than on the right, $lambda(u_l)>lambda(u_r)$ for a scalar law.
For fluids, a shock represents a rapid transition layer from one data state to another.
This layer is usually of the order of the mean free path of the molecules, and thus can be approximated by a discontinuity.
Determining the speed $S$ of this wave is solving the Riemann problem @riemann_problem. This can be solved using the Rankine-Hugoniot conditions @rankine-hugoniot.

#example[
  For Burgers' equation, the Rankine-Hugoniot condition gives $S = 1/2 (u_r + u_l)$, that is the shock travels at the mean velocity between the left and right states.
]

#figure(
  {
    show: char-plot
    lq.diagram(
      height: 5cm,
      lq.plot((0, 2), (0, 1), mark:none, color:black),
      lq.plot((1, 2), (0, 1), mark:none, color:black),
      lq.plot((2, 2), (0, 1), mark:none, color:black),
      lq.plot((2, 4), (1, 3), mark:none, color:black, stroke:(1.5pt) ),
      lq.plot((0, 2), (1, 1), mark:none, color:gray, stroke:(dash:"dashed")),
      lq.place(-0.2,1)[$t_c$]


    )
  },
  caption:[Creation of a shock]
)


Another situation can arise when we consider discontinuities, when the characteristics do not move into, but outwards from a discontinuity. 
For a scalar law, this happens when $lambda(u_l)<lambda(u_r)$ and the problem is again given by the Riemann problem @riemann_problem.

#figure(
{
    show: char-plot
    lq.diagram(
      height: 5cm,
      lq.plot((0, 1), (0, 1), mark:none, color:black, stroke:(1.5pt) ),
      lq.plot((1, 4), (1, 3), mark:none, color:black),
      lq.plot((1, 3), (1, 3), mark:none, color:black),
      lq.plot((1, 2), (1, 3), mark:none, color:black),
      lq.plot((1, 2.5), (1, 3), mark:none, color:black),
      lq.plot((1, 3.5), (1, 3), mark:none, color:black),
    )
  },
  caption:[A discontinuity diverging into multiple characteristics, leaving a region undetermined]
)

A solution would be to continue representing the wave by a discontinuity creating an expansion shock.
This would respect the Rankine-Hugoniot condition @rankine-hugoniot but is not an acceptable solution. 
If we try to determine the value on the rays by tracing it backwards in time, we end up at the discontinuity which doesn't carry any information. We therefore impose an additional admissibility criterion for discontinuities.

#definition([Lax entropy condition])[
  A discontinuity of speed $S$ in the $lambda_i$ field is admissible if
  $
  lambda_i (bold(u)_L) > S > lambda_i (bold(u)_R) .
  $
]

This means that for a discontinuity to be admissible, the characteristics should go into the shock, which isn't the case for the expansion shock.

The correct physical solution called a *rarefaction fan* is constructed using self-similarity, that is it doesn't change depending of the scale of the variables.
#let tu = $tilde(u)$
#let btu = $bold(tilde(u))$

For a scalar law, setting $u(x,t) = tu(xi)$, with $xi=x/t$ gives
$
(lambda(tu(xi)) - xi) tu'(xi) = 0
$
so at each point either $tu' = 0$, giving a constant state, or
$
lambda(tu(xi)) = xi .
$<fan_condition>
and for $lambda$ invertible, we obtain
$
tu(xi) = lambda^(-1)(xi)
$<rarefaction_sol>
A discontinuous jump would have to satisfy @rankine-hugoniot, and the only such jump is the expansion shock rejected above, hence the solution must be continuous.
This forces $tu(xi) -> u_l$ as $xi -> lambda(u_l)$ and $tu(xi) -> u_r$ as $xi -> lambda(u_r)$.
The fan is therefore bounded by the two characteristics of slopes $lambda(u_l)$ and $lambda(u_r)$, called its _head_ and _tail_ respectively, and @rarefaction_sol gives the solution between them.
#example[
  For Burgers' equation $lambda(u) = u$ is the identity, so the fan is simply
  $u = x slash t$, with head $u_l$ and tail $u_r$.
]


=== The Riemann problem for non-linear hyperbolic systems

Now that we looked into non-linear scalar quasi-linear forms, let us now look how systems behave on Riemann problems.
Instead of having only one characteristics, systems have several, ordered by their respective speed $lambda_1 (bold(u)) <= dots <= lambda_i (bold(u)) <= dots <= lambda_(m)(bold(u))$.
Assuming strict hyperbolicity, the eigenvalues are distinct and each family is associated with a single eigenvector direction. Each family then produces one wave separating two constant states.
#figure(
{
    show: char-plot
    lq.diagram(
      height: 5cm,
      lq.plot((2, 5), (0, 3), mark:none, color:black),
      lq.plot((2, 4), (0, 3), mark:none, color:black),
      lq.plot((2, 3), (0, 3), mark:none, color:black),
      lq.plot((2, 2), (0, 3), mark:none, color:black),
      lq.plot((2, 1), (0, 3), mark:none, color:black),
      lq.plot((2, 1), (0, 3), mark:none, color:black),
      lq.place(1,0.5)[$bold(u)_L = bold(u)^((0))$],
      lq.place(4,0.5)[$bold(u)_R = bold(u)^((m))$],
      lq.place(1.6,2.5)[$bold(u)^((1))$],
      lq.place(2.5,2.5)[$bold(u)^((2))$],
      lq.place(3.3,2.5)[$bold(u)^((dots))$],
      lq.place(4.2,2.5)[$bold(u)^((m minus 1))$],
    )
  },
  caption:[Example solution to the Riemann problem for a hyperbolic system of conservation laws]
)
In order to find the solution to @riemann_problem, we have to find each of these $m-1$ intermediate states as we have $m+1$ regions and we know the left and the right one.
In total that represents $m(m-1)$ unknowns as each unknown intermediate state has $m$ components.
// Each wave leaves exactly one freee parameter, and there are m wave so  m scalar unknowns.

Let us look at a specific characteristic family $i$.
It separates two states $bold(u)^((i-1))$ on the left and $bold(u)^((i))$ on the right.
We have $bold(u)^((0)) = u_L$ and $bold(u)^((m)) = u_R$.
As seen for the scalar case, if the associated characteristic speed  on the left is greater than the right one $lambda_i (bold(u)^((i-1))) > lambda_i (bold(u)^((i)))$, we get a shock. The Rankine-Hugoniot conditions give us $m$ equations. If we know one of the two data states, this still leaves one unknown as we have to also determine the shock speed $S$.

If instead the left characteristic is smaller than the right one, we get a rarefaction fan. 
Substituting $bold(u)(x,t) = btu(xi)$ into the quasilinear form @quasilinear gives

$
bold(A)(btu) btu'(xi) = xi btu'(xi)
$<rarefaction_self_similar_sys>

so inside the fan $btu'$ is a right eigenvector of $bold(A)$ with eigenvalue $xi$.
This yields the scalar relation $lambda_i (btu(xi)) = xi$ for one of the $m$ eigenvalues $i$ of $A$ fixed throughout the rarefaction fan.
As for the scalar case, this gives us the solution inside the fan when $lambda_i$ is invertible : $btu(xi) = lambda_(i)^(-1)(xi)$.
But we still want to find a relationship between the left and right data states $bold(u)^((i-1))$ and $bold(u)^((i))$.
Let us use the other information provided by @rarefaction_self_similar_sys, $btu'(xi)$ is an eigenvector of $bold(A)$ associated with the eigenvalue $lambda_i$ as such we have :
$
btu'(xi) in "span"( bold(r)_i (btu)) ,
$
with $bold(r)_i$ the right eigenvector associated with $lambda_i$.
This means there exist $alpha(xi)$ such that $btu'(xi) = alpha(xi) bold(r)_i (btu)$.
Componentwise we get $(btu^('(j))(xi))/(bold(r)_(i)^((j))(btu(xi))) = alpha(xi)$ and multiplying by $dif xi$ gives
$
(btu^('(1)) dif xi)/(bold(r)_(i)^((1))(btu)) = (btu^('(2)) dif xi)/(bold(r)_(i)^((2))(btu)) = dots = (btu^('(m)) dif xi)/(bold(r)_(i)^((m))(btu)) = alpha dif xi
$
We have $btu^('(j)) dif xi = dif btu^((j))$ and dropping $alpha dif xi$ gives $m-1$ relations involving only the state variables called the _Generalized Riemann Invariants_ : 
$
(dif btu^((1)))/(bold(r)_(i)^((1))(btu)) = (dif btu^((2)))/(bold(r)_(i)^((2))(btu)) = dots = (dif btu^((m)))/(bold(r)_(i)^((m))(btu))
$<gri>
Each equality is an ordinary differential equation in two of the state variables.
These $m-1$ relations give us a way to associate $m-1$ variables from the left state to the right state, leaving one undetermined as in the shock case.

Finally in systems we can have one more case, $lambda_i (bold(u)^((i-1))) = lambda_i (bold(u)^((i)))$ with $bold(u)^((i-1)) != bold(u)^((i))$ we get a discontinuity that is neither a shock nor a rarefaction. 
This can't happen in the scalar case when we assume the flux to be convex or concave.
The characteristic speed is given if we know either state across the discontinuity, $S = lambda_i (bold(u)^((i))) = lambda_i (bold(u)^((i-1)))$. 
The Rankine Hugoniot conditions apply as we have a discontinuity and so do the Generalized Riemann Invariants.
For such waves, we have enough equations to determine one state from the other unlike shocks and rarefactions.

=== Characteristic Fields and wave types

Now we have seen 3 types of wave encountered in one dimensional non-linear hyperbolic conservation laws. We have _contact_ waves, _shock_ wave and _rarefaction_ waves.

The one dimensional equations studied during this internship have the property to only have such waves.
Furthermore, we can categorize the waves by the type of the field defined by each eigenvalue $lambda(bold(u))$ in state space.

#definition([Linearly degenerate fields])[ A $lambda_i$-field is called _linearly degenerate_ if 
$
nabla_(bold(u)) lambda_i (bold(u)) dot bold(r)_i (bold(u)) = 0, quad forall bold(u) in RR^m
$
]
This means that $lambda_(i)(bold(u))$ is constant along the direction $bold(r)_i$ which is the only one the state varies across the $i$-th wave. 
Such fields can only produce contact waves.

#definition([Genuinely non-linear fields])[ A $lambda_i$-field is called _genuinely non-linear_ if 
$
nabla_(bold(u)) lambda_i (bold(u)) dot bold(r)_i (bold(u)) != 0, quad forall bold(u) in RR^m
$
]

This means that $lambda_(i)(bold(u))$ is strictly monotone along the direction $"span"(bold(r)_i)$, which again is the only one the state varies across the $i$-th wave. 
Such fields generate either shocks or rarefactions and cannot change their wave type.
This also ensures that $lambda_i$ is invertible and thus make the solutions for rarefactions possible.

One can also have fields that are neither linearly degenerate nor genuinely non-linear creating composite waves. These are encountered in more complicated settings like MHD and do not appear in the classical Euler equations.
For the scalar case, this reduces to assuming that the flux functions are either linear, concave or convex.

=== Analysis of the Euler equations
The conservative formulation of the Euler equations @Euler_cons can be rewritten in quasi-linear form using the Jacobian of the flux $nabla_bold(u) bold(F)(bold(u))$ as 
$
bold(u)_t + nabla_bold(u) bold(F)(bold(u))bold(u)_x = 0
$

Its eigenvalues are
$
lambda_1 = u-a, quad lambda_2 = u " and " lambda_3 = u+a
$
and the associated right eigenvectors are
$
bold(r)_1 = vec(1, u-a, H-u a) , quad bold(r)_2 = vec(1, u, 1/2 u^2) " and " bold(r)_3 = vec(1, u+a, H+u a)
$
where $H = (E+p)/rho$ the total specific enthalpy and $a = sqrt((gamma p)/rho)$ the sound speed.
For the computations see @toro2009Riemann pp. 87-90.

Computing $nabla lambda_i(bold(u)) dot bold(r)_i$, we get that the $lambda_1$ and $lambda_3$ fields are genuinely non-linear.
The waves generated by these two fields are called _accoustic waves_ and $rho, u "and" p$ either change smoothly for a rarefaction or discontinuously for shocks.
On the other hand the $lambda_2$ field is linearly degenerate.
It's generated wave is called the _entropy wave_ and only the density jumps discontinuously while the pressure and the particle velocity stay constant.
As $a>0$, we have that $lambda_1 < lambda_2 < lambda_3$ so the structure of the solution to the Riemann problem will consist of a contact wave surrounded by two non-linear waves (shocks and/or rarefactions).

Until now we only looked at one dimensional systems.
If we introduce additional spatial dimensions, the eigenvalue $u_n$ (velocity normal to the considered direction) is not of multiplicity one any more and we loose strict hyperbolicity.
This linearly degenerate field now carries multiple waves, the entropy wave wich is a density jump and the shear waves wich are tangential velocity jumps.

Finally we will also consider the isentropic euler equations that drops the energy equation which removes the entropy wave.
We are left with the accoustic and the shear waves.
Thus it stays strictly hyperbolic up to two dimensions.


== Numerical methods for hyperbolic systems

Now we will look at the Godunov method, a finite-volumes method designed to solve hyperbolic systems of partial differential equations, like the Euler system.

=== The finite volumes method <finite_volumes>
#let ub = $bold(overline(u))$

The finite volume method results from considering the integral form of conservation laws:
$
dif/(dif t) integral_V bold(u) dif v + integral_Sigma cal(H) dot n dif sigma = 0
$
where $u$ is the conserved quantity vector, $V$ the control volume, $Sigma$ the boundary of $V$, $cal(H)$ the flux tensor and $n$ the outward pointing normal to $Sigma$. 
This is then enforced on each cell of a given discretization of the domain.


Consider a conservation law : 
$
bold(u)_t + sum_(i=1)^m F_(i)(bold(u))_(x_i) = 0
$
the integral form gives : 
$
dif/(dif t) integral_V bold(u) dif v = - sum_(Sigma_k in Sigma) integral_Sigma_k sum_(i=1)^m F_(i)(bold(u)) n_i dif sigma
$

Where $(Sigma_i)$ is a set of faces (or edges in 2D) such that $union Sigma_i = Sigma$ and $Sigma_i inter Sigma_j = emptyset$. 

It is left to choose a numerical approximation for the integral, a timestepping method and most importantly a way to compute the fluxes. 



=== Godunov's method

In 1959, Godunov developed a first-order finite volume method for non-linear hyperbolic conservation laws like the Euler system @godunovFinite. 

The method considers piecewise constant approximation using the cell averages $ub = 1/abs(V) integral_V bold(u) dif V$.
To compute the fluxes it uses an exact one dimensional Riemann solver, that is an algorithm that can find $bold(u)^*(x,t)$, the solution of the Riemann problem @riemann_problem.
The left and right states $bold(u)_L$ and $bold(u)_R$ are the cell averages of the two neighboring cells separated by the face and the numerical flux is given by $F(bold(u)^*(0,t))$, evaluating the solution on the curve $x/t = 0$.

Conceptually we are letting the fluid flow for some small timestep before taking the average of the cell.
As we add and subtract the same amount, the flux, between two cells, this method is conservative.

To determine the timestep we use the criteria that the fastest wave should not travel more than one cell in one timestep. Supposing that no wave acceleration takes place as a consequence of wave interaction, for structured meshes this translates to : 
$
Delta t <= (Delta x_i) / max(abs(lambda))
$

with $lambda$ the eigenvalues of the jacobian of the flux in the direction $x_i$. Then using a coefficient $0<C_"cfl" <1$ we obtain : 
$
Delta t = C_"cfl" (Delta x_i) / max(abs(lambda))
$

Godunov only presented this method for structured meshes, where each direction has its own flux function $F_i$. Under the condition that the system is _rotationally invariant_, which the Euler equations are, this approach is generalized to unstructured meshes. Using the rotation matrix $T(bold(n))$ that transforms vector components aligning them with the normal coordinate axis and leaves the scalar component unchanged, we get : 

$
sum_i F_i (bold(u)) n_i = T(bold(n))^(-1) F_1 (T(bold(n)) bold(u))
$

This in turn makes it possible to approximate the face flux as $integral_Sigma_k sum_i F_(i)(bold(u)) n_i dif sigma approx abs(Sigma_k)  T^(-1) hat(F)$, where $T^(-1) = T(bold(n))^(-1)$. $hat(F) = F_1(bold(u_"rot"^*))$ is the numerical flux obtained using $bold(u_"rot"^*)$ the solution of a one dimensional Riemann solver in the rotated coordinates.


// We finally get : 
// $
// dif/(dif t) integral_V bold(u) dif v = - sum_(Sigma_k in Sigma) abs(Sigma_i)  T^(-1) hat(F)
// $
// We can show that for a cartesian grid and the euler equation in 2D we get the known formula un+1 = un - dt/dx (F--F-) - dt/dy (G--G+)
// how do we compute lengths and volumes?


For a forward in time approximation, the update would then read : 
$
ub_i^(n+1) = ub_i^n - (Delta t)/abs(V_i)  sum_(Sigma_k in Sigma) abs(Sigma_k)  T^(-1) hat(F)
$

// CFL in multi-d ? voir cours chap 3.
// faire les calculs pour euler

==== Riemann solvers for linear systems <riemann_linear>

When the system is linear we have independent transport equations on the characteristic variables $bold(v)_i$ each advected with a speed $lambda_i$. 
The solution to the individual Riemann problems are given by
$
bold(v)_i (x,t) = cases(
  bold(v)_(L,i) "if" x/t > lambda_(i),
  bold(v)_(R,i) "else",
)
$
The full solution is then just a superposition of all these individual solution, $bold(u)(x,t) = sum bold(v)_i (x,t)$.
Finally the flux is computed by $F(bold(u)(0,t))$.

This can be simplified, as we are only interested in $bold(u)(0,t)$ we can only determine $bold(v)_(i)(0,t)$
If $lambda_i>0$ the flow goes to the right then $bold(v)_i = bold(v)_(L,i)$ on the contrary if $lambda_i<0$, $bold(v)_i = v_(R,i)$.


#example("The wave equations")[

The wave equations is a one-dimensional linear system : 

$
cases(
  u_t + 1/rho p_x &= 0,
  p_t + kappa u_x &= 0
)
 
$
It can be rewritten in conservative form : 
$
bold(u)_t + F(bold(u))_x = 0
$

with $bold(u) = vec(u,p)$ and $F(bold(u)) = vec(1/rho p, kappa u)$

Diagonalizing the Jacobian gives : 
$
nabla F(bold(u)) = mat(0, 1/rho; kappa, 0) = mat(1, 1; rho c, -rho c) mat(c, 0; 0, -c) mat(1/2, 1/(2 rho c); 1/2, -1/(2 rho c))
$
with $c = sqrt(kappa/rho)$. We have 
$
bold(u)_t + nabla F(bold(u)) bold(u)_x = 0
$<linear_cons>

The characteristic variables are : 
$
bold(v) = mat(1/2, 1/(2 rho c); 1/2, -1/(2 rho c)) bold(u) = vec(u/2 + p/(2 rho c), u/2 - p/(2 rho c))
$
and we can write :

$
bold(v)_t + mat(c, 0; 0, -c) bold(v)_x = 0
$

We obtain two independent transport equations for each characteristic variable. One is advected to the right with characteristic speed $c$ and the other is advected to the left with a characteristic speed $-c$.
Finally we can use the upwind method and compute the flux.
$
F(u(0,t)) = 
$
]


==== Riemann Solvers for non-linear systems

As we know a more complicated case arises when the characteristic speeds depend on $bold(u)$ resulting in non-linear relationships between the state $bold(u)$ and the flux.
For the scalar case with a convex or concave flux, we have to distinguish between shocks and rarefactions.
We solve them using either the Rankine-Hugoniot conditions or the self similarity solution.

#example("Burgers' equation")[

  The Burgers' equation is given by @burgers is a non linear scalar conservation law with a convex flux.

When $u_L > u_R$ we have a shock. To compute its speed and direction, we use the _Rankine-Hugoniot Condition_ @rankine-hugoniot: 
$
1/2 (u_R^2 - u_L^2) = S (u_R - u_L)
$
and we obtain the shock wave speed as the mean between the left and right state: 
$
S = 1/2 (u_R + u_L)
$
If $S>0$ then the shock is moving to the right and if $S<0$ then the shock is moving to the left.
This gives us the solution : 
$
u(x,t) = cases(
                u_L "if" x-S t <0,
                u_R "if" x-S t >0,
              )
$

Now for the case $u_L < u_R$ we get a rarefaction. The head and tail are given by the flux of the left and right states and we interpolate linearly inbetween them. The solution is given by : 
$
u(x,t) = cases(
                u_L "if" x/t < u_L,
                x/t "if" u_L < x/t < u_R,
                u_R "if" x/t > u_R
              )
$

finally we sample the result at $u(0, t)$ and compute the flux as $f(u(0,t)) = u(0,t)$.
]<exampleburger>

Solving non-linear system on the other hand is more complicated and we have to use more tools.
For more details look at the end of chapter 2 of the book by Toro @toro2009Riemann (pp. 76-85).

==== The exact Riemann Solver for the Euler equations <exact_riemann>

As seen when stuyding the Euler equations, the solution to the Riemann problem has 4 regions, on the left side of the left non-linear wave, we have $bold(u)_L$, then on the right side of the right nonlinear wave, we have $bold(u)_R$.
The region in between the two non linear waves is called the _star region_, it is split in two by a contact wave.
Pressure and particle velocity are constant in the star region but density changes accross the contact wave.
We thus have to determine the 4 unknowns $rho_L^*, rho_R^*, u^* "and" p^* $ as well as the left and right non linear waves (shock speed and the solution inside the rarefaction fans).

The first step of the exact riemann solver is to determine $u^*$ and $p^*$.
For this we distinguish the cases of shocks and rarefactions on the  left and right sides.
For shock waves we start by changing the reference frame to align with the shock speed. This results in the shock speed in the reference frame to be null. We can than apply the Rankine-Hugoniot conditions that are simplified : 
$
bold(F)(bold(hat(u)_(K))) = bold(F)(bold(hat(u))^*) ,
$
with $bold(hat(u))$ the conserved variables in the new frame, where only the particle velocities are adapted and $K$ being either left or right. 
From this we can deduce the relations
$
cases(
  u^* = u_L - f_(L)(p^*, bold(w)_L) &" for the left side case",
  u^* = u_R + f_(R)(p^*, bold(w)_R) &" for the right side case",
)
$<ustar>
where 
$
f_(K)(p^*, bold(w)_K) = (p^* - p_K) sqrt(2/(((gamma+1) rho_K) (p^* + (gamma-1)/(gamma+1) p_K))) , 
$ 
is a function that only depends on $p^*$ and the known state $bold(w)_K$ (left or right).

For the Rarefaction cases, we also obtain the same relations as @ustar but with a different function $f_K$. Using the isentropic law and the generalised Riemann invariants, one can obtain :
$
f_(K)(p^*) = (2 a_K)/(gamma-1) [(p^* /p_K) ^((gamma-1)/(2 gamma)) -1] .
$
We obtain two equation of the form @ustar, one for the left non-linear wave and another for the right. Subtracting them yields 
$
f_(L)(p^*, bold(w_L)) + f_(R)(p^*, bold(w_R)) + u_R - u_L = 0 .
$
This means that $p^*$ is the root of the function : 
$
f(p^*, bold(w_L), bold(w_R)) = f_(L)(p^*, bold(w_L)) + f_(R)(p^*, bold(w_R)) + u_R - u_L 
$
That can be found using Newton's method as it is differentiable and has a simple behaviour. We can then determine $u^*$ using @ustar.
This leaves to determine the densities while determining the non-linear waves. 
The type of non-linear wave is determined by the difference between the pressure in the star region and the outside pressure.
When the pressure in the star region is greater to the exterior one we have a shock. We get the density of the left/right star region from 
$
rho_K^* = rho_K [(p^* /p_K+ (gamma-1)/(gamma+1))/((gamma-1)/(gamma+1) p^* /p_K  +1)]
$
and the shock speed
$
S_L = u_L - a_L sqrt((gamma+1)/(2 gamma) p^* /p_L + (gamma -1)/(2 gamma))
$
or
$
S_R = u_R + a_R sqrt((gamma+1)/(2 gamma) p^* /p_R + (gamma -1)/(2 gamma))
$
with $K$ being either left or right.

When the pressure in the star region is smaller than the exterior one, we have a rarefaction wave. We get the density from the isentropic law 
$
rho^*_K = rho_K (p^* / p_K)^(1/gamma)
$
and the sound speed is
$
a^*_K = a_K (p^* / p_K)^((gamma-1)/(2 gamma))
$

For a left rarefaction wave, the speed of the head and tail are given by 
$
S_(H L) = u_L - a_L &"for the head"\
S_(T L) = u_L^* - a_L^* &"for the tail"\
$
For a right rarefaction wave, the speed of the head and tail are given by 
$
S_(H R) = u_R + a_R &"for the head"\
S_(T R) = u_R^* + a_R^* &"for the tail"\
$

The solution inside a rarefaction fan can be found using the Generalised Riemann Invariants and the slope of the characteristic emanating from the origin to the point $(x,t)$ inside the fan.

Finally we can sample the solution at $(0,t)$ and compute the flux for our finite volume solver.

This solver is exact but may not be computationally friendly for big problems as we have to solve a root finding problem. Furthermore it doesn't allow for vacuum (division by zero) which can happen when removing all the density of a cell.


==== Approximate Riemann Solvers for the Euler system

=== Managing Boundary conditions

=== Computing timestep on unstructured meshes

