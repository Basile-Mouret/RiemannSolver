#import "../packages/theorems/config.typ" : *
#show: thmrules.with(qed-symbol: $square$)

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




= Introduction

This internship is focused on understanding numerical methods for hyperbolic systems of equations, especially for the euler equations. The final objective is to implement a suite of test that trigger a numerical instability at supersonic speeds.

==  The Euler equations

Describing the dynamics of fluids by each particle composing the medium is too difficult. Instead we prefer using macroscopic quantities that result from taking the mean of microscopic values over many particles.

For the dynamics of compressible fluids these macroscopic quantities are :
- the *density* $rho (bold(x),t)$ in $"kg"dot "m"^(-3)$ resulting from taking the mean mass per volume.

- the *velocity vector* $bold(v) (bold(x),t) = mat(u, v, w)^T$ in $"m" dot "s"^(-1)$ resulting from taking the mean velocity of each particle per volume.

- the *pressure* $p(bold(x), t)$ in $"Pa" = "N" dot "m"^(-2)$ resulting from the mean magnitude of force per area.

These are called the _primitive_ variables, but to compute them we use the _conserved variables_ : 

- the *density* $rho (bold(x),t)$
- the *momentum vector* $rho bold(v) (bold(x),t)$
- the *total energy per unit volume* $E(bold(x), t)$

// what is the difference? why do we have multiple representations?
// why mass proportionality is important

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

This is a system of differential equations and thus assumes smooth solutions (partial derivatives exist). In order to handle discontinuous solutions we will need to rewrite it in integral form see @finite_volumes.


// Explaining all the terms
// Equation of state
// Subequations
// Viscous Stresses and NS equations

== Mathematical analysis of the Euler equations
// Todo
== Submodels
=== Barotropic Euler equations
By assuming that the density depends only on pressure, $rho = rho(p)$, the system simplifies a lot 
// what are we loosing?
These are a class of subsystems called the Barotropic Euler equations (from the greek "Baro-" : pressure, "-tropic" : depends on). 
One particular barotropic system arises when we suppose the entropy stays constant, this is called the _isentropic_ Euler equations.
$
dif S = 0 \
dots \
p(rho) = kappa rho^gamma
$
Another Barotropic system comes from the SWE:
explain free surface gravity flows + swe simplifications

$U = vec(phi, phi u)$ and $F(U) = vec(phi u, phi u^2 + 1/2 phi^2)$

we recognize a isentropic euler equation with $kappa=1/2$ and $gamma = 2$.
==== Analysis of the Isentropic Euler equations

$U = vec(rho, rho u)$ and $F(U) = vec(rho u, rho u^2 + 1/2 rho^2)$

The speed of sound is given by 
$
a = sqrt(((partial p) / (partial rho))_s) = sqrt((dif p)/(dif rho)) = sqrt(kappa gamma rho^(gamma-1))
$

We now study the eigenvalues and eigenvectors of the jacobian of the flux.
$
nabla_U F(U) = mat(0, 1; a^2 - u^2, 2u)
$

Its eigenvalues are the solutions to the characteristic polynomial given by
$
lambda^2 - 2 u lambda + u^2 - a^2 = 0\
$
As $Delta = 4 a^2 > 0$, the Jacobian has two real roots
$
lambda_(1,2) = u plus.minus a
$

The corresponding right eigenvectors are
$
K_(1,2) = vec(1, u plus.minus a)
$

Let us now characterize these fields : 
$
&"Let " U in RR^2, &nabla lambda_1 (U) dot K_1(U) &= vec((-u)/rho - ((gamma-1) a) / (2 rho), 1/rho) dot vec(1, u - a) = - ((gamma + 1) a) / (2 rho) != 0\
&"and"             &nabla lambda_2 (U) dot K_2(U) &= vec((-u)/rho + ((gamma-1) a) / (2 rho), 1/rho) dot vec(1, u + a) = ((gamma+1)a)/ (2 rho) != 0\
$

So both field are genuinely non linear.
Now we want to solve the Riemann problem for the isentropic Euler equations.

$
U(x, t=0) = cases(U_l "if" x<=x_m, U_r "if" x>x_m)
$

Applying the Rankine Hugoniot condition across the left wave of speed S_l

$
cases(
  rho_l u_l - rho_* u_* &= S_l (rho_l - rho_*),
  rho_l u_l^2 + kappa rho_l^gamma - rho_* u_*^2 - kappa rho_*^gamma &= S_l (rho_l u_l - rho_* u_*)
)
$
and to the right wave : 
$
cases(
  rho_* u_* - rho_r u_r &= (u-a) (rho_* - rho_r),
  rho_* u_*^2 + kappa rho_*^gamma - rho_r u_r^2 - kappa rho_r^gamma &= (u-a)(rho_* u_* - rho_r u_r)
)
$



= Numerical methods for hyperbolic systems

Following the book on Riemann Solvers by Toro @toro2009Riemann we will implement the Godunov method. It is designed to solve hyperbolic systems of partial differential equations, like the Euler system. We will start by looking into the method as presented in the book, before implementing it in the Julia programming language.

== The finite volumes method <finite_volumes>
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



== Godunov's method

In 1959, Godunov developed a first-order finite volume method for non-linear hyperbolic conservation laws like the Euler system @godunovFinite. 
Hyperbolicity means that the jacobians $nabla bold(F_i)(bold(u))$ have $m$ real eigenvalues and $m$ linearly independent right eigenvectors. It ensure the conservation law is well posed.

The method considers piecewise constant approximation using the cell averages $ub = 1/abs(V) integral_V bold(u) dif V$.
To compute the fluxes it uses an exact one dimensional Riemann solver, that is an algorithm that can find $bold(u)^*(x,t)$ solution of the Riemann problem: 
$
cases(
  bold(u)_t + F(bold(u))_(x) = 0,
  bold(u)(x,0) = cases(bold(u_L) "if" x<0, bold(u_R) "if" x>0)
)
$<riemann_problem>
where the left and right states $bold(u)_L$ and $bold(u)_R$ are the cell averages of the two neighboring cells separated by the face. The flux is then given by $F(bold(u)^*(0,t))$

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

Godunov only presented this method for structured meshes, where each direction has its own flux function $F_i$. Under the condition that the system is _rotationally invariant_, this approach is generalized to unstructured meshes. Using the rotation matrix $T(bold(n))$ that transforms vector components aligning them with the normal coordinate axis and leaves the scalar component unchanged, we get : 

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

== Exact Riemann solvers

Now lets see how to solve the Riemann problem @riemann_problem.
We are considering hyperbolic conservation laws, so the Jacobian $nabla F(bold(u)) = mat((partial u_1)/(partial u_1), dots,(partial u_1)/(partial u_n);, dots.down,;(partial u_n)/(partial u_1), dots, (partial u_n)/(partial u_n))$ is diagonalizable with real eigenvalues. // proof?
Thus there exists an invertible matrix $P(bold(u))$ and a diagonal matrix $D(bold(u))$ such that $nabla F(bold(u)) = P(bold(u)) D(bold(u)) P(bold(u))^(-1)$.

We call _characteristic variables_ the components of the vector 
$
bold(v) = P(bold(u))^(-1) bold(u)
$<characteristics>.
Using them, the conservation law can be rewritten as : 
$
bold(v)_t + D(bold(u)) bold(v)_x = 0
$
which separates the system into $n$ independent equations of the form:
$
(partial v_i)/(partial t) + lambda_(i)(bold(u)) (partial v_i)/(partial x) = 0
$

Where $lambda_(i)(bold(u))$ are the diagonal entries of $D(bold(u))$, i.e. the eigenvalues of $nabla F(bold(u))$. They are called the _characteristic speeds_.



=== Riemann solvers for linear systems <riemann_linear>

When the system is linear, $D$ doesn't depend on $bold(u)$ and we have independent transport equations with an advection speed $d_i$. Now to find $bold(u^*)$, we have to look at the upwind direction, that is the direction in which information is travelling. If $d_i>0$ the flow goes to the right then $u_i^* = u_(L,i)$ on the contrary if $d_i<0$, $u_i^* = u_(R,i)$. Doing this for all the equations we can construct $bold(u)^*$ and compute the flux as $F(bold(u^*))$.


#example("The wave equations")[

The wave equations is a one-dimensional linear system : 

$
cases(
  p_t + 1/rho u_x &= 0,
  u_t + kappa p_x &= 0
)
 
$
It can be rewritten in conservative form : 
$
bold(u)_t + F(bold(u))_x = 0
$

with $bold(u) = vec(p,u)$ and $F(bold(u)) = vec(1/rho u, kappa p)$

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
bold(v) = mat(1/2, 1/(2 rho c); 1/2, -1/(2 rho c)) bold(u) = vec(p/2 + u/(2 rho c), p/2 - u/(2 rho c))
$
and we can write :

$
bold(v)_t + mat(c, 0; 0, -c) bold(v)_x = 0
$

We obtain two independent transport equations for each characteristic variable. One is advected to the right with characteristic speed $c$ and the other is advected to the left with a characteristic speed $-c$.
Finally we can use the upwind method and compute the flux.

]

=== Riemann Solvers for non-linear systems

A more complicated case arises when the characteristic speeds depend on $bold(u)$. This results in non-linear relationships between the state $bold(u)$ and the flux.
// The flow can then have different forms depending on its monotonicity : 
//
// - when $d(u)$ is monotonically increasing, i.e. $d'(u)>0$, the flow is *convex*.
// - when $d(u)$ is monotonically decreasing, i.e. $d'(u)<0$, the flow is *concave*.
// - finally, when $d(u)$ admits and extrema, i.e. $exists u | d'(u) = 0$, the flow is neither convex nor concave.
Still the approach stays similar to the linear case as we look at the orientation of the flow to determine the upwind direction.
But contrary to the linear case, wave distortion can arise from the difference of $d(bold(u))$ between two neighboring points. 
This can lead to the creation of shock and rarefaction waves.

Shock waves are thin layers where properties change very rapidly.
These are generally approximated as mathematical discontinuities.
To find the speed $S$ at which this discontinuity is moving, we use the _Rankine-Hugoniot Condition_ : 
$
Delta bold(f) = S Delta bold(u)
$<rankine-hugoniot>
It is obtained by applying the integral form of the conservation law on a control volume encompassing the shock (see @toro2009Riemann, p. 70).

On the other hand, rarefaction waves are created through a stretching between two states.
To solve it we can't use a discontinuity like in the shock case as this would violate the entropy condition and lead to instabilities.
Instead the rarefaction is modelled by two different waves the head and the tail. The solution inside of the rarefaction fan can then be determined using the self similarity condition and the Generalized Riemann invariants.

#example("Burgers' equation")[

  The Burgers' equation is given by
$
u_t + F(u)_x = 0
$
with $F(u) = 1/2 u^2$ and $(partial F(u))/(partial u) = u$. 

Here the propagation speed depends on the value of $u$ and the flow is non-linear.

include image burgers

To solve the Riemann problem @riemann_problem, we have to distinguish between the shock and the rarefaction case.

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

finally we sample the result at $u^* = u(0, t)$ and compute the flux as $f(u^*)$.
]<exampleburger>

// this part I am not sure
For systems, we have to look at the field of each characteristic@characteristics. 
An eigenvalue $lambda_(i)(bold(u)(bold(x),t))$ characterises the speed of a wave in the physical space ($x-t$). It is associated with a corresponding right eigenvector $bold(K)_(i)(bold(u))$, that defines the paths of the wave in phase space.

A $lambda_i$-characteristic field is called _linearly degenerate_ when 
$
nabla lambda_(i)(bold(u)) dot bold(k)_(i)(bold(u)) = 0, quad forall bold(u) in RR^m
$

this means that the evolution of this wave doesn't affect its speed. Discontinuities on such field are called contact waves. On the other hand, we call a $lambda_i$-characteristic field _genuinely non-linear_ when 
$
nabla lambda_(i)(bold(u)) dot bold(k)_(i)(bold(u)) != 0, quad forall bold(u) in RR^m
$

in that case, the evolution of the wave affects its speed and a discontinuity on such field creates either shocks or rarefactions.

For more details look at the end of chapter 2 of the book by Toro @toro2009Riemann (pp. 76-85).

=== The exact Riemann Solver for the Euler equations <exact_riemann>

Now let's study the one-dimensional euler equation in order to determine the solution to the associated Riemann Problems. 
// The conservative form of the euler equations is 
// $
// bold(u) + bold(F)(bold(u))_x = 0
// $
// with $bold(u) = vec(rho, rho u, E)$ the conserved variables and $bold(F)(bold(u)) = vec(rho u, rho u^2 + p, (E + p) u)$ the flux vector. The total energy per unit volume is given by $E = rho(1/2 u^2 + e)$, with $e$ the specific internal energy given by an equation of state. For ideal gases one has $ e(rho, p) = p/((gamma-1)rho)$, with $gamma = c_p/c_v$ the ratio of specific heats.
//
The conservative formulation @Euler_cons can be rewritten in quasi-linear form using the Jacobian of the flux $nabla_bold(u) bold(F)(bold(u))$ as 
$
bold(u)_t + nabla_bold(u) bold(F)(bold(u))bold(u)_x = 0
$

Its eigenvalues are
$
lambda_1 = u-a, quad lambda_2 = u " and " lambda_3 = u+a
$
and the associated right eigenvectors are
$
bold(k)_1 = vec(1, u-a, H-u a) , quad bold(k)_2 = vec(1, u, 1/2 u^2) " and " bold(k)_3 = vec(1, u+a, H+u a)
$
where $H = (E+p)/rho$ the total specific enthalpy and $a = sqrt((gamma p)/rho)$ the sound speed. For the computations see @toro2009Riemann pp. 87-90.

Computing $nabla lambda_i(bold(u)) dot bold(k)_i$, we get that the $lambda_1$ and $lambda_3$ fields are genuinely non-linear and that the $lambda_2$ field is linearly degenerate.
As $a>0$, we have that $lambda_1 < lambda_2 < lambda_3$ so the structure of the solution to the Riemann problem will consist of a contact wave surrounded by two non-linear waves (shocks and/or rarefactions).

Accross the contact wave, the pressure and particle velocity are constant and the density jumps discontinuously.
Accross rarefaction waves, $rho, u "and" p$ change smoothly.
Finally accross a shock wave, $rho, u "and" p$ change very rapidly.

The solution to the Riemann problem has 4 regions, on the left side of the left non-linear wave, we have $bold(u)_L$, then on the right side of the right nonlinear wave, we have $bold(u)_R$.
The region in between the two non linear waves is called the _star region_, it is split in two by the contact wave.
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

==== Exact Riemann solver for the Barotropic Euler equations

We have $U = vec(rho, rho u)$ and $F(U) = vec(rho u, rho u + kappa rho^gamma)$.



=== Managing Boundary conditions

=== Computing timestep on unstructured meshes

