
= Notions on Hyperbolic Partial Differential Equations

We consider _hyperbolic PDEs_ and _hyperbolic conservations laws_\
Why? Euler equation are mostly hyperbolic (when not considering heat transfers and viscosity) and Riemann solvers are made for hyperbolic PDEs in mind.


== Quasi-Linear Systems

A PDE can be written as : 

$ U_t + A U_x + B = 0 $
where $U$ is the variable vector of size $m$, A is a $m times m$ matrix, and B is a vector of size m that don't depend on $U$.
- If A and B don't depend on $U$, the system is called *linear*.
- If A or B depend on $x$ or $t$, then it becomes a *linear system with variable coefficients*.
- If A or B depend on U (non linearly) : $A(U)$, $B(U)$, then the system is *quasi-linear*
- If $B=0$ then the system is called *homogeneous*.

$=>$ we have to define BC's and IC's



=== Conservation Laws

$ U_t + F(U)_x = 0 $

Using Jacobian Matrix : $ A(U) = (partial F)/(partial U)$, the equation becomes $ U_t + A(U)U_x = 0 $ which is a homogeneous quasi linear PDE.


We can transform the linear advection and Burgers equation into a conservation law :
$ (partial u)/(partial t) + (partial f(u))/(partial t) $
with $f(u) = a u$ for the linear advection and $f(u) = 1/2 u^2$ for Burgers.

=== Hyperbolic Systems

A system is said to be _hyperbolic_ at a point (x,t) if A has m real eigenvalues and a corresponding set of m linearly independent right eigenvectors.

It is *strictly hyperbolic* if the eigenvalues are distinct.

The system is elliptic if non of the eigenvalues are real.

== Linear Advection Equation
Simplest hyperbolic PDE
$
u_t + A nabla u = u
$
_ why is there a source term u on the RHS ?_ 

=== Initial value Problem
$ 
cases(
u_t + a u_x = 0 &quad forall x in RR"," t >0,
u(x,0) = u_0(x) &quad forall x in RR
)
$

$a$ being a constant wave propagation speed.


=== Characteristic Curve

We call characteristic curve the set of points in phase space $t-x$ where $u = "cst"$.

Take a solution $u(x,t)$ of the IVP and parametrize its curve with a variable $v$ : $u(x(s),t(s))$.
Then $(dif u)/(dif s) =  (dif x)/(dif s)(partial u)/(partial x) + (dif t)/(dif s) (partial u)/(partial t)$ 

For $u = "cst"$ we need $(dif u)/(dif s) = 0$. And by identification we obtain : 

$
cases(
  (dif x)/(dif s) = a,
  (dif t)/(dif s) = 0,
)
$

then $t(s) = s$
and $x(s) = x_0 + a s => x(t) = x_0 + a t$

We finally have $u(x,t) = u_0(x_0) = u_0(x-a t)$

// à revoir le cours


=== Riemann Problem

A Riemann Problem is an IVP of the form
$
cases(
  "PDE",
  u_0(x) = cases(u_L "if" x<0, u_R "if" x>0)
)
$
For the linear advection case, the solution is simply given by
$
u(x,t) = u_0(x-a t) = cases(u_L "if" x-a t<0, u_R "if" x - a t>0)
$

== Linear Hyperbolic Systems

$
U_t + A U_x = 0
$

=== Diagonalisation

If $A$ is diagonalisable (as we work with hyperbolic systems it always is), then the system is easy to solve.\
Let $A = K^(-1) Lambda K$, where $Lambda = "diag"(lambda)$ and $K$ the corresponding right eigenvectors.

Then we solve this system using characteristic variables $ W = K^(-1)U$ which corresponds to solving $m$ transport equations.


=== Small perturbation Equations

$
U_x + A U_y = 0
$
with 
$
U = vec(u,v), A = mat(0, -a^2; - 1, 0)
$
$P_A(X) = X^2 - a^2$, so $"sp"(A) = {-a, +a}$

Then let $ K = mat(a, -a; 1, 1)$ and $K^(-1) = 1/(2 a)mat(1, a; -1, a)$

and let $ W = vec(w_1, w_2) = K^(-1)U$

Then $w_1(x,t) = w_1^((0))(x+a t), w_2(x,t) = w_2^((0))(x-a t)$

with $W^((0)) = K^(-1) U^(0) = 1/(2 a)vec(u_1^((0))(x) + a u_2^((0))(x),u_1^((0))(x) - a u_2^((0))(x))$

and then multiplying $W$ by $K$ we find $U$. 


== Conservation Laws
now $lambda_i$ depend on $U$
=== Integral forms of Conservation laws.
To allow for discontinuities, we write conservaiton laws in integral form.

=== Form I

$
dif/(dif t)integral_(x_R)^x_L U(x,t) dif x = F(U(x_L, t)) - F(U(x_L, t))
$

=== Form II

$
integral_(x_L)^x_R U(x,t_1) - U(x,t_2) dif x = integral_(t_1)^(t_2) F(U(x_L, t)) - F(U(x_L, t)) dif t
$

=== Form III

$
integral.cont [U dif x - F(U) dif t] = 
$

=== System of conseration laws
 constructed from linear system using $F(U) = A U$

=== Non linearities an Shock Formation

the conservation law can be rewritten
$
u_t + lambda(u) u_x = 0
$
whith $lambda(u) = (dif f)/(dif u) = f'(u)$ the characteristic speed .

==== Monotonicity of the characteristic speed
the monotonicity is very important:
- $forall u quad lambda'(u)>0 =>$ convex flux 
- $forall u quad lambda'(u)<0 =>$ concave flux 
- $exists u  quad lambda'(u)=0 =>$  non concave and nonconvex flux

e.g Buckley-Leverett equation
$
cases(
  u_t + f(u)_x = 0,
  f(u) = (u^2)/(u^2 + b(1-u)^2)
)
\
$

Wolfram Math alpha tells me $f''$ has two real roots so The buckely Leverett flux is non convex non concave.




