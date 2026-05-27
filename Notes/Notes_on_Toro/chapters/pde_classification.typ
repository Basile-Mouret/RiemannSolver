== Classification of PDEs by Tom-Robin Teschner

This teacher from Cranfield had the same questions and clarified it in the following #link("https://cfd.university/learn/10-key-concepts-everyone-must-understand-in-cfd/what-are-hyperbolic-parabolic-and-elliptic-equations-in-cfd/")[blog]

=== Examples of classical PDE's


Laplace : $Delta u = 0$

Wave : $u_(t t) - c^2 Delta u = 0$

Burger : $u_t + u u_x = 0$ 

Reaction : $u_t + r u = 0$

Transport : $u_t + c dot nabla u =0$

Diffusion : $U_t - nu Delta u = 0$

=== Origin of the names 
Let us consider a second order polynomial $A x^2 + B x y + C x^2 + D x + E y + F = 0$\
We have $Delta = B^2 - 4 A C$
- $Delta >0 => "hyberbola"=> "real roots" $
- $Delta =0 => "parabola" => 0 "is root" $
- $Delta <0 => "ellipse"  => "complex roots"$

=== Classification of a PDE
This is a second order PDE : $a u_(x x) + b u_(x y) + c u_(y y) + d u_x + e u_y + f = 0$\
We only consider the highest order partial derivatives (they drive the equation).\
Let $Delta = B^2 - 4 A C$\
As for the polynomial,\
- $Delta >0 => "hyberbolic" $
- $Delta =0 => "parabolic" $
- $Delta <0 => "elliptic" $

=== Classification of a system

==== Transforming second order PDE into conservtive system
$a u_(x x) + b u_(x y) + c u_(y y) + d u_x + e u_y + f = 0$

Let $v = u_x$ and $w = u_y$

Then it becomes when keeping only the higher orders
$
cases(
  a v_x + b v_y + c w_y = 0,
  v_y = w_x
)
$
which we can write in the conservative form 
$
U_x + A(U)U_y = 0
$
with $U = vec(v, w)$
and $A = 1/a mat(b, c; -a, 0)$

The characteristic polynomial of A is $(b/a-X)(-X) + c/a = 0 <=> a X^2 - b X + c = 0$

$Delta = (-b)^2 + 4 a c = b^2 - 4 a c$

- $Delta > 0 => "2 real distinct eigenvalues" => "hyberbolic" $
- $Delta = 0 => "1 real eigenvalue " =>"parabolic" $
- $Delta < 0 => "two complex conjugate eigenvalues" => "elliptic" $
We get the same classification.

== Considering the time variable

For 1d its the same as for our 2d steady state, just replace x by t and y by x.

Why do elliptic pde's not have a time dependency?\
Why does information travel instantaneously in parabolic systems but not in hyperbolic systems?

== Considering multidimensional case

We can write the system in the form (dropping all partial derivatives that are not of order 1)
$
  A(U)U_x + B(U)U_y =0
$

To get a single matrix we multiply them by a scalar and sum them:
$
T(U) = A(U)n_x + B(U)n_y
$
We then compute the eigenvalues of this matrix as a function of $x = n_x/n_y$ and deduce the type of system.

For a bigger system : 
$
  A U_x + B U_y + C U_z = 0
$

$
T = A n_x + B n_y + C n_z
$

then find the roots of the characteristic funtion (separate into two equation in x,y and x,z using galilean referential) and again deducing the type of the system.

== Adding time variable




== Principal symbol and characteristics

In summary to characterise a PDE we look at the highest order partial derivatives, and find the roots to the characteristic polynomial associated.

We can also set it in conservative form and look at the eigenvalues of the matrix (again computing the roots of the characteristic polynomial)

Elliptic has complex characteristics 
Hyperbolic has real and distincts characteristics => each give a speed and direction of propagation
Parabolic has real characteristics with the same value. 

In hyperbolic systems, information can only travel at finite speeds (the given eigenvalues) => $u(x,t) = u_0(x-lambda t)$
In parabolic systems, information travels instantaneously. as $0 in "spec"(M)$ => there is one direction where $u(x,t) = u_0(x-lambda t) = u_0(x)$ ?? 

Why is this important?

Different solvers for each types, e.g. Riemann solvers are used on hyperbolic systems and multigrid works really well for elliptic problems.



=== Classifications of PDEs (seen in PDE course ensimag)

We classify PDEs into the following categories : 


- *Elliptic* : time-independent, describing smooth equilibrium states
- *Parabolic* : time-dependent and diffusive
- *Hyperbolic* : time dependent and wave like, with finite speed of propagation

The second order coefficient tensor can also be used for second order PDEs, if $Q$ is :

- Definite (i.e $forall x, x^T Q x$ has the same sign and is non zero) $=>$ Elliptic
- Positive or Negative Non-Definite (i.e $forall x, x^T Q_x x$ has the same sign and $exists x != 0, x^T Q x = 0$ ) $=>$ Parabolic 
- Neither Definite nor positive or negative (i.e $exists x,y !=0, x^T Q x <0 and y^T Q y > 0$) $=>$ Hyperbolic


In terms of eigenvalues, this means :

- elliptic => all the eigenvalues have the same sign and are non negative
- parabolic $=> 0 in "sp"(Q)$, the matrix is non invertible
- hyperbolic



== Pending Quesetions

What is special about the time variable?
why only real/imaginary characteristics? what happens when mixed? in the blog he tells about mixed hyperbolic/parabolic systems??
