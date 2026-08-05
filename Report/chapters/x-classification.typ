
/* == Mathematical analysis of the euler equations

To mathematically categorize differential equations we have to look at the principal part, that is the highest order derivative terms.
// Why? highest order derivative define the qualitative features like well-posedness, smoothness and information propagation directions
We first construct the polynomial called _principal symbol_ :
$
P(x, xi) = sum_(abs(alpha) = m) c_(alpha)(x) xi^alpha
$
where $c_alpha$ are the coefficient of the highest order derivative terms in the system.
// done by doing a fourier transform
We then call _characteristics_ the non zero directions $xi$ such that $P(x, xi) = 0$.
// solutions in the fourier space
The number of these characteristics can then be used to classify the system.

#link("https://cfd.university/learn/10-key-concepts-everyone-must-understand-in-cfd/what-are-hyperbolic-parabolic-and-elliptic-equations-in-cfd/")[blog on PDE classification]

Tom-Robin Teschner : "Hyperbolic PDEs have real and distinct characteristics, while parabolic PDEs have real characteristics, which will have the same value. Elliptic PDEs, on the other hand, have only imaginary characteristics." 

// elliptic equations are time independant.
// parabolic equations are time dependent and diffusive
// hyperbolic equations are time dependent problems with finite speed propagation of information.

=== Second order scalar equations
For second order scalar equations,
$
sum_i^n sum_j^n a_(i j)(bold(x)) (partial^2 u)/(partial x_i partial x_j) + sum_(i=1)^n (partial u)/(partial x_i) + c(bold(x)) u = f
$

// keeping only the second order terms and doing a fourier transform yields

the principal symbol is given by the quadratic form 
$
Q(xi) = xi^T bold(A) xi 
$

Wk have three distinct categories, that take their name from geometry : 

- elliptic 
- parabolic 
- hyperbolic


- if Q is definite then $forall xi != 0, p(xi) != 0$ so there are no real characteristics and the equation is elliptic.
- if Q is degenerate/singular then there is at least one real charactersistic but not all so the equation is parabolic.
- if Q is indefinite then all the characteristics are real and the equation is hyperbolic.
=== First order systems of partial differential equations
construct
$ T = sum_k A_k n_k$
and solve for
=== Adding the time variable
=> can't be elliptic
=== 1-dimensional time dependent equations
Now we will distinguish time from space variables. 

For 1-dimensional time dependent first order systems of partial differential equations:
// independent variables x and t
// dependent variables u
$
(partial u_i)/(partial t) + sum_(j) a_(i j)(x, t, u) (partial u_j)/(partial x) + b_(i)(x, t, u) = 0
$
that can be rewritten in matrix form as 
$
u_t + bold(a)(x, t, u) u_x + b(x, t, u) = 0
$


=== 1-dimensional time dependent system of pde's

Now we go from a scalar equation to a system of pde's searching for $n$ unknowns $bold(u)(x,t) = vec(u_(1)(x,t), dots, u_(n)(x,t))$

The system can be written : 
$
(partial u_i)/(partial t) + sum_(j) a_(i j)(bold(x), t, bold(u)) (partial u_j)/(partial x) + b_(i)(bold(x), t, bold(u)) = 0
$
that can be rewritten in matrix form as 
$
bold(u)_t + bold(A)(bold(x), t, bold(u)) bold(u)_x + bold(b)(bold(x), t, bold(u)) = 0
$
The system is *hyperbolic* at $(x, t)$ if the matrix $A(x,t)$ is diagonalizable with real eigenvalues.



=== Multi-dimensional time-dependent systems of PDE's

Finally we generalize the number of spatial coordiantes as our equations generally live in 2 or 3 spatial coordinates.

$
(partial u_i)/(partial t) + sum_(k) sum_(j) (a_k)_(i j)(bold(x), t, bold(u)) (partial u_j)/(partial x_k) + b_(i)(bold(x), t, bold(u)) = 0
$
that can be rewritten in matrix form as 
$
bold(u)_t + sum_k bold(A)_(k)(bold(x), t, bold(u)) nabla bold(u)_k + B(bold(x), t, bold(u)) = 0
$
For 3 spacial dimensions, we have : 
$
bold(I) bold(u)_t + bold(A)(bold(x), t, bold(u)) bold(u)_x + bold(B)(bold(x), t, bold(u)) bold(u)_y + bold(C)(bold(x), t, bold(u)) bold(u)_z + D(bold(x), t, bold(u)) = 0
$

=== Classifying the Euler equation

Euler system in @Euler_cons is a time dependent first order system of PDEs.
*/





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
