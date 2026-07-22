
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
