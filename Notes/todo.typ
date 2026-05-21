Solveur de Riemann pour les équations d'euler en 1D

solve an EVP:

$
cases(
  rho_t + nabla dot (rho U) = 0 &"Conservation of mass",
  (rho U)_t + nabla dot (rho U times.o U + rho I) = 0 &"Conservation of momentum",
  E_t + nabla dot [(E+p)U] = 0 &"Conservation of energy",
  u_0(x) = cases(u_L "when" x<0,u_R "when" x>0) &"Initial condition"
)
$

where $U = vec(u, v, w)$ is the velocity vector (noted $V$ in the book), $rho(x, t)$ is the mass density, $p(x,t)$ is the pressure and $E$ is the total energy by unit volume. 
$E$ is defined as
$
E = rho(1/2 U^2 + e) 
$

This can be rewritten as a divergence form:

$
U_t + nabla dot H = 0
$



