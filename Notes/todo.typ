Solveur de Riemann pour les équations d'euler en 1D

solve the EVP:
$
cases(
  rho_t + (rho u)_x = 0 &"Conservation of mass",
  (rho u)_t + (rho u^2 + p)_x = 0 &"Conservation of momentum",
  E_t + [(E+p)u]_x = 0 &"Conservation of energy",
  u_0(x) = cases(u_L "when" x<0,u_R "when" x>0) &"Initial condition"
)
$
With $E = rho(1/2 u^2 + e)$


Comprendre la différence entre shéma volumes finis, solveur de riemann et modèle de godunov.

- vim : w, b, e, ?, learn right  and "". use around and inner.



