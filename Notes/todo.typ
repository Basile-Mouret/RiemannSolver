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

TODO :
- changer pour utiliser gmsh (extraire donnnées d'un fichier, convertir gmsh -> Mesh2D, et pourquoi pas construire des mesh en julia)
- transport 2D conditions de bord (periodique ou function)

- Waves 1D
- Burgers 1D

- Ecrire rapport sur ça

- Euler 1D (Godunov, Solveur de Riemann exact/approché)

- Comprendre et traduire le solveur de Riemann 1D proposé dans Toro (Résumer chap 3 et surtout 4)

- Comprendre et traduire le modèle de Godunov dans Toro (Résumer chap 5)

VIM :
- avoir le go to definition
- utiliser jk
- utiliser w,b,e
- utiliser inner
- regarder pour commenter rapidement plusieurs lignes
- regarder pour que indentation garde sélection
- regarder pour que (), "", etc garde la sélection

Julia:
- CamelCase pour struct
- snake_case pour variables


Questions :  Utiliser GMSH sur 1D aussi?
