= The Equations of Fluid dynamics

=== Variable types

- primitive/physical variables
- conserved variables
== Euler Equations

$
cases(
  rho_t + nabla dot (rho V) = 0 &"Conservation of mass",
  (rho V)_t + nabla dot (rho V times.o V + rho I) = 0 &"Conservation of momentum",
  E_t + nabla dot [(E+p)V] = 0 &"Conservation of energy",
)
$

where $V = vec(u, v, w)$ is the velocity vector,\
$rho(x, t)$ is the mass density,\
$p(x,t)$ is the pressure,\
$E$ is the total energy by unit volume, defined as
$
E = rho(1/2 V^2 + e) 
$

with $e$ the specific internal energy.

This can be rewritten as a conservation law:
$
  U_t + F(U)_x  + G(U)_y + H(U)_z = 0
$
and in divergence form:
$
U_t + nabla dot H = 0
$

== Thermodynamic Considerations

We need some closure as there are more unknows than equations.
We should define $e$ the specific internal energy (linked to _First Law of Thermodynamics_).
We can also use $s$ the entropy, linked to the _Second Law of Thermodynamics_.
Sometimes, the temperature $T$ can also be necessary.

=== Equation of State (EOS)

Equation of state are approximate statements about the nature of a material.

We use the variables $p$,$v$ and $T$ to describe thermodynamic systems. When it has reached equilibrium it can be completely described by $p$ and $v$. We then have $T = T(p,v)$. 

For example, for _thermally ideal gases_ $T = (p v)/R$

We can define $e$ by the _caloric EOS_ : $e = e(p,v)$ e.g : $e = (p v)/(gamma-1) = p/(rho (gamma - 1))$

$=>$ much more EOS and thermodynamical variables.

$e$ Specific internal enegery\
$Q$ Heat\
$h$ Specific enthalpy\
$f$ Helmholtz free energy\
$alpha, beta$ Volume expansifvity and compressibility\
$c_p, c_v$ Heat capacity at constant pressure and volume\ 
$a$ Speed of sound\

== Viscous Stresses

We augment the Euler equations by adding viscosity.

Stresses in a fluid are due to the effects of pressure $p$ and the viscous stresses. It can be written as the tensor :  
$
S = -p I + Pi
$
$p$ is defined by and EOS\
to define the viscous stress contribution $Pi$ we canuse the _Newtonian approximation_ relating $Pi$ to $V$ via the Deformation tensor.

This introduces new undetermined quantities, _shear viscosity_ and _bulk viscosity_ that can be approximated using Molecular Theory.

== Heat Conduction

Heat $Q$ contributes to the total energy. It can result from multiple sources : temperature gradient in the fluid, diffusion processes in gas mixtures and radiation. The book only considers the first effect.

We can approximate $Q$ using the graidnet of temperature via Fourier's heat conduction law : $Q = -kappa nabla T$

$kappa$ just as $eta$ depends on T and slightly on $p$. We can relate the using Prandl number ($approx$ constant) 
$
  P_r equiv (c_p eta)/kappa
$

When we consider the effects of viscosity and heat conduction to the Euler equations we get the *Navier-Stokes* equations with heat conduction.


== Integral Form of the Equations

The differential forms assume smoothness of the flow variables. to allow for discontinuities we have to work with the integral forms that are based on control volumes and their boundaries.

Another reason comes from the computational point of view, where discrete domain result in finite control volumes/computational cells. Enforcing the fundamental equations on these volumes lead to _Finite Volume_ Methods.

=== Time Derivatives

the _susbstantial_/_material_ derivative:
$
  (D phi.alt)/(D t) = (partial phi.alt)/(partial t) + V dot "grad"phi.alt
$

$=>$ look into Gauss's theorem.

The material derivative and Gauss's theorem are used to find the differential form of the euler equations.


== Submodels

This section looks into simplified versions of the Euler Equations.

=== Flows with Area Variation

$=>$ the area under study can vary in time. For example taking the cross section of a tube that deforms we get $A = A(x,t)$. It depends on the time and the location along the tube $x$.

=== Axi-Symmetric Flows

When flow is symmetric around a an axis, we can reduce the dimensionality by using the coordinates $(r, z)$.

=== Cylindrical and Spherical Symmetry

Here we can reduce the dimensionality even more using only the distance from the center $r$ as a coordinate.

=== Plain One-Dimensional Flow

1D euler : 
$
  U_t + F(U)_x + 0\
$
where

$
  U = vec(rho, rho u, E), F = vec(rho u, rho u^2 + p, u(E+p))
$

Assuming entropy $s$ is constant the EOS becomes $p=p(rho)equiv C rho^gamma$, $C=$constant.
This makes hte Energy equation redundant and we get a $2 times 2$ system.

$
cases(
  U_t + F(U)_x = 0,
  U = vec(rho, rho u)"," F = vec(rho u, rho u^2+rho a^2)
)
$













