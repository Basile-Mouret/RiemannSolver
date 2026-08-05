#import "../packages/theorems/config.typ" : *
#show: thmrules.with(qed-symbol: $square$)

= Numerical Instabilities at supersonic speeds

== Grid aligned instabilities at supersonic speed

- what are these instabilities?
- where do they come from?
- what methods to suprress?
- how to analyze a scheme?

== Test cases
In order to check if a numerical scheme is subject to the carbuncle, one has to go through a lot of tests as described in @ 

=== Shock travelling a duct
Mach 6\
tf = 5.0

=== Odd-even decoupling
Same test as the shock travelling the duct but with a stronger artificial instability.

From Quirk
cartesian mesh, with some sensible instability in the middle :
$
Y_(i, j "mid") = Y_(j "mid") + 10^(-6) "for "i" odd"\
Y_(i, j "mid") = Y_(j "mid") - 10^(-6) "for "i" even"
$


=== Supersonic flow around blunt body (half cylinder)

Mesh : half circle, quadrangle, aligned with the half circle
initial conditions and bc the same
$
rho = 1
"sound" = sqrt(kappa*gamma*rho^(gamma-1))
u = "Mach" times "sound"
v = 0
$

=== Elling test
Sometimes carbuncle like solutions are the observed physical solution to triggered instablities.
Can be a solution sometimes. One of these examples is the Elling test.

==== Elling test for the full Euler equations
Data from Fleischman et al, A low dissipation method to cure the grid aligned shock instability
==== Elling test for the SWE
#let Fr = "Fr"

Let's look at a shock in the SWE.
we have two states:
$
vec(h_l, h_l u_l) "and" vec(h_r, h_r u_r)
$
we set ourselves in the shock reference:
$
hat(u_k) = u_k - S
$
with $S = Fr_r sqrt(g h_r)$ the shock speed
The rankine Hugoniot conditions give
$
& cases(
h_l hat(u_l) = h_r hat(u_r),
h_l hat(u_l)^2 + 1/2 h_l^2 = h_r hat(u_r)^2 + 1/2 h_r^2
)\
$

isolating $hat(u_r)^2$ and using the first relation, we get
$
hat(u_r)^2 = 1/2 h_l / h_r (h_r + h_l)
$
using $hat(u_r) = u_r - S = (Fr_r - Fr_l) sqrt(h)$ we get a quadratic equation for the ratio $h_l / h_r$. We obtain : 
$
h_l / h_r = 1/2 (-1 plus.minus sqrt(1 + 8 Delta Fr ^2))
$.

Both $h_l$ and $h_r$ are positive so 
$
h_l = h_r /2 (sqrt(1+8 Delta Fr^2) - 1)
$
Finally using the conservation of mass : 
$hat(u_l) = hat(u_r) h_r/h_l$, we get 
$
u_l &= Fr_s sqrt(h_r) + (2 (Fr_r - Fr_s) sqrt(h_r)) / (sqrt(1 + 8 Delta Fr ^2) - 1) \
&=  sqrt(h_r) (Fr_s + (1 + sqrt(1 + 8 Delta Fr ^2) ) / (4(Fr_r - Fr_s))) \
$

#example("Elling Test for the Shallow Water Equations")[
  The example from @SWE_Bader_Kemm sets $h_l = 1, Fr_l = 30 "and" S = 0$. Here we are on the opposite case, where the left state is given and we have to determine the right one.
  First we can compute $u_l = Fr_l * sqrt(h_l) = 30$. Then we use our formulas to compute $h_r$ and $u_r$ : 
  $
  h_r &= 1/2(sqrt(1+8 times Fr_l^2) - 1) &approx 41.929 \
  u_r &= h_l/h_r u_l &approx 0.71549
  $
  In order to trigger the instability a horizontal filament ($u=0$) is introduced in the middle of the left part as in the Elling test for the Euler equations.

  
]

=== SWE : Kelvin-Helmholtz instability + small perturbation

$
H(y) = cases(
  -&sin(pi (y+1/4)/omega)  &"if" -&1/4 - omega/2 <= y < -&1/4 + omega/2,
  -&1                      &"if" -&1/4 + omega/2 <= y <  &1/4 - omega/2,
   &sin(pi (y-1/4)/omega)  &"if"  &1/4 - omega/2 <= y <  &1/4 + omega/2,
   &1                      &&"otherwise" 
)\

(rho, u, v) = (gamma + r H(y), M_a H(y), delta M_a sin(3 pi x))
$

$t_"end" = 0.8/M_a$

$M_a = 0.01$


=== Slow moving shock
Models that try to 

=== Expansion shock
shock diffracting around a corner
cartesian mesh, 70x70 with a corner of size 20x30 in the bottom left

=== Kinked mach stem (shock on a ramp)
from Quirk
mesh : ramp, 40x50 aligned mesh, 10 cells before ramp starts
Mach 5.09 shock initialized before the corner, gamma = 1.4 

=== Sedov blast wave problem
=== Forward facing step
=== Noh implosion problem
=== Modified quirk tests : advancing shock wave, reflecting shock wave and steady shock wave
=== Slight perturbation in slowly moving shocks
=== Steady circular hydraulic jump
=== Laminar boundary layer test
