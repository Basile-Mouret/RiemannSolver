#import "../packages/theorems/config.typ" : *
#show: thmrules.with(qed-symbol: $square$)

= Numerical Instabilities at supersonic speeds

== The Carbuncle instability

== Test cases

=== Slow moving shock
=== Kinked mach stem (shock on a ramp)
=== Odd-even decoupling
=== Supersonic flow around blunt body (half circle)
=== Supersonic flow around blunt body (ringleb)
=== Sedov blast wave problem
=== Forward facing step
=== Noh implosion problem
=== Modified quirk tests : advancing shock wave, reflecting shock wave and steady shock wave
=== Slight perturbation in slowly moving shocks
=== SWE : Kelvin-Helmholtz instability + small perturbation
=== Steady circular hydraulic jump
=== Elling test
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

=== Laminar boundary layer test
