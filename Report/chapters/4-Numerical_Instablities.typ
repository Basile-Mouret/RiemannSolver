#import "../packages/theorems/config.typ" : *
#show: thmrules.with(qed-symbol: $square$)

= Numerical Instabilities in Godunov-type methods
// grid aligned instabilities at supersonic speeds

Usually complete Riemann solvers, like exact Riemann solvers, Roe's solver and HLLC are considered better in the sense that they capture shocks more sharply than uncomplete solvers like HLL would.
In 1992, J.Quirk published a paper @quirk1994Contribution raising awareness on flaws observed on Godunov-type methods.
Complete Riemann solvers can create spurious solutions that are not observed in more diffusive solvers. 

- what are these instabilities?
- where do they come from?
- what methods to supress?
// - how to analyze a scheme?

== Test cases
In order to check if a numerical scheme is subject to the carbuncle, one has to go through a lot of tests as described in @reportDelGrosso.
The main test cases were implemented in the Aerosol software for the isentropic Euler equation with $gamma=2$ and $kappa=0.5$.

=== Shock travelling a duct

The first testcase is simply checking the behaviour of a shock travelling through a duct.
The domain is a $800 times 20$ rectangle with a structured mesh.We set a shock wave with Mach 6 at the start of the duct and let it travel to the end of the domain.


#figure(
  grid(
    columns: (auto, 1fr, 1fr),
    rows: (auto, auto, auto, auto),
    gutter: 5pt,
    [], align(center)[Finite Volumes], align(center)[New Method],
    align(horizon)[Roe],  image("assets/results/duct/duct_fv_roe.svg"),  image("assets/results/duct/duct_dgCurlDiv_roe.svg"),
    align(horizon)[HLLC],  image("assets/results/duct/duct_fv_hllc.svg"),  image("assets/results/duct/duct_dgCurlDiv_hllc.svg"),
    align(horizon)[HLL],  image("assets/results/duct/duct_fv_hll.svg"),  image("assets/results/duct/duct_dgCurlDiv_hll.svg"),
  ),
 caption:[Mach 6 shockwave at the end of a $800 times 20$ duct]
)
As expected the HLL scheme smears the shock wave and doesn't devellop any instabilities for both method.
On the other hand using complete Riemann solvers results in the frontal shock being destroyed into a spike like structure for the standart finite volumes method.
This isn't the case for the curl and divergence preserving method where the shock stays sharp.

// Mach 6\
// tf = 5.0
//
// The first let us study how each scheme handles shocks travelling a long tube.
// No instability is introduced but rounding error from the grid can become visible at long time horizons for some schemes.
//


=== Odd-even decoupling
The instabilities seen in the previous case come from rounding errors in the mesh alignement. 
In order to accentuate and control these instabilities  we reproduced the Odd-Even decoupling example introduced by Quirk in @quirk1994Contribution.

The same initial condition and domain are used but we include a more sensible alignement error in the middle points of the mesh, creating a zig-zag: 

$
Y_(i, j "mid") = Y_(j "mid") + 10^(-6) "for "i" odd"\
Y_(i, j "mid") = Y_(j "mid") - 10^(-6) "for "i" even"
$

This makes it possible to watch the previous spikes grow earlier in a more regular and controlled way.
#figure(
  grid(
    columns: (auto, 1fr, 1fr),
    rows: (auto, auto, auto, auto),
    gutter: 5pt,
[], align(center)[Finite Volumes], align(center)[New Method],
    align(horizon)[Roe],  image("assets/results/odd-even-decoupling/480/FV_Roe_contour_black.svg"),  image("assets/results/odd-even-decoupling/480/dgCurlDiv_Roe_contour_black.svg"),
    align(horizon)[HLLC],  image("assets/results/odd-even-decoupling/480/FV_HLLC_contour_black.svg"),  image("assets/results/odd-even-decoupling/480/dgCurlDiv_HLLC_contour_black.svg"),
    align(horizon)[HLL],  image("assets/results/odd-even-decoupling/480/FV_HLL_contour_black.svg"),  image("assets/results/odd-even-decoupling/480/dgCurlDiv_HLL_contour_black.svg"),
  ),
 caption:[Contour lines at $X_s approx 480$ of the momentum magnitude of a Mach 6 shock travelling a duct with odd even decoupled instabilities in the mesh]
)

We can observe the standart finite volume methods behaving like in Quirk's paper @quirk1994Contribution, the complete Riemann solvers leading to instabilities and the uncomplete HLL solver's shock being diffused.
The new method on the other hand keeps a sharp shock for every Riemann solver.



=== Supersonic flow around blunt body

The following test case is the reason the studied instability is called the _Carbuncle instability_.
Placing a blunt body inside of a supersonic flow generates a shock in front of its nose that can destabilises into a sort of carbuncle.
This is especially observed when aligning the grid with the shock front.

The domain consist of the area between two half-circles filled with a regular mesh of quadrangles.
The data is initialized to a given Mach flow perpendicular to the blunt body and the outer circle acts as an inflow.

// $
// rho &= 1\
// a &= sqrt(kappa*gamma*rho^(gamma-1)) &=& 1\
// u &= "Mach" times a &=& "Mach"\
// v &= 0
// $

// #figure(
//   image("assets/results/Carbuncle/carbuncle_mesh.svg", width:100%)
// )

#figure(
  grid(
    columns: (1fr, 1fr),
    rows: (auto, 40%),
    gutter: 5pt,
    align(center)[Finite Volumes], align(center)[New Method],
    image("assets/results/Carbuncle/carbuncle_FV_Roe.svg"),
    image("assets/results/Carbuncle/carbuncle_dgCurlDiv_Roe.svg"),
  ),
 caption:[Density plot of the Carbuncle using the Roe solver]
)

#figure(
  grid(
    columns: (1fr, 1fr),
    rows: (auto, 40%),
    gutter: 5pt,
    align(center)[Finite Volumes], align(center)[New Method],
    image("assets/results/Carbuncle/carbuncle_FV_HLLC.svg"),
    image("assets/results/Carbuncle/carbuncle_dgCurlDiv_HLLC.svg"),
  ),
 caption:[Density plot of the Carbuncle using the HLLC solver]
)

The result using the Roe solver is the case that triggered this internship. 
As expected we see a carbuncle instability in the Finite Volume result. 
The new method doesn't trigger the instability and we get a sharp shock front.
On the other hand using the HLLC solver gets us the opposite results.
The finite volumes doesn't trigger the carbuncle and the new method does.
Weirdly enough using the standard finite volumes with HLLC results in a diffuse shock interface even tho it is a complete Riemann solver.
// test the exact same case with my solver
// redo the experiments

=== Elling test

In Ellings paper, he stated that these kind of instabilities could be weak solutions. 
The Lax-Wendroff theorem tells us that ... solvers converge to weak solutions. Unicity of such weak solutions isn't proved.

In order to illustrate this, quirk designed a test in which the carbuncle can be observed as a weak solution???
did he proved it was a weak solution in this case?
what is his take??

I used the testcases presented in Fleischman et al, A low dissipation method to cure the grid aligned shock instability, 
and @SWE_Bader_Kemm as well as the function from @reportDelGrosso.

Initial condition // show image of IC + mesh
Border conditions
Domain/Mesh
final time

#figure(
  grid(
    columns: (auto, 1fr, 1fr),
    rows: (auto, auto, auto, auto),
    gutter: 5pt,
[], align(center)[Finite Volumes], align(center)[New Method],
    align(horizon)[Roe],  image("assets/results/Elling/FV_Roe.svg"),  image("assets/results/Elling/dgCurlDiv_Roe.svg"),
    align(horizon)[HLLC],  image("assets/results/Elling/FV_HLLC.svg"),  image("assets/results/Elling/dgCurlDiv_HLLC.svg"),
    align(horizon)[HLL],  image("assets/results/Elling/FV_HLL.svg"),  image("assets/results/Elling/dgCurlDiv_HLL.svg"),
  ),
 caption:[Todo]
)


=== Kelvin-Helmholtz instability

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


// === Slow moving shock
// Models that try to 
//
// === Expansion shock
// shock diffracting around a corner
// cartesian mesh, 70x70 with a corner of size 20x30 in the bottom left
//
// === Kinked mach stem (shock on a ramp)
// from Quirk
// mesh : ramp, 40x50 aligned mesh, 10 cells before ramp starts
// Mach 5.09 shock initialized before the corner, gamma = 1.4 
//
// === Sedov blast wave problem
// === Forward facing step
// === Noh implosion problem
// === Modified quirk tests : advancing shock wave, reflecting shock wave and steady shock wave
// === Slight perturbation in slowly moving shocks
// === Steady circular hydraulic jump
// === Laminar boundary layer test
