#import "../packages/theorems/config.typ" : *
#show: thmrules.with(qed-symbol: $square$)

= Numerical Instabilities in Godunov-type methods
// grid aligned instabilities at supersonic speeds

Usually complete Riemann solvers, like exact Riemann solvers, Roe's solver and HLLC are considered better in the sense that they capture shocks more sharply than uncomplete solvers like HLL would.
In 1994, J.Quirk published a paper @quirk1994Contribution raising awareness on flaws observed on Godunov-type methods.
Complete Riemann solvers can create spurious solutions that are not observed in more diffusive solvers.
For a study of the origin of the instability and possible cures see @dumbser2004Matrix. 
// They also propose a Matrix analysis to check the stability of numerical methods.
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
As expected the HLL scheme smears the shock wave and doesn't develop any instabilities for both method.
On the other hand using complete Riemann solvers results in the frontal shock being destroyed into a spike like structure for the standard finite volumes method.
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

We can observe the standard finite volume methods behaving like in Quirk's paper @quirk1994Contribution, the complete Riemann solvers leading to instabilities and the incomplete HLL solver's shock being diffused.
The new method on the other hand keeps a sharp shock for every Riemann solver.



=== Supersonic flow around blunt body

The following test case is the reason the studied instability is called the _Carbuncle instability_.
Placing a blunt body inside of a supersonic flow generates a shock in front of its nose that can destabilizes into a sort of carbuncle.
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

For the Roe solver, we see a carbuncle instability in the Finite Volume result but the new method doesn't trigger the instability and we get a sharp shock front.
On the other hand using the HLLC solver gets us the opposite results.
The finite volumes don't trigger the carbuncle and the new method does.
Weirdly enough using the standard finite volumes with HLLC results in a diffuse shock interface eventho it is a complete Riemann solver.
// test the exact same case with my solver
// redo the experiments

This testcase isn't complete it would be better to do multiple simulation with different Mach values, ranging from transonic to hypersonic speeds.
Even if the Carbuncle didn't appear using the new method and the roe solver at Mach 10 doesn't mean it wouldn't appear at Mach 20.

=== Elling test

The Lax-Wendroff theorem tells us that conservative solvers converge to weak solutions but unicity of such weak solutions is not assured in multiple dimensions.
In @elling2009Carbuncle, Elling stated that these observed instabilities may be admissible entropy solutions and thus are incurable.
In order to illustrate this, Elling designed a test in which the carbuncle can be observed to be "physically correct".
The setting is inspired from the testcases in @fleischmann2020Low and @SWE_Bader_Kemm.
We place a steady Mach 30 shock at $x = 50$ on a $[0,100] times [0,40]$ domain discretized by $800 times 320$.
This gives for the upstream variables $(rho, u, v) = (1.0, 30.0, 0.0)$ and downstream variables $(rho, u, v) = (41.929, 0.71549, 0.0)$.
Then we immitate a filament in the upstream part by setting $u = 0.0$ for $19.75 < x < 20.25$.
The final time is set to $t = 1.0$.

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
 caption:[Density plot of the Elling test for the isentropic Euler Equation ($kappa=0.5, gamma = 2$)]
)

The complete Riemann solver create sharp carbuncle like structure and the uncomplette HLL solver supresses using both methods.
The new curl and divergence preserving method creates spurious oscilation both in the upstream and downstream regions that is not observed in the classical finite volumes.


=== Kelvin-Helmholtz instability

Finally as an additional experience I wanted to test the sensitivity of the new method to the KH instability, originating from a perturbation in a shear wave.
The instability grows into vortexes that are expected to be sharper using the new method as it is preserving the curl.

The domain is $[0,2] times [0,1]$ discretized by a $300 times 150$ regular mesh.
We set the left and right boundaries to be periodic and the top and bottom to a null gradient outflow. 
Then we divide the domain into a left and right moving fluids with a smooth transition layer defined by the function $H$ as described in @reportDelGrosso.
$
H(y) = cases(
  -&sin(pi (y+1/4)/omega)  &"if" -&1/4 - omega/2 <= y < -&1/4 + omega/2,
  -&1                      &"if" -&1/4 + omega/2 <= y <  &1/4 - omega/2,
   &sin(pi (y-1/4)/omega)  &"if"  &1/4 - omega/2 <= y <  &1/4 + omega/2,
   &1                      &&"otherwise" 
)\
$
Furthermore a vertical perturbation is also added to trigger the instability. 
$
(rho, u, v) = (gamma + r H(y), M_a H(y), delta M_a sin(3 pi x))
$
with $M_a = 0.1$. 

#figure(
  image("assets/results/KH/Roe_dgCurlDiv_8.svg"),
  caption:[KH instability after 3 timesteps using the new Method with the Roe solver]
)


Only the complete solvers generate the instability in this setting. The new methods is sensible enough to generate 8 waves in the first few steps whereas the finite volume methods only show bigger whorls that take much more time to show.
These whorls then combine until a single one fills the whole domain for all complete solvers.

#figure(
  image("assets/results/KH/Roe_dgCurlDiv_4.svg"),
  caption:[KH instability after 7 timesteps using the new method with the Roe solver]
)

#figure(
  image("assets/results/KH/Roe_FV_4.svg"),
  caption:[KH instability after 8 timesteps using finite volumes with the Roe solver]
)

Using the HLL solver doesn't trigger the instability in all the cases, even when refining the mesh by $16$ giving a $1200 times 600$ grid.

=== Other test

There are a lot of other test I didn't get to implement.
In Quirk's paper, he also included two other testcases :  an expansion shock that is created when a shock encounters a stepdown, and a shock encountering a ramp that can create a kink.
In @reportDelGrosso they mention further test done in the litterature to check grid-aligned instabilities.

