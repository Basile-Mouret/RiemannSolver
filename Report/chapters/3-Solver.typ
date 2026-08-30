#let eq = $arrow.l.r.double.long$
#let pm = $plus.minus$
#let num(eq) = math.equation(block: true, supplement:none, numbering: it => {numbering("(1.1)", counter(heading).get().first(), it)}, eq)
= Implementation in Julia

In order to truely understand the computational side of Godunov's method, I implemented the algorithm in Julia.
It is a high level programming language designed for scientific computing, ideal for quick implementations that require performant code.

== Structure of the code

The dependencies of the Julia project are defined in the `Project.toml` file and they are versionned in the `Manifest.toml` file.
This makes it easy to use the same packages on other devices.

=== Source code

The source code is divided in multiple folders, each defining a specific part of the solver.
The custom *mesh* structure holds the data of the mesh (points, cells, faces, boundaries) as well as some precomputed values used in the solver loop. 
Some helper functions makes it easy to generate 1D regular meshes and load 2D meshes from a `.msh` file.
The solver uses types to distinguish between different *equations*.
To define an equation one needs to define the number of variables, a numerical flux and a cfl condition.
*Boundary conditions* are also defined in custom structures and are solved using ghost cells.
The solver interface then requires a mesh, an equation type and a dictionnary linking physical borders to boundary conditions.
The solver in itself uses the finite volume method as in @finite_volumes.
The data is written in `.vtu` files during the simulation and can then be visualized with dedicated softwares like Paraview.
The data to be written is defined by a method `output_fields` for each type of equation.

== Models

Some models are already implemented and can be directly used.
We started with linear models like the linear transport and the wave system. 
For such equations we proceed as described in @riemann_linear and compute the numerical flux from the upwind direction.
Then we experienced with the non linear Burgers' equation. To compute the numerical flux we proceeded as in @exampleburger. In comparison to the linear case, we now had to compute the fastest wave dynamically as it is state dependent.
Finally we implemented the full Euler equations for ideal gases using the Godunov method with multiple Riemann solvers.

== Results

=== Linear Advection
The linear advection problem is very simple as the initial and border conditions are simply advected following the vector $bold(c)$.
As such we can compare the result from the finite volume method to the exact solution.

#figure(
  grid(
    columns: 2,
    image("assets/results/advection_2D.svg"),
    image("assets/results/advection_2D-2.svg")
  ),
 caption:[Linear Advection in 2D, a disk being advected from left to right]
)

This example shows a disk defined on the boundary condition being linearly advected to the right.
The exact solution would be a perfect round circle with value $1$. 
We see that the method is diffusive as the circular shape isn't perfectly clean.
Furthermore the front is more diffused than the rear as it has been transported for longer.

// === Burgers' equation
// show discontinuity
// doesn't work anymore :(

=== Euler equations

For the Euler equations we implemented $4$ numerical fluxes, the Godunov flux with an exact Riemann solver as well as the HLL, HLLC and the Roe flux.
To look at something more aerodynamic than a simple circle, we used the NACA0012 airfoil exposed to a Mach 3 flow.
Its contour can be easily defined by a curve as given on #link("https://en.wikipedia.org/wiki/NACA_airfoil#Equation_for_a_cambered_4-digit_NACA_airfoil")[Wikipedia].

#figure(
  grid(
    columns: 2,
    image("assets/results/naca_mach3_isolines.svg", width:100%),
    image("assets/results/naca_mach3_isolines_zoom.svg"),
  ),
 caption:[Density plot with isolines of a NACA0012 airfoil at Mach 3 (exact Riemann solver).]
)

On the front of the airfoil we have a shock as we can see a rapid jump in the state values and packed isolines.
When zooming on the front of the airfoil we can see some irregularities in the density field.
A more severe case can arise with particular conditions. By increasing the flow speed and using an object that is more blunt and a grid that aligns with the shock a carbuncle like disruption is created on the shock.

// do cylinder carbuncle using godunov, hll, hllc and roe.
// compare to aerosol

// border conditions not perfect as we currently do not account for the characteristic that enters the domain.
// Our border conditions aren't perfect, the outlets consider that all of the characteristics go to the right, but this isn't always true.
// We have seen some instabilities arising from flows that are purely tangential to the border that introduce some normal component and destabilize the simulation.



=== Performances

Performance wise the implementation seems to be on par with high performance solvers like Aerosol.
It is still only single threaded, but adding multiprocessing and gpu can be done without changing to much of the code. 
The only part that would need rework is the main finite volume loop inside the solver such that it doesn't break race conditions.

