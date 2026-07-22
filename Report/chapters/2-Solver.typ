#let eq = $arrow.l.r.double.long$
#let pm = $plus.minus$
#let num(eq) = math.equation(block: true, supplement:none, numbering: it => {numbering("(1.1)", counter(heading).get().first(), it)}, eq)

#pagebreak() 

= Implementation in Julia

The first part of this internship was dedicated to understand and implement the Godunov method. I chose to use Julia, a high level programming languge designed for scientific computing.

== Structure of the code

The dependencies of the Julia project are defined in the `Project.toml` file and they are versionned in the `Manifest.toml` file.
This makes it easy to use the same packages on other devices.

=== Source code
The source code is divided in multiple parts, each defining a specific part of the solver.
The custom *mesh* structure holds the data of the mesh (points, cells, faces, boundaries) as well as some precomputed values used in the solver loop. 
Some helper functions makes it easy to generate regular meshes (1D) and load meshes from a `.msh` file.
The solver uses types to distinguish between different *equations*.
To define an equation one needs to define the number of variables, a numerical flux and a cfl condition.
*Boundary conditions* are also defined in structures and are solved using ghost cells.
The solver interface then requires a mesh, an equation type and a dictionnary linking physical borders to boundary conditions.
The solver in itself uses the finite volume method as in @finite_volumes.
The data is written in `.vtu` files during the simulation.
The data to be written is defined by a method `output_fields` for each type of equation.

== Models

Some models are already implemented and can be directly used.
We started with linear models like the linear transport and the wave system. 
For such equations we proceed as described in @riemann_linear and compute the numerical flux from the upwind direction.
Then we experienced with the non linear Burgers' equation. To compute the numerical flux we proceeded as in @exampleburger. In comparison to the linear case, we now had to compute the fastest wave dynamically as the state was influencing it.
Finally we implemented the full Euler equations for ideal gases using the Godunov method with multiple Riemann solvers.
== Tests

