module FiniteVolumes
using StaticArrays, LinearAlgebra

# Includes

# Boundary Conditions
include("boundary_conditions/BoundaryConditions.jl")
include("boundary_conditions/eulerbc.jl")


# Meshes
include("mesh/mesh.jl")
include("mesh/mesh1D.jl")
include("mesh/mesh2D.jl")


# Equations
include("equations/Equations.jl")

# Linear advection equation
include("equations/advection/advection.jl")
include("equations/advection/advection2D.jl")

# Wave equations
include("equations/wave/wave.jl")
include("equations/wave/wave2D.jl")

# Burgers' equation 
include("equations/burgers/burgers.jl")

# Euler equations
include("equations/euler/riemann_solvers/exact_riemann_solver.jl")
include("equations/euler/numerical_fluxes/godunov_flux.jl")
include("equations/euler/numerical_fluxes/hll_flux.jl")
include("equations/euler/numerical_fluxes/hllc_flux.jl")
include("equations/euler/numerical_fluxes/roe_flux.jl")
include("equations/euler/euler.jl")
include("equations/euler/euler2D.jl")


# Solver Interface
include("solver/timestepping.jl")
include("solver/solver.jl")


# Visualization
include("visualization/vtkStreamWriter.jl")

# Exports

#Meshes
export AbstractMesh
export Mesh1D, generate_1DMesh, quadrature_1D, compute_L2_1D
export Mesh2D, load_mesh2D, face_outward_normal

# Equations
export AbstractEquation
export AbstractEquation1D, Advection1D, Wave1D, Burgers1D, Euler1D
export AbstractEquation2D, Advection2D, Wave2D, Euler2D
export num_vars, max_wave_speed, flux, exact_solution!
export OutputField, output_fields
export entropy

# Riemann Solvers
export get_star_values, solve_riemann_exact

# Boundary conditions
export AbstractBC, Outflow, apply_ghost
export AbstractBC1D, Dirichlet, Reflecting, ReflectingEuler1D
export AbstractBC2D, Dirichlet2D, Reflecting2D, ReflectingEuler2D
 
# Solver
export explicit_euler_step!
export solve

# Output
export VTKStreamWriter, write_frame!, maybe_write!, close_writer!

end
