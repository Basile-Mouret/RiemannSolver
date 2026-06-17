module FiniteVolumes

using StaticArrays

include("boundary_conditions/BoundaryConditions.jl")
include("boundary_conditions/eulerbc.jl")


include("mesh/mesh.jl")
include("mesh/mesh1D.jl")
include("mesh/mesh2D.jl")

include("RiemannSolver.jl")

# Equations
include("equations/Equations.jl")
# 1D Equations
include("equations/burgers.jl")
include("equations/advection.jl")
include("equations/wave.jl")
include("equations/euler.jl")
# 2D Equations
include("equations/advection2D.jl")
include("equations/wave2D.jl")
include("equations/euler2D.jl")

include("solver/timestepping.jl")
include("solver/solver.jl")

include("visualization/plot1D.jl")
include("visualization/plot2D.jl")

include("visualization/data_to_vtk.jl")
include("visualization/vtkStreamWriter.jl")

export AbstractMesh
export Mesh1D, generate_1DMesh, quadrature_1D, compute_L2_1D
export Mesh2D, load_mesh2D, face_outward_normal
export solve_riemann_exact
export AbstractEquation
export AbstractEquation1D, Advection1D, Wave1D, Burgers1D, Euler1D
export AbstractEquation2D, Advection2D, Wave2D, Euler2D
export num_vars, max_wave_speed, flux, exact_solution!
export OutputField, output_fields
export entropy
export AbstractBC, Outflow, apply_ghost
export AbstractBC1D, Dirichlet, Reflecting, ReflectingEuler1D
export AbstractBC2D, Dirichlet2D, Reflecting2D, ReflectingEuler2D
 
export explicit_euler_step!
export solve
export plot_cell_values, animate_cell_values, save_animation
export data_to_vtk, VTKStreamWriter, write_frame!, maybe_write!, close_writer!


end
