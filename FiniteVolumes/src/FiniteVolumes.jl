module FiniteVolumes


include("boundary_conditions/BoundaryConditions.jl")

include("mesh/mesh.jl")
include("mesh/mesh1D.jl")
include("mesh/mesh2D.jl")

include("RiemannSolver.jl")

include("equations/Equations.jl")
include("equations/burgers.jl")
include("equations/advection.jl")
include("equations/wave.jl")
include("equations/euler.jl")

include("solver/timestepping.jl")
include("solver/solver.jl")

include("visualization/plot1D.jl")

export AbstractMesh
export Mesh1D, generate_1DMesh, quadrature_1D, compute_L2_1D
export Mesh2D, load_mesh2D, face_outward_normal
export solve_riemann_exact
export AbstractEquation
export AbstractEquation1D, Advection1D, Wave1D, Burgers1D, Euler1D
export num_vars, max_wave_speed, flux, exact_solution!
export entropy
export AbstractBC
export AbstractBC1D, Dirichlet, Reflecting, Outflow
export explicit_euler_step!
export solve
export plot_cell_values, animate_cell_values

end
