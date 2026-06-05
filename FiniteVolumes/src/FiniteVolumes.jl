module FiniteVolumes

include("mesh.jl")
include("RiemannSolver.jl")

include("boundary_conditions/BoundaryConditions.jl")

include("equations/Equations.jl")
include("equations/burgers.jl")
include("equations/advection.jl")
include("equations/wave.jl")
include("equations/euler.jl")

include("solver/timestepping.jl")
include("solver/Solver.jl")

include("visualization/plot1D.jl")

export Mesh1D, generate_1DMesh, cell_centers, cell_width, quadrature_1D, compute_L2_1D
export solve_riemann
export AbstractEquation1D, Advection1D, Wave1D, Burgers1D, Euler1D
export num_vars, max_wave_speed, flux, exact_solution!
export entropy
export AbstractBC1D, Dirichlet, Reflecting, Outflow
export explicit_euler_step!
export solve
export plot1D, animate_1D_solution

end
