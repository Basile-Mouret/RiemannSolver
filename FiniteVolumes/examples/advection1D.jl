using FiniteVolumes
using StaticArrays

x0, x1 = 0.0, 1.0
N = 100

mesh = generate_1DMesh(x0, x1, N, false)
eq = Advection1D(c = 1.5, flux_type = :upwind)
boundary_conditions = Dict(
                           "left" => Dirichlet(t -> SVector(sin(4 * π * t))),
                           "right" => Outflow()
                          )

u0(x) = SVector(0.0)

max_time_steps = 100
final_time = 1.0
CFL = 0.8

solve(mesh, eq, boundary_conditions, u0; max_time_steps = max_time_steps, CFL = CFL, final_time = final_time, output_dir = "out/advection1D")
