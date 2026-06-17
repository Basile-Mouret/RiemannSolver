using FiniteVolumes
using StaticArrays

x0, x1 = 0.0, 1.0
N = 100

mesh = generate_1DMesh(x0, x1, N, false)
eq = Burgers1D()
bcs = Dict("left" => Dirichlet(t -> SVector(abs(cos(4 * π * t)))), "right" => Outflow())

u0(x) = SVector(sin(2 * π * x))

CFL = 0.9
final_time = 3.0
max_time_steps = 1000

solve(mesh, eq, bcs, u0; max_time_steps = max_time_steps, CFL = CFL, final_time = final_time, output_dir = "out/burgers1D")
