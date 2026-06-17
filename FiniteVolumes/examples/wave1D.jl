using FiniteVolumes
using StaticArrays

x0, x1 = 0.0, 1.0
N = 100

kappa, rho = 8.0, 2.0
c = sqrt(kappa / rho)

mesh = generate_1DMesh(x0, x1, N, false)
eq = Wave1D(kappa, rho)
bcs = Dict("left" => Reflecting(), "right" => Reflecting())

p01(x) = 0.2 < x < 0.4 ? 100 * exp(0.1 / ((x - 0.2) * (x - 0.4))) : 0.0
p0(x) = p01(x)
u0_func(x) = 0.0

ic(x) = SVector(p0(x), u0_func(x))

max_time_steps = 10000
final_time = 10.0
CFL = 0.9

solve(mesh, eq, bcs, ic; max_time_steps = max_time_steps, CFL = CFL, final_time = final_time, output_dir = "out/wave1D")
