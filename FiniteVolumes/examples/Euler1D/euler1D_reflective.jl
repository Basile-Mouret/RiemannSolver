using FiniteVolumes
using StaticArrays

x0, x1 = 0.0, 1.0
N = 100

gamma = 1.4

mesh = generate_1DMesh(x0, x1, N, false)
eq = Euler1D(gamma, :Godunov)
bcs = Dict("left" => ReflectingEuler1D(), "right" => ReflectingEuler1D())

# Initial Conditions
rho0 = x -> 1.0
u0 = x -> 0.4 < x < 0.6 ? (exp(0.1 / ((x - 0.4) * (x - 0.6)))) / abs(exp(0.1 / ((0.5 - 0.4) * (0.5 - 0.6)))) : 0.125
p0 = x -> 1.0

function ic(x)
    rho = rho0(x)
    u = u0(x)
    p = p0(x)

    rhou = rho * u
    E = p / (gamma - 1.0) + 0.5 * rho * u^2

    return SVector(rho, rhou, E)
end

max_time_steps = 1000
final_time = 10.0
CFL = 0.9

solve(mesh, eq, bcs, ic; max_time_steps = max_time_steps, CFL = CFL, final_time = final_time, output_dir = "out/euler1D_reflective")
