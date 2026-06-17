using FiniteVolumes
using StaticArrays

x0, x1 = 0.0, 1.0
N = 100

gamma = 1.4

mesh = generate_1DMesh(x0, x1, N, false)
eq = Euler1D(gamma, :exact)
bcs = Dict("left" => Outflow(), "right" => Outflow())

rhoL, uL, pL, rhoR, uR, pR = 1.0, 0.0, 1000.0, 1.0, 0.0, 0.01
xm = 0.5

# Initial Conditions (Left state for x <= 0.5, Right state for x > 0.5)
rho0(x) = x <= xm ? rhoL : rhoR
u0(x)   = x <= xm ? uL : uR
p0(x)   = x <= xm ? pL : pR

function ic(x)
    rho = rho0(x)
    u = u0(x)
    p = p0(x)

    rhou = rho * u
    E = p / (gamma - 1.0) + 0.5 * rho * u^2

    return SVector(rho, rhou, E)
end

max_time_steps = 1000
final_time = 0.012
CFL = 0.9

solve(mesh, eq, bcs, ic; max_time_steps = max_time_steps, CFL = CFL, final_time = final_time, output_dir = "out/euler_godunov_3")
