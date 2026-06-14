using CairoMakie
using FiniteVolumes

x0, x1 = 0.0, 1.0
N = 100 

gamma = 1.4

mesh = generate_1DMesh(x0, x1, N, false)
eq = Euler1D(gamma, :exact) 
bcs = Dict("left" => ReflectingEuler1D(), "right" => ReflectingEuler1D())

# Initial Conditions
rho0 = x -> 1.0
u0 = x -> 0.4 < x < 0.6 ? (exp(0.1 / ((x - 0.4) * (x - 0.6))))/abs(exp(0.1 / ((0.5 - 0.4) * (0.5 - 0.6)))) : 0.125
p0 = x -> 1.0

function ic(x)
    rho = rho0(x)
    u = u0(x)
    p = p0(x)
    
    rhou = rho * u
    E = p / (gamma - 1.0) + 0.5 * rho * u^2
    
    return [rho, rhou, E]
end

max_time_steps = 1000
final_time = 10.0
CFL = 0.9

U_hist, dt_hist = solve(mesh, eq, bcs, ic;
                                      max_time_steps = max_time_steps,
                                      CFL = CFL,
                                      final_time = final_time)

xmid = mesh.cell_centers

rho_hist = [mat[:, 1] for mat in U_hist]
u_hist = [mat[:, 2]./mat[:,1] for mat in U_hist]
p_hist = [(gamma - 1.0) .* (mat[:,3] .- 0.5 .* mat[:,2].^2 ./ mat[:, 1]) for mat in U_hist]
e_hist = [mat[:,3]./mat[:,1] .- 0.5 .* (mat[:,2]./mat[:,1]).^2 for mat in U_hist]

folder = "media/euler1D/example1_godunov/"
mkpath(folder)

animate_cell_values(mesh, u_hist; dt_hist = dt_hist)
