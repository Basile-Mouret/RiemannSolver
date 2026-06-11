using CairoMakie
using FiniteVolumes

x0, x1 = 0.0, 1.0
N = 100 

gamma = 1.4

mesh = generate_1DMesh(x0, x1, N, false)
eq = Euler1D(gamma, :exact) 
bcs = Dict("left" => Outflow(), "right" => Outflow())


rhoL, uL, pL, rhoR, uR, pR = 1.0, -19.59745, 1000.0, 1.0, -19.59745, 0.01
xm = 0.8

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
    
    return [rho, rhou, E]
end

max_time_steps = 1000
final_time = 0.012
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

folder = "media/euler1D/example5_godunov/"
mkpath(folder)

animate_cell_values(mesh, rho_hist; dt_hist = dt_hist)
animate_cell_values(mesh, u_hist; dt_hist = dt_hist)
animate_cell_values(mesh, p_hist; dt_hist = dt_hist)
animate_cell_values(mesh, e_hist; dt_hist = dt_hist)

fig = Figure(size = (1000, 1000))

ax1 = Axis(fig[1, 1], xlabel = "x", ylabel = "ρ")
ax2 = Axis(fig[1, 2], xlabel = "x", ylabel = "u")
ax3 = Axis(fig[2, 1], xlabel = "x", ylabel = "p")
ax4 = Axis(fig[2, 2], xlabel = "x", ylabel = "e")

stairs!(ax1, xmid, rho_hist[end], step = :center)
stairs!(ax2, xmid, u_hist[end], step = :center)
stairs!(ax3, xmid, p_hist[end], step = :center)
stairs!(ax4, xmid, e_hist[end], step = :center)

display(fig)

save(folder*"end_results.png", fig)
println("saved final result")


