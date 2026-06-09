using CairoMakie
using FiniteVolumes

x0, x1 = 0.0, 1.0
N = 100 

gamma = 1.4

mesh = generate_1DMesh(x0, x1, N, false)
eq = Euler1D(gamma, :exact) 
bcs = Dict(:left => Outflow(), :right => Outflow())


# Initial Conditions (Left state for x <= 0.5, Right state for x > 0.5)
rho0(x) = x <= 0.5 ? 1.0 : 1.0
u0(x)   = x <= 0.5 ? -2.0 : 2.0
p0(x)   = x <= 0.5 ? 0.4 : 0.4 

function ic(x)
    rho = rho0(x)
    u = u0(x)
    p = p0(x)
    
    rhou = rho * u
    E = p / (gamma - 1.0) + 0.5 * rho * u^2
    
    return [rho, rhou, E]
end

xmid = cell_centers(mesh)

max_time_steps = 1000
final_time = 0.15
CFL = 0.8

xmid, U_hist, U_exact_hist, dt_hist = solve(mesh, eq, bcs, ic; 
                                            max_time_steps = max_time_steps, 
                                            CFL = CFL, 
                                            final_time = final_time)

rho_hist = [mat[:, 1] for mat in U_hist]
u_hist = [mat[:, 2]./mat[:,1] for mat in U_hist] 
p_hist = [(gamma - 1.0) .* (mat[:,3] .- 0.5 .* mat[:,2].^2 ./ mat[:, 1]) for mat in U_hist] 
e_hist = [mat[:,3]./mat[:,1] .- 0.5 .* (mat[:,2]./mat[:,1]).^2 for mat in U_hist]

folder = "media/euler1D/ex2/"

animate_1D_solution(xmid, rho_hist, folder*"rho.mp4"; dt_hist = dt_hist)
animate_1D_solution(xmid, u_hist, folder*"u.mp4"; dt_hist = dt_hist)
animate_1D_solution(xmid, p_hist, folder*"p.mp4"; dt_hist = dt_hist)
animate_1D_solution(xmid, e_hist, folder*"e.mp4"; dt_hist = dt_hist)

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
