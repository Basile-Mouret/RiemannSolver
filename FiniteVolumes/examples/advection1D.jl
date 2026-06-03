using CairoMakie
using FiniteVolumes

x0, x1 = 0.0, 1.0
N = 100

mesh = generate_1DMesh(x0, x1, N, false)
eq = Advection1D(c = 1.5, flux_type = :upwind)
bcs = Dict(:left => Dirichlet(t -> sin(4 * π * t)), :right => Outflow())
u0(x) = [0.0]

max_time_steps = 100
final_time = 1.0
CFL = 0.8

xmid, U_hist, U_exact_hist = solve(mesh, eq, bcs, u0; max_time_steps = max_time_steps, CFL = CFL, final_time=1.0)

u0_vals = [u0(x)[1] for x in xmid]
display(plot1D(xmid, u0_vals; title = "Initial condition"))

anim_file = "media/advection_1d.mp4"
U_scalar = [mat[:, 1] for mat in U_hist]
U_exact_scalar = [mat[:, 1] for mat in U_exact_hist]
animate_1D_solution(xmid, U_scalar, anim_file; U_exact_hist = U_exact_scalar)

run(`xdg-open $(anim_file)`)
