using CairoMakie
using FiniteVolumes

x0, x1 = 0.0, 1.0
N = 100

mesh = generate_1DMesh(x0, x1, N, false)
eq = Euler1D(1.0)
bcs = Dict(:left => Outflow(), :right => Outflow())

rho0(x) = 1.0+0.0*x 
u0(x) = sin(2*pi*x)
p0(x) = 1.0 + 0.0*x

gamma = 5/3

ic(x) = [rho0(x), u0(x)*rho0(x), rho0(x)*(0.5*u0(x)*u0(x) + 1/((gamma-1)*rho0(x)) * p0(x))]

xmid = cell_centers(mesh)

# display(plot1D(xmid, u0.(xmid); title = "Initial condition"))

max_time_steps = 100
final_time = 1.0
CFL = 0.4

xmid, U_hist, U_exact_hist, dt_hist = solve(mesh, eq, bcs, ic; max_time_steps = max_time_steps, CFL = CFL, final_time=1.0)


rho_hist = [mat[:, 1] for mat in U_hist]

animate_1D_solution(xmid, rho_hist, "media/euler_rho_1d.mp4"; dt_hist = dt_hist)
run(`xdg-open media/euler_rho_1d.mp4`)

