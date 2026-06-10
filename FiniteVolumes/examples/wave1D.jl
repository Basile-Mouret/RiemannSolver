using CairoMakie
using FiniteVolumes

x0, x1 = 0.0, 1.0
N = 100

kappa, rho = 8.0, 2.0
c = sqrt(kappa / rho)

mesh = generate_1DMesh(x0, x1, N, false)
eq = Wave1D(kappa, rho)
bcs = Dict("left" => Reflecting(), "right" => Outflow())

p0(x) = 0.2 < x < 0.4 ? 100 * exp(0.1 / ((x - 0.2) * (x - 0.4))) : 0.0
p0(x) = sin(2*pi*x)
u0_func(x) = 0.0

ic(x) = [p0(x), u0_func(x)]


xmid = mesh.cell_centers
display(plot_cell_values(mesh, p0.(xmid); title = "Initial condition (p)"))
display(plot_cell_values(mesh, u0_func.(xmid); title = "Initial condition (u)"))

v0(x) = 0.5 * p0(x) - u0_func(x) * 0.5 / (rho * c)
w0(x) = 0.5 * p0(x) + u0_func(x) * 0.5 / (rho * c)
display(plot_cell_values(mesh, v0.(xmid); title = "Initial values for a characteristic variable"))

max_time_steps = 100
final_time = 1.0
CFL = 0.8

U_hist, dt_hist = solve(mesh, eq, bcs, ic; max_time_steps = max_time_steps, CFL = CFL, final_time=1.0)

entropy_hist = [entropy(eq, mat, mesh.cells_length[1]) for mat in U_hist]
display(plot(entropy_hist))

p_hist = [mat[:, 1] for mat in U_hist]
u_hist = [mat[:, 2] for mat in U_hist]

animate_cell_values(mesh, u_hist, "media/wave_u_1d.mp4"; dt_hist = dt_hist)
run(`xdg-open media/wave_u_1d.mp4`)

animate_cell_values(mesh, p_hist, "media/wave_p_1d.mp4"; dt_hist = dt_hist)
run(`xdg-open media/wave_p_1d.mp4`)
