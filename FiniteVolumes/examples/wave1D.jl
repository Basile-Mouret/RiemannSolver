using CairoMakie
using FiniteVolumes

x0, x1 = 0.0, 1.0
N = 100

kappa, rho = 8.0, 2.0
c = sqrt(kappa / rho)

mesh = generate_1DMesh(x0, x1, N, false)
eq = Wave1D(kappa, rho)
bcs = Dict("left" => Reflecting(), "right" => Reflecting())

p01(x) = 0.2 < x < 0.4 ? 100 * exp(0.1 / ((x - 0.2) * (x - 0.4))) : 0.0
#p02(x) = 0.5 < x < 0.7 ? 100 * exp(0.1 / ((x - 0.5) * (x - 0.7))) : 0.0
p0(x) = p01(x)#+ p02(x)
u0_func(x) = 0.0

ic(x) = [p0(x), u0_func(x)]
xmid = mesh.cell_centers

v0(x) = 0.5 * p0(x) - u0_func(x) * 0.5 / (rho * c)
w0(x) = 0.5 * p0(x) + u0_func(x) * 0.5 / (rho * c)

max_time_steps = 10000
final_time = 1.0
CFL = 0.9
final_time = 10.0

U_hist, dt_hist = solve(mesh, eq, bcs, ic; max_time_steps = max_time_steps, CFL = CFL, final_time=final_time)

entropy_hist = [entropy(eq, mat, mesh.cell_measure[1]) for mat in U_hist]

p_hist = [mat[:, 1] for mat in U_hist]
u_hist = [mat[:, 2] for mat in U_hist]

#display(plot_cell_values(mesh, p0.(xmid); title = "Initial condition (p)"))
#display(plot_cell_values(mesh, u0_func.(xmid); title = "Initial condition (u)"))
#display(plot_cell_values(mesh, v0.(xmid); title = "Initial values for a characteristic variable"))
#display(plot(entropy_hist))

animate_cell_values(mesh, u_hist; dt_hist = dt_hist)
animate_cell_values(mesh, p_hist; dt_hist = dt_hist)

