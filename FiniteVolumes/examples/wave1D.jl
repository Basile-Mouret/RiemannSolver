using CairoMakie
using FiniteVolumes

x0, x1 = 0.0, 1.0
N = 100

kappa, rho = 8.0, 2.0
c = sqrt(kappa / rho)

mesh = generate_1DMesh(x0, x1, N, false)
eq = Wave1D(kappa, rho)
bcs = Dict(:left => Reflecting(), :right => Outflow())

p0(x) = 0.2 < x < 0.4 ? 100 * exp(0.1 / ((x - 0.2) * (x - 0.4))) : 0.0
p0(x) = sin(2*pi*x)
u0_func(x) = 0.0

ic(x) = [p0(x), u0_func(x)]

xmid = cell_centers(mesh)

display(plot1D(xmid, p0.(xmid), u0_func.(xmid); title = "Initial condition"))

v0(x) = 0.5 * p0(x) - u0_func(x) * 0.5 / (rho * c)
w0(x) = 0.5 * p0(x) + u0_func(x) * 0.5 / (rho * c)
display(plot1D(xmid, v0.(xmid), w0.(xmid); title = "Initial characteristic variables"))

max_time_steps = 100
final_time = 1.0
CFL = 0.8

xmid, U_hist, U_exact_hist = solve(mesh, eq, bcs, ic; max_time_steps = max_time_steps, CFL = CFL, final_time=1.0)

entropy_hist = [entropy(eq, mat, cell_width(mesh)) for mat in U_hist]
display(plot(entropy_hist))

p_hist = [mat[:, 1] for mat in U_hist]
u_hist = [mat[:, 2] for mat in U_hist]
p_exact_hist = [mat[:, 1] for mat in U_exact_hist]
u_exact_hist = [mat[:, 2] for mat in U_exact_hist]

animate_1D_solution(xmid, u_hist, "media/wave_u_1d.mp4"; U_exact_hist = u_exact_hist)
run(`xdg-open media/wave_u_1d.mp4`)

animate_1D_solution(xmid, p_hist, "media/wave_p_1d.mp4"; U_exact_hist = p_exact_hist)
run(`xdg-open media/wave_p_1d.mp4`)
