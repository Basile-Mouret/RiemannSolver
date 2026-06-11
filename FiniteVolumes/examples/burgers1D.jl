using CairoMakie
using FiniteVolumes

x0, x1 = 0.0, 1.0
N = 100

mesh = generate_1DMesh(x0, x1, N, false)
eq = Burgers1D()
bcs = Dict("left" => Dirichlet(t -> [abs(cos(4 * π * t))]), "right" => Outflow())

function u0(x)
    return [sin(2*π*x)]
end

CFL = 0.9
final_time = 3.0
max_time_steps = 1000

U_hist, dt_hist = solve(mesh, eq, bcs, u0; max_time_steps = max_time_steps, CFL = CFL, final_time=final_time)

u0_vals = [u0(x)[1] for x in mesh.cell_centers]
display(plot_cell_values(mesh, u0_vals; title = "Initial condition"))

U_scalar = [mat[:, 1] for mat in U_hist]
animate_cell_values(mesh, U_scalar; dt_hist = dt_hist)

