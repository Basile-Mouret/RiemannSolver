using FiniteVolumes

x0, x1 = 0.0, 1.0
N = 100

mesh = generate_1DMesh(x0, x1, N, false)
eq = Advection1D(c = 1.5, flux_type = :upwind)
boundary_conditions = Dict(
                           "left" => Dirichlet(t->[sin(4*π*t)]),
                           "right" => Outflow()
                          )

u0(x) = [0.0]

max_time_steps = 100
final_time = 1.0
CFL = 0.8

U_hist, dt_hist = solve(mesh, eq, boundary_conditions, u0; max_time_steps = max_time_steps, CFL = CFL, final_time=final_time)

xmid = mesh.cell_centers
u0_vals = [u0(x)[1] for x in xmid]
display(plot_cell_values(mesh, u0_vals; title = "Initial condition"))

U_scalar = [mat[:, 1] for mat in U_hist]
animate_cell_values(mesh, U_scalar; dt_hist = dt_hist)

