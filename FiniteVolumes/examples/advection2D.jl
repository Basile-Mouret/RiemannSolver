using FiniteVolumes

N = 100
mesh = load_mesh2D("meshes/triangle_conv2.msh")

eq = Advection2D(c = (1.0, 0.0))

u0 = x-> x[1]

x0 = 0.5
t0 = 0.7
r = 0.3

circle_dirichlet(x,t) = (t-t0)^2 + (x[2] - x0)^2 <= r*r ? [1.0] : [0.0]



boundary_conditions = Dict(
                           "boundary_1" => Outflow(),
                           "boundary_2" => Outflow(),
                           "boundary_3" => Outflow(),
                           "boundary_4" => Dirichlet2D(circle_dirichlet),
                          )

max_time_steps = 1000
final_time = 3.0
CFL = 0.9

U_hist, dt_hist = solve(mesh, eq, boundary_conditions, u0; max_time_steps = max_time_steps, CFL = CFL, final_time=final_time)

display(plot_cell_values(mesh, u0.(mesh.cell_centers), title = "initial conditions"))

U_scalar = [mat[:, 1] for mat in U_hist]
animate_cell_values(mesh, U_scalar; dt_hist = dt_hist)
# save_animation(mesh, U_scalar, "output.mp4"; dt_hist = dt_hist)

