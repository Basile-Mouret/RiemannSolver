using FiniteVolumes

N = 100
#mesh = load_mesh2D("meshes/triangle_conv2.msh")
mesh = load_mesh2D("meshes/triangle_conv1.msh")

eq = Advection2D(c = (-0.30, 0.20))

u0(x) = sin(x[1])*sin(x[2]) 

x0 = 0.5
t0 = 0.7
r = 0.1

circle_dirichlet(x,t) = t0 - sqrt(r - (x-x0)*(x-x0)) <= t <= t0 - sqrt(r + (x-x0)*(x-x0)) ? 1.0 : 0.0



boundary_conditions = Dict(
                           "boundary_1" => Dirichlet2D(t->[1.0, 1.0]),
                           "boundary_2" => Outflow(),
                           "boundary_3" => Outflow(),
                           "boundary_4" => Outflow(),
                          )

max_time_steps = 1000
final_time = 3.0
CFL = 0.9

U_hist, dt_hist = solve(mesh, eq, boundary_conditions, u0; max_time_steps = max_time_steps, CFL = CFL, final_time=final_time)

display(plot_cell_values(mesh, u0.(mesh.cell_centers), title = "initial conditions"))

anim_file = "media/advection_2d.mp4"
U_scalar = [mat[:, 1] for mat in U_hist]
animate_cell_values(mesh, U_scalar, anim_file; dt_hist = dt_hist)

run(`xdg-open $(anim_file)`)
