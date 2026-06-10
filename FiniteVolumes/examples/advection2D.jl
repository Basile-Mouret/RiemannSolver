using FiniteVolumes

N = 100
#mesh = load_mesh2D("meshes/triangle_conv2.msh")
#mesh = load_mesh2D("meshes/triangle_conv0.msh")
mesh = load_mesh2D("meshes/quad_conv0.msh")
println(mesh)

eq = Advection2D(c = (1.0, 2.0))

u0(x) = sin(x[1])*sin(x[2]) 

display(show_heatmap(mesh, u0.(mesh.cell_centers), "initial conditions"))

boundary_conditions = Dict(
                           "boundary_100" => Outflow(),
                           "boundary_200" => Outflow(),
                           "boundary_300" => Outflow(),
                           "boundary_400" => Outflow(),
                          )

max_time_steps = 10
final_time = 1.0
CFL = 0.8

U_hist, dt_hist = solve(mesh, eq, boundary_conditions, u0; max_time_steps = max_time_steps, CFL = CFL, final_time=1.0)


display(show_heatmap(mesh, U_hist[end-1][:, 1], "just before last"))
display(show_heatmap(mesh, U_hist[end][:, 1], "just last"))

