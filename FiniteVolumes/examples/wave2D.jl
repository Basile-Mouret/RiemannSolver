using FiniteVolumes

mesh = load_mesh2D("meshes/triangle_conv3.msh")

kappa, rho = 8.0, 2.0
c = sqrt(kappa / rho)

eq = Wave2D(kappa=kappa, rho=rho)

x0 = 0.5
y0 = 0.5
r = 0.2

# initial conditions
p0 = x -> ((x[1]-x0)^2 + (x[2]-y0)^2)<= r^2 ? 10*(exp(-1.0/(1.0-(1/r) * ((x[1]-x0)^2 + (x[2]-y0)^2)))) : 0.0
u0 = x->0.0
v0 = x->0.0

ic(x) = [p0(x), u0(x), v0(x)]

boundary_conditions = Dict(
                           "boundary_1" => Reflecting2D(),
                           "boundary_2" => Reflecting2D(),
                           "boundary_3" => Reflecting2D(),
                           "boundary_4" => Reflecting2D(),
                          )

max_time_steps = 1000
final_time = 3.0
CFL = 0.9

U_hist, dt_hist = solve(mesh, eq, boundary_conditions, ic; max_time_steps = max_time_steps, CFL = CFL, final_time=final_time)

p_scalar = [mat[:, 1] for mat in U_hist]
animate_cell_values(mesh, p_scalar; dt_hist = dt_hist)
# save_animation(mesh, U_scalar, "output.mp4"; dt_hist = dt_hist)

