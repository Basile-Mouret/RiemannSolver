using FiniteVolumes
using StaticArrays

mesh = load_mesh2D("meshes/Triangles/triangle_conv3.msh")

kappa, rho = 8.0, 2.0
c = sqrt(kappa / rho)

eq = Wave2D(kappa=kappa, rho=rho)

x0 = 0.5
y0 = 0.5
r = 0.2

# initial conditions
p0 = x -> ((x[1]-x0)^2 + (x[2]-y0)^2) <= r^2 ? 10*(exp(-1.0/(1.0-(1/r) * ((x[1]-x0)^2 + (x[2]-y0)^2)))) : 0.0
u0 = x -> 0.0
v0 = x -> 0.0

ic(x) = SVector(p0(x), u0(x), v0(x))

boundary_conditions = Dict(
                           "boundary_1" => Reflecting2D(),
                           "boundary_2" => Reflecting2D(),
                           "boundary_3" => Reflecting2D(),
                           "boundary_4" => Reflecting2D(),
                          )

max_time_steps = 1000
final_time = 3.0
CFL = 0.9

solve(mesh, eq, boundary_conditions, ic; max_time_steps = max_time_steps, CFL = CFL, final_time = final_time, output_dir = "out/wave2D")
