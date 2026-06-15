using FiniteVolumes

mesh = load_mesh2D("meshes/triangle_conv3.msh")

gamma = 1.4
eq = Euler2D(gamma, :exact)

x0, y0 = 0.5, 0.5

function ic(x)
    r2 = (x[1]-x0)^2 + (x[2]-y0)^2
    rho = 1.0
    u, v = 0.0, 0.0
    p = r2 <= 0.3^2 ? 1.0 : 0.1   # high pressure circle without vacuum
    E = p/(gamma-1.0) + 0.5*rho*(u^2+v^2)
    return [rho, rho*u, rho*v, E]
end

boundary_conditions = Dict(
                           "boundary_1" => ReflectingEuler2D(),
                           "boundary_2" => ReflectingEuler2D(),
                           "boundary_3" => ReflectingEuler2D(),
                           "boundary_4" => ReflectingEuler2D(),
                          )

max_time_steps = 1000
final_time = 5.0
CFL = 0.9


U_hist, dt_hist = solve(mesh, eq, boundary_conditions, ic; max_time_steps = max_time_steps, CFL = CFL, final_time=final_time)

u_scalar = [mat[:, 2].^2 + mat[:, 3].^2 for mat in U_hist]
animate_cell_values(mesh, u_scalar; dt_hist = dt_hist)

