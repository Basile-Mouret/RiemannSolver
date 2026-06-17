using FiniteVolumes
using StaticArrays

mesh = load_mesh2D("meshes/mesh_hole.msh")
println(keys(mesh.boundary_tags))


γ = 1.4
eq = Euler2D(γ, :exact)


M = 3.3  # inflow Mach number
ρ  = 1.0
u  = M
v  = 0.0
p  = 1.0 / γ
E  = p / (γ - 1) + 0.5 * ρ * (u^2 + v^2)

function ic(x)
    rho = 1.0
    u, v = 0.0, 0.0
    p = 1.0/γ
    E = p/(γ-1.0) + 0.5*rho*(u^2+v^2)
    return SVector(rho, rho*u, rho*v, E)
end

W_inflow = SVector(ρ, ρ*u, ρ*v, E)


boundary_conditions = Dict{String, Union{Outflow, Dirichlet2D, ReflectingEuler2D}}(
                           "Bottom"         => Outflow(),
                           "Top"            => Outflow(),
                           # let-bind the captured value so the closure is type-stable
                           # (capturing a non-const global would box it and allocate)
                           "Left"           => let W = W_inflow; Dirichlet2D((x,t) -> W) end,
                           "Right"          => Outflow(),
                           "HoleBoundary"   => ReflectingEuler2D(),
                          )

max_time_steps = 10000
final_time = 3.0
CFL = 0.9


U_hist, dt_hist = solve(mesh, eq, boundary_conditions, ic; max_time_steps = max_time_steps, CFL = CFL, final_time=final_time)

ρ_hist = [mat[:, 1] for mat in U_hist]
u_hist = [mat[:, 2].^2 + mat[:, 3].^2 for mat in U_hist]
animate_cell_values(mesh, ρ_hist; dt_hist = dt_hist)

