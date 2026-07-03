using FiniteVolumes
using StaticArrays

mesh = load_mesh2D("meshes/Euler2D/cylinder2D.msh")
println(keys(mesh.boundary_tags))


const γ = 1.4
eq = Euler2D(γ, :exact)


function ic(x)
    rho  = 1.2;
    p  = 1e5;
    u  = 0.0;
    v  = 0.0;
    E = p/(γ-1.0) + 0.5*rho*(u^2+v^2)
    return SVector(rho, rho*u, rho*v, E)
end;

const ρ  = 1.2;
const p  = 1e5;
const a = sqrt(γ*p/ρ)
const u  = 4.0 * a;
const v  = 0.0;
const E  = p / (γ - 1) + 0.5 * ρ * (u^2 + v^2);

const W_inflow = SVector(ρ, ρ*u, ρ*v, E)
dirichlet_inflow = Dirichlet2D((x,t) -> W_inflow)

boundary_conditions = Dict{String, Union{Outflow, typeof(dirichlet_inflow), ReflectingEuler2D}}(
    "Bottom"       => Outflow(),
    "Top"          => Outflow(),
    "Left"         => dirichlet_inflow,
    "Right"        => Outflow(),
    "Cylinder"     => ReflectingEuler2D(),
)

max_time_steps = 100000
final_time = 0.1
CFL = 0.9


solve(mesh, eq, boundary_conditions, ic; max_time_steps = max_time_steps, CFL = CFL, final_time=final_time, output_dir="out/euler2D")

