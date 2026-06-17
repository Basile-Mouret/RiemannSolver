using FiniteVolumes
using BenchmarkTools
using StaticArrays
using InteractiveUtils
using Cthulhu

mesh = load_mesh2D("meshes/mesh_hole.msh")
println(keys(mesh.boundary_tags))


const γ = 1.4
eq = Euler2D(γ, :exact)


function ic(x)
    gamma = 1.4
    rho = 1.0
    u, v = 0.0, 0.0
    p = 1.0
    E = p/(gamma-1.0) + 0.5*rho*(u^2+v^2)
    return SVector(rho, rho*u, rho*v, E)
end;

const ρ  = 1.0;
const u  = 3.3;
const v  = 0.0;
const p  = 1.0 / γ;
const E  = p / (γ - 1) + 0.5 * ρ * (u^2 + v^2);

const W_inflow = SVector(ρ, ρ*u, ρ*v, E)
dirichlet_inflow = Dirichlet2D((x,t) -> W_inflow)

boundary_conditions = Dict{String, Union{Outflow, typeof(dirichlet_inflow), ReflectingEuler2D}}(
    "Bottom"       => Outflow(),
    "Top"          => Outflow(),
    "Left"         => dirichlet_inflow,
    "Right"        => Outflow(),
    "HoleBoundary" => ReflectingEuler2D(),
)

max_time_steps = 10000
final_time = 3.0
CFL = 0.9

solve(mesh, eq, boundary_conditions, ic; max_time_steps = max_time_steps, CFL = CFL, final_time=final_time);

#@code_warntype solve(mesh, eq, boundary_conditions, ic; max_time_steps = max_time_steps, CFL = CFL, final_time=final_time)
#@descend solve(mesh, eq, boundary_conditions, ic; max_time_steps = max_time_steps, CFL = CFL, final_time=final_time)
@benchmark solve(mesh, eq, boundary_conditions, ic; max_time_steps = max_time_steps, CFL = CFL, final_time=final_time)

using PProf
using Profile
Profile.Allocs.clear();
results = Profile.Allocs.fetch();
Profile.Allocs.@profile sample_rate=1.0 solve(mesh, eq, boundary_conditions, ic; max_time_steps=300, CFL=0.9, final_time=3.0);
results = Profile.Allocs.fetch();
PProf.Allocs.pprof(results)
