using FiniteVolumes
using StaticArrays

mesh = load_mesh2D("meshes/Euler2D/cylinder2D.msh")
println(keys(mesh.boundary_tags))


const γ = 1.4
eq = Euler2D(γ, :exact)
const Mach = 4.0


const ρ_inf = 1.2
const p_inf = 1e5
const a_inf = sqrt(γ*p_inf/ρ_inf)
const u_inf = Mach * a_inf
const v_inf = 0.0
const E_inf = p_inf/(γ-1.0) + 0.5*ρ_inf*(u_inf^2+v_inf^2)

function ic(x)
    return SVector(ρ_inf, ρ_inf*u_inf, ρ_inf*v_inf, E_inf)
end;

const W_inflow = SVector(ρ_inf, ρ_inf*u_inf, ρ_inf*v_inf, E_inf)
dirichlet_inflow = Dirichlet2D((x,t) -> W_inflow)

boundary_conditions = Dict{String, Union{Outflow, typeof(dirichlet_inflow), ReflectingEuler2D}}(
    "Bottom"       => Outflow(),
    "Top"          => Outflow(),
    "Left"         => dirichlet_inflow,
    "Right"        => Outflow(),
    "Cylinder"     => ReflectingEuler2D(),
)

max_time_steps = 100000
final_time = 0.01
CFL = 0.8


solve(mesh, eq, boundary_conditions, ic; max_time_steps = max_time_steps, CFL = CFL, final_time=final_time, output_dir="out/cylinder_Mach$(replace(string(Mach),"."=>"_"))")

