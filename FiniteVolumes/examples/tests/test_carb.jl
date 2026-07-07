using FiniteVolumes
using StaticArrays
using BenchmarkTools

mesh = load_mesh2D("meshes/aerosol/carbuncleQuad.msh")
println(keys(mesh.boundary_tags))


const γ = 1.4
eq = Euler2D(γ, :Godunov)
const Mach = 0.8


const ρ = 1.0
const u = 1.0 
const v = 0.0
const a = u/Mach
const p = ρ  * a * a / γ
const E = p / (γ - 1.0) + 0.5 * ρ * (u^2 + v^2)

function ic(x)
    return SVector(ρ, ρ*u, ρ*v, E)
end;

const W_inflow = SVector(ρ, ρ*u, ρ*v, E)
dirichlet_inflow = Dirichlet2D((x,t) -> W_inflow)

boundary_conditions = Dict{String, Union{Outflow, typeof(dirichlet_inflow), ReflectingEuler2D}}(
    "Outlet"       => Outflow(),
    "Inlet"        => dirichlet_inflow,
    "Cylinder"     => ReflectingEuler2D(),
)

max_time_steps = 100000
final_time = 5.0
CFL = 0.8



@benchmark solve(mesh,
      eq,
      boundary_conditions,
      ic;
      max_time_steps = max_time_steps,
      CFL = CFL,
      final_time=final_time,
      dt_max = 0.05,
      output_dir="out/test"
     )


