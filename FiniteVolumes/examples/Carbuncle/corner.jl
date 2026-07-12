using FiniteVolumes
using StaticArrays

mesh = load_mesh2D("meshes/Euler2D/Carbuncle/corner.msh")
println(keys(mesh.boundary_tags))

const γ = 1.4
const Mr = 0.0
const Ms = 5.09

# before the shock
const ρr = 1.2
const pr = 1e5
const ar = sqrt(γ*pr/ρr)
const ur = Mr*ar
const vr = 0.0
const Er = pr/(γ-1.0) + 0.5*ρr*(ur^2+vr^2)

# after the shock
const ρs = ρr * ((γ+1)*(Mr-Ms)^2) / ((γ-1)*(Mr-Ms)^2 + 2)
const ps = pr * (2*γ*(Mr-Ms)^2 - (γ-1)) / (γ+1)
const S = Ms*ar
const us = (1- ρr/ρs) * S + ur * ρr/ρs
const vs = vr
const Es = ps/(γ-1.0) + 0.5*ρs*(us^2+vs^2)

num_flux = IdealGasHLLC()
eq = Euler2D(γ, num_flux)

const xs = 0.5 # position of the shock
function ic(x)
    if x[1] < xs
        return SVector(ρs, ρs*us, ρs*vs, Es)
    else
        return SVector(ρr, ρr*ur, ρr*vr, Er)
    end
end

const W_inflow = SVector(ρs, ρs*us, ρs*vs, Es)
dirichlet_inflow = Dirichlet2D((x,t) -> W_inflow)

boundary_conditions = Dict{String, Union{Outflow, typeof(dirichlet_inflow), ReflectingEuler2D}}(
    "Top"        => Outflow(),
    "Bottom"     => Outflow(),
    "Right"      => Outflow(),
    "Left"       => dirichlet_inflow,
    "Object"     => ReflectingEuler2D(),
)
println(boundary_conditions)

max_time_steps = 1
final_time = 0.01
CFL = 0.5

output_dir="out/Carbuncle/Corner/HLLC"
n_info = 1


solve(mesh,
      eq,
      boundary_conditions,
      ic;
      max_time_steps = max_time_steps,
      CFL = CFL,
      final_time=final_time,
      output_dir=output_dir,
      n_info = n_info
     )


