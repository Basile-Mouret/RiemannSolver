abstract type AbstractEquation end

num_vars(::AbstractEquation)::Int = error("num_vars not implemented")

"""
    OutputField(name, kind, extract)

Describes one field an equation wants written to output (e.g. VTK).

- `name`    : label shown in the viewer (ParaView, ...).
- `kind`    : `:scalar` or `:vector`.
- `extract` : maps a single cell's conserved state `U::SVector` to the field
              value — a `Real` for `:scalar`, an `SVector` of components for
              `:vector`.
"""
struct OutputField
    name    :: String
    kind    :: Symbol
    extract :: Function
end

"""
    output_fields(eq) -> Vector{OutputField}

Fields written to output for equation `eq`. The default dumps each conserved
variable as a named scalar; equations override this to expose physically
meaningful (possibly derived) quantities such as velocity or pressure.
"""
output_fields(eq::AbstractEquation) =
    [OutputField("u$i", :scalar, U -> U[i]) for i in 1:num_vars(eq)]

max_wave_speed(::AbstractMesh, ::AbstractEquation, cell_values::Matrix{Float64})::Float64 = error("max_wave_speed not implemented")

compute_dt(::AbstractMesh, eq::AbstractEquation, cell_values::Matrix{Float64}, CFL::Float64)::Float64 = error("compute_dt not implemented for $(typeof(eq))")

flux(eq::AbstractEquation, uL::AbstractVector{Float64}, uR::AbstractVector{Float64}, normal) = error("flux not implemented for $(typeof(eq))")

# TODO
exact_solution!(utrue::Matrix{Float64}, eq::AbstractEquation, xmid::Vector{Float64}, ic::Function, bcs::Dict, x0::Float64, x1::Float64, t::Float64) = error("exact_solution! not implemented for $(typeof(eq))")

# 1D
abstract type AbstractEquation1D <:AbstractEquation end

flux(eq::AbstractEquation1D, uL::AbstractVector{Float64}, uR::AbstractVector{Float64}, normal) = flux(eq, uL, uR)

# 2D
abstract type AbstractEquation2D <: AbstractEquation end

