abstract type AbstractEquation end

"""
number of variables of the system
"""
num_vars(::AbstractEquation)::Int = error("num_vars not implemented")

"""
returns the maximum wave speed for a specific equation
"""
function max_wave_speed(::AbstractMesh, ::AbstractEquation, cell_values::Matrix{T})::T where {T<:Real}
    error("max_wave_speed not implemented")
end

"""
computes the timestep from data and the given Courant number 
"""
function compute_dt(::AbstractMesh, eq::AbstractEquation, cell_values::Matrix{T}, CFL::T)::T where {T<:Real}
    error("compute_dt not implemented for $(typeof(eq))")
end

"""
computes the flux between two neighbooring cells by solving the corresponding Riemann problem
"""
flux(eq::AbstractEquation, uL::AbstractVector{T}, uR::AbstractVector{T}, normal) where {T<:Real} = error("flux not implemented for $(typeof(eq))")

# 1D
abstract type AbstractEquation1D <:AbstractEquation end

flux(eq::AbstractEquation1D, uL::AbstractVector{T}, uR::AbstractVector{T}, normal) where {T<:Real} = flux(eq, uL, uR)


# 2D
abstract type AbstractEquation2D <: AbstractEquation end

"""
defines how to save the data

- `name`    : label shown in the viewer
- `kind`    : `:scalar` or `:vector`
- `extract` : function defining how to compute the wanted values from the data
"""
struct OutputField
    name    :: String
    kind    :: Symbol
    extract :: Function
end

"""
defines the specific output fields for each type of equation
by default just gives the values returned by the solver
"""
output_fields(eq::AbstractEquation) = [OutputField("u$i", :scalar, U -> U[i]) for i in 1:num_vars(eq)]
