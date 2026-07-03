abstract type AbstractEquation end

"""
number of variables of the system
"""
num_vars(::AbstractEquation)::Int = error("num_vars not implemented")

"""
returns the maximum wave speed for a specific equation
"""
max_wave_speed(::AbstractMesh, ::AbstractEquation, cell_values::Matrix{Float64})::Float64 = error("max_wave_speed not implemented")

"""
computes the timestep from data and the given Courant number 
"""
compute_dt(::AbstractMesh, eq::AbstractEquation, cell_values::Matrix{Float64}, CFL::Float64)::Float64 = error("compute_dt not implemented for $(typeof(eq))")

"""
computes the flux between two neighbooring cells by solving the corresponding Riemann problem
"""
flux(eq::AbstractEquation, uL::AbstractVector{Float64}, uR::AbstractVector{Float64}, normal) = error("flux not implemented for $(typeof(eq))")

# # TODO
# """
# computes the exact solution of a simulation
# """
# exact_solution!(utrue::Matrix{Float64}, eq::AbstractEquation, xmid::Vector{Float64}, ic::Function, bcs::Dict, x0::Float64, x1::Float64, t::Float64) = error("exact_solution! not implemented for $(typeof(eq))")

# 1D
abstract type AbstractEquation1D <:AbstractEquation end

flux(eq::AbstractEquation1D, uL::AbstractVector{Float64}, uR::AbstractVector{Float64}, normal) = flux(eq, uL, uR)

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
