abstract type AbstractEquation end

num_vars(::AbstractEquation)::Int = error("num_vars not implemented")

max_wave_speed(::AbstractMesh, ::AbstractEquation, cell_values::Matrix{Float64})::Float64 = error("max_wave_speed not implemented")

compute_dt(::AbstractMesh, eq::AbstractEquation, cell_values::Matrix{Float64}, CFL::Float64)::Float64 = error("compute_dt not implemented for $(typeof(eq))")

flux(eq::AbstractEquation, uL::Vector{Float64}, uR::Vector{Float64}, normal::Vector{Float64})::Vector{Float64} = error("flux not implemented for $(typeof(eq))")

# TODO
exact_solution!(utrue::Matrix{Float64}, eq::AbstractEquation, xmid::Vector{Float64}, ic::Function, bcs::Dict, x0::Float64, x1::Float64, t::Float64) = error("exact_solution! not implemented for $(typeof(eq))")

# 1D
abstract type AbstractEquation1D <:AbstractEquation end

flux(eq::AbstractEquation1D, uL::Vector{Float64}, uR::Vector{Float64}, normal::Vector{Float64})::Vector{Float64} = flux(eq, uL, uR)

# 2D
abstract type AbstractEquation2D <: AbstractEquation end

