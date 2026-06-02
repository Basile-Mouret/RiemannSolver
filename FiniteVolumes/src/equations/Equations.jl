abstract type AbstractEquation1D end

num_vars(::AbstractEquation1D)::Int = error("num_vars not implemented")
max_wave_speed(::AbstractEquation1D)::Float64 = error("max_wave_speed not implemented")

function flux(eq::AbstractEquation1D, uL::Vector{Float64}, uR::Vector{Float64})::Vector{Float64}
    error("flux not implemented for $(typeof(eq))")
end

function exact_solution!(utrue::Matrix{Float64}, eq::AbstractEquation1D, xmid::Vector{Float64},
                         ic::Function, bcs::Dict, x0::Float64, x1::Float64, t::Float64)
    error("exact_solution! not implemented for $(typeof(eq))")
end