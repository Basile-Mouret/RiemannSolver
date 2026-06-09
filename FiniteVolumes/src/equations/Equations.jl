# 1D
abstract type AbstractEquation1D end

num_vars(::AbstractEquation1D)::Int = error("num_vars not implemented")
max_wave_speed(::AbstractEquation1D, values::Matrix{Float64}, Mesh::Mesh1D)::Float64 = error("max_wave_speed not implemented")
flux(eq::AbstractEquation1D, uL::Vector{Float64}, uR::Vector{Float64})::Vector{Float64} = error("flux not implemented for $(typeof(eq))")
exact_solution!(utrue::Matrix{Float64}, eq::AbstractEquation1D, xmid::Vector{Float64}, ic::Function, bcs::Dict, x0::Float64, x1::Float64, t::Float64) = error("exact_solution! not implemented for $(typeof(eq))")
compute_dt(mesh::Mesh1D, eq::AbstractEquation1D, CFL::Float64, values::Matrix{Float64})::Float64 = error("compute_cfl not implemented for $(typeof(eq))")




# 2D
abstract type AbstractEquation2D end

num_vars(::AbstractEquation2D)::Int = error("num_vars not implemented")
max_wave_speed(::AbstractEquation2D, values::Matrix{Float64}, Mesh::Mesh2D)::Float64 = error("max_wave_speed not implemented")
compute_dt(mesh::Mesh2D, eq::AbstractEquation2D, CFL::Float64, values::Matrix{Float64})::Float64 = error("compute_cfl not implemented for $(typeof(eq))")

