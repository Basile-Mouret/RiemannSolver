struct Advection2D <: AbstractEquation1D
    c::NTuple{2, Float64}
end

Advection2D(; c::Float64) = Advection1D(c)

num_vars(::Advection2D) = 2

max_wave_speed(eq::Advection2D, ::Matrix{Float64}, ::Mesh1D ) = abs.(eq.c)

function flux(eq::Advection1D, UL::Vector{Float64}, UR::Vector{Float64}, normal::NTuple{2, Float64})
    cn = eq.c ⋅ normal
    if cn > 0
        return cn*UL
    else
        return cn*UR
    end
end

function exact_solution!(utrue::Matrix{Float64}, eq::Advection2D, xmid::Vector{Float64},
                         ic::Function, bcs::Dict, x0::Float64, x1::Float64, t::Float64)
    #TODO
end
