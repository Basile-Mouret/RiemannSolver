struct Advection2D <: AbstractEquation1D
    c::NTuple{2, Float64}
end

Advection2D(; c::Float64) = Advection1D(c)

num_vars(::Advection2D) = 2


function flux(eq::Advection1D, UL::Vector{Float64}, UR::Vector{Float64}, normal::NTuple{2, Float64})
    cn = eq.c ⋅ normal
    if cn > 0
        return cn*UL
    else
        return cn*UR
    end
end

