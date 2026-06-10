using LinearAlgebra

struct Advection2D <: AbstractEquation2D
    c::NTuple{2, Float64}
end

Advection2D(; c::NTuple{2, Float64}) = Advection2D(c)

num_vars(::Advection2D) = 2

function flux(eq::Advection2D, UL::Vector{Float64}, UR::Vector{Float64}, normal::NTuple{2, Float64})
    cn = eq.c ⋅ normal
    if cn > 0
        return cn*UL
    else
        return cn*UR
    end
end

function compute_dt(mesh::Mesh2D, eq::Advection2D, values::Matrix{Float64}, CFL::Float64)::Float64
    # this is wrong but I want to see if it runs before fixing
    return CFL * minimum(mesh.cell_measure) / abs(maximum(eq.c))
end


