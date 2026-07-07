struct ReflectingEuler1D <: AbstractBC1D end

function apply_ghost(::ReflectingEuler1D, u_interior::AbstractVector{Float64}, x, t::Float64, normal::AbstractVector)
    return SVector(u_interior[1], -u_interior[2], u_interior[3])
end

struct ReflectingEuler2D <: AbstractBC2D end

function apply_ghost(::ReflectingEuler2D, u_interior::AbstractVector{Float64}, x, t::Float64, normal::AbstractVector)
    rhou_n = normal[1] * u_interior[2] + normal[2] * u_interior[3]
    return SVector(u_interior[1], u_interior[2] - 2*rhou_n*normal[1], u_interior[3] - 2*rhou_n*normal[2], u_interior[4])
end

