struct ReflectingEuler1D <: AbstractBC1D end

function apply_ghost(::ReflectingEuler1D, u_interior::Vector{Float64}, x, t::Float64)
    return [u_interior[1], -u_interior[2], u_interior[3]]
end

struct ReflectingEuler2D <: AbstractBC2D end

function apply_ghost(::ReflectingEuler2D, u_interior::Vector{Float64}, x, ::Float64)
    return [u_interior[1], -u_interior[2], -u_interior[3], u_interior[4]]
end

