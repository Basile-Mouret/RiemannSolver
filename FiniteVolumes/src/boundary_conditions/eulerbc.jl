struct ReflectingEuler1D <: AbstractBC1D end

function apply_ghost(bc::ReflectingEuler1D, u_interior::Vector{Float64}, x, t::Float64)
    return [u_interior[1], -u_interior[2], u_interior[3]]
end
