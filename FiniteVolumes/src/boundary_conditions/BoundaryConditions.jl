abstract type AbstractBC1D end

struct Dirichlet <: AbstractBC1D
    value::Function
end

struct Reflecting <: AbstractBC1D end

struct Outflow <: AbstractBC1D end

function apply_ghost(bc::AbstractBC1D, u_interior::Vector{Float64}, t::Float64)::Vector{Float64}
    error("apply_ghost not implemented for $(typeof(bc))")
end

function apply_ghost(bc::Outflow, u_interior::Vector{Float64}, ::Float64)
    return copy(u_interior)
end

function apply_ghost(bc::Dirichlet, u_interior::Vector{Float64}, t::Float64)
    return [bc.value(t)]
end

abstract type AbstractBC2D end

struct Dirichlet2D <: AbstractBC2D
    value::Function
end

struct Reflecting2D <: AbstractBC2D end

struct Outflow2D <: AbstractBC2D end

function apply_ghost(bc::AbstractBC2D, u_interior::Vector{Float64}, t::Float64)::Vector{Float64}
    error("apply_ghost not implemented for $(typeof(bc))")
end

